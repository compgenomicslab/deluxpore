#!/usr/bin/env python3
"""
Post-demultiplexing Illumina-adapter trimming and optional chimera splitting.

Runs per-sample on the already-demultiplexed FASTA. Terminal adapter
trimming matches the behaviour of bin/06 (73bp window). When --remove_chimeras
is true, internal BLAST hits that cover >= --min_chimera_coverage of the
adapter template are treated as genuine chimeric junctions and the read is
split into fragments at each one.
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import csv
import os
import sys
from collections import defaultdict


def check_arg(args=None):
    parser = argparse.ArgumentParser(prog='09-trim_and_remove_chimeras.py')

    parser.add_argument('--fasta_reads', '-i', required=True,
                        help='Per-sample FASTA produced by demultiplexing')
    parser.add_argument('--blast_output', '-b', required=True,
                        help='BLAST -outfmt 6 output (sample reads vs complete index DB)')
    parser.add_argument('--complete_indexes_fna', '-c', required=True,
                        help='Complete adapter+barcode FASTA used to build the BLAST DB '
                             '(needed to look up template lengths for coverage calculation)')
    parser.add_argument('--remove_chimeras', type=lambda x: x.lower() == 'true',
                        default=False,
                        help='Split reads at confident internal adapter occurrences '
                             '[default: False]')
    parser.add_argument('--min_chimera_coverage', type=float, default=0.7,
                        help='Minimum fraction of the adapter template an internal '
                             'alignment must cover to be treated as a genuine chimeric '
                             'junction rather than noise [default: 0.7]')
    parser.add_argument('--output', '-o', required=True,
                        help='Output FASTA path')
    parser.add_argument('--report', '-r', required=True,
                        help='Output TSV report (one row per split read)')

    return parser.parse_args()


def load_template_lengths(complete_indexes_fna):
    """Return {index_id: template_length_bp} from the complete-adapter FASTA."""
    lengths = {}
    for record in SeqIO.parse(complete_indexes_fna, 'fasta'):
        lengths[record.id] = len(record.seq)
    return lengths


def parse_blast_hits(blast_file):
    """Group BLAST hits by query read ID.

    Returns {read_id: [(qstart_1based, qend_1based, sseqid, sstart_1based, send_1based), ...]}
    All coordinates are 1-based inclusive (native BLAST tabular format).
    """
    hits = defaultdict(list)
    with open(blast_file) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            read_id = parts[0]
            sseqid  = parts[1]
            qstart  = int(parts[6])
            qend    = int(parts[7])
            sstart  = int(parts[8])
            send    = int(parts[9])
            hits[read_id].append((qstart, qend, sseqid, sstart, send))
    return hits


def chimera_coverage(sstart, send, template_len):
    """Subject alignment span / template length (both 1-based inclusive coords)."""
    if template_len == 0:
        return 0.0
    return (abs(send - sstart) + 1) / template_len


def trim_and_split(read_seq, hits, template_lens, remove_chimeras, min_chimera_coverage):
    """
    Apply terminal trimming and optional chimera splitting to a single read.

    All BLAST coordinates are 1-based inclusive (as BLAST outputs them).

    Chimera junctions are always detected (for reporting); splitting only happens
    when remove_chimeras=True.

    Returns:
        (fragments, chimera_junctions)
        fragments: list of (fragment_seq, orig_start_0based, orig_end_0based_excl)
                   Single-element list when no split occurred.
        chimera_junctions: list of (qstart_1based, qend_1based, coverage) for
                           every qualifying internal hit, regardless of remove_chimeras.
    """
    read_len = len(read_seq)

    start_trim_qends = []   # 1-based qend of start-terminal hits
    end_trim_qstarts = []   # 1-based qstart of end-terminal hits
    chimera_junctions = []  # (qstart_1based, qend_1based, coverage) — always detected

    for qstart, qend, sseqid, sstart, send in hits:
        if qstart < 73 and qend <= 73:
            start_trim_qends.append(qend)
        elif qend > read_len - 73:
            end_trim_qstarts.append(qstart)
        else:
            tlen = template_lens.get(sseqid, 66)
            cov = chimera_coverage(sstart, send, tlen)
            if cov >= min_chimera_coverage:
                chimera_junctions.append((qstart, qend, cov))

    # Terminal trimming. BLAST coords are 1-based; Python slicing uses 0-based.
    # seq[qend:]      keeps everything after the start-terminal adapter (qend is
    #                 1-based inclusive last adapter base → 0-based start of kept region)
    # seq[:qstart-1]  keeps everything before the end-terminal adapter (qstart-1
    #                 is 0-based exclusive slice end = last kept genomic position + 1)
    T_start = max(start_trim_qends) if start_trim_qends else 0  # 0-based start of trimmed seq
    T_end   = min(end_trim_qstarts) - 1 if end_trim_qstarts else read_len  # 0-based exclusive end

    trimmed = read_seq[T_start:T_end]
    if len(trimmed) == 0:
        trimmed = read_seq  # safety: don't produce an empty record
        T_start = 0

    if not chimera_junctions or not remove_chimeras:
        return [(trimmed, T_start, T_start + len(trimmed))], chimera_junctions

    # Chimera splitting. Convert junction coords to trimmed-sequence space.
    # j_start_trimmed = (qstart_1based - 1) - T_start  (0-based in trimmed)
    # j_end_trimmed   = qend_1based - T_start           (0-based excl in trimmed)
    chimera_junctions.sort(key=lambda x: x[0])
    fragments = []
    cursor = 0
    for qstart, qend, cov in chimera_junctions:
        j_start = (qstart - 1) - T_start
        j_end   = qend - T_start
        if j_start <= 0 or j_end > len(trimmed):
            continue  # junction entirely in a trimmed region; skip
        frag = trimmed[cursor:j_start]
        if len(frag) > 0:
            orig_s = cursor + T_start
            orig_e = j_start + T_start
            fragments.append((frag, orig_s, orig_e))
        cursor = max(cursor, j_end)

    tail = trimmed[cursor:]
    if len(tail) > 0:
        orig_s = cursor + T_start
        orig_e = orig_s + len(tail)
        fragments.append((tail, orig_s, orig_e))

    # If all junctions were in trimmed regions, fall back to single record
    return (fragments if fragments else [(trimmed, T_start, T_start + len(trimmed))]), chimera_junctions


if __name__ == '__main__':
    args = check_arg()

    template_lens = load_template_lengths(args.complete_indexes_fna)
    blast_hits = parse_blast_hits(args.blast_output)
    reads = SeqIO.index(args.fasta_reads, 'fasta')

    report_rows = []

    with open(args.output, 'w') as out_fasta:
        for read_id, record in reads.items():
            hits = blast_hits.get(read_id, [])

            if not hits:
                # No BLAST hits for this read — write it unchanged
                SeqIO.write(SeqRecord(record.seq, id=read_id, description=""), out_fasta, 'fasta')
                continue

            fragments, chimera_junctions = trim_and_split(
                record.seq, hits, template_lens,
                args.remove_chimeras, args.min_chimera_coverage
            )

            if not chimera_junctions:
                # No internal adapter hits — write trimmed read as-is
                frag_seq, _, _ = fragments[0]
                SeqIO.write(SeqRecord(frag_seq, id=read_id, description=""), out_fasta, 'fasta')
            elif len(fragments) == 1:
                # Chimera detected but not split (remove_chimeras=False)
                frag_seq, orig_s, orig_e = fragments[0]
                SeqIO.write(SeqRecord(frag_seq, id=read_id, description=""), out_fasta, 'fasta')
                best_cov = max(cov for _, _, cov in chimera_junctions)
                report_rows.append({
                    'original_read_id': read_id,
                    'num_fragments': 1,
                    'fragment_id': read_id,
                    'fragment_start': orig_s,
                    'fragment_end': orig_e,
                    'junction_coverage': f"{best_cov:.3f}",
                    'chimera_split': 'no',
                })
            else:
                # Chimera detected and split
                best_cov = max(cov for _, _, cov in chimera_junctions)
                for idx, (frag_seq, orig_s, orig_e) in enumerate(fragments, start=1):
                    frag_id = f"{read_id}_frag{idx}"
                    SeqIO.write(SeqRecord(frag_seq, id=frag_id, description=""), out_fasta, 'fasta')
                    report_rows.append({
                        'original_read_id': read_id,
                        'num_fragments': len(fragments),
                        'fragment_id': frag_id,
                        'fragment_start': orig_s,
                        'fragment_end': orig_e,
                        'junction_coverage': f"{best_cov:.3f}",
                        'chimera_split': 'yes',
                    })

    with open(args.report, 'w', newline='') as report_f:
        writer = csv.DictWriter(report_f, fieldnames=[
            'original_read_id', 'num_fragments', 'fragment_id',
            'fragment_start', 'fragment_end', 'junction_coverage', 'chimera_split'
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(report_rows)

    n_detected = len({r['original_read_id'] for r in report_rows})
    n_split = len({r['original_read_id'] for r in report_rows if r['chimera_split'] == 'yes'})
    print(f"Processed {len(reads)} reads. Detected {n_detected} chimeric reads "
          f"({n_split} split into fragments).",
          file=sys.stderr)