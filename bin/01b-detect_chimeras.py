#!/usr/bin/env python3
"""
Pre-demultiplexing chimera detection and splitting.

Runs on raw (post-Nanopore-filter, pre-demultiplexing) reads, per chunk.
Detects internal (non-terminal) BLAST hits against the complete Illumina
index database that cover >= --min_chimera_coverage of the adapter
template -- these mark genuine chimeric junctions (two ligated fragments,
each carrying its own i5/i7 pair), as opposed to the single terminal i5 and
terminal i7 adapter every normal read has.

Detection always runs and is always reported. Splitting only happens when
--remove_chimeras is true; otherwise chimeric reads are passed through
unchanged (still flagged in the report) so demultiplexing sees the original
read.

Splitting happens BEFORE demultiplexing (unlike terminal adapter trimming,
which happens after -- see bin/09-trim_illumina_indexes.py) so each fragment
carries only its own sample's i5/i7 pair into extraction, distance matching
and assignment. Left unsplit, a chimeric read's i5 (from one ligated
fragment) and i7 (from another) get paired into an invalid cross-sample
combination and the whole read is excluded.
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import csv
import sys
from collections import defaultdict


def check_arg(args=None):
    parser = argparse.ArgumentParser(prog='01b-detect_chimeras.py')

    parser.add_argument('--fasta_reads', '-i', required=True,
                        help='Raw per-chunk FASTA (pre-demultiplexing)')
    parser.add_argument('--blast_output', '-b', required=True,
                        help='BLAST -outfmt 6 output (raw reads vs complete index DB)')
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
    parser.add_argument('--min_fragment_length', '-ml', type=int, default=100,
                        help='Minimum length for a split fragment to be kept '
                             '[default: 100]')
    parser.add_argument('--output', '-o', required=True,
                        help='Output FASTA path (chimera-checked reads: split into '
                             'fragments where a chimera was detected and removal is on, '
                             'unchanged otherwise)')
    parser.add_argument('--report', '-r', required=True,
                        help='Output TSV report (one row per detected chimera fragment)')

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


def find_junctions(read_len, hits, template_lens, min_chimera_coverage):
    """Return internal (non-terminal) hits that qualify as chimera junctions.

    A hit counts as terminal -- the normal i5/i7 adapter every read has -- if
    it starts within the first 73bp or ends within the last 73bp of the read
    (adapter templates are ~66-71bp, so a hit starting/ending near a
    boundary can never extend far past it). Everything else is a candidate
    internal junction, gated by --min_chimera_coverage to reject short
    coincidental matches.

    Returns list of (qstart_1based, qend_1based, coverage), sorted by qstart.
    """
    junctions = []
    for qstart, qend, sseqid, sstart, send in hits:
        if qstart < 73 or qend > read_len - 73:
            continue  # terminal adapter -- not a junction
        tlen = template_lens.get(sseqid, 66) #66 is the fallback
        cov = chimera_coverage(sstart, send, tlen)
        if cov >= min_chimera_coverage:
            junctions.append((qstart, qend, cov))
    junctions.sort(key=lambda x: x[0])
    return junctions


def split_at_junctions(read_seq, junctions):
    """Split a read into fragments at each junction.

    No terminal trimming happens here -- that is handled separately, after
    demultiplexing, by bin/09-trim_illumina_indexes.py.

    Returns list of (fragment_seq, orig_start_0based, orig_end_0based_excl).
    """
    read_len = len(read_seq)
    fragments = []
    cursor = 0
    for qstart, qend, cov in junctions:
        j_start = qstart - 1  # 0-based
        j_end   = qend        # 0-based exclusive
        if j_start <= cursor or j_end > read_len:
            continue  # junction outside remaining sequence; skip
        frag = read_seq[cursor:j_start] #first iteration would be 0-first junction qstart (jstart in 1-based)
        if len(frag) > 0:
            fragments.append((frag, cursor, j_start))
        cursor = max(cursor, j_end) #the evaluation goes to j_end at the end of the first junction

    tail = read_seq[cursor:]
    if len(tail) > 0:
        fragments.append((tail, cursor, cursor + len(tail)))

    return fragments if fragments else [(read_seq, 0, read_len)]


if __name__ == '__main__':
    args = check_arg()

    template_lens = load_template_lengths(args.complete_indexes_fna)
    blast_hits = parse_blast_hits(args.blast_output)
    reads = SeqIO.index(args.fasta_reads, 'fasta')

    report_rows = []
    dropped_short = 0

    with open(args.output, 'w') as out_fasta:
        for read_id in reads:
            record = reads[read_id]
            hits = blast_hits.get(read_id, [])

            junctions = find_junctions(
                len(record.seq), hits, template_lens, args.min_chimera_coverage
            ) if hits else []

            if not junctions:
                # No chimera junction -- pass the read through unchanged.
                SeqIO.write(SeqRecord(record.seq, id=read_id, description=""), out_fasta, 'fasta')
                continue

            if not args.remove_chimeras:
                # Detected but removal is off -- report it, but still pass the
                # original read through unchanged so demultiplexing sees it.
                best_cov = max(cov for _, _, cov in junctions)
                report_rows.append({
                    'original_read_id': read_id,
                    'num_fragments': 1,
                    'fragment_id': read_id,
                    'fragment_start': 0,
                    'fragment_end': len(record.seq),
                    'junction_coverage': f"{best_cov:.3f}",
                    'chimera_split': 'no',
                })
                SeqIO.write(SeqRecord(record.seq, id=read_id, description=""), out_fasta, 'fasta')
                continue

            fragments = split_at_junctions(record.seq, junctions)
            kept = [f for f in fragments if len(f[0]) >= args.min_fragment_length]
            dropped_short += len(fragments) - len(kept)

            if len(fragments) == 1:
                # All junctions fell outside the usable sequence -- nothing split.
                if kept:
                    frag_seq, _, _ = kept[0]
                    SeqIO.write(SeqRecord(frag_seq, id=read_id, description=""), out_fasta, 'fasta')
                continue

            best_cov = max(cov for _, _, cov in junctions)
            for idx, (frag_seq, orig_s, orig_e) in enumerate(kept, start=1):
                frag_id = f"{read_id}_frag{idx}"
                SeqIO.write(SeqRecord(frag_seq, id=frag_id, description=""), out_fasta, 'fasta')
                report_rows.append({
                    'original_read_id': read_id,
                    'num_fragments': len(kept),
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
          f"({n_split} split into fragments). Dropped {dropped_short} fragments "
          f"below --min_fragment_length={args.min_fragment_length}.",
          file=sys.stderr)
