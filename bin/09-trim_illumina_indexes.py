#!/usr/bin/env python3
"""
Post-demultiplexing terminal Illumina-adapter trimming.

Runs per-sample on the already-demultiplexed FASTA. Trims the terminal i5/i7
adapter (and everything before/after it) using the same 73bp terminal-window
logic as bin/06 (production's pre-demux trimming). Chimera detection and
splitting happens separately and earlier, in bin/01b-detect_chimeras.py,
before demultiplexing -- this script only removes leftover adapter sequence
from the ends of each already-correctly-assigned read.
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import sys


def check_arg(args=None):
    parser = argparse.ArgumentParser(prog='09-trim_illumina_indexes.py')

    parser.add_argument('--fasta_reads', '-i', required=True,
                        help='Per-sample FASTA produced by demultiplexing')
    parser.add_argument('--blast_output', '-b', required=True,
                        help='BLAST -outfmt 6 output (sample reads vs complete index DB)')
    parser.add_argument('--min_fragment_length', '-ml', type=int, default=100,
                        help='Minimum length for a terminally-trimmed read to be written '
                             'to the output FASTA [default: 100]')
    parser.add_argument('--output', '-o', required=True,
                        help='Output FASTA path')

    return parser.parse_args()


def parse_blast_hits(blast_file):
    """Group BLAST hits by query read ID. Returns {read_id: [(qstart_1based, qend_1based), ...]}"""
    hits = {}
    with open(blast_file) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            read_id = parts[0]
            qstart  = int(parts[6])
            qend    = int(parts[7])
            hits.setdefault(read_id, []).append((qstart, qend))
    return hits


def terminal_trim(read_seq, hits):
    """Trim the terminal i5/i7 adapter (and everything before/after it).

    A hit counts as terminal if it starts within the first 73bp or ends
    within the last 73bp of the read -- adapter templates are ~66-71bp, so a
    hit starting/ending near a boundary can never extend far past it.
    """
    read_len = len(read_seq)
    start_trim_qends = []
    end_trim_qstarts = []

    for qstart, qend in hits:
        if qstart < 73:
            start_trim_qends.append(qend)
        elif qend > read_len - 73:
            end_trim_qstarts.append(qstart)

    T_start = max(start_trim_qends) if start_trim_qends else 0
    T_end   = min(end_trim_qstarts) - 1 if end_trim_qstarts else read_len

    trimmed = read_seq[T_start:T_end]
    if len(trimmed) == 0:
        trimmed = read_seq  # safety: don't produce an empty record
    return trimmed


if __name__ == '__main__':
    args = check_arg()

    blast_hits = parse_blast_hits(args.blast_output)
    reads = SeqIO.index(args.fasta_reads, 'fasta')

    dropped_short = 0

    with open(args.output, 'w') as out_fasta:
        for read_id in reads:
            record = reads[read_id]
            hits = blast_hits.get(read_id, [])

            trimmed = terminal_trim(record.seq, hits) if hits else record.seq

            if len(trimmed) >= args.min_fragment_length:
                SeqIO.write(SeqRecord(trimmed, id=read_id, description=""), out_fasta, 'fasta')
            else:
                dropped_short += 1

    print(f"Processed {len(reads)} reads. Dropped {dropped_short} reads "
          f"below --min_fragment_length={args.min_fragment_length}.",
          file=sys.stderr)
