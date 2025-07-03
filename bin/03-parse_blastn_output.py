#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
import re
import argparse


def check_arg(args=None):
    '''
        Description:
        Function to collect arguments from command line using argparse
    Input:
        args # command line arguments
    Constant:
        None
    Variables
        parser
    Return
        parser.parse_args() # Parsed arguments
    '''
    parser = argparse.ArgumentParser(prog='03-parse_blastn_output.py', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='03-parse_blastn_output.py parses blastn reads to index database results to extract' \
                                     'complete index sequences from trimmed reads and their corresponding coordinates in the read')
    
    parser.add_argument('--input', '-i', required=True,
                        help='Path to BLAST output file')

    parser.add_argument('--reads', '-r', required=True,
                        help='Path to trimmed reads fasta file')

    parser.add_argument('--output', '-o', required=True,
                        help='Output fasta file name to write exact index sequence from trimmed reads')

    return parser.parse_args()

#################
### FUNCTIONS ###
#################

# Parse BLASTN output file, read trimmed reads and modify headers
def parse_blast_output_write_fasta(blast_out, reads_dict, reads_out):
    with open(blast_out, 'r') as blast_file, open(reads_out, 'w') as out_file:
        writer = FastaWriter(out_file, wrap=None)
        for line in blast_file:
            parts = line.strip().split('\t')
            read_id = parts[0]
            index_id = parts[1]
            qstart = int(parts[6]) - 1
            qend = int(parts[7])

            full_seq = reads_dict[read_id]

            extracted_seq = SeqRecord(
                seq=full_seq.seq[qstart:qend],
                id=f"{full_seq.id}.{index_id}.{qstart+1}.{qend}",
                description=""
            )

            writer.write_record(extracted_seq) #write record to output file

if __name__ == '__main__':
    args = check_arg()

    # Read input sequences to dictionary
    with open(args.reads, 'r') as read_seqs:
        reads_dict = SeqIO.to_dict(SeqIO.parse(read_seqs, "fasta"))

    # Parse BLAST output and write new fasta file containing index sequences for each read
    parse_blast_output_write_fasta(args.input, reads_dict, args.output)
