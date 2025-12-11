#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import re
import json
import argparse
from collections import defaultdict, namedtuple


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
    parser = argparse.ArgumentParser(prog='06-extract_best_index.py', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='06-extract_best_index.py ')

    parser.add_argument('--fasta_reads', '-i', required=True,
                        help='Path to reads file in fasta format')

    parser.add_argument('--complete_index_reads', '-ir', required=True,
                        help='Path to complete index query sequences in fasta format')

    parser.add_argument('--output', '-o', required=True,
                        help='Output path for files')

    return parser.parse_args()

#################
### FUNCTIONS ###
#################
def parse_complete_index_read_seqs(comp_index_reads):
    parsed_comp_index_reads = defaultdict(lambda: {})

    for read_id in comp_index_reads:
        read_name = read_id.split('.')[0]  # Get the read name without any additional info
        read_positions = '.'.join(read_id.split('.')[2:4])
        
        parsed_comp_index_reads[read_name][read_positions] = str(comp_index_reads[read_id].seq)
        
    return parsed_comp_index_reads

def trim_illumina_indexes_short(parsed_comp_index_reads, reads_dict, chunkID, output_path):
    trimmed_reads = {}
    for read_id in parsed_comp_index_reads.keys():

        if read_id not in reads_dict:
                continue  # skip reads not found in reads_dict
            
        seq_to_trim = reads_dict[read_id].seq

        index_seqs_per_read = parsed_comp_index_reads[read_id]
            
        # Collect trim positions
        start_trim_positions = []
        end_trim_positions = []

        for positions, seq in index_seqs_per_read.items():
            start, end = map(int, positions.split('.'))

            if start < 73 and end <= 73:
                # Trimming at the start
                start_trim_positions.append(end)
            elif start >= len(seq_to_trim) - 73 and end <= len(seq_to_trim):
                # Trimming at the end
                end_trim_positions.append(start)

            # If in the middle, ignore

        trimmed_seq = seq_to_trim

        # Trim at the end first if needed
        if end_trim_positions:
            trim_end_pos = min(end_trim_positions)
            trimmed_seq = trimmed_seq[:trim_end_pos]

        # Then trim at the start
        if start_trim_positions:
            trim_start_pos = max(start_trim_positions)
            trimmed_seq = trimmed_seq[trim_start_pos:]

        # Safety check for empty sequences
        if trimmed_seq is None or len(trimmed_seq) == 0:
            trimmed_seq = seq_to_trim

        trimmed_reads[read_id] = trimmed_seq
    
    # Write trimmed reads to output file
    output_file = os.path.join(output_path, f"{chunkID}.NI_trimm.fna")
    with open(output_file, 'w') as out_f:
        for read_id, trimmed_seq in trimmed_reads.items():
            seq_record = SeqRecord(trimmed_seq, id=read_id, description="")
            SeqIO.write(seq_record, out_f, 'fasta')


if __name__ == "__main__":
    args = check_arg()

    nanopore_trimmed_reads = args.fasta_reads #fasta file with nanopore trimmed reads
    complete_index_reads = args.complete_index_reads #fasta file with complete illumina indexes for each read (all possible indexes based on the mappings)
    # chunkID = (os.path.splitext(os.path.basename(nanopore_trimmed_reads))[0]).split(".")[0]
    chunkID = os.path.basename(nanopore_trimmed_reads).replace(".nano_trimmed_filtered.fasta", "")
    
    # # load nanopore trimmed reads into dictionary
    # with open(nanopore_trimmed_reads, 'r') as reads_file:
    #     reads_dict = SeqIO.to_dict(SeqIO.parse(reads_file, 'fasta'))

    # # load complete index sequences for each read into dictionary
    # with open(complete_index_reads, 'r') as comp_index_file:
    #     comp_index_reads = SeqIO.to_dict(SeqIO.parse(comp_index_file, 'fasta'))
    
    reads_dict = SeqIO.index(nanopore_trimmed_reads, 'fasta')
    comp_index_reads = SeqIO.index(complete_index_reads, 'fasta')

    parsed_complete_index_read_dicts = parse_complete_index_read_seqs(comp_index_reads)

    trim_illumina_indexes_short(parsed_complete_index_read_dicts, reads_dict, chunkID, args.output)