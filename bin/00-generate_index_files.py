#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import namedtuple

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
    parser = argparse.ArgumentParser(prog='00-generate_index_files.py', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='00-generate_index_files.py reads index ids from experimental design file and generates two' \
                                     'files: complete_index_ids.faa and unique_index.faa. These will be specific of the experimental design')
    
    parser.add_argument('--input', '-i', required=True,
                        help='Path to experimental design file in tsv format')

    parser.add_argument('--indexKit', '-ik', required=True,
                        help='index sequences kit name used in library construction')
    
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory to write the results')

    return parser.parse_args()

#################
### FUNCTIONS ###
#################

element = namedtuple('seq_id', ['i5', 'i7'])

def read_index_seqs_into_dict(indexkit):
    with open(f"/home/acatalina/PHD/long_nanopore/run_nextflow_demultiplex/nextflow_demultiplex/assets/{indexkit}.complete_indexes.fna", "r") as complete_indexes:
        with open(f"/home/acatalina/PHD/long_nanopore/run_nextflow_demultiplex/nextflow_demultiplex/assets/{indexkit}.unique_indexes.fna", "r") as unique_indexes:
            complete_dict = SeqIO.to_dict(SeqIO.parse(complete_indexes, "fasta"))
            unique_dict = SeqIO.to_dict(SeqIO.parse(unique_indexes, "fasta"))

    return complete_dict, unique_dict

def read_experimental_design(input_file):
    seq_indexes = {}
    
    with open(input_file, "r") as file:
        lines = file.readlines()

        for line in lines:
            info = line.strip().split(" ")
            sample = info[0]
            i5 = str(info[1])
            i7 = str(info[2])

            seq_indexes[sample] = element(i5=i5, i7=i7)

    return seq_indexes

def create_file_handlers(output):
    # Create file handlers using context manager
    files = {
        "complete": open(f'{output}/complete_index_seqs.fna', 'w'),
        "unique": open(f'{output}/unique_index_seqs.fna', 'w'),
        "rc": open(f'{output}/unique_index_seqs_rc.fna', 'w')
        }
    return files


def write_index_files(seq_indexes, complete_dict, unique_dict, output_dir):

    file_handlers = create_file_handlers(output_dir)

    for seq in seq_indexes:
        print(seq)
        i5 = seq_indexes[seq].i5
        i7 = seq_indexes[seq].i7

        i5_seq_complete = SeqRecord(complete_dict[i5].seq, id=complete_dict[i5].id, description="")
        i7_seq_complete = SeqRecord(complete_dict[i7].seq, id=complete_dict[i7].id, description="")

        i5_seq_unique = SeqRecord(unique_dict[i5].seq, id=unique_dict[i5].id, description="")
        i7_seq_unique = SeqRecord(unique_dict[i7].seq, id=unique_dict[i7].id, description="")

        i5_seq_unique_rc = SeqRecord((unique_dict[i5].seq).reverse_complement(), id=unique_dict[i5].id, description="")
        i7_seq_unique_rc = SeqRecord((unique_dict[i7].seq).reverse_complement(), id=unique_dict[i7].id, description="")


        SeqIO.write(i5_seq_complete, file_handlers["complete"], 'fasta')
        SeqIO.write(i7_seq_complete, file_handlers["complete"], 'fasta')

        SeqIO.write(i5_seq_unique, file_handlers["unique"], 'fasta')
        SeqIO.write(i7_seq_unique, file_handlers["unique"], 'fasta')

        SeqIO.write(i5_seq_unique_rc, file_handlers["rc"], 'fasta')
        SeqIO.write(i7_seq_unique_rc, file_handlers["rc"], 'fasta')

if __name__ == '__main__':
    args = check_arg()

    #Parse index seq files and store info into ditionary
    complete_dict, unique_dict = read_index_seqs_into_dict(args.indexKit)

    my_seq_indexes = read_experimental_design(args.input)
    # print(my_seq_indexes)
    write_index_files(my_seq_indexes, complete_dict, unique_dict, args.output)


