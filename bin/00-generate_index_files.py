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

    parser.add_argument('--completeIndexes', '-ci', required=True,
                        help='Path to complete index sequences FASTA file (adapter + barcode)')

    parser.add_argument('--uniqueIndexes', '-ui', required=True,
                        help='Path to unique barcode-only index sequences FASTA file')

    parser.add_argument('--output', '-o', required=True,
                        help='Output directory to write the results')

    return parser.parse_args()

#################
### FUNCTIONS ###
#################

element = namedtuple('seq_id', ['i5', 'i7'])

def read_index_seqs_into_dict(complete_path, unique_path):
    with open(complete_path, "r") as complete_indexes:
        with open(unique_path, "r") as unique_indexes:
            complete_dict = SeqIO.to_dict(SeqIO.parse(complete_indexes, "fasta"))
            unique_dict = SeqIO.to_dict(SeqIO.parse(unique_indexes, "fasta"))
    return complete_dict, unique_dict

def read_experimental_design(input_file):
    # seq_indexes = {}
    i5_indexes = set()
    i7_indexes = set()
    with open(input_file, "r") as file:
        lines = file.readlines()

        for line in lines:
            info = line.strip().split("\t")
            sample = info[0]
            i5 = str(info[1])
            i7 = str(info[2])

            i5_indexes.add(i5)
            i7_indexes.add(i7)
            # seq_indexes[sample] = element(i5=i5, i7=i7)

    return i5_indexes, i7_indexes

def create_file_handlers(output):
    files = {
        "complete": open(f'{output}/complete_index_seqs.fna', 'w'),
        "unique": open(f'{output}/unique_index_seqs.fna', 'w'),
        "rc": open(f'{output}/unique_index_seqs_rc.fna', 'w')
        }
    return files


def write_index_files(i5_indexes, i7_indexes, complete_dict, unique_dict, output_dir):

    file_handlers = create_file_handlers(output_dir)

    for i5 in i5_indexes:
        i5_seq_complete = SeqRecord(complete_dict[i5].seq, id=complete_dict[i5].id, description="")
        i5_seq_unique = SeqRecord(unique_dict[i5].seq, id=unique_dict[i5].id, description="")
        i5_seq_unique_rc = SeqRecord((unique_dict[i5].seq).reverse_complement(), id=unique_dict[i5].id, description="")

        SeqIO.write(i5_seq_complete, file_handlers["complete"], 'fasta')
        SeqIO.write(i5_seq_unique, file_handlers["unique"], 'fasta')
        SeqIO.write(i5_seq_unique_rc, file_handlers["rc"], 'fasta')
    
    for i7 in i7_indexes:
        i7_seq_complete = SeqRecord(complete_dict[i7].seq, id=complete_dict[i7].id, description="")
        i7_seq_unique = SeqRecord(unique_dict[i7].seq, id=unique_dict[i7].id, description="")
        i7_seq_unique_rc = SeqRecord((unique_dict[i7].seq).reverse_complement(), id=unique_dict[i7].id, description="")

        SeqIO.write(i7_seq_complete, file_handlers["complete"], 'fasta')
        SeqIO.write(i7_seq_unique, file_handlers["unique"], 'fasta')
        SeqIO.write(i7_seq_unique_rc, file_handlers["rc"], 'fasta')
    
    # for seq in seq_indexes:
    #     print(seq)
    #     i5 = seq_indexes[seq].i5
    #     i7 = seq_indexes[seq].i7

    #     i5_seq_complete = SeqRecord(complete_dict[i5].seq, id=complete_dict[i5].id, description="")
    #     i7_seq_complete = SeqRecord(complete_dict[i7].seq, id=complete_dict[i7].id, description="")

    #     i5_seq_unique = SeqRecord(unique_dict[i5].seq, id=unique_dict[i5].id, description="")
    #     i7_seq_unique = SeqRecord(unique_dict[i7].seq, id=unique_dict[i7].id, description="")

    #     i5_seq_unique_rc = SeqRecord((unique_dict[i5].seq).reverse_complement(), id=unique_dict[i5].id, description="")
    #     i7_seq_unique_rc = SeqRecord((unique_dict[i7].seq).reverse_complement(), id=unique_dict[i7].id, description="")


    #     SeqIO.write(i5_seq_complete, file_handlers["complete"], 'fasta')
    #     SeqIO.write(i7_seq_complete, file_handlers["complete"], 'fasta')

    #     SeqIO.write(i5_seq_unique, file_handlers["unique"], 'fasta')
    #     SeqIO.write(i7_seq_unique, file_handlers["unique"], 'fasta')

    #     SeqIO.write(i5_seq_unique_rc, file_handlers["rc"], 'fasta')
    #     SeqIO.write(i7_seq_unique_rc, file_handlers["rc"], 'fasta')

if __name__ == '__main__':
    args = check_arg()

    #Parse index seq files and store info into dictionary
    complete_dict, unique_dict = read_index_seqs_into_dict(args.completeIndexes, args.uniqueIndexes)

    my_i5_indexes, my_i7_indexes = read_experimental_design(args.input)
    # print(my_seq_indexes)
    write_index_files(my_i5_indexes, my_i7_indexes, complete_dict, unique_dict, args.output)


