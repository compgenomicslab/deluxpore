#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Levenshtein import distance


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
    parser = argparse.ArgumentParser(prog='05-calculate_uniq_index_text_similarity.py', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='05-calculate_uniq_index_tesxt_similarity.py uses levenshtein distance to calculate the similarity'
                                     'between sequences read unique index sequences and illumina index unique sequences. It compares each sequence '
                                     'against a set of unique index sequences and their reverse complements, printing the minimum distance for each '
                                     'sequence.')
    
    parser.add_argument('--read_unique_indexes', '-ru', required=True,
                        help='Path to read unique index sequences in FASTA format')

    parser.add_argument('--project_unique_indexes', '-iu', required=True,
                        help='Path to project specific unique index sequences in FASTA format')

    parser.add_argument('--project_unique_indexes_rc', '-iur', required=True,
                        help='Path to project specific unique index reverse complement sequences in FASTA format')

    parser.add_argument('--output', '-o', required=True,
                        help='Output tsv file with distance matrix')

    return parser.parse_args()

#################
### FUNCTIONS ###
#################

# def iterative_levenshtein(s, t):
#     """ 
#         iterative_levenshtein(s, t) -> ldist
#         ldist is the Levenshtein distance between the strings 
#         s and t.
#         For all i and j, dist[i,j] will contain the Levenshtein 
#         distance between the first i characters of s and the 
#         first j characters of t
#     """

#     rows = len(s)+1
#     cols = len(t)+1
#     dist = [[0 for x in range(cols)] for x in range(rows)]

#     # source prefixes can be transformed into empty strings 
#     # by deletions:
#     for i in range(1, rows):
#         dist[i][0] = i

#     # target prefixes can be created from an empty source string
#     # by inserting the characters
#     for i in range(1, cols):
#         dist[0][i] = i

#     for col in range(1, cols):
#         for row in range(1, rows):
#             if s[row-1] == t[col-1]:
#                 cost = 0
#             else:
#                 cost = 1
#             dist[row][col] = min(dist[row-1][col] + 1,      # deletion
#                                  dist[row][col-1] + 1,      # insertion
#                                  dist[row-1][col-1] + cost) # substitution

#     return dist[row][col]

# Function to parse a FASTA file and return a dictionary mapping headers to sequences
def parse_fasta(fasta_file):
    head2seqs = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        head2seqs[record.id] = str(record.seq)  # Storing each sequence in the dictionary using its header as the key
    return head2seqs


# Function to calculate and print distances of a sequence to each index
def index_distances(record, uniq_index, uniq_index_rc, index_order):
    distances = {}
    read_name = record.id
    read_seq = str(record.seq)
    distances[read_name] = {}

    for index in index_order:
        dist = distance(read_seq, uniq_index[index])  # Distance to forward index
        dist_rc = distance(read_seq, uniq_index_rc[index])  # Distance to reverse complement index
        distances[read_name][index] = min(dist, dist_rc)  # Store the minimum of the two distances

    # Print the distances for the given read_name, tab-separated
    distances = (read_name + "\t" + "\t".join([str(distances[read_name][i]) for i in index_order]) + "\n")

    return distances


if __name__ == '__main__':
    args = check_arg()

    # Parsing the index files to create dictionaries of sequences
    uniq_index = parse_fasta(args.project_unique_indexes)        # Dictionary of normal index sequences
    uniq_index_rc = parse_fasta(args.project_unique_indexes_rc)  # Dictionary of reverse complement index sequences

    # Creating a list of index headers in the order they appear
    index_order = []
    for index in uniq_index:
        index_order.append(index)
    
    with open(args.output, 'w') as out_file:
        out_file.write("query_id\t" + "\t".join([i for i in index_order]) + "\n")  # Write header line with index names
        
        # Iterate over each record in the fasta file containing adapters and calculate distances
        for record in SeqIO.parse(args.read_unique_indexes, 'fasta'):
            distances = index_distances(record, uniq_index, uniq_index_rc, index_order)  # Calculate distances for each sequence
            out_file.write(distances)

