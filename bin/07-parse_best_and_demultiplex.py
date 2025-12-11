#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
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
    
    parser.add_argument('--index_distance_table', '-id', required=True,
                        help='Path to read index to illumina index distance table in tsv format')

    parser.add_argument('--experimental_design', '-ed', required=True,
                        help='Path to experimental design file in tsv format')

    parser.add_argument('--output', '-o', required=True,
                        help='Output path for files')

    return parser.parse_args()

#################
### FUNCTIONS ###
#################


def parse_exp_design(file_path):
    exp_des_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            info = line.strip().split('\t')
            i5 = info[1]
            i7 = info[2]
            exp_des_dict[info[0]] = [i5, i7]
    return exp_des_dict

def validate_index_pairs(exp_des_dict, index_pair):
    return index_pair in exp_des_dict.values()

def create_best_distance_dict(distance_file):
    best = defaultdict(list)
    with open(distance_file, 'r') as tsv:
        for line in tsv:
            if line.strip().startswith("query_id"):
                header = line
                header_list = header.strip().split('\t')
                print(header_list)
            else:
                info = line.strip().split('\t')
                query_id = info[0]
                distances = list(map(int, info[1:]))
                min_value_index = distances.index(min(distances))
                min_value = distances[min_value_index]

                best_match = header_list[min_value_index + 1]
                if min_value <= 3:
                    best[query_id] = [min_value, best_match]
        return best


def parse_best_dictionary_should_update(best_dict, exp_des_dict):
    GroupData = namedtuple('GroupData', ['i5', 'i5_dist', 'i5_qstart', 'i7', 'i7_dist', 'i7_qstart'])
    grouped_data = {}
    for key, (min_value, index_name) in best_dict.items():
        print(key)
        read_name, _, align_info = key.split('.')  # split key name into read name and the potential best index it was mapped to
        qstart = int(align_info.split('-')[0].replace('q', ''))

        if read_name not in grouped_data:  # check if the read was already checked
            # initialize read if not seen before
            grouped_data[read_name] = GroupData(
                i5='', i5_dist=float('inf'), i5_qstart=float('inf'), 
                i7='', i7_dist=float('inf'), i7_qstart=float('inf'))
        current = grouped_data[read_name]
        
        if index_name.startswith("i5"):
            should_update = False

            # Priority: 1) lower distance, 2) alignment closer to start, 3) valid combo
            if min_value < current.i5_dist:
                should_update = True
            elif min_value == current.i5_dist:
                # Same distance → prefer alignment closer to start (lower qstart)
                if qstart < current.i5_qstart:
                    should_update = True
                    print(f"Same distance, but closer to start ({qstart} < {current.i5_qstart}), updating i5")
                elif qstart == current.i5_qstart and index_name != current.i5:
                    # Same position, different index → check valid combo
                    print(read_name)
                    print(f"min value '{min_value}' is the same as the current min value '{current.i5_dist}' " 
                        f"and new i5 '{index_name}' is different to current i5 '{current.i5}'")
                    # Handle tie case - check if this combination is valid
                    current_combo = [current.i5, current.i7] if current.i5 and current.i7 else None
                    new_combo = [index_name, current.i7] if current.i7 else None
                    # Prefer valid combinations over invalid ones
                    current_is_valid = validate_index_pairs(exp_des_dict, current_combo)
                    new_is_valid = validate_index_pairs(exp_des_dict, new_combo)

                    if new_is_valid and not current_is_valid:
                        should_update = True
                        print(f"update, new combination '{new_combo}' is valid and current '{current_combo}' is not")
                        print("updated i5 value for read:", read_name)

                    elif current_is_valid and not new_is_valid:
                        # Current is valid
                        should_update = False
                        print(f"do not update, combination '{current_combo}' is valid and '{new_combo}' is not")

                    elif new_is_valid and current_is_valid:
                        grouped_data[read_name] = GroupData(i5='', i5_dist=float('inf'), i5_qstart=float('inf'), i7='', i7_dist=float('inf'), i7_qstart=float('inf'))
                        print(read_name, " could be attributed to two different samples, reverting to empty record")


            if should_update:
                grouped_data[read_name] = GroupData(
                    i5=index_name,
                    i5_dist=min_value,
                    i5_qstart=qstart,
                    i7=current.i7,
                    i7_dist=current.i7_dist,
                    i7_qstart=current.i7_qstart
                )

        elif index_name.startswith("i7"):
            should_update = False
            if min_value < current.i7_dist:
                should_update = True
            elif min_value == current.i7_dist:
                # Same distance → prefer alignment closer to end (higher qstart)
                if qstart > current.i7_qstart:
                    should_update = True
                    print(f"Same distance, but closer to end ({qstart} > {current.i7_qstart}), updating i7")
                elif qstart == current.i7_qstart and index_name != current.i7:
                    # Same position, different index → check valid combo
                    print(read_name)
                    print(f"min value '{min_value}' is the same as the current min value '{current.i7_dist}' " 
                        f"and new i7 '{index_name}' is different to current i7 '{current.i7}'")
                    # Handle tie case - check if this combination is valid
                    current_combo = [current.i5, current.i7] if current.i5 and current.i7 else None
                    new_combo = [current.i5, index_name] if current.i5 else None
                    print(current_combo, new_combo)
                    # Prefer valid combinations over invalid ones
                    current_is_valid = validate_index_pairs(exp_des_dict, current_combo)
                    new_is_valid = validate_index_pairs(exp_des_dict, new_combo)
                    print(current_is_valid, new_is_valid)
                    if new_is_valid and not current_is_valid:
                        should_update = True
                        print(f"update, new combination '{new_combo}' is valid and current '{current_combo}' is not")
                        print("updated i7 value for read:", read_name)
                        
                    elif current_is_valid and not new_is_valid:
                        # Current is valid
                        should_update = False
                        print(f"do not update, combination '{current_combo}' is valid and '{new_combo}' is not")

                    elif new_is_valid and current_is_valid:
                        grouped_data[read_name] = GroupData(i5='', i5_dist=float('inf'), i5_qstart=float('inf'), i7='', i7_dist=float('inf'), i7_qstart=float('inf'))
                        print(read_name, "could be attributed to two different samples, reverting to empty record")

            if should_update:
                grouped_data[read_name] = GroupData(
                    i5=current.i5,
                    i5_dist=current.i5_dist,
                    i5_qstart=current.i5_qstart,
                    i7=index_name,
                    i7_dist=min_value,
                    i7_qstart=qstart
                )
    return grouped_data

def count_inconclusive(grouped_data):
    count = 0
    for key in grouped_data:
        if grouped_data[key].i5_dist == float('inf') and grouped_data[key].i7_dist == float('inf'):
            print(key)
            count = count + 1
    print(f"A total of {count} reads result in inconclusive mapping")

def write_info_into_file(grouped_data, chunkID, output_path):
    #Write namedtuple dict to json file for tracability purposes
    serializable_grouped_data = {
    key: value._asdict() for key, value in grouped_data.items()
        }

    # Write to JSON
    with open(f'{output_path}/grouped_data.{chunkID}.json', 'w') as f:
        json.dump(serializable_grouped_data, f, indent=4)

def write_fasta_files_per_sample(grouped_data, chunkID, exp_des_dict, reads, output_path):

    i5_to_sample = {}
    i7_to_sample = {}

    for sample, indexes in exp_des_dict.items():
        i5, i7 = indexes
        if i5 not in i5_to_sample:
            i5_to_sample[i5] = []
        i5_to_sample[indexes[0]].append(sample)
        if i7 not in i7_to_sample:
            i7_to_sample[i7] = []
        i7_to_sample[i7].append(sample)

    per_sample_chunk_output_file = { 
        sample:open(f'{output_path}/{sample}.{chunkID}.fna', 'w')
        for sample in exp_des_dict.keys()
    }
    
    try:
        for read in grouped_data:
            i5 = grouped_data[read].i5
            i7 = grouped_data[read].i7

            assigned_sample = None

            #Case 1: both indexes present
            if i5 != "" and i7 != "":
                index_pair = [i5, i7]
                for sample_name, sample_indexes in exp_des_dict.items():
                    if sample_indexes == index_pair:
                        assigned_sample = sample_name
                
            #Case 2: only i5 present
            elif i5 != "" and i7 == "":
                if len(i5_to_sample[i5]) == 1: #only one sample linked to the index
                    assigned_sample = i5_to_sample[i5][0]
                else:
                    #i5 appears in more that one sample, therefore inconclusive
                    pass

            #Case 3: only i7 present
            elif i5 == "" and i7 != "":
                if len(i7_to_sample[i7]) == 1: #only one sample linked to the index
                    assigned_sample = i7_to_sample[i7][0]
                else:
                    #i7 appears in more that one sample, therefore inconclusive
                    pass
            
            if assigned_sample:
                new_seq = SeqRecord(
                    reads[read].seq, 
                    id=reads[read].id, 
                    description=""
                )
                SeqIO.write(new_seq, per_sample_chunk_output_file[assigned_sample], 'fasta')
                print(f"'{read}' was sucessfully demultiplexed")
                
    finally:
        for file_handle in per_sample_chunk_output_file.values():
            file_handle.close()



if __name__ == "__main__":
    args = check_arg()

    distance_table = args.index_distance_table
    clean_reads_file = args.fasta_reads
    # chunkID = (os.path.splitext(os.path.basename(distance_table))[0]).split(".")[0]
    chunkID = os.path.basename(distance_table).replace(".distance_matrix.tsv", "")

    # with open(clean_reads_file, 'r') as reads_file:
    #     reads_dict = SeqIO.to_dict(SeqIO.parse(reads_file, 'fasta'))
    reads_dict = SeqIO.index(clean_reads_file, 'fasta')

    exp_des_dict = parse_exp_design(args.experimental_design)
    best_dict = create_best_distance_dict(distance_table)
    mapped_data = parse_best_dictionary_should_update(best_dict, exp_des_dict)
    count_inconclusive(mapped_data)
    write_info_into_file(mapped_data, chunkID, args.output)

    write_fasta_files_per_sample(mapped_data, chunkID, exp_des_dict, reads_dict, args.output)

