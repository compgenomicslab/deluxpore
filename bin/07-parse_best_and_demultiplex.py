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

    parser.add_argument('--rc_collision_events', '-rc', required=True,
                        help='Path to RC collision events TSV from bin/05')

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
                min_val = min(distances)

                if min_val > 3:
                    continue

                # Collect all columns that tied at the minimum distance.
                # With slot-restricted matching (bin/05), cross-slot columns are
                # set to 999, so ties here only occur for same-slot candidates
                # under normal operation. The cross-slot check below catches the
                # fallback case where slot could not be determined from the record
                # ID and all columns were compared (slot = None in bin/05).
                min_indices = [i for i, d in enumerate(distances) if d == min_val]
                min_names = [header_list[i + 1] for i in min_indices]

                if len(min_names) > 1:
                    has_i5 = any(n.startswith('i5') for n in min_names)
                    has_i7 = any(n.startswith('i7') for n in min_names)
                    if has_i5 and has_i7:
                        # Cross-slot tie: the slot could not be determined and
                        # two different slots match equally well. Slot is
                        # unresolvable from sequence alone; skip this extraction.
                        print(f"Cross-slot tie for {query_id}: {min_names} all at distance {min_val}, skipping")
                        continue

                # Single winner or same-slot tie: take the first minimum.
                # Same-slot ties within a read are resolved downstream by
                # parse_best_dictionary_should_update using qstart and valid-combo checks.
                min_value_index = distances.index(min_val)
                best_match = header_list[min_value_index + 1]
                best[query_id] = [min_val, best_match]
    return best


def parse_best_dictionary_should_update(best_dict, exp_des_dict, ambiguous_events):
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
                        current_sample = next((s for s, idxs in exp_des_dict.items() if idxs == current_combo), "unknown")
                        new_sample = next((s for s, idxs in exp_des_dict.items() if idxs == new_combo), "unknown")
                        ambiguous_events.append((
                            read_name, "tie_both_valid",
                            f"i5={current.i5}+i7={current.i7} OR i5={index_name}+i7={current.i7}",
                            f"{current_sample}|{new_sample}"
                        ))


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
                        current_sample = next((s for s, idxs in exp_des_dict.items() if idxs == current_combo), "unknown")
                        new_sample = next((s for s, idxs in exp_des_dict.items() if idxs == new_combo), "unknown")
                        ambiguous_events.append((
                            read_name, "tie_both_valid",
                            f"i5={current.i5}+i7={current.i7} OR i5={current.i5}+i7={index_name}",
                            f"{current_sample}|{new_sample}"
                        ))

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

def write_fasta_files_per_sample(grouped_data, chunkID, exp_des_dict, reads, output_path, ambiguous_events):

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
                    #i5 appears in more than one sample, therefore inconclusive
                    ambiguous_events.append((
                        read, "single_barcode_multi_sample",
                        f"i5={i5} (no i7 found)",
                        "|".join(i5_to_sample[i5])
                    ))

            #Case 3: only i7 present
            elif i5 == "" and i7 != "":
                if len(i7_to_sample[i7]) == 1: #only one sample linked to the index
                    assigned_sample = i7_to_sample[i7][0]
                else:
                    #i7 appears in more than one sample, therefore inconclusive
                    ambiguous_events.append((
                        read, "single_barcode_multi_sample",
                        f"i7={i7} (no i5 found)",
                        "|".join(i7_to_sample[i7])
                    ))
            
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



def load_rc_collision_events(rc_collision_file):
    events = []
    with open(rc_collision_file, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 4:
                events.append(parts)  # [read_id, extraction_slot, colliding_index, actual_distance]
    return events

def append_rc_collision_ambiguous_events(rc_collision_events, grouped_data, exp_des_dict, ambiguous_events):
    # For each RC collision event where the read was actually assigned, add a
    # reporting entry. The read is NOT excluded — it already went to its correct
    # sample. This entry exists only so users can inspect RC collision cases.
    idx_to_samples = {}
    for sample, (i5, i7) in exp_des_dict.items():
        for idx in (i5, i7):
            idx_to_samples.setdefault(idx, []).append(sample)

    seen = set()
    for read_id, extraction_slot, colliding_index, actual_distance in rc_collision_events:
        if read_id in seen:
            continue
        seen.add(read_id)
        if read_id not in grouped_data:
            continue
        data = grouped_data[read_id]
        assigned_i5, assigned_i7 = data.i5, data.i7
        if not assigned_i5 and not assigned_i7:
            continue  # read was not assigned at all; skip
        colliding_samples = idx_to_samples.get(colliding_index, ['unknown'])
        ambiguous_events.append((
            read_id, "rc_collision",
            f"assigned={assigned_i5}+{assigned_i7} | rc_match={colliding_index}(dist={actual_distance},slot={extraction_slot})",
            "|".join(colliding_samples)
        ))

def write_ambiguous_fasta_files(ambiguous_events, reads, chunkID, output_path):
    ambiguous_dir = f'{output_path}/ambiguous'
    os.makedirs(ambiguous_dir, exist_ok=True)

    handlers = {
        'tie_both_valid':              open(f'{ambiguous_dir}/tie_both_valid.{chunkID}.fna', 'w'),
        'single_barcode_multi_sample': open(f'{ambiguous_dir}/single_barcode_multi_sample.{chunkID}.fna', 'w'),
        'rc_collision':                open(f'{ambiguous_dir}/rc_collision.{chunkID}.fna', 'w'),
    }
    try:
        for read_id, ambiguity_type, _, _ in ambiguous_events:
            if ambiguity_type in handlers and read_id in reads:
                new_seq = SeqRecord(reads[read_id].seq, id=reads[read_id].id, description="")
                SeqIO.write(new_seq, handlers[ambiguity_type], 'fasta')
    finally:
        for fh in handlers.values():
            fh.close()


def write_ambiguous_report(ambiguous_events, chunkID, output_path):
    report_path = f'{output_path}/ambiguous_reads.{chunkID}.tsv'
    with open(report_path, 'w') as f:
        f.write("read_id\tambiguity_type\tbarcode_info\tpossible_samples\n")
        for event in ambiguous_events:
            f.write("\t".join(str(x) for x in event) + "\n")
    print(f"Ambiguous reads report written to: {report_path} ({len(ambiguous_events)} events)")


if __name__ == "__main__":
    args = check_arg()

    distance_table = args.index_distance_table
    clean_reads_file = args.fasta_reads
    chunkID = os.path.basename(distance_table).replace(".distance_matrix.tsv", "")

    reads_dict = SeqIO.index(clean_reads_file, 'fasta')

    exp_des_dict = parse_exp_design(args.experimental_design)
    best_dict = create_best_distance_dict(distance_table)
    ambiguous_events = []
    mapped_data = parse_best_dictionary_should_update(best_dict, exp_des_dict, ambiguous_events)
    count_inconclusive(mapped_data)
    write_info_into_file(mapped_data, chunkID, args.output)

    write_fasta_files_per_sample(mapped_data, chunkID, exp_des_dict, reads_dict, args.output, ambiguous_events)

    rc_events = load_rc_collision_events(args.rc_collision_events)
    append_rc_collision_ambiguous_events(rc_events, mapped_data, exp_des_dict, ambiguous_events)

    write_ambiguous_fasta_files(ambiguous_events, reads_dict, chunkID, args.output)
    write_ambiguous_report(ambiguous_events, chunkID, args.output)

