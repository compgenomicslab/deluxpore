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
    parser = argparse.ArgumentParser(prog='07-parse_best_and_demultiplex.py', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='07-parse_best_and_demultiplex.py ')

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

    parser.add_argument('--rc_collision_withhold_dist', '-rcw', type=int, default=0,
                        help='Reads whose best RC collision distance is <= this value are '
                             'excluded from their sample FASTA and reported as '
                             'rc_collision_excluded in ambiguous_reads.tsv. '
                             '0 (default) = exclude only exact RC mirrors (dist=0); '
                             '-1 = flag only, never exclude.')

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

# Maximum accepted Levenshtein distance for a same-slot barcode match (8bp
# barcodes). At <=3, ~44% of matched reads turned out to be invalid_index_pair
# noise -- borderline dist-2/3 matches on one or both slots, not confident
# barcode reads. <=2 cuts that to ~14% while only losing ~4% of total reads
# (measured on real data); <=1 goes further (~4% invalid pairs) at the cost of
# losing about half the reads. <=2 is the current best noise-vs-coverage
# trade-off found.
MAX_BARCODE_MATCH_DIST = 2


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

                if min_val > MAX_BARCODE_MATCH_DIST:
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
                            f"{current_sample}|{new_sample}",
                            "excluded"
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
                            f"{current_sample}|{new_sample}",
                            "excluded"
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

def write_fasta_files_per_sample(grouped_data, chunkID, exp_des_dict, reads, output_path, ambiguous_events, rc_suppressed=None, assignment_counts=None):
    """Write per-sample FASTA files.

    rc_suppressed: optional set of read names that have an RC collision but no
    unambiguous sample assignment.  These reads are flagged in ambiguous_reads.tsv
    (by append_rc_collision_ambiguous_events, called before this function) but must
    not be written to any sample FASTA — they have no quorum.

    assignment_counts: optional dict mutated in place with counts of how reads
    matched -- 'both' (i5+i7 confirmed), 'i5_only', 'i7_only', or 'unassigned'
    (no valid/unique index found at all, including invalid i5+i7 pairs and
    single-barcode multi-sample ties). This reflects the barcode-matching
    outcome only, independent of RC-collision exclusion -- a read excluded by
    rc_suppressed still counts under whichever path it matched via, since RC
    collision is an orthogonal downstream QC filter on top of a valid match.
    """
    if assignment_counts is None:
        assignment_counts = defaultdict(int)

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
            assignment_path = "unassigned"  # overwritten below on a successful match

            #Case 1: both indexes present
            if i5 != "" and i7 != "":
                index_pair = [i5, i7]
                for sample_name, sample_indexes in exp_des_dict.items():
                    if sample_indexes == index_pair:
                        assigned_sample = sample_name
                        assignment_path = "both"

                if not assigned_sample:
                    # Both barcodes matched confidently, but not to any real sample
                    # pair -- e.g. index hopping or a chimeric read. Previously this
                    # fell through completely unreported; now flagged so it's visible
                    # in ambiguous_reads.tsv instead of silently vanishing into
                    # "unassigned" with no trace.
                    ambiguous_events.append((
                        read, "invalid_index_pair",
                        f"i5={i5}+i7={i7} (no sample uses this combination)",
                        "unassigned",
                        "excluded"
                    ))

            #Case 2: only i5 present
            elif i5 != "" and i7 == "":
                if len(i5_to_sample[i5]) == 1: #only one sample linked to the index
                    assigned_sample = i5_to_sample[i5][0]
                    assignment_path = "i5_only"
                else:
                    #i5 appears in more than one sample, therefore inconclusive
                    ambiguous_events.append((
                        read, "single_barcode_multi_sample",
                        f"i5={i5} (no i7 found)",
                        "|".join(i5_to_sample[i5]),
                        "excluded"
                    ))

            #Case 3: only i7 present
            elif i5 == "" and i7 != "":
                if len(i7_to_sample[i7]) == 1: #only one sample linked to the index
                    assigned_sample = i7_to_sample[i7][0]
                    assignment_path = "i7_only"
                else:
                    #i7 appears in more than one sample, therefore inconclusive
                    ambiguous_events.append((
                        read, "single_barcode_multi_sample",
                        f"i7={i7} (no i5 found)",
                        "|".join(i7_to_sample[i7]),
                        "excluded"
                    ))

            assignment_counts[assignment_path] += 1

            if assigned_sample:
                # If this read's only reason to be assigned is an RC-colliding barcode
                # combination that doesn't unambiguously point to a sample, skip the
                # FASTA — it is already flagged in ambiguous_reads.tsv.
                if rc_suppressed and read in rc_suppressed:
                    print(f"'{read}' suppressed from FASTA: RC collision with no quorum")
                    continue
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

    return assignment_counts



def load_rc_collision_events(rc_collision_file):
    """Return list of (read_name, extraction_slot, colliding_index, actual_distance).

    The query_id in the TSV is the full bin/04 record ID:
    '{read_name}.{subject_id}.{alignment_info}' — we extract just the read_name
    (first dot-separated field) so it matches the keys in grouped_data.
    """
    events = []
    with open(rc_collision_file, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 4:
                query_id, extraction_slot, colliding_index, actual_distance = parts
                read_name = query_id.split('.')[0]
                events.append((read_name, extraction_slot, colliding_index, actual_distance))
    return events


def append_rc_collision_ambiguous_events(rc_collision_events, grouped_data, exp_des_dict,
                                         ambiguous_events, withhold_dist=-1):
    """Flag RC collision reads in ambiguous_events; return the set of excluded read names.

    Always reported as ambiguity_type "rc_collision"; the decision column
    ("included"/"excluded") carries the outcome, and barcode_info notes
    "dual_confirmed" when i5+i7 independently confirm a valid pair (in which
    case it's always "included", regardless of withhold_dist).

    Withholding only ever applies to single-barcode assignments (i5_only /
    i7_only). A read whose i5 AND i7 independently confirm a valid,
    self-consistent sample pair is never withheld, even if one of its
    barcodes also collides with an unrelated sample's barcode in RC
    orientation: the slot itself (i5 vs i7) is established by BLAST-aligning
    ~60bp of adapter backbone, not the 8bp barcode, so an 8bp cross-slot
    coincidence doesn't actually undermine a match corroborated by a second,
    independently-extracted barcode. Single-barcode assignments have no such
    corroboration, so those remain eligible for withholding.

    Deduplication: bin/05 records one row per BLAST alignment × colliding index.
    We keep the best (lowest) distance per (read_name, slot, colliding_index), then
    reduce further to the globally minimum distance per read. Events above
    RC_COLLISION_MAX_DIST are dropped; only the minimum-distance event(s) are reported.
    """
    # Deduplicate: keep best (lowest) distance per (read_name, extraction_slot, colliding_index)
    best = {}
    for read_name, extraction_slot, colliding_index, actual_distance in rc_collision_events:
        key = (read_name, extraction_slot, colliding_index)
        dist = int(actual_distance)
        if key not in best or dist < best[key]:
            best[key] = dist

    # Find the global minimum distance per read; only report at that minimum and within cap.
    # bin/05 emits RC collisions with dist <= 3. For reporting we require <= RC_COLLISION_MAX_DIST
    # (exact or near-exact RC mirror); anything above is not a genuine collision risk.
    RC_COLLISION_MAX_DIST = 1

    read_min_dist = {}
    for (read_name, extraction_slot, colliding_index), dist in best.items():
        if read_name not in read_min_dist or dist < read_min_dist[read_name]:
            read_min_dist[read_name] = dist

    suppressed = set()

    for (read_name, extraction_slot, colliding_index), dist in best.items():
        if dist > read_min_dist[read_name]:
            continue  # not the closest RC match for this read
        if dist > RC_COLLISION_MAX_DIST:
            continue  # closest match still too far to be a genuine RC collision
        if read_name not in grouped_data:
            continue
        data = grouped_data[read_name]
        assigned_i5, assigned_i7 = data.i5, data.i7
        if not assigned_i5 and not assigned_i7:
            continue  # read was not assigned at all

        dual_confirmed = False  # both i5 and i7 independently confirm a valid, real sample pair

        if assigned_i5 and assigned_i7:
            assigned_sample = next(
                (s for s, idxs in exp_des_dict.items() if idxs == [assigned_i5, assigned_i7]),
                "unassigned"
            )
            dual_confirmed = (assigned_sample != "unassigned")
        elif assigned_i5:
            matches = [s for s, idxs in exp_des_dict.items() if idxs[0] == assigned_i5]
            assigned_sample = matches[0] if len(matches) == 1 else "unassigned"
        else:
            matches = [s for s, idxs in exp_des_dict.items() if idxs[1] == assigned_i7]
            assigned_sample = matches[0] if len(matches) == 1 else "unassigned"

        excluded = (
            not dual_confirmed
            and withhold_dist >= 0
            and read_min_dist[read_name] <= withhold_dist
        )
        if excluded:
            suppressed.add(read_name)

        # ambiguity_type is always "rc_collision" -- the decision column carries
        # the outcome:
        #   included — i5+i7 independently confirm a valid pair (dual_confirmed),
        #              so it's kept in its sample FASTA regardless of withhold_dist;
        #              or a single-barcode collision below withhold_dist.
        #   excluded — single-barcode assignment with no second barcode to
        #              corroborate, and dist <= withhold_dist.
        decision = "excluded" if excluded else "included"
        ambiguous_events.append((
            read_name, "rc_collision",
            f"{extraction_slot}_extracted; RC_matches_{colliding_index} (dist={dist})"
            f"{'; dual_confirmed' if dual_confirmed else ''}",
            assigned_sample,
            decision
        ))

    return suppressed

def write_ambiguous_fasta_files(ambiguous_events, reads, chunkID, output_path):
    ambiguous_dir = f'{output_path}/ambiguous'
    os.makedirs(ambiguous_dir, exist_ok=True)

    handlers = {
        'tie_both_valid':              open(f'{ambiguous_dir}/tie_both_valid.{chunkID}.fna', 'w'),
        'single_barcode_multi_sample': open(f'{ambiguous_dir}/single_barcode_multi_sample.{chunkID}.fna', 'w'),
        'invalid_index_pair':          open(f'{ambiguous_dir}/invalid_index_pair.{chunkID}.fna', 'w'),
    }
    try:
        for read_id, ambiguity_type, _, _, _ in ambiguous_events:
            if ambiguity_type in handlers and read_id in reads:
                new_seq = SeqRecord(reads[read_id].seq, id=reads[read_id].id, description="")
                SeqIO.write(new_seq, handlers[ambiguity_type], 'fasta')
    finally:
        for fh in handlers.values():
            fh.close()


def write_ambiguous_report(ambiguous_events, chunkID, output_path):
    report_path = f'{output_path}/ambiguous_reads.{chunkID}.tsv'
    with open(report_path, 'w') as f:
        f.write("read_id\tambiguity_type\tbarcode_info\tpossible_samples\tdecision\n")
        for event in ambiguous_events:
            f.write("\t".join(str(x) for x in event) + "\n")
    print(f"Ambiguous reads report written to: {report_path} ({len(ambiguous_events)} events)")


def write_index_assignment_summary(assignment_counts, chunkID, output_path):
    """Write a one-row-per-chunk summary of how reads were matched: via a
    confirmed i5+i7 pair ('both'), a single unambiguous barcode ('i5_only' /
    'i7_only'), or not at all ('unassigned'). This reflects the barcode
    matching outcome only -- independent of RC-collision exclusion, which is
    a separate downstream QC filter (see write_fasta_files_per_sample)."""
    both = assignment_counts.get("both", 0)
    i5_only = assignment_counts.get("i5_only", 0)
    i7_only = assignment_counts.get("i7_only", 0)
    unassigned = assignment_counts.get("unassigned", 0)
    total = both + i5_only + i7_only + unassigned

    report_path = f'{output_path}/index_assignment_summary.{chunkID}.tsv'
    with open(report_path, 'w') as f:
        f.write("chunkID\tboth\ti5_only\ti7_only\tunassigned\ttotal_reads\n")
        f.write(f"{chunkID}\t{both}\t{i5_only}\t{i7_only}\t{unassigned}\t{total}\n")
    print(f"Index assignment summary written to: {report_path} "
          f"(both={both}, i5_only={i5_only}, i7_only={i7_only}, unassigned={unassigned})")


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

    # RC collision events must be processed BEFORE writing FASTA files so we can
    # build the suppressed-read set and pass it to write_fasta_files_per_sample.
    rc_events = load_rc_collision_events(args.rc_collision_events)
    withhold_dist = args.rc_collision_withhold_dist
    rc_suppressed = append_rc_collision_ambiguous_events(
        rc_events, mapped_data, exp_des_dict, ambiguous_events, withhold_dist
    )
    if rc_suppressed:
        print(f"{len(rc_suppressed)} RC collision reads excluded from sample FASTA "
              f"(best collision dist <= {withhold_dist})")

    assignment_counts = write_fasta_files_per_sample(mapped_data, chunkID, exp_des_dict, reads_dict, args.output, ambiguous_events, rc_suppressed)

    write_ambiguous_fasta_files(ambiguous_events, reads_dict, chunkID, args.output)
    write_ambiguous_report(ambiguous_events, chunkID, args.output)
    write_index_assignment_summary(assignment_counts, chunkID, args.output)

