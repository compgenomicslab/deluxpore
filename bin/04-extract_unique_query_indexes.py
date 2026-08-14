#!/usr/bin/env python3
"""
Efficient extraction of unique index sequences from BLAST alignments
Handles both i5 and i7 indexes with forward and reverse strand alignments
"""

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
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
    parser = argparse.ArgumentParser(prog='04-extract_uniq_index_seqs.py', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='04-extract_uniq_index_seqs.py extract unique index sequences from the previously parsed' \
                                     'query index sequences. As the position of the unique index sequence in the complete index sequence'
                                     'is fixed, a simple logical step is used to extract the unique index sequences.')
    
    parser.add_argument('--input', '-i', required=True,
                        help='Path to BLAST output file')

    parser.add_argument('--reads', '-r', required=True,
                        help='Path to complete query sequences in FASTA format')

    parser.add_argument('--complete_indexes_fna', '-ic', required=True,
                        help='Path to project complete index sequences (adapter + barcode) '
                             'in FASTA format (used, together with --unique_indexes_fna, to '
                             'locate the barcode position within the adapter)')

    parser.add_argument('--unique_indexes_fna', '-iu', required=True,
                        help='Path to project unique index sequences in FASTA format '
                             '(used, together with --complete_indexes_fna, to locate the '
                             'barcode position within the adapter)')

    parser.add_argument('--output', '-o', required=True,
                        help='Output fasta file name to write exact unique index sequence from complete query index sequences')

    return parser.parse_args()

#################
### FUNCTIONS ###
#################

# Populated at runtime by determine_index_positions() -- see __main__ below.
# {'i7': {'start': int, 'end': int}, 'i5': {'start': int, 'end': int}}
# start/end are 0-based positions of the unique barcode within the complete
# adapter+barcode sequence (exclusive end), in subject/template coordinates.
INDEX_POSITIONS = {}


def determine_index_positions(complete_indexes_fna, unique_indexes_fna):
    """Locate where each slot's unique barcode sits within its complete
    adapter+barcode sequence, for this project's index kit.

    Adapter chemistry places the barcode at a fixed offset for every index
    within a given slot (i5/i7) of a kit, but that offset, and the barcode
    length, differs between kits (NEBNext vs NEXTERA vs custom). Rather
    than hardcoding a single kit's positions, this locates the unique
    barcode as an exact substring of its complete sequence (matched by
    index ID) for every index, and uses the (start, length) pair shared by
    the largest number of indexes per slot.

    Returns {'i5': {'start': int, 'end': int}, 'i7': {'start': int, 'end': int}}.
    """
    complete = SeqIO.to_dict(SeqIO.parse(complete_indexes_fna, 'fasta'))

    slot_windows = {'i5': {}, 'i7': {}}  # (start, end) -> number of indexes agreeing
    for record in SeqIO.parse(unique_indexes_fna, 'fasta'):
        index_id = record.id
        slot = get_index_type(index_id)
        if slot is None or index_id not in complete:
            continue

        full_seq = str(complete[index_id].seq)
        unique_seq = str(record.seq)
        start = full_seq.find(unique_seq)
        if start == -1:
            print(f"Warning: unique index '{index_id}' was not found as a substring "
                  f"of its complete index sequence -- skipping it for position detection",
                  file=sys.stderr)
            continue

        window = (start, start + len(unique_seq))
        slot_windows[slot][window] = slot_windows[slot].get(window, 0) + 1

    positions = {}
    for slot, windows in slot_windows.items():
        if not windows:
            raise ValueError(
                f"Could not determine the barcode position for slot '{slot}': no unique "
                f"index sequence was found within its complete index sequence. Check that "
                f"{complete_indexes_fna} and {unique_indexes_fna} use matching index IDs "
                f"and that the barcode is a substring of the complete sequence."
            )
        if len(windows) > 1:
            print(f"Warning: inconsistent {slot} barcode positions detected across indexes "
                  f"in this kit: {windows} (window -> index count) -- using the most common one",
                  file=sys.stderr)
        best_window = max(windows.items(), key=lambda kv: kv[1])[0]
        positions[slot] = {'start': best_window[0], 'end': best_window[1]}

    return positions


def parse_blast_line(line):
    """Parse BLAST output line and return relevant fields in 0-based indexing when appropriate"""
    fields = line.strip().split('\t')
    return {
        'qseqid': fields[0],    # Query sequence ID
        'sseqid': fields[1],    # Subject sequence ID
        'qstart': int(fields[6])-1, # Query alignment start
        'qend': int(fields[7])-1,   # Query alignment end
        'sstart': int(fields[8])-1, # Subject alignment start
        'send': int(fields[9])-1    # Subject alignment end
    }

def get_index_type(subject_id):
    """Determine index type from subject sequence ID"""
    if subject_id.startswith('i7'):
        return 'i7'
    elif subject_id.startswith('i5'):
        return 'i5'
    return None


def retrieve_unique_positions_subject(index_type):
    """Calculate the positions of unique sequence in subject coordinates"""
    if index_type not in INDEX_POSITIONS:
        return None, None
    
    pos_info = INDEX_POSITIONS[index_type]
    unique_start = pos_info['start']
    unique_end = pos_info['end']
    
    return unique_start, unique_end


def is_unique_region_covered(sstart, send, unique_start, unique_end):
    """Check if the unique sequence region is fully covered by alignment"""
    align_start = min(sstart, send)
    align_end = max(sstart, send)
    
    return align_start <= unique_start and align_end >= unique_end

def extract_unique_sequence(query_seq, blast_data, index_type):
    """
    Extract unique index sequence from query based on BLAST alignment

    Args:
        query_seq: Full query sequence string
        blast_data: Parsed BLAST line data
        index_type: 'i5' or 'i7'

    Returns:
        Extracted unique sequence or None if extraction fails
    """

    # Get alignment coordinates already in 0-based indexing
    qstart = blast_data['qstart']
    qend = blast_data['qend']
    sstart = blast_data['sstart']
    send = blast_data['send']

    # Determine unique sequence positions in subject
    unique_start_subj, unique_end_subj = retrieve_unique_positions_subject(index_type)
    if unique_start_subj is None:
        return None
    # Slice width comes directly from this slot's detected window, rather than
    # a separate length constant that could silently drift out of sync with it.
    index_length = unique_end_subj - unique_start_subj

    # Check if unique region is covered by alignment
    if not is_unique_region_covered(sstart, send, unique_start_subj, unique_end_subj):
        return None

    # Check minimum subject alignment length
    if abs(send - sstart) < 20:
        return None

    subject_reverse = sstart > send
    if subject_reverse:
        offset_from_align_end = sstart - unique_end_subj
        unique_start_query = qstart + offset_from_align_end
        unique_end_query = unique_start_query + index_length
        extracted_seq = query_seq[unique_start_query:unique_end_query]
    else:
        offset_from_align_start = unique_start_subj - sstart
        unique_start_query = qstart + offset_from_align_start
        unique_end_query = unique_start_query + index_length
        extracted_seq = query_seq[unique_start_query:unique_end_query]

    # Validate extracted sequence length
    if len(extracted_seq) != index_length:
        return None

    return extracted_seq

def process_blast_output(blast_file, fasta_file, output_file):
    """
    Process BLAST output and extract unique index sequences

    Args:
        blast_file: Path to BLAST output file
        fasta_file: Path to query FASTA file
        output_file: Optional output file path
    """

    # Load query sequences
    query_sequences = {}
    # for record in SeqIO.parse(fasta_file, 'fasta'):
    #     query_sequences[record.id] = str(record.seq)
    query_sequences = SeqIO.index(fasta_file, 'fasta')

    # Process BLAST output
    results = []

    with open(blast_file, 'r') as f:
        for line in f:
            if line.strip():
                try:
                    blast_data = parse_blast_line(line) #returns a dictionary with parsed fields per line in 0-based indexing when appropriate
                    index_type = get_index_type(blast_data['sseqid'])
                    
                    if index_type and blast_data['qseqid'] in query_sequences:
                        # query_seq = query_sequences[blast_data['qseqid']]
                        query_seq = str(query_sequences[blast_data['qseqid']].seq)
                        unique_seq = extract_unique_sequence(query_seq, blast_data, index_type)
                        
                        if unique_seq:
                            results.append({
                                'query_id': blast_data['qseqid'],
                                'subject_id': blast_data['sseqid'], 
                                'index_type': index_type,
                                'unique_sequence': unique_seq,
                                'alignment_info': f"q{blast_data['qstart']}-{blast_data['qend']}_s{blast_data['sstart']}-{blast_data['send']}"
                            })
                            
                except Exception as e:
                    print(f"Warning: Failed to process line: {line.strip()}", file=sys.stderr)
                    print(f"Error: {e}", file=sys.stderr)
                    continue
    
    with open(output_file, 'w') as out_file:
        writer = FastaWriter(out_file, wrap=None)
        for result in results:
            final_seq = SeqRecord(
                seq=result["unique_sequence"],
                id=f"{result['query_id']}.{result['subject_id']}.{result['alignment_info']}",
                description=""
            )

            writer.write_record(final_seq)

    return results


if __name__ == "__main__":
    args = check_arg()

    # Barcode position and length are derived from the actual complete/unique
    # index sequences in use, rather than a single hardcoded window -- works
    # for NEBNext, NEXTERA, and custom kits alike, each of which places the
    # barcode at a different offset within its complete adapter sequence.
    INDEX_POSITIONS = determine_index_positions(args.complete_indexes_fna, args.unique_indexes_fna)
    print(f"Detected index positions: {INDEX_POSITIONS}", file=sys.stderr)

    results = process_blast_output(args.input, args.reads, args.output)
    print(f"Processed {len(results)} successful extractions", file=sys.stderr)
