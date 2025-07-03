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

    parser.add_argument('--output', '-o', required=True,
                        help='Output fasta file name to write exact unique index sequence from complete query index sequences')

    return parser.parse_args()

#################
### FUNCTIONS ###
#################

# Index position constants
INDEX_POSITIONS = {
    'i7': {'start': 24, 'end': 32},  # 8-base unique sequence at positions 24-32
    'i5': {'start': 29, 'end': 37}   # 8-base unique sequence at positions 29-37
}

INDEX_LENGTH = 8


def parse_blast_line(line):
    """Parse BLAST output line and return relevant fields"""
    fields = line.strip().split('\t')
    return {
        'qseqid': fields[0],    # Query sequence ID
        'sseqid': fields[1],    # Subject sequence ID  
        'qstart': int(fields[6]), # Query alignment start
        'qend': int(fields[7]),   # Query alignment end
        'sstart': int(fields[8]), # Subject alignment start
        'send': int(fields[9])    # Subject alignment end
    }

def get_index_type(subject_id):
    """Determine index type from subject sequence ID"""
    if subject_id.startswith('i7'):
        return 'i7'
    elif subject_id.startswith('i5'):
        return 'i5'
    return None


def calculate_unique_positions_subject(index_type):
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

    # Get alignment coordinates
    qstart = blast_data['qstart']
    qend = blast_data['qend']
    sstart = blast_data['sstart']
    send = blast_data['send']

    # Determine unique sequence positions in subject
    unique_start_subj, unique_end_subj = calculate_unique_positions_subject(index_type)
    if unique_start_subj is None:
        return None

    # Check if unique region is covered by alignment
    if not is_unique_region_covered(sstart, send, unique_start_subj, unique_end_subj):
        return None

    subject_reverse = sstart > send
    if subject_reverse:
        offset_from_align_end = sstart - unique_end_subj
        unique_start_query = qstart + offset_from_align_end
        unique_end_query = unique_start_query + INDEX_LENGTH 
        extracted_seq = query_seq[unique_start_query-1:unique_end_query-1] # Convert to 0-based
    else:
        offset_from_align_start = unique_start_subj - sstart
        unique_start_query = qstart + offset_from_align_start
        unique_end_query = unique_start_query + INDEX_LENGTH 
        extracted_seq = query_seq[unique_start_query-1:unique_end_query-1] # Convert to 0-based

    # Validate extracted sequence length
    if len(extracted_seq) != INDEX_LENGTH:
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
    for record in SeqIO.parse(fasta_file, 'fasta'):
        query_sequences[record.id] = str(record.seq)

    # Process BLAST output
    results = []

    with open(blast_file, 'r') as f:
        for line in f:
            if line.strip():
                try:
                    blast_data = parse_blast_line(line) #returns a dictionary with parsed fields per line
                    index_type = get_index_type(blast_data['sseqid'])
                    
                    if index_type and blast_data['qseqid'] in query_sequences:
                        query_seq = query_sequences[blast_data['qseqid']]
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
                id=f"{result['query_id']}.{result['subject_id']}",
                description=""
            )

            writer.write_record(final_seq) #write record to output file

    return results


if __name__ == "__main__":
    args = check_arg()

    results = process_blast_output(args.input, args.reads, args.output)
    print(f"Processed {len(results)} successful extractions", file=sys.stderr)