#!/usr/bin/env python3

import os
import sys
import argparse
import glob

def parse_breadth_file(filepath, species_ids, coverage_rate_threshold, num_covered_nt_threshold):
    sequences = set()
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            seq_id = fields[0]
            try:
                num_covered_nt = float(fields[1])
                orig_length = float(fields[2])
                coverage_rate = float(fields[3])
            except ValueError:
                continue
            # Apply filtering criteria
            if seq_id in species_ids and (coverage_rate >= coverage_rate_threshold and num_covered_nt >= num_covered_nt_threshold):
                sequences.add(seq_id)
    return sequences

def read_species_ids(species_id_file):
    species_ids = set()
    with open(species_id_file, 'r') as f:
        for line in f:
            species_ids.add(line.strip())
    return species_ids

def main():
    parser = argparse.ArgumentParser(description='Generate a presence/absence table for sequences across samples.')
    parser.add_argument('-i', '--input', required=True, help='Input breadth files pattern (e.g., "work/Aliarcobacter_skirrowii/mapping/*.breadth")')
    parser.add_argument('-s', '--species_ids', required=True, help='File containing species IDs to filter (one per line)')
    parser.add_argument('-o', '--output', required=True, help='Output table file')
    parser.add_argument('--coverage_rate_threshold', type=float, default=0.90, help='Minimum coverage rate threshold (default: 0.90)')
    parser.add_argument('--num_covered_nt_threshold', type=float, default=1000, help='Minimum number of covered nucleotides threshold (default: 1000)')

    args = parser.parse_args()

    # Find all breadth files matching the input pattern
    breadth_files = glob.glob(args.input)
    if not breadth_files:
        print('No breadth files found with pattern:', args.input)
        sys.exit(1)

    # Read species IDs to filter
    species_ids = read_species_ids(args.species_ids)

    all_sequences = set()
    sample_sequences = {}  # key: sample_name, value: set of sequences that pass the filters

    for filepath in breadth_files:
        sample_name = os.path.basename(filepath).replace('.breadth', '')
        sequences = parse_breadth_file(filepath, species_ids, args.coverage_rate_threshold, args.num_covered_nt_threshold)
        all_sequences.update(sequences)
        sample_sequences[sample_name] = sequences

    all_sequences = sorted(all_sequences)

    # Write the presence/absence table
    with open(args.output, 'w') as out_f:
        # Write header
        header = ['Sequence_ID'] + sorted(sample_sequences.keys())
        out_f.write('\t'.join(header) + '\n')
        for seq_id in all_sequences:
            row = [seq_id]
            for sample_name in sorted(sample_sequences.keys()):
                if seq_id in sample_sequences[sample_name]:
                    row.append('1')
                else:
                    row.append('0')
            out_f.write('\t'.join(row) + '\n')

if __name__ == '__main__':
    main()

