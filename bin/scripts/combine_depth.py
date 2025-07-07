#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(description='Combine .depth files by merging the first two columns and removing duplicates.')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Input .depth files (supports wildcard expansion).')
    parser.add_argument('-o', '--output', required=True, help='Output combined .depth file.')

    args = parser.parse_args()

    unique_entries = set()

    for depth_file in args.input:
        with open(depth_file, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                ref_id = parts[0]
                position = parts[1]
                unique_entries.add((ref_id, position))

    # Write the unique entries to the output file, sorted for consistency
    with open(args.output, 'w') as out_f:
        for ref_id, position in sorted(unique_entries):
            out_f.write(f'{ref_id}\t{position}\n')

if __name__ == '__main__':
    main()

