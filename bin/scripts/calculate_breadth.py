#!/usr/bin/env python3
"""
This script calculates the breadth (coverage rate) from a depth file.
It outputs four columns:
1. Reference ID
2. Number of covered nucleotides
3. Original length of the reference
4. Coverage rate (breadth)

Usage:
    python calculate_breadth.py <depth_file>

The depth file should be in the format:
<reference_id>\t<position>\t<depth>

Example line from the depth file:
scDNA_Inf_10374_NODE_2_length_2799_cov_4.639213    2784    10
"""

import sys
import re

def extract_length(ref_name):
    # Extracts the length information from the reference name using regex
    pattern = r'length_(\d+)'
    match = re.search(pattern, ref_name)
    if match:
        return int(match.group(1))
    else:
        print(f"Cannot extract length from {ref_name}")
        return None

def main(depth_file):
    ref_data = {}

    with open(depth_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Parse each line of the depth file
            ref_name, pos_str, depth_str = line.split('\t')
            position = int(pos_str)
            depth = int(depth_str)
            if ref_name not in ref_data:
                length = extract_length(ref_name)
                if length is None:
                    continue
                ref_data[ref_name] = {'length': length, 'positions': set()}
            ref_data[ref_name]['positions'].add(position)

    # Output the results
    for ref_name, data in ref_data.items():
        length = data['length']
        num_positions_covered = len(data['positions'])
        breadth = num_positions_covered / length
        print(f"{ref_name}\t{num_positions_covered}\t{length}\t{breadth:.4f}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python calculate_breadth.py <depth_file>")
    else:
        depth_file = sys.argv[1]
        main(depth_file)

