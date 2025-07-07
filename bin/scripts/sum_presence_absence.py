#!/usr/bin/env python3

import sys
import os
import pandas as pd

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Process presence_absence_table.txt to sum counts.')
    parser.add_argument('-i', '--input', required=True, help='Input presence_absence_table.txt')
    parser.add_argument('-o', '--output', required=True, help='Output two-column file')
    parser.add_argument('-s', '--species', required=False, help='Species name to use as header. If not provided, derived from directory name.')
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    species_name = args.species

    if species_name is None:
        # Try to derive species name from directory name
        # Assume the input file is in work/<species>/presence_absence_table.txt
        input_dir = os.path.dirname(os.path.abspath(input_file))
        parent_dir = os.path.basename(input_dir)
        species_name = parent_dir

    # Read the presence_absence_table.txt
    df = pd.read_csv(input_file, sep='\t')

    # Sum across rows (excluding the first column)
    df[species_name] = df.iloc[:, 1:].sum(axis=1)

    # Keep only Sequence_ID and the new column
    df_out = df[['Sequence_ID', species_name]]

    # Write the output
    df_out.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    main()

