#!/usr/bin/env python3

import sys
import pandas as pd
import glob

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Merge multiple two-column tables on Sequence_ID.')
    parser.add_argument('-i', '--input_pattern', required=True, help='Input file pattern (e.g., "work/*/summed_table.txt")')
    parser.add_argument('-o', '--output', required=True, help='Output merged table')
    args = parser.parse_args()

    input_files = glob.glob(args.input_pattern)
    if not input_files:
        print('No input files found matching pattern:', args.input_pattern)
        sys.exit(1)

    merged_df = None
    for input_file in input_files:
        df = pd.read_csv(input_file, sep='\t')
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on='Sequence_ID', how='outer')

    # Fill NaN with 0
    merged_df = merged_df.fillna(0)
    # Convert counts to integers
    for col in merged_df.columns[1:]:
        merged_df[col] = merged_df[col].astype(int)
    # Write the merged table
    merged_df.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()

