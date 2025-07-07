#!/usr/bin/env python3
import os
import argparse
import json

def generate_json_files(directory):
    # List all files in the directory
    files = os.listdir(directory)
    # Filter for group_*_cells.txt files
    group_files = [f for f in files if f.endswith('_cells.txt') and ('_group_0_cells.txt' in f or '_group_1_cells.txt' in f)]

    for filename in group_files:
        filepath = os.path.join(directory, filename)
        # Extract the basename and group number
        if '_group_0_cells.txt' in filename:
            basename = filename.replace('_group_0_cells.txt', '')
            group_number = '0'
        elif '_group_1_cells.txt' in filename:
            basename = filename.replace('_group_1_cells.txt', '')
            group_number = '1'
        else:
            continue  # Skip any files that don't match the expected pattern

        # Read the cell IDs from the file
        with open(filepath, 'r') as f:
            cell_ids = [line.strip() for line in f if line.strip()]

        # Create the JSON content
        key_name = f"{basename}_{group_number}"
        json_content = {key_name: cell_ids}

        # Define the output JSON filename
        json_filename = f"{basename}_{group_number}.json"
        json_filepath = os.path.join(directory, json_filename)

        # Write the JSON file
        with open(json_filepath, 'w') as json_file:
            json.dump(json_content, json_file, indent=2)

        print(f"Generated JSON file: {json_filepath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate JSON files from group_*_cells.txt files.')
    parser.add_argument('--directory', type=str, required=True, help='Directory containing group_*_cells.txt files.')
    args = parser.parse_args()

    generate_json_files(args.directory)
