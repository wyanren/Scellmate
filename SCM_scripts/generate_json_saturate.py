#!/usr/bin/env python3
import os
import argparse
import json

def generate_json_files(directory):
    # List all files in the directory
    files = os.listdir(directory)
    # Filter for files that do NOT contain 'group' and end with '_cells.txt'
    non_group_files = [f for f in files if f.endswith('_cells.txt') and 'group' not in f]

    for filename in non_group_files:
        filepath = os.path.join(directory, filename)

        # Read the cell IDs from the file
        with open(filepath, 'r') as f:
            cell_ids = [line.strip() for line in f if line.strip()]

        # Create the JSON content
        key_name = filename.replace('_cells.txt', '')
        json_content = {key_name: cell_ids}

        # Define the output JSON filename
        json_filename = f"{key_name}.json"
        json_filepath = os.path.join(directory, json_filename)

        # Write the JSON file
        with open(json_filepath, 'w') as json_file:
            json.dump(json_content, json_file, indent=2)

        print(f"Generated JSON file: {json_filepath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate JSON files from non-group *_cells.txt files.')
    parser.add_argument('--directory', type=str, required=True, help='Directory containing *_cells.txt files.')
    args = parser.parse_args()

    generate_json_files(args.directory)

