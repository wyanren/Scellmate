#!/usr/bin/env python3

import os
import json
import re
import argparse

def ordinal(n):
    """
    Convert an integer n to its ordinal representation:
    1 -> '1st', 2 -> '2nd', 3 -> '3rd', etc.
    """
    if 10 <= n % 100 <= 20:
        suffix = 'th'
    else:
        suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')
    return f"{n}{suffix}"

def get_raw_cluster_file(base_dir, current_round):
    """
    Given the current round, return the path to the raw cluster data file.
    """
    if current_round == 'round_2_inconsistent':
        # For round_2_inconsistent, use specific file names
        raw_cluster_file = os.path.join(base_dir, 'round_2_inconsistent', 'cluster_members_2nd.raw.json')
    else:
        # Extract the round number, e.g., 'round_5' -> 5
        match = re.match(r'round_(\d+)', current_round)
        if not match:
            print(f"Unable to extract round number from {current_round}")
            return None
        round_number = int(match.group(1))
        if round_number < 1:
            print(f"Invalid round number {round_number} in {current_round}")
            return None
        # Get the previous round number
        prev_round_number = round_number - 1
        if prev_round_number < 1:
            print(f"{current_round} has no previous round")
            return None
        # Get ordinal suffix for raw_cluster_file
        ordinal_suffix = ordinal(prev_round_number)
        raw_cluster_file = os.path.join(base_dir, current_round, f'cluster_members_{ordinal_suffix}.raw.json')
    return raw_cluster_file

def get_final_cluster_file(base_dir, round_name):
    """
    Given a round name, return the path to the final cluster data file.
    """
    if round_name == 'round_2_inconsistent':
        # For round_2_inconsistent, use specific file names
        final_cluster_file = os.path.join(base_dir, 'round_2_inconsistent', 'cluster_members_round_2_inconsistent.json')
    else:
        # Extract the round number, e.g., 'round_4' -> 4
        match = re.match(r'round_(\d+)', round_name)
        if not match:
            print(f"Unable to extract round number from {round_name}")
            return None
        round_number = int(match.group(1))
        if round_number < 1:
            print(f"Invalid round number {round_number} in {round_name}")
            return None
        final_cluster_file = os.path.join(base_dir, round_name, f'cluster_members_round_{round_number}.final.json')
    return final_cluster_file

def get_parent_clusters(base_dir, current_round, cluster_id):
    """
    Retrieve the list of parent cluster IDs for a given cluster_id in current_round.
    """
    raw_cluster_file = get_raw_cluster_file(base_dir, current_round)
    if not raw_cluster_file:
        print(f"Unable to determine raw cluster file path for {current_round}")
        return []
    if not os.path.exists(raw_cluster_file):
        print(f"Raw cluster data file not found: {raw_cluster_file}")
        return []
    try:
        with open(raw_cluster_file, 'r') as f:
            raw_cluster_data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from {raw_cluster_file}: {e}")
        return []
    # Assuming cluster IDs are strings in JSON
    parent_clusters = raw_cluster_data.get(str(cluster_id), [])
    print(f"Retrieved parent clusters for cluster {cluster_id} in {current_round}: {parent_clusters}")
    return parent_clusters

def get_barcodes(base_dir, prev_round, parent_cluster_id):
    """
    Retrieve the list of barcodes for a given parent_cluster_id in prev_round.
    """
    final_cluster_file = get_final_cluster_file(base_dir, prev_round)
    if not final_cluster_file:
        print(f"Unable to determine final cluster file path for {prev_round}")
        return []
    if not os.path.exists(final_cluster_file):
        print(f"Final cluster data file not found: {final_cluster_file}")
        return []
    try:
        with open(final_cluster_file, 'r') as f:
            final_cluster_data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from {final_cluster_file}: {e}")
        return []
    # Assuming cluster IDs are strings in JSON
    barcodes = final_cluster_data.get(str(parent_cluster_id), [])
    if not barcodes:
        print(f"No barcodes found for cluster {parent_cluster_id} in {prev_round}")
    return barcodes

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Process overfit clusters and generate JSON files.')
    parser.add_argument('--base_dir', type=str, required=True, help='Base directory path.')
    parser.add_argument('--current_dir', type=str, required=True, help='Current directory path (e.g., 1st_align).')
    args = parser.parse_args()

    # Assign arguments to variables
    base_dir = args.base_dir
    current_dir = args.current_dir

    # Validate base_dir
    if not os.path.exists(base_dir):
        print(f"Base directory does not exist: {base_dir}")
        exit(1)

    # Validate current_dir
    full_current_dir = os.path.join(base_dir, current_dir)
    if not os.path.exists(full_current_dir):
        print(f"Current directory does not exist: {full_current_dir}")
        exit(1)

    # Set the split_output_overfit directory based on current_dir
    OVERFIT_DIR = os.path.join(full_current_dir, "split_output_overfit")

    # Check if OVERFIT_DIR exists
    if not os.path.exists(OVERFIT_DIR):
        print(f"Overfit directory not found: {OVERFIT_DIR}")
        exit(1)

    # Collect unique cluster IDs to process
    unique_clusters = set()

    # Iterate over the overfitting cluster files in split_output_overfit/
    for filename in os.listdir(OVERFIT_DIR):
        if not filename.endswith('.txt'):
            continue  # Skip non-text files
        # Example filename: 'round_5_87_group_0_cells.txt'
        match = re.match(r'(round_\d+(_inconsistent)?)_(\d+)_group_\d+_cells\.txt', filename)
        if not match:
            print(f"Filename does not match expected pattern, skipping: {filename}")
            continue  # Skip files that don't match the expected pattern
        current_round = match.group(1)  # e.g., 'round_5' or 'round_2_inconsistent'
        cluster_id = match.group(3)      # e.g., '87'
        unique_clusters.add((current_round, cluster_id))

    # Process each unique cluster
    for current_round, cluster_id in unique_clusters:
        print(f'Processing overfitting cluster {cluster_id} in {current_round}')
        
        # Get parent clusters from the current round's raw cluster data
        parent_clusters = get_parent_clusters(base_dir, current_round, cluster_id)
        if not parent_clusters:
            print(f'No parent clusters found for cluster {cluster_id} in {current_round}')
            print('')
            continue
        
        # Determine the previous round
        if current_round == 'round_2_inconsistent':
            prev_round = 'round_1'  # Assuming 'round_1' is the previous round
            print(f"Current round is 'round_2_inconsistent', setting previous round to 'round_1'")
        else:
            # Extract round number
            match_round = re.match(r'round_(\d+)', current_round)
            if not match_round:
                print(f'Unable to extract round number from {current_round}')
                print('')
                continue
            round_number = int(match_round.group(1))
            prev_round_number = round_number - 1
            prev_round = f'round_{prev_round_number}'
            # Special case: if previous round is 2 and 'round_2_inconsistent' exists, use it
            if prev_round_number == 2:
                inconsistent_final_cluster_file = os.path.join(base_dir, 'round_2_inconsistent', 'cluster_members_round_2_inconsistent.json')
                if os.path.exists(inconsistent_final_cluster_file):
                    prev_round = 'round_2_inconsistent'
                    print(f"Previous round number is 2 and 'round_2_inconsistent' exists, setting previous round to 'round_2_inconsistent'")
                else:
                    print(f"Previous round number is 2, setting previous round to 'round_2'")
            print(f'Previous round: {prev_round}')
        
        # Iterate through each parent cluster and generate JSON files
        for parent_cluster_id in parent_clusters:
            barcodes = get_barcodes(base_dir, prev_round, parent_cluster_id)
            if barcodes:
                # Add prefix to the cluster ID
                prefixed_cluster_id = f'{prev_round}_{parent_cluster_id}'
                # Prepare the output JSON data
                output_data = {prefixed_cluster_id: barcodes}
                # Output JSON file name
                output_filename = f'{prefixed_cluster_id}.json'
                # Set output path to split_output_overfit directory
                output_filepath = os.path.join(OVERFIT_DIR, output_filename)
                # Write the JSON data
                try:
                    with open(output_filepath, 'w') as outfile:
                        json.dump(output_data, outfile, indent=4)
                    print(f'Generated JSON file: {output_filename}')
                except IOError as e:
                    print(f"Error writing to file {output_filename}: {e}")
            else:
                print(f'No barcodes found for cluster {parent_cluster_id} in {prev_round}')
        print('')  # Add a newline for readability

if __name__ == '__main__':
    main()

