#!/usr/bin/env python3

import pandas as pd
import json
import os
import argparse

parser = argparse.ArgumentParser(description='Merge JSON data based on cluster assignments.')
parser.add_argument('-c', '--cluster', type=str, required=True, help='Input CSV file containing cluster assignments')
parser.add_argument('-j', '--json', type=str, required=True, help='Input file of original JSON data')
parser.add_argument('-o', '--output', type=str, required=True, help='Output JSON file for merged data')
parser.add_argument('-x', '--suffix', type=str, default='.500.fasta', help='Suffix to remove from sequence paths (default: .500.fasta)')
parser.add_argument('-t', '--txt', type=str, default='clustered_sequences.txt', help='Output TXT file listing clustered sequences')
args = parser.parse_args()

# Read cluster assignments
cluster_df = pd.read_csv(args.cluster)

# Find clusters containing multiple sequences
cluster_counts = cluster_df['cluster'].value_counts()
clusters_to_merge = cluster_counts[cluster_counts > 1].index.tolist()

# Create a mapping from original cluster ID to new cluster ID starting from 1
cluster_mapping = {original_id: new_id+1 for new_id, original_id in enumerate(clusters_to_merge)}

# Filter cluster_df to only include clusters to be merged
filtered_cluster_df = cluster_df[cluster_df['cluster'].isin(clusters_to_merge)]

# Read original JSON data
with open(args.json, 'r') as json_file:
    round_data = json.load(json_file)

# Create a dictionary to store merged elements for each cluster
cluster_dict = {}

# Record merged sequences
merged_sequences = {}

# Record used round_ids
used_round_ids = set()

# Iterate over filtered cluster assignments and merge JSON elements
for _, row in filtered_cluster_df.iterrows():
    sequence_path = row['sequence']
    original_cluster_id = row['cluster']
    new_cluster_id = cluster_mapping[original_cluster_id]
    cluster_id = f"cluster_{new_cluster_id}"

    # Extract round ID from sequence path using the specified suffix
    round_id = os.path.basename(sequence_path).replace(args.suffix, '')

    # Add this round's elements to the corresponding cluster
    if cluster_id not in cluster_dict:
        cluster_dict[cluster_id] = []

    if round_id in round_data:
        cluster_dict[cluster_id].extend(round_data[round_id])
        used_round_ids.add(round_id)

        # Record which sequences are merged
        if cluster_id not in merged_sequences:
            merged_sequences[cluster_id] = []
        merged_sequences[cluster_id].append(round_id)

# Remove used round_ids from the original data
for round_id in used_round_ids:
    del round_data[round_id]

# Save merged data to new JSON file (only merged clusters)
with open(args.output, 'w') as output_file:
    json.dump(cluster_dict, output_file, indent=4)

# Save the list of clustered sequences to TXT file
with open(args.txt, 'w') as txt_file:
    for cluster_id, sequences in merged_sequences.items():
        txt_file.write(f"{cluster_id}: {', '.join(sequences)}\n")

print(f"Merging complete. Results saved to '{args.output}'.")
print(f"List of clustered sequences saved to '{args.txt}'.")

