import pandas as pd
import json
import os
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Merge JSON data based on cluster assignments.')
parser.add_argument('-c', '--cluster', type=str, required=True, help='Input CSV file with cluster assignments')
parser.add_argument('-j', '--json', type=str, required=True, help='Input JSON file with original data')
parser.add_argument('-o', '--output', type=str, required=True, help='Output JSON file for merged data')
parser.add_argument('-x', '--suffix', type=str, default='.500.fasta', help='Suffix to remove from sequence paths (default: .500.fasta)')
args = parser.parse_args()

# Load the cluster assignment
cluster_df = pd.read_csv(args.cluster)

# Load the original JSON data
with open(args.json, 'r') as json_file:
    round_data = json.load(json_file)

# Create a dictionary to hold the merged elements for each cluster
cluster_dict = {}

# Loop through the cluster assignments and merge JSON elements
for _, row in cluster_df.iterrows():
    sequence_path = row['sequence']
    cluster_id = f"cluster_{row['cluster']}"

    # Extract the round ID from the sequence path using the specified suffix
    round_id = os.path.basename(sequence_path).replace(args.suffix, '')

    # Add the elements of this round to the corresponding cluster
    if cluster_id not in cluster_dict:
        cluster_dict[cluster_id] = []

    if round_id in round_data:
        cluster_dict[cluster_id].extend(round_data[round_id])

# Save the merged data into a new JSON file
with open(args.output, 'w') as output_file:
    json.dump(cluster_dict, output_file, indent=4)

print(f"Merging complete. Results saved to '{args.output}'.")
