import pandas as pd
import networkx as nx
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Cluster sequences based on ANI values.')
parser.add_argument('-i', '--input', type=str, required=True, help='Input file (fastANI data)')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file (clustering results)')
args = parser.parse_args()

# Load fastANI data
df = pd.read_csv(args.input, sep='\t', header=None, names=['seq1', 'seq2', 'ani', 'len1', 'len2'])

# Filter pairs with ANI greater than 97
filtered_df = df[df['ani'] > 97]

# Create a graph where nodes are sequences; edges represent ANI > 97 relationships
G = nx.Graph()
for _, row in filtered_df.iterrows():
    G.add_edge(row['seq1'], row['seq2'])

# Find connected components (clusters)
clusters = list(nx.connected_components(G))

# Assign a cluster ID to each sequence
cluster_dict = {}
for cluster_id, cluster in enumerate(clusters, start=1):
    for seq in cluster:
        cluster_dict[seq] = cluster_id

# Output the clustering results
cluster_df = pd.DataFrame(list(cluster_dict.items()), columns=['sequence', 'cluster'])
cluster_df.to_csv(args.output, index=False)

print("Clustering complete. Results saved to '{}'.".format(args.output))
