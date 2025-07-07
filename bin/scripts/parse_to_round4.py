#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import scipy.spatial.distance as sdist
import os
import subprocess
import time
import json
from scipy.cluster.hierarchy import linkage, fcluster
import argparse

parser = argparse.ArgumentParser(description='Process some directories.')
parser.add_argument('--dirt-base', required=True, help='Base directory path')
parser.add_argument('--dirt-fastq', required=True, help='Directory path for fastq files')
parser.add_argument('--qc-file',    required=True, help='QC table with genus annotation')
args = parser.parse_args()

# Convert relative paths to absolute paths and ensure they end with '/'
dirt_base = os.path.abspath(args.dirt_base) + '/'
dirt_fastq = os.path.abspath(args.dirt_fastq) + '/'
qc_file = os.path.abspath(args.qc_file) # ../QC_non_mock-genus_known.tsv

dirt_csv_dirt = dirt_base + "round_3/"
out_put_dirt = dirt_base + "round_4/"
temp_work_dirt = dirt_base + "round_4/temp/"
file_name = '3rd_round_cmp.csv'

if not os.path.exists(out_put_dirt):
    os.makedirs(out_put_dirt)
if not os.path.exists(temp_work_dirt):
    os.makedirs(temp_work_dirt)

mat_file = pd.read_csv(dirt_csv_dirt+file_name)
mat_file = 1.0-mat_file
for i in range(len(mat_file)):
        mat_file.iloc[i,i] = 0


square_form = sdist.squareform(mat_file)
assignments = fcluster(linkage(square_form, method = 'complete'), 0.97)
print(max(assignments))

cluster_result = {}
for j in range(max(assignments)):
    a_temp_list = []
    for i in range(len(assignments)):
        if assignments[i] == j + 1:
            # Remove the prefix and only keep the file name
            file_path = mat_file.columns.values[i]
            file_base_name = os.path.splitext(os.path.basename(file_path))[0]  # Get the file name without extension
            a_temp_list.append(file_base_name)
    cluster_result[str(j)] = a_temp_list

raw_json_path = os.path.join(out_put_dirt, f'cluster_members_{file_name[:3]}.raw.json')

with open(raw_json_path, 'w') as out:
    json.dump(cluster_result, out, indent=4)

print(f"Raw clustering files for Round 4 have been generated at {dirt_base}")

# Load previous round mappings
with open(dirt_base + "round_2_inconsistent/cluster_members_round_2_inconsistent.json") as input_1st:
    cluster_result_readin_1st = json.load(input_1st)
with open(dirt_base + "round_3/cluster_members_round_3.final.json") as input_2nd:
    cluster_result_readin_2nd = json.load(input_2nd)

# Define a function to map Round 4 clusters back to Round 1 based on previous rounds
def generate_final_mapping_to_barcodes(current_round_result, previous_rounds, output_dir, round_name):
    final_mapping = {}
    for cluster_id, members in current_round_result.items():
        barcodes = []
        for member in members:
            if member in previous_rounds["round_3"]:
                barcodes.extend(previous_rounds["round_3"][member])
        final_mapping[cluster_id] = barcodes

    # Save the mapping as JSON and TSV
    json_path = os.path.join(output_dir, f'cluster_members_{round_name}.raw.original.json')
    tsv_path = os.path.join(output_dir, f'cluster_members_{round_name}.raw.original.tsv')

    with open(json_path, 'w') as json_file:
        json.dump(final_mapping, json_file, indent=4)

    with open(tsv_path, 'w') as tsv_file:
        for cluster_id, barcodes in final_mapping.items():
            for barcode in barcodes:
                tsv_file.write(f"{cluster_id}\t{barcode}\n")

    print(f"TSV and JSON files for {round_name} have been generated at {output_dir}")

# Create the final mapping based on Round 4 clusters
previous_rounds = {
    "round_2": cluster_result_readin_1st,
    "round_3": cluster_result_readin_2nd  # Ensure "round_3" is correctly added here
}
generate_final_mapping_to_barcodes(cluster_result, previous_rounds, out_put_dirt, "round_4")


import pandas as pd
import json
import os

def merge_tsv_files(qc_file, cluster_file, output_file):
    qc_df = pd.read_csv(qc_file, sep='\t')
    cluster_df = pd.read_csv(cluster_file, sep='\t', header=None, names=['cluster_id', 'variable'])
    
    merged_df = pd.merge(cluster_df, qc_df, on='variable', how='left')
    merged_df.to_csv(output_file, sep='\t', index=False)
    return merged_df

def evaluate_clusters(merged_df, output_file):
    cluster_evaluation = []
    
    for cluster_id in merged_df['cluster_id'].unique():
        cluster_data = merged_df[merged_df['cluster_id'] == cluster_id]
        if len(cluster_data) < 4:
            cluster_evaluation.append([cluster_id, 'N/A'])
        else:
            most_common_genus = cluster_data['rank.1_taxa'].mode()
            if not most_common_genus.empty:
                most_common_genus = most_common_genus[0]
                if (cluster_data['rank.1_taxa'] == most_common_genus).mean() >= 0.75:
                    cluster_evaluation.append([cluster_id, most_common_genus])
                else:
                    cluster_evaluation.append([cluster_id, 'N/A'])
            else:
                cluster_evaluation.append([cluster_id, 'N/A'])

    evaluation_df = pd.DataFrame(cluster_evaluation, columns=['cluster_id', 'evaluation'])
    evaluation_df.to_csv(output_file, sep='\t', index=False)
    return evaluation_df



def main():
    cluster_file = dirt_base + 'round_4/cluster_members_round_4.raw.original.tsv'
    merged_output = dirt_base + 'round_4/merged_output_round_4.tsv'
    evaluation_output = dirt_base + 'round_4/cluster_evaluation_round_4.tsv'
    filtered_tsv = dirt_base + 'round_4/cluster_members_round_4.tsv'
    filtered_samples_file = dirt_base + 'round_4/contaminated_samples.tsv'
    
    merged_df = merge_tsv_files(qc_file, cluster_file, merged_output)
    
    evaluation_df = evaluate_clusters(merged_df, evaluation_output)
if __name__ == '__main__':
    main()


import pandas as pd
import json

def compare_annotations_and_save_results(dirt_base, round_3_file, round_4_file, json_mapping, merged_output_file, original_json_file, final_json_file, log_file):
    # Load evaluation data
    round_3_eval = pd.read_csv(dirt_base + round_3_file, sep='\t')
    round_4_eval = pd.read_csv(dirt_base + round_4_file, sep='\t')
    merged_output = pd.read_csv(dirt_base + merged_output_file, sep='\t')

    # Replace NaN with 'N/A'
    round_3_eval['evaluation'] = round_3_eval['evaluation'].fillna('N/A')
    round_4_eval['evaluation'] = round_4_eval['evaluation'].fillna('N/A')

    # Load cluster mapping JSON
    with open(dirt_base + json_mapping, 'r') as file:
        cluster_mapping = json.load(file)

    # Load original JSON data
    with open(dirt_base + original_json_file, 'r') as file:
        original_json_data = json.load(file)

    # Mapping from round 3
    round_3_mapping = round_3_eval.set_index('cluster_id')['evaluation'].to_dict()

    # Collect clusters to be removed
    clusters_to_delete = []

    # Compare evaluations
    for cluster_id, members in cluster_mapping.items():
        round_4_evaluation = round_4_eval.loc[round_4_eval['cluster_id'] == int(cluster_id), 'evaluation'].values[0]
        round_3_evaluation_set = {round_3_mapping.get(int(member), 'N/A') for member in members}

        # Check taxa data
        taxa_data = merged_output[merged_output['cluster_id'] == int(cluster_id)]['rank.1_taxa'].dropna().unique()
        if round_4_evaluation == 'N/A' and all(taxa == 'N/A' or taxa in round_3_evaluation_set for taxa in taxa_data):
            continue

        # Record inconsistent annotations
        for member in members:
            round_3_evaluation = round_3_mapping.get(int(member), 'N/A')
            if round_3_evaluation != 'N/A' and (round_4_evaluation == 'N/A' or round_3_evaluation != round_4_evaluation):
                print(f'Round 4 Cluster {cluster_id}: {round_4_evaluation}, Round 3 Cluster {member}: {round_3_evaluation}')
                clusters_to_delete.append(str(cluster_id))
                with open(dirt_base + log_file, 'a') as log:
                    log.write(f'Round 4 Cluster {cluster_id}: {round_4_evaluation}, Round 3 Cluster {member}: {round_3_evaluation}\n')

    # Delete specified clusters from original JSON data
    for cluster_id in clusters_to_delete:
        if cluster_id in original_json_data:
            del original_json_data[cluster_id]

    # Save the modified JSON data
    with open(dirt_base + final_json_file, 'w') as file:
        json.dump(original_json_data, file, indent=4)


def main():
    round_3_file = 'round_3/cluster_evaluation_round_3.tsv'
    round_4_file = 'round_4/cluster_evaluation_round_4.tsv'
    json_mapping = 'round_4/cluster_members_3rd.raw.json'
    merged_output_file = 'round_4/merged_output_round_4.tsv'
    original_json_file = 'round_4/cluster_members_round_4.raw.original.json'
    final_json_file = 'round_4/cluster_members_round_4.de_overfit.json'
    log_file = 'round_4/overfitted.log'
    
    compare_annotations_and_save_results(
        dirt_base,
        round_3_file,
        round_4_file,
        json_mapping,
        merged_output_file,
        original_json_file,
        final_json_file,
        log_file
    )

if __name__ == '__main__':
    main()


import json
import pandas as pd

def filter_tsv_by_json(json_file, tsv_file, output_file):
    with open(json_file, 'r') as f:
        json_data = json.load(f)

    json_elements = set()
    for key, values in json_data.items():
        json_elements.update(values)

    print(f"JSON contains {len(json_elements)} unique elements.")

    tsv_data = pd.read_csv(tsv_file, sep='\t')

    print("TSV file columns:", tsv_data.columns)
    print("First few rows of TSV file:")
    print(tsv_data.head())

    filtered_data = tsv_data[tsv_data['variable'].isin(json_elements)]

    print(f"Filtered data contains {len(filtered_data)} rows.")

    filtered_data.to_csv(output_file, sep='\t', index=False)

def main():

    json_file = dirt_base + 'round_4/cluster_members_round_4.de_overfit.json'
    tsv_file = dirt_base + 'round_4/merged_output_round_4.tsv'
    output_file = dirt_base + 'round_4/merged_output_round_4-de_overfit.tsv'

    filter_tsv_by_json(json_file, tsv_file, output_file)

if __name__ == '__main__':
    main()


def filter_and_generate_output(merged_df, evaluation_df, filtered_tsv, filtered_json, filtered_samples_file):
    final_filtered_rows = []
    filtered_samples = []

    if 'cluster_id' not in evaluation_df.columns or 'evaluation' not in evaluation_df.columns:
        raise ValueError("Evaluation DataFrame must contain 'cluster_id' and 'evaluation' columns.")
    
    print(f"Evaluation DataFrame head:\n{evaluation_df.head()}")

    for cluster_id, evaluation in evaluation_df[['cluster_id', 'evaluation']].itertuples(index=False):
        cluster_data = merged_df[merged_df['cluster_id'] == cluster_id]

        if pd.isna(evaluation) or evaluation == 'N/A':
            final_filtered_rows.append(cluster_data)
        else:
            filtered_data = cluster_data[(cluster_data['rank.1_taxa'] == evaluation) | (cluster_data['rank.1_taxa'].isnull())]
            removed_data = cluster_data[~cluster_data['variable'].isin(filtered_data['variable'])]
            if not removed_data.empty:
                removed_data = removed_data.copy()
                removed_data['annotation'] = evaluation
                filtered_samples.append(removed_data)
            final_filtered_rows.append(filtered_data)

    final_filtered_df = pd.concat(final_filtered_rows, ignore_index=True)
    final_filtered_df.to_csv(filtered_tsv, sep='\t', index=False)

    filtered_json_data = {}
    for cluster_id, group in final_filtered_df.groupby('cluster_id'):
        filtered_json_data[str(cluster_id)] = group['variable'].tolist()

    with open(filtered_json, 'w') as json_file:
        json.dump(filtered_json_data, json_file, indent=4)

    if filtered_samples:
        filtered_samples_df = pd.concat(filtered_samples, ignore_index=True)
        filtered_samples_df.to_csv(filtered_samples_file, sep='\t', index=False)

    print(f"Filtering completed. Filtered samples file: {filtered_samples_file}")



def main():
    merged_output_file = dirt_base + 'round_4/merged_output_round_4-de_overfit.tsv'
    evaluation_output_file = dirt_base + 'round_4/cluster_evaluation_round_4.tsv'

    filtered_tsv = dirt_base + 'round_4/cluster_members_round_4.tsv'
    filtered_json = dirt_base + 'round_4/cluster_members_round_4.final.json'
    filtered_samples_file = dirt_base + 'round_4/contaminated_samples.tsv'
    
    filter_and_generate_output(pd.read_csv(merged_output_file, sep='\t'), pd.read_csv(evaluation_output_file, sep='\t'), filtered_tsv, filtered_json, filtered_samples_file)

if __name__ == '__main__':
    main()


final_json_file = dirt_base + 'round_4/cluster_members_round_4.final.json'
with open(final_json_file, 'r') as json_file:
    filtered_cluster_result = json.load(json_file)

# Process fastq files based on filtered JSON
for cluster_id, barcodes in filtered_cluster_result.items():
    a_temp_list = barcodes
    for BC_temp_with_fastq in a_temp_list:
        BC_temp = os.path.splitext(os.path.basename(BC_temp_with_fastq))[0]  # Get the basename
        subprocess.call(f"cp {dirt_fastq}{BC_temp}_R1_paired.fastq {temp_work_dirt}", shell=True)
        subprocess.call(f"cp {dirt_fastq}{BC_temp}_R2_paired.fastq {temp_work_dirt}", shell=True)
    
    subprocess.call(f"cat {temp_work_dirt}*_R1_paired.fastq > {out_put_dirt}{cluster_id}_R1.fastq", shell=True)
    subprocess.call(f"cat {temp_work_dirt}*_R2_paired.fastq > {out_put_dirt}{cluster_id}_R2.fastq", shell=True)
    subprocess.call(f"rm {temp_work_dirt}*.fastq", shell=True)
    subprocess.call(f"echo process cluster {cluster_id}", shell=True)

print("FASTQ file processing based on filtered JSON completed.")
