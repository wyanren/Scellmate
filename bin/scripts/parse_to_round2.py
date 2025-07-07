#!/usr/bin/env python
# coding: utf-8
"""
Introduction:
1.  Read the round-1 similarity matrix, cluster SAGs (97 % cutoff), and write raw cluster membership JSON/TSV for “round_2_inconsistent”.
2.  Use 1st-QC annotations as supervision, evaluate genus‐level consistency, filter clusters, and emit cleaned cluster files + list of removed (contaminated) SAGs.
3.  For each retained cluster, combine its paired FASTQ files into cluster-level R1/R2 FASTQs.
"""

import pandas as pd
import scipy.spatial.distance as sdist
import os, subprocess, json
from scipy.cluster.hierarchy import linkage, fcluster
import argparse

parser = argparse.ArgumentParser(description='Process some directories.')
parser.add_argument('--dirt-base', required=True, help='Base directory path')
parser.add_argument('--dirt-fastq', required=True, help='Directory path for fastq files')
parser.add_argument('--qc-file',    required=True, help='QC table with genus annotation')
args = parser.parse_args()

# Convert to absolute paths and append trailing slash (end with '/')
dirt_base = os.path.abspath(args.dirt_base) + '/'
dirt_fastq = os.path.abspath(args.dirt_fastq) + '/'
qc_file = os.path.abspath(args.qc_file) # ../QC_non_mock-genus_known.tsv

dirt_csv_dirt = dirt_base + "round_1/"
out_put_dirt = dirt_base + "round_2_inconsistent/" # dirt of output dirt

if not os.path.exists(out_put_dirt):
    os.makedirs(out_put_dirt)
temp_work_dirt = dirt_base + "round_2_inconsistent/temp/" # temp dirt used to cat fastq files
if not os.path.exists(temp_work_dirt):
    os.makedirs(temp_work_dirt)
file_name = '1st_round_cmp.csv' # sourmash produce

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
        if assignments[i] == j+1:
            file_base_name = os.path.splitext(os.path.basename(mat_file.columns.values[i]))[0]
            a_temp_list.append(file_base_name)
            # a_temp_list += [mat_file.columns.values[i][:-len('.fasta')]]
    cluster_result[str(j)]=a_temp_list
with open(out_put_dirt+'cluster_members_'+file_name[:3]+'.raw.json', 'w') as out:
    json.dump(cluster_result, out)
with open(out_put_dirt+'cluster_members_'+file_name[:3]+'.raw.json') as input:
    cluster_result_readin = json.load(input)

# raw JSON + TSV
round_2_inconsistent_raw_tsv_path = os.path.join(out_put_dirt, 'cluster_members_round_2_inconsistent.raw.tsv')
round_2_inconsistent_raw_json_path = os.path.join(out_put_dirt, 'cluster_members_round_2_inconsistent.raw.json')

with open(round_2_inconsistent_raw_json_path, 'w') as json_file:
    json.dump(cluster_result, json_file, indent=4)

with open(round_2_inconsistent_raw_tsv_path, 'w') as tsv_file:
    for cluster_id, barcodes in cluster_result.items():
        for barcode in barcodes:
            tsv_file.write(f"{cluster_id}\t{barcode}\n")

print(f"Raw clustering files for round_2_inconsistent have been generated at {dirt_base}")
print(f"Start decontamination based on genus-level annotation by marker gene")

import time
# part 2: supervision + cluster filtering
def merge_tsv_files(qc_file, cluster_file, output_file):
    qc_df = pd.read_csv(qc_file, sep='\t')
    cluster_df = pd.read_csv(cluster_file, sep='\t', header=None, names=['cluster_id', 'variable'])
    merged_df = pd.merge(cluster_df, qc_df, on='variable', how='left')
    merged_df.to_csv(output_file, sep='\t', index=False)
    return merged_df

# genus-level annotation of temporary cluster
def evaluate_clusters(merged_df, output_file):
    cluster_evaluation = []
    
    for cluster_id in merged_df['cluster_id'].unique():
        cluster_data = merged_df[merged_df['cluster_id'] == cluster_id]
        if len(cluster_data) < 4:
            cluster_evaluation.append([cluster_id, 'N/A'])
        else:
            # if more than 75% SAGs of the genus in one cluster is consistent —— than give cluster annotation
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

# part 3: combine FASTQ per cluster
def filter_and_generate_output(merged_df, evaluation_df, filtered_tsv, filtered_json, filtered_samples_file):
    final_filtered_rows = []
    filtered_samples = []

    for cluster_id, evaluation in evaluation_df.values:
        cluster_data = merged_df[merged_df['cluster_id'] == cluster_id]
        if evaluation != 'N/A':
            filtered_data = cluster_data[(cluster_data['rank.1_taxa'] == evaluation) | (cluster_data['rank.1_taxa'].isnull())]
            removed_data = cluster_data[~cluster_data['variable'].isin(filtered_data['variable'])]
            if not removed_data.empty:
                removed_data = removed_data.copy()
                removed_data['annotation'] = evaluation
                filtered_samples.append(removed_data)
            final_filtered_rows.append(filtered_data)
        else:
            final_filtered_rows.append(cluster_data)

    final_filtered_df = pd.concat(final_filtered_rows)
    final_filtered_df.to_csv(filtered_tsv, sep='\t', index=False)

    filtered_json_data = {}
    for cluster_id, group in final_filtered_df.groupby('cluster_id'):
        filtered_json_data[str(cluster_id)] = group['variable'].tolist()

    with open(filtered_json, 'w') as json_file:
        json.dump(filtered_json_data, json_file, indent=4)

    if filtered_samples:
        filtered_samples_df = pd.concat(filtered_samples)
        filtered_samples_df.to_csv(filtered_samples_file, sep='\t', index=False)

def main():
    cluster_file = round_2_inconsistent_raw_tsv_path
    merged_output = dirt_base + 'round_2_inconsistent/merged_output_round_2_inconsistent.tsv'
    evaluation_output = dirt_base + 'round_2_inconsistent/cluster_evaluation_round_2_inconsistent.tsv'
    filtered_tsv = dirt_base + 'round_2_inconsistent/cluster_members_round_2_inconsistent.tsv'
    filtered_json = dirt_base + 'round_2_inconsistent/cluster_members_round_2_inconsistent.json'
    filtered_samples_file = dirt_base + 'round_2_inconsistent/contaminated_samples.tsv'
    
    merged_df = merge_tsv_files(qc_file, cluster_file, merged_output)
    evaluation_df = evaluate_clusters(merged_df, evaluation_output)
    filter_and_generate_output(merged_df, evaluation_df, filtered_tsv, filtered_json, filtered_samples_file)

    print(f"All tasks completed. Files generated:\n{merged_output}\n{evaluation_output}\n{filtered_tsv}\n{filtered_json}\n{filtered_samples_file}")

if __name__ == '__main__':
    main()

filtered_json_path = dirt_base + 'round_2_inconsistent/cluster_members_round_2_inconsistent.json'
with open(filtered_json_path, 'r') as json_file:
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

