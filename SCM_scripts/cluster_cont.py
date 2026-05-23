import pandas as pd
import scipy.spatial.distance as sdist
import os
import subprocess
import time
import json
from scipy.cluster.hierarchy import linkage, fcluster
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Process directories and JSON file for clustering.')
parser.add_argument('--dirt_base', required=True, help='Base directory path')
parser.add_argument('--dirt_fastq', required=True, help='Directory path for FASTQ files')
parser.add_argument('-j', '--json', required=True, help='Input JSON file with cluster data')
args = parser.parse_args()

# Convert paths to absolute paths and ensure they end with '/'
dirt_base = os.path.abspath(args.dirt_base)
if not dirt_base.endswith('/'):
    dirt_base += '/'

dirt_fastq = os.path.abspath(args.dirt_fastq)
if not dirt_fastq.endswith('/'):
    dirt_fastq += '/'

# Set up output and temporary directories
out_put_dirt = dirt_base + "cluster_cont/"  # Output directory for files
temp_work_dirt = dirt_base + "cluster_cont/temp/"  # Temporary directory for FASTQ files

if not os.path.exists(out_put_dirt):
    os.makedirs(out_put_dirt)
if not os.path.exists(temp_work_dirt):
    os.makedirs(temp_work_dirt)

# Load the JSON file with cluster data
final_json_file = args.json
with open(final_json_file, 'r') as json_file:
    filtered_cluster_result = json.load(json_file)

# Process FASTQ files based on the filtered JSON
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
