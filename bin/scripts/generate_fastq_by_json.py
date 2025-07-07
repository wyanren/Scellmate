import os
import subprocess
import json
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Process directories and JSON file for clustering.')
parser.add_argument('--dirt_base', required=True, help='Base directory path')
parser.add_argument('--dirt_fastq', required=True, help='Directory path for FASTQ files')
parser.add_argument('-j', '--json', required=True, help='Input JSON file with cluster data')
parser.add_argument('--output_dir', required=True, help='Output directory path')
args = parser.parse_args()

# Convert paths to absolute paths and ensure they end with '/'
dirt_base = os.path.abspath(args.dirt_base)
if not dirt_base.endswith('/'):
    dirt_base += '/'

dirt_fastq = os.path.abspath(args.dirt_fastq)
if not dirt_fastq.endswith('/'):
    dirt_fastq += '/'

# Set up output and temporary directories
out_put_dirt = os.path.abspath(args.output_dir)
if not out_put_dirt.endswith('/'):
    out_put_dirt += '/'

temp_work_dirt = os.path.join(out_put_dirt, "temp/")
if not temp_work_dirt.endswith('/'):
    temp_work_dirt += '/'

# Ensure output and temporary directories exist
if not os.path.exists(out_put_dirt):
    os.makedirs(out_put_dirt)
if not os.path.exists(temp_work_dirt):
    os.makedirs(temp_work_dirt)

# Load the JSON file with cluster data
final_json_file = args.json
with open(final_json_file, 'r') as json_file:
    filtered_cluster_result = json.load(json_file)

# Process FASTQ files based on the filtered JSON
for cluster_dict in filtered_cluster_result:
    for cluster_id, barcodes in cluster_dict.items():
        a_temp_list = barcodes
        for barcode in a_temp_list:
            # Since your barcodes are like "scDNA_Inf_00628", we assume the FASTQ files are named accordingly
            BC_temp = barcode  # Use the barcode directly
            # Copy the FASTQ files to the temporary directory
            subprocess.call(f"cp {dirt_fastq}{BC_temp}_R1_paired.fastq {temp_work_dirt}", shell=True)
            subprocess.call(f"cp {dirt_fastq}{BC_temp}_R2_paired.fastq {temp_work_dirt}", shell=True)

        # Concatenate the FASTQ files into the output directory with the desired names
        subprocess.call(f"cat {temp_work_dirt}*_R1_paired.fastq > {out_put_dirt}{cluster_id}_R1.fastq", shell=True)
        subprocess.call(f"cat {temp_work_dirt}*_R2_paired.fastq > {out_put_dirt}{cluster_id}_R2.fastq", shell=True)
        # Remove the temporary FASTQ files
        subprocess.call(f"rm {temp_work_dirt}*.fastq", shell=True)
        print(f"Processed cluster {cluster_id}")

