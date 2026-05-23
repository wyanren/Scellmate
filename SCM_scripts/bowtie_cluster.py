#!/usr/bin/env python3
import os
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import shutil
import glob

def get_contig_lengths(fasta_file):
    contig_lengths = {}
    with open(fasta_file, 'r') as f:
        contig_name = ''
        sequence = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if contig_name:
                    contig_lengths[contig_name] = len(sequence)
                contig_name = line[1:].split()[0]
                sequence = ''
            else:
                sequence += line
        # Don't forget the last contig
        if contig_name:
            contig_lengths[contig_name] = len(sequence)
    return contig_lengths

def process_cluster(cluster_number):
    print(f"Processing cluster {cluster_number} with barcodes: {', '.join(cluster_info[cluster_number])}")
    temp_output_dirt = os.path.join(align_base, cluster_number + '/')
    os.makedirs(temp_output_dirt, exist_ok=True)
    
    # Copy the cluster fasta file to temp_output_dirt
    cluster_fasta_src = os.path.join(dirt_bins, f"{cluster_number}.fasta")
    cluster_fasta_dst = os.path.join(temp_output_dirt, f"{cluster_number}.fasta")
    shutil.copy(cluster_fasta_src, cluster_fasta_dst)
    
    # Copy the FASTQ files for each barcode
    for barcode_temp in cluster_info[cluster_number]:
        fastq_r1_src = os.path.join(dirt_fastq, f"{barcode_temp}_R1_paired.fastq")
        fastq_r2_src = os.path.join(dirt_fastq, f"{barcode_temp}_R2_paired.fastq")
        fastq_r1_dst = os.path.join(temp_output_dirt, f"{barcode_temp}_R1_paired.fastq")
        fastq_r2_dst = os.path.join(temp_output_dirt, f"{barcode_temp}_R2_paired.fastq")
        shutil.copy(fastq_r1_src, fastq_r1_dst)
        shutil.copy(fastq_r2_src, fastq_r2_dst)
    
    # Build bowtie2 index
    cluster_fasta = cluster_fasta_dst
    reference_prefix = os.path.join(temp_output_dirt, f"{cluster_number}_reference")
    bowtie2_build_command = f"bowtie2-build {cluster_fasta} {reference_prefix}"
    with open(log_file_path, 'a') as log_file:
        subprocess.call(bowtie2_build_command, shell=True, stdout=log_file, stderr=log_file)
    
    # Get the list of samples
    fastq_r1_files = glob.glob(os.path.join(temp_output_dirt, "*_R1_paired.fastq"))
    
    # For each sample
    for fastq_r1_file in fastq_r1_files:
        sample_prefix = fastq_r1_file.replace('_R1_paired.fastq', '')
        fastq_r1 = sample_prefix + '_R1_paired.fastq'
        fastq_r2 = sample_prefix + '_R2_paired.fastq'
        sam_file = sample_prefix + '.sam'
        bam_file = sample_prefix + '.bam'
        tsv_file = sample_prefix + '.tsv'
        idxstats_file = sample_prefix + '.read_coverage'
        
        # Run bowtie2 alignment
        bowtie2_command = f"bowtie2 -x {reference_prefix} -1 {fastq_r1} -2 {fastq_r2} -S {sam_file}"
        with open(log_file_path, 'a') as log_file:
            subprocess.call(bowtie2_command, shell=True, stdout=log_file, stderr=log_file)
        
        # Convert SAM to sorted BAM
        samtools_sort_command = f"samtools sort {sam_file} -o {bam_file}"
        with open(log_file_path, 'a') as log_file:
            subprocess.call(samtools_sort_command, shell=True, stdout=log_file, stderr=log_file)
        # Remove the SAM file
        os.remove(sam_file)
        
        # Run samtools depth and process depth
        samtools_depth_command = f"samtools depth {bam_file}"
        with open(log_file_path, 'a') as log_file:
            depth_process = subprocess.Popen(samtools_depth_command, shell=True, stdout=subprocess.PIPE, stderr=log_file, text=True)
        # Get contig lengths
        contig_lengths = get_contig_lengths(cluster_fasta)
        depth_per_contig = {}
        for line in depth_process.stdout:
            if line.strip():
                contig, position, depth = line.strip().split('\t')
                depth = int(depth)
                depth_per_contig[contig] = depth_per_contig.get(contig, 0) + depth
        depth_process.wait()
        
        # Write depth info to TSV file
        with open(tsv_file, 'w') as tsv_out:
            for contig in depth_per_contig:
                total_depth = depth_per_contig[contig]
                length = contig_lengths.get(contig, 0)
                if length > 0:
                    average_depth = total_depth / length
                else:
                    average_depth = 0
                tsv_out.write(f"{contig}\t{total_depth}\t{average_depth}\n")
        
        # Index the BAM file
        samtools_index_command = f"samtools index {bam_file}"
        with open(log_file_path, 'a') as log_file:
            subprocess.call(samtools_index_command, shell=True, stdout=log_file, stderr=log_file)
        
        # Get read coverage stats
        samtools_idxstats_command = f"samtools idxstats {bam_file} > {idxstats_file}"
        with open(log_file_path, 'a') as log_file:
            subprocess.call(samtools_idxstats_command, shell=True, stdout=log_file, stderr=log_file)
            
    print(f"Cluster {cluster_number} processing completed.")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Process clusters and align reads to references.')
    parser.add_argument('--dirt_base', required=True, help='Base directory path for alignment outputs')
    parser.add_argument('--dirt_fastq', required=True, help='Directory path for FASTQ files')
    parser.add_argument('--dirt_bins', required=True, help='Directory path for bins (spades_output)')
    parser.add_argument('-j', '--json', required=True, help='Input JSON file with cluster data')
    parser.add_argument('--threads', type=int, default=15, help='Number of threads to use (default: 15)')
    args = parser.parse_args()
    
    # Convert paths to absolute paths and ensure they end with '/'
    dirt_base = os.path.abspath(args.dirt_base)
    if not dirt_base.endswith('/'):
        dirt_base += '/'
    
    dirt_fastq = os.path.abspath(args.dirt_fastq)
    if not dirt_fastq.endswith('/'):
        dirt_fastq += '/'
    
    dirt_bins = os.path.abspath(args.dirt_bins)
    if not dirt_bins.endswith('/'):
        dirt_bins += '/'
    
    # Assign the JSON file
    json_file = args.json
    
    # Set up the alignment base directory (align_base)
    align_base = dirt_base
    
    # Create the align_base directory if it doesn't exist
    if not os.path.exists(align_base):
        os.makedirs(align_base)
    
    # Set up the log file path
    log_file_path = os.path.join(dirt_base, 'pipeline.log')
    
    # Load the cluster information from JSON
    with open(json_file, 'r') as input_file:
        cluster_info = json.load(input_file)
    
    # Use ThreadPoolExecutor to process clusters in parallel
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_cluster = {executor.submit(process_cluster, cluster_number): cluster_number for cluster_number in cluster_info.keys()}
        for future in as_completed(future_to_cluster):
            cluster_number = future_to_cluster[future]
            try:
                future.result()
            except Exception as exc:
                print(f'{cluster_number} generated an exception: {exc}')
                with open(log_file_path, 'a') as log_file:
                    log_file.write(f'{cluster_number} generated an exception: {exc}\n')


