#!/usr/bin/env python3
import json
import subprocess
import os
import pandas as pd
import numpy as np
import math
import random
import sys
from scipy.cluster.hierarchy import linkage, cut_tree
import argparse
import shutil

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Process clustering and splitting based on contig coverage.')
parser.add_argument('folder_name', type=str, help='Folder name to process.')
parser.add_argument('--dirt_base', type=str, required=True, help='Base directory for initial clustering folders.')
parser.add_argument('--output_fastq', type=str, required=True, help='Directory path to output final fasta files after splitting.')
parser.add_argument('--output_condition', type=str, required=True, help='Directory path to output condition files.')
parser.add_argument('--cutoff_of_contig', type=int, default=1000, help='Contigs shorter than this length will be removed.')
parser.add_argument('--minimum_cell_number', type=int, default=10, help='Cell clusters with no more than this number of barcodes will be removed.')
parser.add_argument('--minimum_cell_removal', type=int, default=4, help='Minimum number of cells to consider when removing cells from a bin.')
parser.add_argument('--threshold', type=float, default=0.95, help='Threshold used for selecting clean cells in clustering.')
parser.add_argument('--cell_kept_ratio_cutoff', type=float, default=0.6, help='Proportion of cells to be kept after splitting.')
parser.add_argument('--overfit', action='store_true', help='Enable overfit logic')
parser.add_argument('--spades_output_dir', type=str, required=False, help='Directory path to spades_output.')
parser.add_argument('--spades_checkm', type=str, required=False, help='Path to the qa.txt file for CheckM analysis.')  # New argument

args = parser.parse_args()

# Get the arguments
folder_name = args.folder_name
dirt_base = args.dirt_base
output_fastq = args.output_fastq
output_condition = args.output_condition
cutoff_of_contig = args.cutoff_of_contig
minimum_cell_number = args.minimum_cell_number
minimum_cell_removal = args.minimum_cell_removal
threshold = args.threshold
cell_kept_ratio_cutoff = args.cell_kept_ratio_cutoff
spades_output_dir = args.spades_output_dir  # New argument
spades_checkm_path = args.spades_checkm  # New argument

# Define working directory based on folder_name
if '_' in folder_name:
    working_dir = os.path.abspath(os.path.join(dirt_base, 'split_output_dirt', '..', folder_name))
else:
    working_dir = os.path.abspath(os.path.join(dirt_base, folder_name))
print(f"Working directory: {working_dir}")

# Define function read_coverage_file_parse
def read_coverage_file_parse(temp_working_dirt, cutoff_of_contig, output_fastq, spades_output_dir):
    import os
    import pandas as pd
    import subprocess

    file_list = os.listdir(temp_working_dirt)
    list_of_barcodes = []
    assembly_file_name = None

    for file_temp in file_list:
        if file_temp.endswith('.fasta'):
            assembly_file_name = file_temp
        if file_temp.endswith('.bam'):
            list_of_barcodes.append(file_temp.split('.bam')[0])

    if assembly_file_name is None:
        raise ValueError("Assembly file not found in the directory.")

    # Read assembly file and extract contig information
    contig_index_list = []
    contig_length_list = []
    contig_coverage_list = []
    contig_name_list = []

    with open(os.path.join(temp_working_dirt, assembly_file_name)) as finput:
        while True:
            line_read = finput.readline()
            if len(line_read) == 0:
                break
            if line_read.startswith(">"):
                line_read = line_read.strip()[1:]  # Remove '>' and newline
                line_split_temp = line_read.split("_")
                contig_name_list.append(line_read)
                contig_index_list.append(int(line_split_temp[1]))
                contig_length_list.append(int(line_split_temp[3]))
                contig_coverage_list.append(float(line_split_temp[5]))

    # Create DataFrame with contig information
    self_mapping_df = pd.DataFrame({
        'contig_name': contig_name_list,
        'contig_index': contig_index_list,
        'contig_length': contig_length_list,
        'contig_coverage': contig_coverage_list
    })

    # Initialize dictionary for barcode data
    barcode_data = {}

    # Read read_coverage files for each barcode
    for i in list_of_barcodes:
        BC_temp = str(i)
        column_name_temp = "BC_" + BC_temp
        bam_read_count = pd.read_csv(
            os.path.join(temp_working_dirt, BC_temp + ".read_coverage"),
            sep='\t',
            header=None
        )
        bam_read_count.columns = ['contig', 'length', 'read_count', 'unpaired']
        barcode_data[column_name_temp] = bam_read_count['read_count'][:-1].values

    # Convert barcode data to DataFrame
    barcode_df = pd.DataFrame(barcode_data)

    # Merge self_mapping_df and barcode_df
    self_mapping_df = pd.concat([self_mapping_df.reset_index(drop=True), barcode_df.reset_index(drop=True)], axis=1)

    # Remove contigs based on length
    to_remove_list_length = self_mapping_df[
        self_mapping_df['contig_length'] < cutoff_of_contig
    ].index.tolist()

    # Remove contigs listed in *_results.filter.id file
    if spades_output_dir is None:
        spades_output_dir = os.path.abspath(os.path.join(temp_working_dirt, '..', 'spades_output'))
    else:
        spades_output_dir = os.path.abspath(spades_output_dir)


    assembly_base_name = os.path.splitext(assembly_file_name)[0]
    filter_id_file_name = assembly_base_name + '_results.filter.id'
    filter_id_file = os.path.join(spades_output_dir, filter_id_file_name)

    if os.path.exists(filter_id_file):
        with open(filter_id_file) as f:
            MGE_contig_list = [line.strip() for line in f if line.strip()]
    else:
        MGE_contig_list = []

    to_remove_list_MGE = self_mapping_df[
        self_mapping_df['contig_name'].isin(MGE_contig_list)
    ].index.tolist()

    # Combine contig lists to remove
    to_remove_list = list(set(to_remove_list_length + to_remove_list_MGE))

    # Remove specified contigs from DataFrame
    self_mapping_df = self_mapping_df.drop(to_remove_list).reset_index(drop=True)

    # Check if DataFrame is empty after removal
    if len(self_mapping_df) <= 1:
        subprocess.call(
            "cp " + os.path.join(temp_working_dirt, "*.fasta") + " " + output_fastq,
            shell=True
        )
        return 0

    # Save DataFrame for further processing
    self_mapping_df.to_csv(os.path.join(temp_working_dirt, 'self_mapping_df.csv'), index=False)
    return 1

# Define function cluster_contig_cleanup
def cluster_contig_cleanup(temp_working_dirt, output_condition, threshold, folder_name, cell_kept_ratio_cutoff):
    import pandas as pd
    import numpy as np
    import os

    # Read the CSV file without setting an index
    self_mapping_df = pd.read_csv(os.path.join(temp_working_dirt, 'self_mapping_df.csv'), index_col=False)

    # Make a copy of the DataFrame
    self_mapping_df_temp = self_mapping_df.copy()
    self_mapping_df_temp.reset_index(drop=True, inplace=True)

    # Ensure 'contig_name' is not used as index
    if 'contig_name' in self_mapping_df_temp.columns:
        self_mapping_df_temp.reset_index(drop=True, inplace=True)

    # Identify non-numeric columns
    non_numeric_cols = self_mapping_df_temp.select_dtypes(include=['object']).columns.tolist()

    # Columns to drop
    columns_to_drop = ['contig_name', 'contig_index', 'contig_length', 'contig_coverage'] + non_numeric_cols

    # Drop non-numeric columns
    self_mapping_df_temp.drop(columns=columns_to_drop, inplace=True, errors='ignore')

    # Convert all remaining columns to numeric type
    self_mapping_df_temp = self_mapping_df_temp.apply(pd.to_numeric, errors='coerce')

    # Calculate median
    median_list_temp = self_mapping_df_temp.median(axis=1, skipna=True)

    # Transform values based on median
    for i in range(len(self_mapping_df_temp)):
        self_mapping_df_temp.iloc[i, :] = np.where(self_mapping_df_temp.iloc[i, :] > median_list_temp.iloc[i], 1, 0)

    # Hierarchical clustering
    linked = linkage(self_mapping_df_temp, 'ward')
    cluster_number = 2
    cutree_contigs = cut_tree(linked, n_clusters=cluster_number).flatten()

    # Calculate length of each group
    group_lengths = [
        self_mapping_df.loc[cutree_contigs == i, 'contig_length'].sum() for i in range(cluster_number)
    ]
    total_length = sum(group_lengths)
    print(f"Total contig length is {total_length}, group 0 length is {group_lengths[0]}, group 1 length is {group_lengths[1]}.\n")

    # Use clustering results to select cells
    self_mapping_df_temp = self_mapping_df_temp.T
    self_mapping_df_temp['total_of_each_cell'] = self_mapping_df_temp.sum(axis=1)

    for i in range(cluster_number):
        index_of_this_cluster = np.where(cutree_contigs == i)[0]
        self_mapping_df_temp[f'hierarchy_cluster_{i}'] = self_mapping_df_temp.iloc[:, index_of_this_cluster].sum(axis=1)

    # Determine clean cells
    BC_list_0_temp = self_mapping_df_temp[
        self_mapping_df_temp[f'hierarchy_cluster_0'] / self_mapping_df_temp['total_of_each_cell'] > threshold
    ].index.tolist()
    BC_list_1_temp = self_mapping_df_temp[
        self_mapping_df_temp[f'hierarchy_cluster_1'] / self_mapping_df_temp['total_of_each_cell'] > threshold
    ].index.tolist()

    # Strip 'BC_' prefix from the cell lists for accurate comparison
    BC_list_0_temp = [cell.replace('BC_', '') for cell in BC_list_0_temp]
    BC_list_1_temp = [cell.replace('BC_', '') for cell in BC_list_1_temp]

    all_cells = [cell.replace('BC_', '') for cell in self_mapping_df_temp.index.tolist()]

    # Calculate number of unclassified cells
    BC_list_not_clean_temp = len(all_cells) - len(BC_list_0_temp) - len(BC_list_1_temp)

    print(f"Total number of cells: {len(all_cells)}, group 0 cell count: {len(BC_list_0_temp)}, group 1 cell count: {len(BC_list_1_temp)}, unclassified cells: {BC_list_not_clean_temp}.\n")

    # Print details of each group
    folder_prefix = folder_name + '_'
    os.makedirs(output_condition, exist_ok=True)

    print("Group 0 details:")
    print(f"Cell count: {len(BC_list_0_temp)}")
    print("Cell list:", BC_list_0_temp[:10], "...")  # Print first 10 for brevity

    print("Group 1 details:")
    print(f"Cell count: {len(BC_list_1_temp)}")
    print("Cell list:", BC_list_1_temp[:10], "...")  # Print first 10 for brevity

    print("Unclassified cells:")
    print(f"Cell count: {BC_list_not_clean_temp}")

    # Save details to files
    with open(os.path.join(output_condition, f'{folder_prefix}group_0_cells.txt'), 'w') as f:
        for cell in BC_list_0_temp:
            f.write(f"{cell}\n")

    with open(os.path.join(output_condition, f'{folder_prefix}group_1_cells.txt'), 'w') as f:
        for cell in BC_list_1_temp:
            f.write(f"{cell}\n")

    with open(os.path.join(output_condition, f'{folder_prefix}unclassified_cells.txt'), 'w') as f:
        for cell in all_cells:
            if cell not in BC_list_0_temp and cell not in BC_list_1_temp:
                f.write(f"{cell}\n")

    # Check if too many cells were removed
    if len(BC_list_0_temp) + len(BC_list_1_temp) < cell_kept_ratio_cutoff * len(all_cells):
        print("Too many cells were removed, be careful.")

# Define function process_overfit_logic
def process_overfit_logic(output_condition, overfit_folder, folder_prefix, spades_checkm_path, folder_name):
    import pandas as pd

    # Get list of files in output_condition
    files = os.listdir(output_condition)
    # Filter to get only *_cells.txt files related to the current folder
    cell_files = [f for f in files if f.startswith(folder_prefix) and f.endswith('_cells.txt')]
    # Create a dictionary to group files by basename
    basename_dict = {}
    for filename in cell_files:
        # Remove suffixes to get the basename
        if filename.endswith('_group_0_cells.txt'):
            basename = filename[:-len('_group_0_cells.txt')]
        elif filename.endswith('_group_1_cells.txt'):
            basename = filename[:-len('_group_1_cells.txt')]
        elif filename.endswith('_unclassified_cells.txt'):
            basename = filename[:-len('_unclassified_cells.txt')]
        else:
            continue  # Skip files that do not match expected patterns
        # Add filename to the group
        if basename not in basename_dict:
            basename_dict[basename] = []
        basename_dict[basename].append(filename)

    # Read the CheckM qa.txt file if path is provided
    contamination = None
    genome_size = None
    if spades_checkm_path:
        try:
            qa_df = pd.read_csv(spades_checkm_path, sep='\t')
            # Construct the Bin Id to match
            bin_id = f"{folder_name}.500"
            # Locate the row corresponding to the bin_id
            bin_row = qa_df[qa_df['Bin Id'] == bin_id]
            if not bin_row.empty:
                contamination = bin_row.iloc[0]['Contamination']
                genome_size = bin_row.iloc[0]['Genome size (bp)']
                print(f"Contamination for {bin_id}: {contamination}")
                print(f"Genome size for {bin_id}: {genome_size}")
            else:
                print(f"No matching Bin Id '{bin_id}' found in CheckM file.")
        except Exception as e:
            print(f"Error reading CheckM file: {e}")

    # Initialize list to keep track of overfit basenames
    overfit_basenames = []

    # Now process each group
    for basename, filenames in basename_dict.items():
        total_lines = 0
        unclassified_lines = 0
        # For each file in filenames, count the number of lines
        for filename in filenames:
            filepath = os.path.join(output_condition, filename)
            with open(filepath, 'r') as f:
                lines = f.readlines()
                num_lines = len(lines)
                total_lines += num_lines
                if filename.endswith('_unclassified_cells.txt'):
                    unclassified_lines += num_lines

        # Determine overfit based on existing conditions
        overfit_existing = (total_lines > 30) and (unclassified_lines / total_lines > 0.6)

        # Determine overfit based on CheckM conditions
        overfit_checkm = False
        if contamination is not None and genome_size is not None:
            if (contamination > 50) or (genome_size > 8000000):
                overfit_checkm = True

        # Combine conditions with OR
        if overfit_existing or overfit_checkm:
            print(f"Marking basename '{basename}' as overfit (Existing: {overfit_existing}, CheckM: {overfit_checkm})")
            # Copy the files to overfit_folder instead of moving
            for filename in filenames:
                src = os.path.join(output_condition, filename)
                dst = os.path.join(overfit_folder, filename)
                shutil.copy(src, dst)
                print(f"Copied '{filename}' to overfit folder.")
            # Add basename to overfit_basenames list
            overfit_basenames.append(basename)

    return overfit_basenames  # Return the list of overfit basenames

# Define function process_filter_logic
def process_filter_logic(
    output_condition,
    filter_folder,
    minimum_cell_number,
    minimum_cell_removal,
    folder_prefix,
    contamination_file,
    overfit_basenames,
    cell_kept_ratio_cutoff  # Ensure this parameter is included
):
    import os

    # Get list of files in output_condition
    files = os.listdir(output_condition)
    # Filter to get only *_cells.txt files related to the current folder
    cell_files = [f for f in files if f.startswith(folder_prefix) and f.endswith('_cells.txt')]
    # Create a dictionary to group files by basename
    basename_dict = {}
    for filename in cell_files:
        # Remove suffixes to get the basename
        if filename.endswith('_group_0_cells.txt'):
            basename = filename[:-len('_group_0_cells.txt')]
        elif filename.endswith('_group_1_cells.txt'):
            basename = filename[:-len('_group_1_cells.txt')]
        elif filename.endswith('_unclassified_cells.txt'):
            basename = filename[:-len('_unclassified_cells.txt')]
        else:
            continue  # Skip files that do not match expected patterns
        # Add filename to the group
        if basename not in basename_dict:
            basename_dict[basename] = []
        basename_dict[basename].append(filename)

    # Now process each group
    for basename, filenames in basename_dict.items():
        # Skip processing if basename is in overfit_basenames
        if basename in overfit_basenames:
            print(f"Skipping basename '{basename}' as it is marked as overfit.")
            continue  # Skip this basename as it's already classified as overfit

        total_lines = 0
        group0_lines = 0
        group1_lines = 0
        unclassified_lines = 0
        # Keep track of filenames
        group0_filename = None
        group1_filename = None
        unclassified_filename = None
        for filename in filenames:
            filepath = os.path.join(output_condition, filename)
            with open(filepath, 'r') as f:
                lines = f.readlines()
                num_lines = len(lines)
                total_lines += num_lines
                if filename.endswith('_group_0_cells.txt'):
                    group0_lines = num_lines
                    group0_filename = filename
                elif filename.endswith('_group_1_cells.txt'):
                    group1_lines = num_lines
                    group1_filename = filename
                elif filename.endswith('_unclassified_cells.txt'):
                    unclassified_lines = num_lines
                    unclassified_filename = filename

        # Calculate the ratio of unclassified cells
        if total_lines > 0:
            unclassified_ratio = unclassified_lines / total_lines
        else:
            unclassified_ratio = 0

        # Initialize a flag to determine if merging is needed
        should_merge = False

        # Determine if merging is required based on unclassified_ratio
        if unclassified_ratio > cell_kept_ratio_cutoff:
            print(f"Unclassified ratio {unclassified_ratio:.2f} exceeds cutoff {cell_kept_ratio_cutoff}. Merging all cells.")
            should_merge = True

        # Determine if merging is required based on total_lines
        if total_lines < minimum_cell_number:
            print(f"Total number of cells {total_lines} is less than the minimum required {minimum_cell_number}. Merging all cells.")
            should_merge = True

        if should_merge:
            # Define merged filename as *_cells.txt
            merged_filename = basename + '_cells.txt'
            merged_filepath = os.path.join(filter_folder, merged_filename)
            with open(merged_filepath, 'w') as outfile:
                for filename in filenames:
                    src = os.path.join(output_condition, filename)
                    with open(src, 'r') as infile:
                        content = infile.read()
                        outfile.write(content)
                        print(f"Copied '{filename}' into merged file '{merged_filename}'.")
            continue  # Skip further processing for this basename

        # Existing Logic: Apply filter conditions
        # Check if group0 or group1 have sample counts less than minimum_cell_removal
        if (group0_filename is not None and group0_lines < minimum_cell_removal):
            # Collect filtered-out cells into contamination.txt
            with open(contamination_file, 'a') as contam_file:
                src = os.path.join(output_condition, group0_filename)
                with open(src, 'r') as infile:
                    content = infile.read()
                    contam_file.write(content)
                    print(f"Added '{group0_filename}' to contamination file.")
            group0_filename = None  # Do not copy group0
        if (group1_filename is not None and group1_lines < minimum_cell_removal):
            # Collect filtered-out cells into contamination.txt
            with open(contamination_file, 'a') as contam_file:
                src = os.path.join(output_condition, group1_filename)
                with open(src, 'r') as infile:
                    content = infile.read()
                    contam_file.write(content)
                    print(f"Added '{group1_filename}' to contamination file.")
            group1_filename = None  # Do not copy group1
        if unclassified_filename is not None:
            # Collect unclassified cells into contamination.txt
            with open(contamination_file, 'a') as contam_file:
                src = os.path.join(output_condition, unclassified_filename)
                with open(src, 'r') as infile:
                    content = infile.read()
                    contam_file.write(content)
                    print(f"Added '{unclassified_filename}' to contamination file.")

        # Copy valid group files into filter_folder instead of moving
        for filename in [group0_filename, group1_filename]:
            if filename is not None:
                src = os.path.join(output_condition, filename)
                dst = os.path.join(filter_folder, filename)
                shutil.copy(src, dst)
                print(f"Copied '{filename}' to filter folder.")



# Main execution
if __name__ == "__main__":
    # Call function and handle the result
    flag_if_no_large_cluster = read_coverage_file_parse(
        working_dir,
        cutoff_of_contig,
        output_fastq,
        spades_output_dir
    )
    if flag_if_no_large_cluster == 0:
        print("Assembly too small, fewer than 2 contigs longer than the specified length.")
    else:
        cluster_contig_cleanup(
            working_dir,
            output_condition,
            threshold,
            folder_name,
            cell_kept_ratio_cutoff
        )
        # Create overfit and filter folders
        overfit_folder = output_condition + '_overfit'
        filter_folder = output_condition + '_filter'
        os.makedirs(overfit_folder, exist_ok=True)
        os.makedirs(filter_folder, exist_ok=True)
        folder_prefix = folder_name + '_'
        # Define contamination file path
        contamination_file = os.path.join(dirt_base, 'contamination.txt')
        # Initialize overfit_basenames list
        overfit_basenames = []
        if args.overfit:
            overfit_basenames = process_overfit_logic(
                output_condition,
                overfit_folder,
                folder_prefix,
                spades_checkm_path,
                folder_name
            )
        # Now process filter logic, excluding overfit basenames
        process_filter_logic(
            output_condition,
            filter_folder,
            minimum_cell_number,
            minimum_cell_removal,
            folder_prefix,
            contamination_file,
            overfit_basenames,
            cell_kept_ratio_cutoff  # Ensure this argument is passed
        )



