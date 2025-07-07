#!/usr/bin/env python3

import csv
from collections import defaultdict

# Define the input and output file paths
input_file = 'gtdbtk_2.summary.tsv'
output_file = 'annotated_bins.tsv'
checkm_file = 'qa_2.txt'  # Adjust the path as necessary

# Define the taxonomic ranks in order
ranks = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
rank_names = {
    'd__': 'Domain',
    'p__': 'Phylum',
    'c__': 'Class',
    'o__': 'Order',
    'f__': 'Family',
    'g__': 'Genus',
    's__': 'Species'
}

# Function to parse taxonomy string into a dictionary
def parse_taxonomy(taxonomy_str):
    taxonomy = {}
    for taxon in taxonomy_str.strip().split(';'):
        if '__' in taxon:
            prefix, name = taxon.split('__', 1)
            taxonomy[prefix + '__'] = name
    return taxonomy

# Function to read CheckM results and return a dictionary
def read_checkm_results(checkm_file):
    checkm_data = {}
    with open(checkm_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        # Get the indices of Bin Id, Completeness, Contamination
        bin_id_idx = headers.index('Bin Id')
        completeness_idx = headers.index('Completeness')
        contamination_idx = headers.index('Contamination')
        for row in reader:
            bin_id = row[bin_id_idx]
            completeness = float(row[completeness_idx])
            contamination = float(row[contamination_idx])
            checkm_data[bin_id] = {
                'completeness': completeness,
                'contamination': contamination
            }
    return checkm_data

# Read the CheckM results
checkm_data = read_checkm_results(checkm_file)

# First pass: Read the input file and count occurrences of taxa at each rank
taxon_counts = {rank: defaultdict(int) for rank in ranks}

with open(input_file, 'r') as infile:
    reader = csv.reader(infile, delimiter='\t')
    next(reader)  # Skip header
    for row in reader:
        taxonomy_str = row[1]
        taxonomy = parse_taxonomy(taxonomy_str)
        for rank in ranks:
            taxon = taxonomy.get(rank, '')
            if taxon:
                taxon_counts[rank][taxon] += 1

# Initialize counts per taxonomic level for 'sp. N' annotations
taxon_sp_counts = {rank: defaultdict(int) for rank in ranks}

# Second pass: Annotate each bin using the counts and new logic
with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')

    # Write header for the output file
    writer.writerow(['bin_name', 'taxonomy', 'annotate', 'level'])

    # Skip header line in input file
    next(reader)

    for row in reader:
        bin_name = row[0]
        taxonomy_str = row[1]
        taxonomy = parse_taxonomy(taxonomy_str)

        annotate = ''
        level = ''

        # Get completeness and contamination from CheckM data
        completeness = checkm_data.get(bin_name, {}).get('completeness', 0)
        contamination = checkm_data.get(bin_name, {}).get('contamination', 100)

        # If species-level annotation exists, use it directly
        if 's__' in taxonomy and taxonomy['s__']:
            annotate = taxonomy['s__']
            level = 'Species'
        else:
            # New higher-priority logic: Check for bins meeting CheckM criteria
            if completeness > 50 and contamination < 10:
                # Find the lowest available taxonomic rank for annotation
                for rank in ['g__', 'f__', 'o__', 'c__', 'p__']:
                    taxon = taxonomy.get(rank, '')
                    if taxon:
                        # Increment taxon count for sp. N annotation
                        taxon_sp_counts[rank][taxon] += 1
                        sp_number = taxon_sp_counts[rank][taxon]
                        annotate = f"{taxon} sp. {sp_number}"
                        level = 'Species'
                        break
                else:
                    # No higher-level taxon assigned
                    taxon_sp_counts['Unclassified']['Unclassified'] += 1
                    sp_number = taxon_sp_counts['Unclassified']['Unclassified']
                    annotate = f"Unclassified sp. {sp_number}"
                    level = 'Unclassified'
            else:
                # Existing logic for bins not meeting CheckM criteria
                # Start from species level and move up
                for rank in reversed(ranks):  # Start from 's__' and move up
                    taxon = taxonomy.get(rank, '')
                    if taxon:
                        count = taxon_counts[rank][taxon]
                        rank_name = rank_names[rank]
                        if count == 1:
                            # Unique annotation
                            annotate = taxon
                            level = rank_name
                        else:
                            annotate = f"{taxon} spp."
                            level = rank_name
                        break
                else:
                    # If all ranks are empty
                    annotate = 'Unclassified'
                    level = 'Unclassified'

        # Write the output row
        writer.writerow([bin_name, taxonomy_str, annotate, level])

