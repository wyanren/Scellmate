import os
import shutil
import json
import re

def safe_name(name):
    """Clean species name to create valid directory name, preserve '-' character"""
    return re.sub(r'[^\w\.\s\-]', '_', name)

# Read species_to_nocluster_deduplicated.tsv, create species to ID mapping
species_to_ids = {}
with open('species_to_nocluster_deduplicated.tsv', 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            id_, species = parts[0], ' '.join(parts[1:])
            species = safe_name(species)
            species_to_ids.setdefault(species, []).append(id_)
        else:
            print(f"Warning: line '{line}' has incorrect format")

# Read final_updated3.species.json, create species to scDNA ID mapping
with open('final.reformat.add_cluster.edit_cont.annotate.species.json', 'r') as f:
    species_to_scDNA = json.load(f)

def copy_file(src, dst):
    """Copy file if source file exists"""
    if os.path.exists(src):
        shutil.copy(src, dst)
    else:
        print(f"Warning: {src} does not exist")

# Define source directories
chromosome_dir = 'chromosome/chromosome_masking/'
fastq_dir = '../01_trim_SAGs//'

# Work directory prefix
work_dir = 'work/'

# Iterate through each species, create directories and copy files
for species in species_to_ids:
    species_dir = os.path.join(work_dir, safe_name(species))
    chromosome_output_dir = os.path.join(species_dir, 'chromosome')
    fastq_output_dir = os.path.join(species_dir, 'fastq')
    os.makedirs(chromosome_output_dir, exist_ok=True)
    os.makedirs(fastq_output_dir, exist_ok=True)
    
    # Copy fasta files to chromosome directory
    ids = species_to_ids[species]
    for id_ in ids:
        src_file = os.path.join(chromosome_dir, id_ + '.fasta')
        dst_file = os.path.join(chromosome_output_dir, id_ + '.fasta')
        copy_file(src_file, dst_file)
    
    # Copy fastq files to fastq directory
    if species in species_to_scDNA:
        scDNA_ids = species_to_scDNA[species]
        for scDNA_id in scDNA_ids:
            for read in ['1', '2']:
                src_file = os.path.join(fastq_dir, f"{scDNA_id}_R{read}_paired.fastq")
                dst_file = os.path.join(fastq_output_dir, f"{scDNA_id}_R{read}_paired.fastq")
                copy_file(src_file, dst_file)
    else:
        print(f"Warning: {species} not found in final_updated3.species.json")

