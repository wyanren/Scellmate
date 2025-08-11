# Scellmate

Scellmate is a contamination-aware, open-reference Snakemake pipeline for large-scale microbial single-cell genomics. It performs two-stage quality control of Single-Amplified Genomes (SAGs), builds context-specific reference databases, assembles and clusters SAGs into high-quality composite SAGs (CoSAGs), reconstructs non-redundant mobile genetic element (eMGE) catalogs, and finally links SAGs to their resident eMGEs.

---

## Table of Contents

1. [Installation](#installation)  
2. [Usage](#usage)  

---

## Installation

1. **Create environment via mamba**  
   ```bash
   mamba create -n scellmate -c wyanren scellmate

2. **Download example databases (including default and testing database)**
   ```bash
   zenodo_get -d 10.5281/zenodo.16593955

   for f in database_GTDB_related_for_test.tar.gz database.tar.gz sub_for_test.tar.gz; do
      tar -xvzf "$f"
   done

3. **Set up database path**
   ```bash
   conda activate scellmate
   scellmate set_default_db database/
   scellmate set_GTDB_db database_GTDB_related_for_test/
   scellmate set_tmp /dev/shm/
   conda activate scellmate
   
## Usage

1. **preprocess – Organize input SAGs**
   ```bash
   scellmate preprocess -i fastq/ -o <path/to/workdir/> --prefix <prefix> -t <threads>

2. **first_qc – Reference-based curation of SAGs**
   ```bash
   scellmate first_qc --workdir <path/to/workdir> -o <path/to/workdir/> -t <threads>

3. **second_qc — Co-assembly-based curation & CoSAG generation**
   ```bash
   scellmate second_qc --workdir <path/to/workdir> -t <threads>

4. **link_eMGE — eMGE identification & SAG–mobilome linkage**
   ```bash
   scellmate link_eMGE --workdir <path/to/workdir> -t <threads>
