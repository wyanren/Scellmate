# Scellmate
# Scellmate

Scellmate is a contamination-aware, open-reference Snakemake pipeline for large-scale microbial single-cell genomics. It performs two-stage quality control of Single-Amplified Genomes (SAGs), builds context-specific reference databases, assembles and clusters SAGs into high-quality composite SAGs (CoSAGs), reconstructs non-redundant mobile genetic element (eMGE) catalogs, and finally links SAGs to their resident eMGEs.

---

## Table of Contents

1. [Features](#features)  
2. [Installation](#installation)  
3. [Requirements](#requirements)  
4. [Usage](#usage)  
5. [Pipeline Modules & Methods](#pipeline-modules--methods)  
   1. [1. Construct Sample-Specific Reference Database](#1-construct-sample–specific-reference-database)  
   2. [2. Filter Unknown SAGs (First-Round QC)](#2-filter-unknown-sags-first-round-qc)  
   3. [3. First-Round QC via Marker-Gene Profiles](#3-first-round-qc-via-marker-gene-profiles)  
   4. [4. Second-Round QC via Supervised Co-Assembly](#4-second-round-qc-via-supervised-co-assembly)  
   5. [5. Curation & Species Assignment of CoSAGs](#5-curation--species-assignment-of-cosags)  
   6. [6. Reconstruct eMGE Catalogue](#6-reconstruct-emge-catalogue)  
   7. [7. SAG–eMGE Linkage Analysis](#7-sag–emge-linkage-analysis)  
6. [Citation](#citation)  
7. [License](#license)  

---

## Features

- **Modular**: each step can be run independently or as part of an end-to-end workflow  
- **Contamination-aware**: dual QC stages (reference-based and co-assembly) to remove contaminated or chimeric SAGs  
- **Adaptive subsetting**: builds a context-specific Kraken2 database from GTDB plus de-novo eMGEs  
- **Iterative co-assembly**: clusters SAGs into provisional CoSAGs at 97% ANI, flags overfitting, and refines clusters  
- **Comprehensive eMGE recovery**: de-novo detection, clustering, extension with CoSAG contigs, and non-redundant catalog building  
- **Host–mobilome linkage**: stringent two-step filtering to generate high-confidence presence/absence matrices  

---

## Installation

1. **Clone the repository**  
   ```bash
   git clone https://github.com/YourOrg/Scellmate.git
   cd Scellmate
