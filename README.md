# Scellmate

Scellmate is an open-reference, contamination-aware pipeline that performs two-stage QC (also include identification) against SAGs and eMGE linkage for **large-scale microbial single-cell genomic data**.

---

## Code availability

During the peer review process, Scellmate is distributed to reviewers and
collaborators via a restricted conda channel.

Complete installation are provided through the private distribution link supplied to reviewers.

A fully documented public release will be made available upon acceptance.

---

## System requirements

Operating system:
Scellmate has been tested and validated on Linux operating systems. 

Software dependencies:
Scellmate integrates the following open-source software packages for data analysis:
BBMap (v39.52), Trimmomatic (v0.40), Bowtie2 (v2.5.4), SAMtools (v1.22.1),
SPAdes (v3.13.1), CheckM (v1.2.4), FastANI (v1.34), Kraken2 (v2.17.1),
GTDB-Tk (v2.6.1), Prodigal (v2.6.3), sourmash (v4.8.14), and geNomad (v1.11.2).

All dependencies are managed and installed automatically through the provided
conda environment.

Hardware requirements:
The pipeline has been tested on a Linux workstation running Ubuntu 18.04 LTS,
equipped with dual Intel Xeon Silver 4214R CPUs (48 CPU threads total) and 754 GB of RAM.
No non-standard hardware (e.g. GPU or accelerator) is required.
Performance scales with available CPU cores.


---

## Installation

1. **Create environment**  ﻿
   ```bash
   mamba create -n Scellmate_env \
     -c <LINK_CONDA_CHANNEL> \
     -c conda-forge \
     scellmate jq

   # <LINK_CONDA_CHANNEL> should be replaced with the channel name provided to reviewers.

   conda activate Scellmate_env﻿

Installation time is approximately 15–25 minutes on a linux workstation. 
Please be patient while dependencies are being resolved and installed.

2. **Download small example datasets (including single-cell data and database)**
   ```bash
   zenodo_get -d 10.5281/zenodo.16593955

   for f in database_GTDB_related_for_test.tar.gz database.tar.gz sub_for_test.tar.gz; do
      tar -xvzf "$f"
   done
   ```

3. **Set up path for scellmate**
   ```bash
   conda activate Scellmate_env
   
   scellmate set_default_db database/                     # Set the default Scellmate database
   scellmate set_GTDB_db database_GTDB_related_for_test/  # Set the GTDB-related database (subset here for quick testing)
   scellmate set_tmp /dev/shm/                            # Set temporary directory to speed up Kraken2 processing
   
   conda deactivate && conda activate Scellmate_env
   ```

---

## Usage
After installation and path setting up, run `scellmate -h` for usage information.
```bash
scellmate -h
Usage: scellmate [module] [options]

Modules:
  end_to_end        Run full pipeline from preprocessing to linkage
  preprocess        Organize and QC input SAGs
  first_qc          Reference-based curation of SAGs
  second_qc         Co-assembly-based curation of SAGs and generation of CoSAGs
  link_eMGE         eMGE identification and SAG-eMGE linkage

General Options:
  -h, --help        Show this help message
  -v, --version     Show version

Database Path Setup:
  scellmate set_default_db /path/to/db        Set default Scellmate database path (env: SCELLMATE_DB)
  scellmate set_GTDB_db /path/to/gtdb         Set GTDB-Tk database path (env: SCELLMATE_GTDB_PATH)
  scellmate set_tmp /path/to/tmp              Set tmp filesystem path (env: TMP_FILESYSTEM_PATH)
```


**end-to-end**

Scellmate provides the `end_to_end` module to execute the complete pipeline in a single command:

   ```bash
   scellmate end_to_end -i <path/to/SAG_fastq> -o <path/to/workdir> --prefix <prefix> -t <num_threads>
   ```

Running the end-to-end pipeline on the provided test dataset with 48 threads
takes approximately **1 hour and 30 minutes** on a Linux workstation.

**module-by-module**

Scellmate also supports running each module independently for flexible control over the workflow:

1. preprocess – Organize input SAGs

   ```bash
   scellmate preprocess -i fastq/ -o <path/to/workdir/> --prefix <prefix> -t <threads>
   ```

2. first\_qc – Reference-based curation of SAGs

   ```bash
   scellmate first_qc --workdir <path/to/workdir> -o <path/to/workdir/> -t <threads>
   ```

3. second\_qc — Co-assembly-based curation & CoSAG generation

   ```bash
   scellmate second_qc --workdir <path/to/workdir> -t <threads>
   ```

4. link\_eMGE — eMGE identification & SAG–mobilome linkage

   ```bash
   scellmate link_eMGE --workdir <path/to/workdir> -t <threads>
   ```


---

## Understanding Output

### 1. `id.gtdb_rep.txt`

Located in:

```
<workdir>/03_reference_db/chromosome/id.gtdb_rep.txt
```

This file lists GTDB representative genomes that are **present in the default database and potentially represented in your SAG reads**.
It is generated during the **`first_qc`** module by mapping SAG reads to the GTDB representative genome database.

Example:

```
GCF_000012825.1
GCF_000012845.1
GCF_000428125.1
...
```


### 2. `QC_1st.tsv`

Located in:

```
<workdir>/05_first_QC/QC_1st.tsv
```

This is the main output table of the **`first_qc`** module.
It summarizes **taxonomic annotation** and **purity metrics** for each SAG after reference-based QC.

Example:

```
SAG                |Genus annotation   |Status         |rank.1 taxa ratio   |classified_read_count   |total_mapped_percentage   |marker_gene_count
scDNA_test_00001   |N/A                |Unknown        |0.515               |3.37                    |0.93                      |33
scDNA_test_00008   |Bacteroides_I      |pass-1st-QC    |0.997               |80.3                    |76.5                      |14993
scDNA_test_00068   |N/A                |fail-1st-QC    |0.857               |14.25                   |72.1                      |259
...
```

**Column descriptions:**
* Genus annotation – Genus containing the most abundant single-copy marker genes.
* Status
  * `Unknown`: too many unclassified reads or insufficient marker gene reads.
  * `fail-1st-QC`: purity (rank.1 taxa ratio) below threshold.
  * `pass-1st-QC`: passes purity threshold.
* rank.1 taxa ratio – Proportion of marker genes belonging to the most abundant genus.
* classified\_read\_count – Reads classified by Kraken2.
* total\_mapped\_percentage – Mapping ratio from reference mapping.
* marker\_gene\_count – Read number of single-copy marker genes detected.


### 3. `final.reformat.add_cluster.edit_cont.annotate.species.txt`

Located in:

```
<workdir>/06_CoSAG_assembly/round_ending/round_all/combined_full/final.reformat.add_cluster.edit_cont.annotate.species.txt
```

This tabular file is generated by the **`second_qc`** module. It links each SAG to its corresponding **CoSAG genome** after co-assembly-based curation and provides genus and species annotations.

**Column descriptions:**

* **CoSAG\_ID** – Identifier of the CoSAG genome to which the SAG belongs. This corresponds to the cluster name assigned in the CoSAG assembly process (e.g., `round_2_inconsistent_16`, `round_3_0`).
* **SAG\_ID** – Identifier of the individual SAG within the CoSAG genome.
* **Genus\_annotation** – Genus-level taxonomic annotation from the **`first_qc`** module, based on single-copy marker gene content of the SAG before co-assembly.
* **Species\_annotation** – Species-level annotation assigned to the CoSAG genome using GTDB-Tk after co-assembly.

**Example:**

```
round_2_inconsistent_16   scDNA_test_00019   Megamonas      Megamonas_funiformis
round_3_1                 scDNA_test_00078   Bacteroides_I  Bacteroides_I_graminisolvens
```

---


### 4. `presence_absence_table.txt`

Located in:

```
<workdir>/07_eMGE_linkage/work/<species_name>/presence_absence_table.txt
```

This file is generated by the **`link_eMGE`** module for each species, summarizing **SAG–eMGE linkage** within that species.

**Column descriptions:**

* **Sequence\_ID** – Identifier of the eMGE sequence. Format follows `<CoSAG_ID>_length_<bp>_cov_<coverage>`, or `<SAG_ID>_NODE_<n>_length_<bp>_cov_<coverage>` if derived from a specific SAG assembly.
* **SAG\_ID columns** – One column per SAG belonging to the species.

  * Value `1` – eMGE is present in the SAG.
  * Value `0` – eMGE is absent in the SAG.

**Example:**

```
Sequence_ID                           scDNA_test_00008  scDNA_test_00013  scDNA_test_00020
round_3_1_NODE_148_length_4819_cov_7.6        0                  0                  1
scDNA_test_00013_NODE_34_length_2452_cov_5.3  0                  1                  0
```

This table allows you to:

* Identify which SAGs within a species contain a specific eMGE.
* Explore intra-species variability in eMGE presence/absence patterns.

---

### 5. `merged_presence_absence_table.txt`

Located in:

```
<workdir>/07_eMGE_linkage/work/merged_presence_absence_table.txt
```

This file merges all **species-level presence/absence tables** into a single **species–eMGE linkage** table.

**Column descriptions:**

* **Sequence\_ID** – Identifier of the eMGE sequence (same format as above).
* **Species columns** – One column per species detected in the dataset.

  * Value `1` – eMGE is present in at least one SAG belonging to the species.
  * Value `0` – eMGE is absent from all SAGs of the species.

**Example:**

```
Sequence_ID                           Bacteroides_I_graminisolvens  Megamonas_funiformis
round_3_1_NODE_148_length_4819_cov_7.6               1                         0
scDNA_test_00013_NODE_34_length_2452_cov_5.3         1                         0
```

This table enables:

* A species-level overview of eMGE distribution across the dataset.
* Direct comparison of eMGE host range across different species.
