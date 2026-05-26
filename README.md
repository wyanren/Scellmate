# Scellmate

Scellmate is an open-reference, contamination-aware pipeline that performs two-stage QC (also include identification) against SAGs and eMGE linkage for **large-scale microbial single-cell genomic data**.

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
   mamba create -y -n Scellmate_env \
     -c wyanren -c conda-forge -c bioconda \
     scellmate==0.1.1 openjdk=23 jq==1.8.1 zenodo_get==1.6.1 genomad==1.11.2 gtdbtk=2.6.1 mmseqs2==18.8cc5c setuptools=80.9.0 "dendropy>=4.1,<5"

   conda activate Scellmate_env
   conda env config vars set _JAVA_OPTIONS="-Xms32m -Xmx1g -Xss1m"
   conda env config vars set JAVA_TOOL_OPTIONS="-Xms32m -Xmx1g -Xss1m"

   sudo apt-get update
   sudo apt-get install -y bc parallel
﻿
Installation time is approximately 15–25 minutes on a linux workstation. 
Please be patient while dependencies are being resolved and installed.

2. **Download small example datasets (including single-cell data and database)**
   ```bash
   zenodo_get -r 18092394 # https://doi.org/10.5281/zenodo.18092394

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
   # Example: scellmate end_to_end -i fastq/ -o scellmate_demo/ --prefix demo --score_cutoff 0.5 -t 48
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
   scellmate first_qc --workdir <path/to/workdir> -t <threads>
   ```

3. second\_qc — Co-assembly-based curation & CoSAG generation

   ```bash
   scellmate second_qc --workdir <path/to/workdir> -t <threads>  --max_CoSAG_cont <value>
   ```

4. link\_eMGE — eMGE identification & SAG–mobilome linkage

   ```bash
   scellmate link_eMGE --workdir <path/to/workdir> -t <threads> --score_cutoff <value>
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


### 4. `candidate_linkages_<score_cutoff>score.tsv`

Located in:

```
<workdir>/07_eMGE_linkage/work/candidate_linkages_<score_cutoff>score.tsv
```

This file is generated by the **`link_eMGE`** module for each species, summarizing **SAG–eMGE linkage** within that species.

**Column descriptions:**

* **Sequence\_ID** – Identifier of the eMGE sequence. Format follows `<SAG_ID>_NODE_<n>_length_<bp>_cov_<coverage>` derived from a specific SAG assembly.
* **SAG\_ID columns** – One column per SAG belonging to the species.

  * Value `1` – eMGE is present in the SAG.
  * Value `0` – eMGE is absent in the SAG.

**Example:**

```
Sample             scDNA_demo_00003_NODE_17_length_7510_cov_6.336821  scDNA_demo_00003_NODE_51_length_3118_cov_4.636957  scDNA_demo_00003_NODE_60_length_2542_cov_5.565742  scDNA_demo_00005_NODE_28_length_3162_cov_4.369810
scDNA_demo_00003   1                                                   1                                                   1                                                   0
scDNA_demo_00005   0                                                   0                                                   0                                                   1
scDNA_demo_00008   0                                                   0                                                   0                                                   0
scDNA_demo_00009   1                                                   1                                                   0                                                   0
```

This table allows you to:

* Identify candidate linkage that eMGEs linked to each SAG across the whole dataset.
* Examine global SAG–eMGE presence/absence patterns after applying a length-weighted candidate score cutoff.

---
---

### 5. `candidate_linkages_<score_cutoff>score.refinement_result.tsv`

Located in:

```
<workdir>/07_eMGE_linkage/work/candidate_linkages_<score_cutoff>score.refinement_result.tsv
```

This file is generated by the `link_eMGE` module after applying species-enrichment refinement via Fisher's exact test. It summarizes the statistical evaluation of each candidate eMGE–species pair and records whether each linkage is retained in the final output.

**Column descriptions:**

* **eMGE** – Identifier of the eMGE sequence.
* **Focal species** – The species being tested for enrichment of candidate linkages.
* **SAG number in focal species** – Total number of SAGs assigned to the focal species.
* **Self positive (linkages)** – Number of SAGs in the focal species with a candidate linkage to the eMGE.
* **Self negative (linkages)** – Number of SAGs in the focal species without a candidate linkage.
* **Other positive (linkages)** – Number of SAGs outside the focal species with a candidate linkage.
* **Other negative (linkages)** – Number of SAGs outside the focal species without a candidate linkage.
* **P value** – One-sided p-value from Fisher's exact test evaluating species-level enrichment.
* **Fisher passed** – Whether the eMGE–species pair passed the significance threshold.
* **Singleton removed** – Whether the eMGE was flagged as a singleton (detected in only one SAG across the species) and excluded.
* **Final (kept)** – Whether the linkage is retained after all filters.
* **Result** – Summary outcome: `Pass` indicates a retained species-enriched linkage; `Fail` indicates rejection, with the reason noted in parentheses (e.g., `Fisher failed`, `Fisher failed; singleton`).

**Example:**

```
eMGE                                                        Focal species            SAG number  Self+  Self-  Other+  Other-  P value   Fisher  Singleton  Final  Result
scDNA_demo_00003_NODE_17_length_7510_cov_6.336821           Wocania_sp009498295      23          6      17     0       18      0.0225    True    False      True   Pass
scDNA_demo_00003_NODE_51_length_3118_cov_4.636957           Wocania_sp009498295      23          2      21     0       18      0.3085    False   False      False  Fail (Fisher failed)
scDNA_demo_00003_NODE_60_length_2542_cov_5.565742           Wocania_sp009498295      23          1      22     0       18      0.5610    False   False      False  Fail (Fisher failed; singleton)
```

This table enables:

* Statistical basis for each retained or rejected SAG–eMGE linkage.
* Identification of SAG–eMGE linkage excluded due to insufficient species enrichment or singleton status.
* Tracing the refinement logic from candidate co-capture signals to final species-enriched linkages.

---
