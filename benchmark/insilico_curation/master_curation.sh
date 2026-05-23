# ====================================================================================================
# 1. Extract 12 genomes from GTDB r220
# ====================================================================================================
python pick_4_genus_triplets.py \
  --taxonomy_tsv /mnt/md0/DB/GTDB/database_GTDB_related_full_v220/taxonomy_r220.processed_with_headers.tsv \
  --rep_dir /mnt/md0/DB/GTDB/database_GTDB_related_full_v220/gtdb_rep/ \
  --seed 62 \
  --output_tsv selected_4_genus_triplets.no_placeholder__62.tsv

# ====================================================================================================
# 2. Unzip the genome files
# ====================================================================================================
TSV="selected_4_genus_triplets.no_placeholder__62.tsv"
OUTDIR="genome_dir_62"
mkdir -p "${OUTDIR}"
awk -F'\t' 'NR>1 {print $6}' "${TSV}" | while read -r f; do
    cp -v "$f" "${OUTDIR}/"
done
cd "${OUTDIR}" || exit 1
gunzip -v *.gz
cd ../

# ====================================================================================================
# (Optional) Check the ANI matrix
# ====================================================================================================
./fastANI_check.sh

# ====================================================================================================
# 3. Droplet generation, MDA simulation, Sequencing reads simulation, and Sequencing reads subsetting
# ====================================================================================================
bash simulation_run.sh 62 run

# ====================================================================================================
# 4. Prepare the input SAG for Scellmate
# ====================================================================================================
OUTDIR=fastq_renamed_62
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

for f in *-sub_1.fq *-sub_2.fq; do
    [ -e "$f" ] || continue
    new=$(echo "$f" | sed -E 's/-sub_1\.fq$/-R1.fastq/; s/-sub_2\.fq$/-R2.fastq/')
    cp -v "$f" "$OUTDIR/$new"
done
for f in *-sub_1.fq *-sub_2.fq; do
    [ -e "$f" ] || continue
    new=$(echo "$f" | sed -E 's/-sub_1\.fq$/-R1.fastq/; s/-sub_2\.fq$/-R2.fastq/')
    cp -v "$f" "$OUTDIR/$new"
done

# ====================================================================================================
# Backup the intermediate files
# ====================================================================================================
mkdir MDA_Amplicon-intermediate-files-62
mv GC*_* MDA_Amplicon-intermediate-files-62/
mv all_tasks.tsv MDA_Amplicon-intermediate-files-62/
mv doublet_pairs.list MDA_Amplicon-intermediate-files-62/
mv doublet_tasks.tsv MDA_Amplicon-intermediate-files-62/
mv singlet_genomes.list MDA_Amplicon-intermediate-files-62/
mv singlet_tasks.tsv MDA_Amplicon-intermediate-files-62/

# ====================================================================================================
# 5. Scellmate Stage-1 curation with fixed purity-gating cutoff
# ====================================================================================================
conda activate Scellmate_env3
scellmate preprocess -i fastq_renamed_62/ -o scellmate_workdir__seed_62_fix --prefix insilico_62 -t 30
scellmate first_qc --workdir scellmate_workdir__seed_62_fix -t 40
# MGE is used for eMGE masking during Stage-2
cp scellmate_workdir__seed_62_fix/06_CoSAG_assembly/MGE* scellmate_workdir__seed_62_fix/
rm -rf scellmate_workdir__seed_62_fix/06_*
mkdir scellmate_workdir__seed_62_fix/06_CoSAG_assembly
cp scellmate_workdir__seed_62_fix/MGE* scellmate_workdir__seed_62_fix/06_CoSAG_assembly/
rm scellmate_workdir__seed_62_fix/05_first_QC/QC_*
Rscript /mnt/md0/wangyanren/first_curation_fixed090.R --workdir ./scellmate_workdir__seed_62_fix/ --skip-unk-filter

# ====================================================================================================
# 6. Scellmate Complete curation workflow (Stage-1 + Stage-2)
# ====================================================================================================
scellmate second_qc --workdir scellmate_workdir__seed_62_fix -t 40

# ====================================================================================================
# 7. Scellmate Stage-1-only baseline
# ====================================================================================================
cd /mnt/md0/wangyanren/workdir/202603_MDA/scellmate_workdir__seed_62_fix
mv 06_CoSAG_assembly 06_CoSAG_assembly_complete
cd 05_first_QC
cp QC_1st.tsv QC_1st.tsv.backup
cp QC_non_mock-genus.tsv QC_non_mock-genus.tsv.backup
# Mark all SAGs as "reference-undetectable" SAGs
awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} {$2="N/A"; $3="pass-1st-QC"; print}' QC_1st.tsv > QC_1st.tsv.tmp && mv QC_1st.tsv.tmp QC_1st.tsv
cd /mnt/md0/wangyanren/workdir/202603_MDA/
mkdir scellmate_workdir__seed_62_fix/06_CoSAG_assembly
cp scellmate_workdir__seed_62_fix/MGE* scellmate_workdir__seed_62_fix/06_CoSAG_assembly/
scellmate second_qc --workdir scellmate_workdir__seed_62_fix -t 40
