#!/bin/bash
set -euo pipefail
trap 'echo "[ERROR] Failed at line $?: $BASH_COMMAND" >&2' ERR

# Default value for threads
threads=20

# Parameter parsing
while [[ $# -gt 0 ]]; do
  case $1 in
    --workdir)
      workdir=$(realpath "$2")
      shift 2
      ;;
    -t|--threads)
      threads="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 --workdir <workdir> [-t <threads>]"
      exit 0
      ;;
    *)
      echo "Unknown parameter: $1"
      echo "Usage: $0 --workdir <workdir> [-t <threads>]"
      exit 1
      ;;
  esac
done

# Check required parameters
if [[ -z "$workdir" ]]; then
  echo "Error: --workdir is required"
  echo "Usage: $0 --workdir <workdir> [-t <threads>]"
  exit 1
fi

half_threads=$((threads / 2))
export workdir threads half_threads

mkdir -p "$workdir/05_first_QC"

logfile="$workdir/05_first_QC/mapping_steps.log"

log_step() {
  echo "$1" >> "$logfile"
}

step_done() {
  grep -Fxq "$1" "$logfile" 2>/dev/null
}

# step 1: reference index
if step_done "reference index"; then
  echo "[INFO] Skipping reference index"
else
  echo "[INFO] Start reference index"
  # source /home/wangyanren/conda/bin/activate microbe-seq
  # provisionally generate reference.fna
  cp "$workdir/06_CoSAG_assembly/MGE.fna" "$workdir/03_reference_db/MGE.fna"
  cat "$workdir/03_reference_db/MGE.fna" "$workdir/03_reference_db/chromosome/combine_reference_DB.fasta" > "$workdir/03_reference_db/reference.fna"
  bowtie2-build --threads "$threads" "$workdir/03_reference_db/reference.fna" "$workdir/03_reference_db/reference" > /dev/null 2>&1
  log_step "reference index"
fi

# step 2: mapping
if step_done "mapping"; then
  echo "[INFO] Skipping mapping"
else
  echo "[INFO] Start mapping"
  # source /home/wangyanren/conda/bin/activate microbe-seq
  mkdir -p "$workdir/05_first_QC/bowtie2_sam" "$workdir/05_first_QC/bowtie2_flagstat" "$workdir/05_first_QC/bowtie2_bam"
  parallel -j "$half_threads" --bar '
    sample=$(basename {} _R2_paired.fastq)
    # echo ${sample}
    bowtie2 -p 2 -x "$workdir/03_reference_db/reference" --maxins 5000 -1 "$workdir/01_trim_SAGs/${sample}_R1_paired.fastq" -2 "$workdir/01_trim_SAGs/${sample}_R2_paired.fastq" -S "$workdir/05_first_QC/bowtie2_sam/bowtie2-${sample}-Chromosome_VR_PR.sam" > /dev/null 2>&1
    samtools flagstat "$workdir/05_first_QC/bowtie2_sam/bowtie2-${sample}-Chromosome_VR_PR.sam" > "$workdir/05_first_QC/bowtie2_flagstat/bowtie2-${sample}-Chromosome_VR_PR.flagstat"
    samtools view -@ 3 -bS "$workdir/05_first_QC/bowtie2_sam/bowtie2-${sample}-Chromosome_VR_PR.sam" > "$workdir/05_first_QC/bowtie2_bam/bowtie2-${sample}-Chromosome_VR_PR.bam"
    samtools sort -@ 3 "$workdir/05_first_QC/bowtie2_bam/bowtie2-${sample}-Chromosome_VR_PR.bam" -o "$workdir/05_first_QC/bowtie2_bam/bowtie2-${sample}-Chromosome_VR_PR.sort.bam" > /dev/null 2>&1
    rm -f "$workdir/05_first_QC/bowtie2_sam/bowtie2-${sample}-Chromosome_VR_PR.sam"
    # echo "${sample} finish mapping in $(date)."
  ' ::: $(ls "$workdir/01_trim_SAGs"/*_R2_paired.fastq)
  log_step "mapping"
fi

# step 3: featureCounts
if step_done "featureCounts"; then
  echo "[INFO] Skipping featureCounts"
else
  echo "[INFO] Start featureCounts"
  mkdir -p "$workdir/05_first_QC/featureCounts"
  featureCounts \
      -p -T 30 -t gene -f -F SAF \
      -a "$workdir/03_reference_db/chromosome/metadata-SCG.txt" \
      -o "$workdir/05_first_QC/featureCounts/featureCounts.txt" \
      $(ls "$workdir/05_first_QC/bowtie2_bam"/bowtie2-*-Chromosome_VR_PR.sort.bam) > /dev/null 2>&1
  log_step "featureCounts"
fi

# step 4: featurecount filtering
awk 'BEGIN { FS="\t"; OFS="\t" }
     NR <= 2 { print $0 }
     NR > 2 {
         sum = 0;
         for (i = 7; i <= NF; i++) {
             sum += $i;
         }
         if (sum > 1) print $0;
     }' "$workdir/05_first_QC/featureCounts/featureCounts.txt" > "$workdir/05_first_QC/featureCounts/featureCounts_filtered.txt"

# step 5: mapped stat
output_file="$workdir/05_first_QC/mapped_stat.txt"
> "$output_file"
echo "sample|total_reads|total_mapped|total_mapped_percentage" > "$output_file"
process_file() {
    i=$1
    output_file=$2

    sample=$(basename "$i" .flagstat)
    SAG_id=$(echo "$sample" | cut -d'-' -f2)

    total_reads=$(grep -m 1 ' in total ' "$i" | awk '{print $1}')

    mapped_line=$(grep -m 1 ' mapped (' "$i")
    total_mapped=$(echo "$mapped_line" | awk '{print $1}')
    total_mapped_percentage=$(echo "$mapped_line" | awk -F'[()%]' '{print $2}')

    echo "$SAG_id|$total_reads|$total_mapped|$total_mapped_percentage" >> "$output_file"
}
export -f process_file
find "$workdir/05_first_QC/bowtie2_flagstat/" -name "bowtie2-*.flagstat" | parallel -j "$threads" process_file {} "$output_file"

# step 6: kraken2 stat
output_file="$workdir/05_first_QC/kraken2_stat.txt"
> "$output_file"
echo "sample|unclassified_read_count" > "$output_file"
process_report() {
    report_file=$1
    output_file=$2

    sample=$(basename "$report_file" .report)

    unclassified_read_count=$(awk '$4 == "U" {print $1}' "$report_file")

    echo "$sample|$unclassified_read_count" >> "$output_file"
}
export -f process_report
find "$workdir/02_kraken2_output/" -name "*.report" | parallel -j "$threads" process_report {} "$output_file"

echo -e "\033[0;32m[Scellmate] Mapping step in first QC have been completed.\033[0m"

