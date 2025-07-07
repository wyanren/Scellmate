#!/usr/bin/env bash
set -euo pipefail

# --------------------- Scellmate: CoSAG Round Ending + Multi-Round Overfit Processing ---------------------

help_message() {
cat <<EOF

Usage: $(basename "$0") -w <workdir> -s <script_dir> [-t <threads>] [--overwrite]

Required:
    -w, --workdir STR        Working directory
    -s, --script STR         Script directory (should contain bowtie_cluster.py, etc.)
Optional:
    -t, --threads INT        Number of parallel threads (default: 20)
    --overwrite              Force rerun all steps and overwrite log
    -h, --help               Show this help message

EOF
}

# Default parameters
threads=20
workdir=""
script=""
overwrite=false

# Parse arguments
OPTS=$(getopt -o w:s:t:h --long workdir:,script:,threads:,overwrite,help -- "$@")
if [ $? -ne 0 ]; then
  help_message >&2
  exit 1
fi
eval set -- "$OPTS"

while true; do
  case "$1" in
    -w|--workdir) workdir="$2"; shift 2 ;;
    -s|--script)  script="$2"; shift 2 ;;
    -t|--threads) threads="$2"; shift 2 ;;
    --overwrite)  overwrite=true; shift ;;
    -h|--help)    help_message; exit 0 ;;
    --)           shift; break ;;
    *)            break ;;
  esac
done

# Check required arguments
if [[ -z "$workdir" || -z "$script" ]]; then
  echo "[ERROR] -w/--workdir and -s/--script are required." >&2
  help_message >&2
  exit 1
fi

# Convert workdir and script to absolute paths
workdir=$(realpath "$workdir")
script=$(realpath "$script")

# ========== Global Parameters ==========
cutoff_of_contig=1000
minimum_cell_number=10
minimum_cell_removal=4
gene_threshold_file="$workdir/05_first_QC/marker_gene_thresold"
if [[ -f "$gene_threshold_file" ]]; then
  threshold=$(cat "$gene_threshold_file" | head -n1 | tr -d '\r\n ')
  if [[ -z "$threshold" ]]; then
    echo "[ERROR] marker_gene_thresold file is empty: $gene_threshold_file" >&2
    exit 1
  fi
else
  echo "[ERROR] marker_gene_thresold file not found: $gene_threshold_file" >&2
  exit 1
fi
cell_kept_ratio_cutoff=0.6

half_threads=$((threads / 2))
RUN_SPADES="$script/spades_traversal.py"
MGE_DB="$workdir/06_CoSAG_assembly/MGE"
FILTER_AWK="$workdir/06_CoSAG_assembly/filter_high_identity.awk"
export workdir threads cutoff_of_contig minimum_cell_number minimum_cell_removal threshold cell_kept_ratio_cutoff script half_threads RUN_SPADES MGE_DB FILTER_AWK

# ========== Stepwise Logging Setup ==========
logfile="$workdir/06_CoSAG_assembly/round_ending/Record-splitting_CoSAG.log"
if [ "$overwrite" = true ]; then
  rm -f "$logfile"
fi

log_step() {
  local step_name="$1"
  local now_date=$(date '+%Y-%m-%d %H:%M:%S')
  echo -e "$step_name\t$now_date" >> "$logfile"
}

step_done() {
  local step_name="$1"
  grep -qx "$step_name" <(cut -f1 "$logfile" 2>/dev/null)
}

# ========== Display Parameters ==========
echo "=== Parameters Used ==="
echo "workdir: $workdir"
echo "script: $script"
echo "threads: $threads"
echo "half_threads: $half_threads"
echo "cutoff_of_contig: $cutoff_of_contig"
echo "minimum_cell_number: $minimum_cell_number"
echo "minimum_cell_removal: $minimum_cell_removal"
echo "threshold: $threshold"
echo "cell_kept_ratio_cutoff: $cell_kept_ratio_cutoff"
echo "overwrite: $overwrite"
echo "logfile: $logfile"
echo "========================"

# ========== Step 1: Initial Alignment and JSON Generation ==========
if step_done "Initial Alignment and JSON Generation"; then
  echo "Skipping Initial Alignment and JSON Generation"
else
  echo "Initial Alignment and JSON Generation"
  (cd "$workdir/06_CoSAG_assembly/round_ending" && \
    python "$script/bowtie_cluster.py" \
      --dirt_base "round_all/1st_align" \
      --dirt_fastq "$workdir/01_trim_SAGs/" \
      --dirt_bins "round_all/spades_output" \
      -j "round_all/round_all.json" \
      --threads "$half_threads")
  log_step "Initial Alignment and JSON Generation"
fi

# ========== Step 2: Filter Assembled Contigs and Run CheckM ==========
if step_done "Filter Assembled Contigs and Run CheckM"; then
  echo "Skipping Filter Assembled Contigs and Run CheckM"
else
  echo "Filter Assembled Contigs and Run CheckM"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output"
  mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/1000" "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500" "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/rRNA"
  parallel -j "$threads" '
    FILE_TEMP={};
    BASENAME=$(basename $FILE_TEMP .fasta);
    echo $BASENAME;
    reformat.sh in="$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/${BASENAME}.fasta" out="$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/1000/${BASENAME}.1000.fasta" minlength=1000;
    reformat.sh in="$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/${BASENAME}.fasta" out="$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/${BASENAME}.500.fasta" minlength=500;
    # barrnap -i "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/${BASENAME}.fasta" -o "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/rRNA/${BASENAME}.rRNA.fna"
  ' ::: "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output"/*.fasta

  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500"
  checkm tree -x fasta "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500" -t "$threads" "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/checkm"
  checkm lineage_set "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/checkm" "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/checkm/marker_file"
  checkm analyze "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/checkm/marker_file" -x fasta "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500" "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/checkm/analyze" -t "$threads"
  checkm tree_qa "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/checkm" -o 2 --tab_table > "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/checkm/tree_qa.txt"
  checkm qa "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/checkm/marker_file" "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/checkm/analyze" -o 2 -t "$threads" -q -f "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/checkm/qa.txt" --tab_table
  log_step "Filter Assembled Contigs and Run CheckM"
fi

# ========== Step 3: Split Overfitted CoSAGs (1st Round) ==========
if step_done "Split Overfitted CoSAGs (1st Round)"; then
  echo "Skipping Split Overfitted CoSAGs (1st Round)"
else
  echo "Split Overfitted CoSAGs (1st Round)"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align"
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/" -maxdepth 1 -type d -name 'round*' | sed 's/\.\///g' | \
    parallel -j "$threads" '
      sample=$(basename {});
      python '"$script"'/split_cluster9.py "$sample" \
        --dirt_base '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align"' \
        --output_fastq '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/output_fastq"' \
        --output_condition '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output"' \
        --spades_output_dir '"$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output"' \
        --cutoff_of_contig '"$cutoff_of_contig"' \
        --minimum_cell_number '"$minimum_cell_number"' \
        --minimum_cell_removal '"$minimum_cell_removal"' \
        --threshold '"$threshold"' \
        --cell_kept_ratio_cutoff '"$cell_kept_ratio_cutoff"' \
        --overfit \
        --spades_checkm '"$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/500/checkm/qa.txt"'
    '
  log_step "Split Overfitted CoSAGs (1st Round)"
fi

# ========== Step 4: Process 1st-Round Overfitted Results ==========
if step_done "Process 1st-Round Overfitted Results"; then
  echo "Skipping Process 1st-Round Overfitted Results"
else
  echo "Process 1st-Round Overfitted Results"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align"
  python "$script/process_overfit.py" \
    --base_dir "$workdir/06_CoSAG_assembly/" \
    --current_dir "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align"

  jq -s . "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output_overfit"/*.json > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output_overfit.json"
  jq 'reduce .[] as $item ({}; . + $item)' "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output_overfit.json" > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output_overfit.reformat.json"

  python "$script/generate_fastq_by_json.py" \
    --dirt_base "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align" \
    --output_dir "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/fastq_split_output_overfit" \
    -j "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output_overfit.json" \
    --dirt_fastq "$workdir/01_trim_SAGs/"

  base_dir="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align"
  output_dir="$base_dir/overfit_results/spades_output"
  mkdir -p "$output_dir"

  grep 'round' "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output_overfit.json" | while read line; do
    cluster_id=$(echo "$line" | grep -oP '(round_\d+_inconsistent_\d+|round_3_\d+|round_4_\d+)')
    round_dir=$(echo "$cluster_id" | grep -oP 'round_\d+_inconsistent|round_3|round_4')
    file_num=$(echo "$cluster_id" | grep -oP '\d+$')
    fasta_file="$base_dir/../../../$round_dir/spades_output/$file_num.fasta"
    if [ -f "$fasta_file" ]; then
      cp "$fasta_file" "$output_dir/${cluster_id}.fasta"
    else
      echo "File not found: $fasta_file"
    fi
  done

  echo "[INFO] 1st-round overfit extraction complete."
  log_step "Process 1st-Round Overfitted Results"
fi

# ========== Step 5: Align SAGs to Overfitted CoSAGs (2nd Round) ==========
if step_done "Align SAGs to Overfitted CoSAGs (2nd Round)"; then
  echo "Skipping Align SAGs to Overfitted CoSAGs (2nd Round)"
else
  echo "Align SAGs to Overfitted CoSAGs (2nd Round)"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align"
  python "$script/bowtie_cluster.py" \
    --dirt_base "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align" \
    --dirt_fastq "$workdir/01_trim_SAGs/" \
    --dirt_bins "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output" \
    -j "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output_overfit.reformat.json" \
    --threads "$half_threads"
  log_step "Align SAGs to Overfitted CoSAGs (2nd Round)"
fi

# ========== Step 6: CheckM and Barrnap on 2nd-Round CoSAGs ==========
if step_done "CheckM and Barrnap on 2nd-Round CoSAGs"; then
  echo "Skipping CheckM and Barrnap on 2nd-Round CoSAGs"
else
  echo "CheckM and Barrnap on 2nd-Round CoSAGs"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output"
  mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/1000" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/rRNA"
  parallel -j 20 '
    FILE_TEMP={};
    BASENAME=$(basename $FILE_TEMP .fasta);
    echo $BASENAME;
    reformat.sh in="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/${BASENAME}.fasta" out="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/1000/${BASENAME}.1000.fasta" minlength=1000 ow=t;
    reformat.sh in="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/${BASENAME}.fasta" out="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/${BASENAME}.500.fasta" minlength=500 ow=t;
    # barrnap -i "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/${BASENAME}.fasta" -o "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/rRNA/${BASENAME}.rRNA.fna";
    blastn -query "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/${BASENAME}.fasta" -db "$MGE_DB" -outfmt 6 -out "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/${BASENAME}_results.out";
    awk -f "$FILTER_AWK" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/${BASENAME}_results.out" | sort | uniq > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/${BASENAME}_results.filter.id"
  ' ::: "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output"/*.fasta

  # cd 500
  checkm tree -x fasta "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500" -t 35 "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/checkm"
  checkm lineage_set "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/checkm" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/checkm/marker_file"
  checkm analyze "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/checkm/marker_file" -x fasta "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/checkm/analyze" -t 35
  checkm tree_qa "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/checkm" -o 2 --tab_table > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/checkm/tree_qa.txt"
  checkm qa "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/checkm/marker_file" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/checkm/analyze" -o 2 -t 30 -q -f "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/checkm/qa.txt" --tab_table
  log_step "CheckM and Barrnap on 2nd-Round CoSAGs"
fi

# ========== Step 7: Split Overfitted CoSAGs (2nd Round) ==========
if step_done "Split Overfitted CoSAGs (2nd Round)"; then
  echo "Skipping Split Overfitted CoSAGs (2nd Round)"
else
  echo "Split Overfitted CoSAGs (2nd Round)"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align"
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/" -maxdepth 1 -type d -name 'round*' | sed 's/\.\///g' | \
    parallel -j "$threads" '
      sample=$(basename {});
      python '"$script"'/split_cluster9.py "$sample" \
        --dirt_base '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align"' \
        --output_fastq '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/output_fastq"' \
        --output_condition '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/split_output"' \
        --spades_output_dir '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output"' \
        --cutoff_of_contig '"$cutoff_of_contig"' \
        --minimum_cell_number '"$minimum_cell_number"' \
        --minimum_cell_removal '"$minimum_cell_removal"' \
        --threshold '"$threshold"' \
        --cell_kept_ratio_cutoff '"$cell_kept_ratio_cutoff"' \
        --overfit \
        --spades_checkm '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/spades_output/500/checkm/qa.txt"'
    '
  log_step "Split Overfitted CoSAGs (2nd Round)"
fi

# ========== Step 8: Process 2nd-Round Overfitted Results ==========
if step_done "Process 2nd-Round Overfitted Results"; then
  echo "Skipping Process 2nd-Round Overfitted Results"
else
  echo "Process 2nd-Round Overfitted Results"
  # cd $workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align
  python "$script/process_overfit.py" \
    --base_dir "$workdir/06_CoSAG_assembly/" \
    --current_dir "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/"

  jq -s . "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/split_output_overfit"/*.json > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/split_output_overfit.json"
  jq 'reduce .[] as $item ({}; . + $item)' "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/split_output_overfit.json"  > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/split_output_overfit.reformat.json"

  base_dir="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align"
  output_dir="$base_dir/overfit_results/spades_output"
  mkdir -p "$output_dir"
  grep 'round' "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/split_output_overfit.json" | while read line; do
    cluster_id=$(echo "$line" | grep -oP '(round_\d+_inconsistent_\d+|round_3_\d+|round_4_\d+)')
    round_dir=$(echo "$cluster_id" | grep -oP 'round_\d+_inconsistent|round_3|round_4')
    file_num=$(echo "$cluster_id" | grep -oP '\d+$')
    fasta_file="$workdir/06_CoSAG_assembly/$round_dir/spades_output/$file_num.fasta"
    if [ -f "$fasta_file" ]; then
      cp "$fasta_file" "$output_dir/${cluster_id}.fasta"
    else
      echo "File not found: $fasta_file"
    fi
  done

  echo "[INFO] 2nd-round overfit extraction complete."
  log_step "Process 2nd-Round Overfitted Results"
fi

# ========== Step 9: Align SAGs to Overfitted CoSAGs (3rd Round) ==========
if step_done "Align SAGs to Overfitted CoSAGs (3rd Round)"; then
  echo "Skipping Align SAGs to Overfitted CoSAGs (3rd Round)"
else
  echo "Align SAGs to Overfitted CoSAGs (3rd Round)"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align"

  python "$script/bowtie_cluster.py" \
    --dirt_base "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align" \
    --dirt_fastq "$workdir/01_trim_SAGs/" \
    --dirt_bins "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output" \
    -j "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/split_output_overfit.reformat.json" \
    --threads "$half_threads"
  log_step "Align SAGs to Overfitted CoSAGs (3rd Round)"
fi

# ========== Step 10: CheckM and Barrnap on 3rd-Round CoSAGs ==========
if step_done "CheckM and Barrnap on 3rd-Round CoSAGs"; then
  echo "Skipping CheckM and Barrnap on 3rd-Round CoSAGs"
else
  echo "CheckM and Barrnap on 3rd-Round CoSAGs"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output"
  mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/1000" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/rRNA"
  parallel -j 20 '
    FILE_TEMP={};
    BASENAME=$(basename $FILE_TEMP .fasta);
    echo $BASENAME;
    reformat.sh in="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/${BASENAME}.fasta" out="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/1000/${BASENAME}.1000.fasta" minlength=1000;
    reformat.sh in="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/${BASENAME}.fasta" out="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/${BASENAME}.500.fasta" minlength=500;
    # barrnap -i "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/${BASENAME}.fasta" -o "$workdir/0
    blastn -query "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/${BASENAME}.fasta" -db "$MGE_DB" -outfmt 6 -out "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/${BASENAME}_results.out";
    awk -f "$FILTER_AWK" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/${BASENAME}_results.out" | sort | uniq > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/${BASENAME}_results.filter.id"
  ' ::: "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output"/*.fasta

  # cd 500
  checkm tree -x fasta "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500" -t 40 "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/checkm"
  checkm lineage_set "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/checkm" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/checkm/marker_file"
  checkm analyze "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/checkm/marker_file" -x fasta "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/checkm/analyze" -t 40
  checkm tree_qa "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/checkm" -o 2 --tab_table > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/checkm/tree_qa.txt"
  checkm qa "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/checkm/marker_file" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/checkm/analyze" -o 2 -t 30 -q -f "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/checkm/qa.txt" --tab_table
  log_step "CheckM and Barrnap on 3rd-Round CoSAGs"
fi

# ========== Step 11: Split Overfitted CoSAGs (3rd Round) ==========
if step_done "Split Overfitted CoSAGs (3rd Round)"; then
  echo "Skipping Split Overfitted CoSAGs (3rd Round)"
else
  echo "Split Overfitted CoSAGs (3rd Round)"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/"
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/" -maxdepth 1 -type d -name 'round*' | sed 's/\.\///g' | \
    parallel -j "$threads" '
      sample=$(basename {});
      python '"$script"'/split_cluster9.py "$sample" \
        --dirt_base '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align"' \
        --output_fastq '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/output_fastq"' \
        --output_condition '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output"' \
        --spades_output_dir '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output"' \
        --cutoff_of_contig '"$cutoff_of_contig"' \
        --minimum_cell_number '"$minimum_cell_number"' \
        --minimum_cell_removal '"$minimum_cell_removal"' \
        --threshold '"$threshold"' \
        --cell_kept_ratio_cutoff '"$cell_kept_ratio_cutoff"' \
        --overfit \
        --spades_checkm '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output/500/checkm/qa.txt"'
    '
  log_step "Split Overfitted CoSAGs (3rd Round)"
fi

# ========== Step 12: Process 3rd-Round Overfitted Results ==========
if step_done "Process 3rd-Round Overfitted Results"; then
  echo "Skipping Process 3rd-Round Overfitted Results"
else
  echo "Process 3rd-Round Overfitted Results"
  # cd $workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/
  python "$script/process_overfit.py" \
    --base_dir "$workdir/06_CoSAG_assembly/" \
    --current_dir "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/"

  jq -s . "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output_overfit"/*.json > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output_overfit.json"
  jq 'reduce .[] as $item ({}; . + $item)' "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output_overfit.json" > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output_overfit.reformat.json"

  base_dir="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align"
  output_dir="$base_dir/overfit_results/spades_output"
  mkdir -p "$output_dir"
  grep 'round' "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output_overfit.json" | while read line; do
    cluster_id=$(echo "$line" | grep -oP '(round_\d+_inconsistent_\d+|round_3_\d+|round_4_\d+)')
    round_dir=$(echo "$cluster_id" | grep -oP 'round_\d+_inconsistent|round_3|round_4')
    file_num=$(echo "$cluster_id" | grep -oP '\d+$')
    fasta_file="$workdir/06_CoSAG_assembly/$round_dir/spades_output/$file_num.fasta"
    if [ -f "$fasta_file" ]; then
      cp "$fasta_file" "$output_dir/${cluster_id}.fasta"
    else
      echo "File not found: $fasta_file"
    fi
  done

  echo "[INFO] 3rd-round overfit extraction complete."
  log_step "Process 3rd-Round Overfitted Results"
fi

# ========== Step 13: Align SAGs to Overfitted CoSAGs (4th Round) ==========
if step_done "Align SAGs to Overfitted CoSAGs (4th Round)"; then
  echo "Skipping Align SAGs to Overfitted CoSAGs (4th Round)"
else
  echo "Align SAGs to Overfitted CoSAGs (4th Round)"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align"
  python "$script/bowtie_cluster.py" \
    --dirt_base "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align" \
    --dirt_fastq "$workdir/01_trim_SAGs/" \
    --dirt_bins "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output" \
    -j "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output_overfit.reformat.json" \
    --threads 20
  log_step "Align SAGs to Overfitted CoSAGs (4th Round)"
fi

# ========== Step 14: CheckM and Barrnap on 4th-Round CoSAGs ==========
if step_done "CheckM and Barrnap on 4th-Round CoSAGs"; then
  echo "Skipping CheckM and Barrnap on 4th-Round CoSAGs"
else
  echo "CheckM and Barrnap on 4th-Round CoSAGs"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/"
  mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/1000" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/rRNA"
  parallel -j 20 '
    FILE_TEMP={};
    BASENAME=$(basename $FILE_TEMP .fasta);
    echo $BASENAME;
    reformat.sh in="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/${BASENAME}.fasta" out="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/1000/${BASENAME}.1000.fasta" minlength=1000;
    reformat.sh in="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/${BASENAME}.fasta" out="$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500/${BASENAME}.500.fasta" minlength=500;
    # barrnap -i "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/${BASENAME}.fasta" -o "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/rRNA/${BASENAME}.rRNA.fna"
    blastn -query "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/${BASENAME}.fasta" -db "$MGE_DB" -outfmt 6 -out "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/${BASENAME}_results.out";
    awk -f "$FILTER_AWK" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/${BASENAME}_results.out" | sort | uniq > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/${BASENAME}_results.filter.id"
  ' ::: "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/"/*.fasta

  # cd 500
  # checkm tree -x fasta "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500" -t 40 "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500/checkm"
  # checkm lineage_set "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500/checkm" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500/checkm/marker_file"
  # checkm analyze "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500/checkm/marker_file" -x fasta "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500/checkm/analyze" -t 40
  # checkm tree_qa "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500/checkm" -o 2 --tab_table > "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500/checkm/tree_qa.txt"
  # checkm qa "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500/checkm/marker_file" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500/checkm/analyze" -o 2 -t 30 -q -f "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output/500/checkm/qa.txt" --tab_table
  log_step "CheckM and Barrnap on 4th-Round CoSAGs"
fi

# ========== Step 15: Split Overfitted CoSAGs (4th Round) ==========
if step_done "Split Overfitted CoSAGs (4th Round)"; then
  echo "Skipping Split Overfitted CoSAGs (4th Round)"
else
  echo "Split Overfitted CoSAGs (4th Round)"
  # cd "${workdir}/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align"

  # no overfit
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align/" -maxdepth 1 -type d -name 'round*' | sed 's/\.\///g' | \
    parallel -j "$threads" '
      sample=$(basename {});
      python '"$script"'/split_cluster9.py "$sample" \
        --dirt_base "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align/" \
        --output_fastq "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align/output_fastq" \
        --output_condition "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align/split_output" \
        --spades_output_dir '"$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output"' \
        --cutoff_of_contig '"$cutoff_of_contig"' \
        --minimum_cell_number '"$minimum_cell_number"' \
        --minimum_cell_removal '"$minimum_cell_removal"' \
        --threshold '"$threshold"' \
        --cell_kept_ratio_cutoff '"$cell_kept_ratio_cutoff"' \
    '
  log_step "Split Overfitted CoSAGs (4th Round)"
fi

# ========== Step 16: End overfit checking ==========
if step_done "End overfit checking"; then
  echo "Skipping End overfit checking"
else
  echo "End overfit checking"
  (cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align" && \
    rm -f split_output_filter/*.json 2>/dev/null || true && \
    python "$script/generate_json.py" --directory split_output_filter/ && \
    jq -s . split_output_filter/*.json > split_output_filter.json && \
    echo "Splitting 3rd Overfitting Collecting")

  (cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align" && \
    rm -f split_output_filter/*.json 2>/dev/null || true && \
    python "$script/generate_json.py" --directory split_output_filter/ && \
    jq -s . split_output_filter/*.json > split_output_filter.json && \
    echo "Splitting 2nd Overfitting Collecting")

  (cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align" && \
    rm -f split_output_filter/*.json 2>/dev/null || true && \
    python "$script/generate_json.py" --directory split_output_filter/ && \
    jq -s . split_output_filter/*.json > split_output_filter.json && \
    echo "Splitting 1st Overfitting Collecting")

  (cd "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align" && \
    rm -f split_output_filter/*.json 2>/dev/null || true && \
    python "$script/generate_json.py" --directory split_output_filter/ && \
    jq -s . split_output_filter/*.json > split_output_filter.json && \
    echo "1st Splitting Collecting")

  (cd "$workdir/06_CoSAG_assembly/round_ending/round_all" && \
    mkdir 2nd_align 2>/dev/null || true && \
    cd 2nd_align && \
    jq -s . "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output_filter/"*.json \
            "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/split_output_filter/"*.json \
            "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output_filter/"*.json \
            "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align/split_output_filter/"*.json \
            > split_output_filter.json && \
    python "$script/generate_fastq_by_json.py" --dirt_base . --output_dir fastq_split_output_filter -j split_output_filter.json --dirt_fastq "$workdir/01_trim_SAGs/")
  log_step "End overfit checking"
fi



# ========== Step 17: Start re-co-assembly - from round 1 to round 2 ==========
if step_done "Start re-co-assembly - from round 1 to round 2"; then
  echo "Skipping Start re-co-assembly - from round 1 to round 2"
else
  echo "Start re-co-assembly - from round 1 to round 2"
  # cd "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align"
  mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align/spades_output"
  $RUN_SPADES \
    -i "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align/fastq_split_output_filter" \
    -o "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align/spades_output" \
    --threads "$threads" \
    --per-job 5 \
    --suffix_R1 _R1.fastq
  log_step "Start re-co-assembly - from round 1 to round 2"
fi

# ========== Step 18: splitting prepare - round 2 ==========
if step_done "Start splitting prepare - round 2"; then
  echo "Skipping Start splitting prepare - round 2"
else
  echo "Start splitting prepare - round 2"
  cd ${workdir}/06_CoSAG_assembly/round_ending/round_all/2nd_align
  mkdir -p ./spades_output_filter
  parallel -j "$threads" '
    FILE_TEMP={};
    BASENAME=$(basename $FILE_TEMP .fasta);
    echo $BASENAME;
    cd ./spades_output
    blastn -query ${BASENAME}.fasta -db /mnt/md0/wangyanren/workdir/2024_March/scDNA_Influent/12-scassembly-v2/MGE -outfmt 6 -out ${BASENAME}_results.out
    awk -f /mnt/md0/wangyanren/workdir/2024_March/scDNA_Influent/12-scassembly-v2/filter_high_identity.awk ${BASENAME}_results.out | sort | uniq > ${BASENAME}_results.filter.id
    filterbyname.sh substring=name ow=t include=f minlen=0 names=${BASENAME}_results.filter.id in=${BASENAME}.fasta out=../spades_output_filter/${BASENAME}.fasta
' ::: ./spades_output/*.fasta

  cd ${workdir}/06_CoSAG_assembly/round_ending/round_all/2nd_align
  jq 'reduce .[] as $item ({}; . + $item)' split_output_filter.json > split_output_filter.reformat.json
  python "$script/bowtie_cluster.py" \
      --dirt_base ../2nd_align_result \
      --dirt_fastq "$workdir/01_trim_SAGs/" \
      --dirt_bins ./spades_output \
      -j split_output_filter.reformat.json \
      --threads "$half_threads"

  log_step "Start splitting prepare - round 2"
fi

# ========== Step 19: splitting - round 2 ==========
if step_done "Start splitting - round 2"; then
  echo "Skipping Start splitting - round 2"
else
  echo "Start splitting - round 2"
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/" -maxdepth 1 -type d -name 'round*' | sed 's/\.\///g' | \
    parallel -j "$threads" '
      sample=$(basename {});
      python '"$script"'/split_cluster9.py "$sample" \
        --dirt_base '"$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/"' \
        --output_fastq '"$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/output_fastq"' \
        --output_condition '"$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/split_output"' \
        --cutoff_of_contig '"$cutoff_of_contig"' \
        --minimum_cell_number '"$minimum_cell_number"' \
        --minimum_cell_removal '"$minimum_cell_removal"' \
        --threshold '"$threshold"' \
        --cell_kept_ratio_cutoff '"$cell_kept_ratio_cutoff"' \
    '
    # --spades_output_dir ./spades_output \
    cd ${workdir}/06_CoSAG_assembly/round_ending/round_all/2nd_align_result
    rm split_output_filter/*json 2>/dev/null || true
    python "$script/generate_json.py" --directory split_output_filter/
    jq -s . split_output_filter/*.json > split_output_filter.json
    python "$script/generate_fastq_by_json.py" --dirt_base . --output_dir fastq_split_output_filter -j split_output_filter.json --dirt_fastq "$workdir/01_trim_SAGs/"
  log_step "Start splitting - round 2"
fi



# ========== Step 20: Start re-co-assembly - from round 2 to round 3 ==========
if step_done "Start re-co-assembly - from round 2 to round 3"; then
  echo "Skipping Start re-co-assembly - from round 2 to round 3"
else
  echo "Start re-co-assembly - from round 2 to round 3"
  mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/spades_output"
  $RUN_SPADES \
    -i "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/fastq_split_output_filter" \
    -o "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/spades_output" \
    --threads "$threads" \
    --per-job 5 \
    --suffix_R1 _R1.fastq
  log_step "Start re-co-assembly - from round 2 to round 3"
fi


# ========== Step 21: splitting prepare - round 3 ==========
if step_done "Start splitting prepare - round 3"; then
  echo "Skipping Start splitting prepare - round 3"
else
  echo "Start splitting prepare - round 3"
  cd ${workdir}/06_CoSAG_assembly/round_ending/round_all/2nd_align_result
  mkdir -p ./spades_output_filter
  parallel -j "$threads" '
    FILE_TEMP={};
    BASENAME=$(basename $FILE_TEMP .fasta);
    echo $BASENAME;
    cd ./spades_output
    blastn -query ${BASENAME}.fasta -db "$MGE_DB" -outfmt 6 -out ${BASENAME}_results.out
    awk -f "$FILTER_AWK" ${BASENAME}_results.out | sort | uniq > ${BASENAME}_results.filter.id
    filterbyname.sh substring=name ow=t include=f minlen=0 names=${BASENAME}_results.filter.id in=${BASENAME}.fasta out=../spades_output_filter/${BASENAME}.fasta
' ::: ./spades_output/*.fasta

  cd ${workdir}/06_CoSAG_assembly/round_ending/round_all/2nd_align_result
  jq 'reduce .[] as $item ({}; . + $item)' split_output_filter.json > split_output_filter.reformat.json
  python "$script/bowtie_cluster.py" \
      --dirt_base ../3rd_align \
      --dirt_fastq "$workdir/01_trim_SAGs/" \
      --dirt_bins ./spades_output \
      -j split_output_filter.reformat.json \
      --threads "$half_threads"

  log_step "Start splitting prepare - round 3"
fi


# ========== Step 22: splitting - round 3 ==========
if step_done "Start splitting - round 3"; then
  echo "Skipping Start splitting - round 3"
else
  echo "Start splitting - round 3"
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/" -maxdepth 1 -type d -name 'round*' | sed 's/\.\///g' | \
    parallel -j "$threads" '
      sample=$(basename {});
      python '"$script"'/split_cluster9.py "$sample" \
        --dirt_base '"$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/"' \
        --output_fastq '"$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/output_fastq"' \
        --output_condition '"$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/split_output"' \
        --cutoff_of_contig '"$cutoff_of_contig"' \
        --minimum_cell_number '"$minimum_cell_number"' \
        --minimum_cell_removal '"$minimum_cell_removal"' \
        --threshold '"$threshold"' \
        --cell_kept_ratio_cutoff '"$cell_kept_ratio_cutoff"' \
    '
    
    cd ${workdir}/06_CoSAG_assembly/round_ending/round_all/3rd_align
    rm split_output_filter/*json 2>/dev/null || true
    python "$script/generate_json.py" --directory split_output_filter/
    jq -s . split_output_filter/*.json > split_output_filter.json
    python "$script/generate_fastq_by_json.py" --dirt_base . --output_dir fastq_split_output_filter -j split_output_filter.json --dirt_fastq "$workdir/01_trim_SAGs/"
  log_step "Start splitting - round 3"
fi


# ========== Step 23: Start re-co-assembly - from round 3 to end ==========
if step_done "Start re-co-assembly - from round 3 to end"; then
  echo "Skipping Start re-co-assembly - from round 3 to end"
else
  echo "Start re-co-assembly - from round 3 to end"
  mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/spades_output"
  $RUN_SPADES \
    -i "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/fastq_split_output_filter" \
    -o "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/spades_output" \
    --threads "$threads" \
    --per-job 5 \
    --suffix_R1 _R1.fastq
  log_step "Start re-co-assembly - from round 3 to end"
fi
