#!/usr/bin/env bash
set -euo pipefail

help_message() {
cat <<'USAGE'
Usage: splitting_CoSAG_clean.sh -w <workdir> -s <script_dir> [-t <threads>] [--overwrite]

Required:
  -w, --workdir STR        Scellmate workdir
  -s, --script STR         Script directory (contains bowtie_cluster.py, split_cluster8.py, process_overfit.py, generate_json.py, generate_fastq_by_json.py)
Optional:
  -t, --threads INT        Threads (default: 20)
  --overwrite              Remove splitting log and rerun
  -h, --help               Show help

Scope:
  This clean version is meant to replace splitting_CoSAG.sh.
  It follows the original bash control flow as closely as possible,
  and prepares outputs for downstream final_CoSAG_wash.sh.
USAGE
}

threads=20
workdir=""
script=""
overwrite=false

OPTS=$(getopt -o w:s:t:h --long workdir:,script:,threads:,overwrite,help -- "$@")
if [[ $? -ne 0 ]]; then
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
    --) shift; break ;;
    *) break ;;
  esac
done

if [[ -z "$workdir" || -z "$script" ]]; then
  echo "[ERROR] -w/--workdir and -s/--script are required." >&2
  exit 1
fi

workdir=$(realpath "$workdir")
script=$(realpath "$script")

# -----------------------------
# Parameters
# -----------------------------
cutoff_of_contig=1000
minimum_cell_number=10
minimum_cell_removal=4
gene_threshold_file="$workdir/05_first_QC/marker_gene_threshold"
if [[ ! -f "$gene_threshold_file" ]]; then
  echo "[ERROR] marker_gene_threshold not found: $gene_threshold_file" >&2
  exit 1
fi
threshold=$(head -n1 "$gene_threshold_file" | tr -d '\r\n ')
if [[ -z "$threshold" ]]; then
  echo "[ERROR] marker_gene_threshold is empty: $gene_threshold_file" >&2
  exit 1
fi
cell_kept_ratio_cutoff=0.6

COSAG_DIR="$workdir/06_CoSAG_assembly"
ROUND_ENDING="$COSAG_DIR/round_ending"
ROUND_ALL="$ROUND_ENDING/round_all"
TRIM_DIR="$workdir/01_trim_SAGs"
MGE_DB="$COSAG_DIR/MGE"
FILTER_AWK="$COSAG_DIR/filter_high_identity.awk"
LOGFILE="$ROUND_ENDING/Record-splitting_CoSAG.log"
CLUSTER_LOG="$COSAG_DIR/Record-clustering_CoSAG.log"
SPLIT_CLUSTER="$script/split_cluster8.py"
BOWTIE_CLUSTER="$script/bowtie_cluster.py"
PROCESS_OVERFIT="$script/process_overfit.py"
GENERATE_JSON="$script/generate_json.py"
GENERATE_FASTQ_BY_JSON="$script/generate_fastq_by_json.py"

for f in "$SPLIT_CLUSTER" "$BOWTIE_CLUSTER" "$PROCESS_OVERFIT" "$GENERATE_JSON" "$GENERATE_FASTQ_BY_JSON"; do
  [[ -f "$f" ]] || { echo "[ERROR] Missing script: $f" >&2; exit 1; }
done
[[ -f "$FILTER_AWK" ]] || { echo "[ERROR] Missing filter awk: $FILTER_AWK" >&2; exit 1; }

if [[ "$overwrite" == true ]]; then
  rm -f "$LOGFILE"
fi
mkdir -p "$ROUND_ENDING"

log_step() {
  printf '%s\t%s\n' "$1" "$(date '+%Y-%m-%d %H:%M:%S')" >> "$LOGFILE"
}
step_done() {
  local step_name="$1"
  [[ -f "$LOGFILE" ]] && cut -f1 "$LOGFILE" | grep -Fxq "$step_name"
}

msg() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

require_stop_step() {
  local stop_line stop_step
  stop_line=$(grep -n '^stop_clustering' "$CLUSTER_LOG" | head -n1 | cut -d: -f1 || true)
  [[ -n "$stop_line" ]] || { echo "[ERROR] stop_clustering not found in $CLUSTER_LOG" >&2; exit 1; }
  stop_step=$(sed -n "$((stop_line-1))p" "$CLUSTER_LOG" | cut -f1)
  [[ -n "$stop_step" ]] || { echo "[ERROR] Failed to parse stop step from $CLUSTER_LOG" >&2; exit 1; }
  echo "$stop_step"
}

STOP_STEP=$(require_stop_step)
STOP_NUM=$(echo "$STOP_STEP" | grep -oP '\d+' | head -n1)

case "$STOP_NUM" in
  3) MAX_SPLIT_ROUND=1 ;;
  4) MAX_SPLIT_ROUND=2 ;;
  5) MAX_SPLIT_ROUND=3 ;;
  *) MAX_SPLIT_ROUND=3 ;;
esac

is_round_ge() {
  local target="$1"
  [[ "$STOP_NUM" -ge "$target" ]]
}

count_group_txt() {
  local dir="$1"
  [[ -d "$dir" ]] || { echo 0; return; }
  find "$dir" -maxdepth 1 -type f -name '*group*_cells.txt' | wc -l
}

count_json() {
  local dir="$1"
  [[ -d "$dir" ]] || { echo 0; return; }
  find "$dir" -maxdepth 1 -type f -name '*.json' | wc -l
}

count_fasta() {
  local dir="$1"
  [[ -d "$dir" ]] || { echo 0; return; }
  find "$dir" -maxdepth 1 -type f -name '*.fasta' | wc -l
}

merge_jsons_in_dir() {
  local dir="$1"
  local out_json="$2"
  local out_reformat="$3"
  local n
  n=$(count_json "$dir")
  if (( n == 0 )); then
    return 1
  fi
  jq -s . "$dir"/*.json > "$out_json"
  jq 'reduce .[] as $item ({}; . + $item)' "$out_json" > "$out_reformat"
}

copy_overfit_fastas_from_json() {
  local json_file="$1"
  local output_dir="$2"
  local base_dir="$3"
  mkdir -p "$output_dir"
  grep 'round' "$json_file" | while read -r line; do
    local cluster_id round_dir file_num fasta_file
    cluster_id=$(echo "$line" | grep -oP '(round_\d+_inconsistent_\d+|round_3_\d+|round_4_\d+|round_5_\d+)' || true)
    [[ -n "$cluster_id" ]] || continue
    round_dir=$(echo "$cluster_id" | grep -oP 'round_\d+_inconsistent|round_3|round_4|round_5')
    file_num=$(echo "$cluster_id" | grep -oP '\d+$')
    fasta_file="$base_dir/$round_dir/spades_output/$file_num.fasta"
    if [[ -f "$fasta_file" ]]; then
      cp "$fasta_file" "$output_dir/${cluster_id}.fasta"
    else
      echo "[WARN] File not found: $fasta_file" >&2
    fi
  done
}

mask_mge() {
  local in_dir="$1"
  local out_dir="$2"
  local jobs="${3:-20}"
  mkdir -p "$out_dir"
  local n
  n=$(count_fasta "$in_dir")
  if (( n == 0 )); then
    msg "No fasta found in $in_dir, skipping MGE masking"
    return 0
  fi
  parallel -j "$jobs" '
    FILE_TEMP={};
    BASENAME=$(basename "$FILE_TEMP" .fasta);
    cd "'"$in_dir"'";
    blastn -query "${BASENAME}.fasta" -db "'"$MGE_DB"'" -outfmt 6 -out "${BASENAME}_results.out";
    awk -f "'"$FILTER_AWK"'" "${BASENAME}_results.out" | sort | uniq > "${BASENAME}_results.filter.id";
    filterbyname.sh substring=name ow=t include=f minlen=0 names="${BASENAME}_results.filter.id" in="${BASENAME}.fasta" out="'"$out_dir"'/${BASENAME}.fasta" >/dev/null 2>&1
  ' ::: "$in_dir"/*.fasta
}

build_1000_500_and_checkm() {
  local base="$1"
  local jobs="${2:-20}"
  local checkm_threads="${3:-35}"
  local n
  n=$(count_fasta "$base")
  if (( n == 0 )); then
    msg "No fasta found in $base, skipping reformat/checkm"
    return 0
  fi
  mkdir -p "$base/1000" "$base/500" "$base/rRNA"
  parallel -j "$jobs" '
    FILE_TEMP={};
    BASENAME=$(basename "$FILE_TEMP" .fasta);
    cd "'"$base"'";
    reformat.sh in="${BASENAME}.fasta" out="1000/${BASENAME}.1000.fasta" minlength=1000 ow=t;
    reformat.sh in="${BASENAME}.fasta" out="500/${BASENAME}.500.fasta" minlength=500 ow=t;
    barrnap -i "${BASENAME}.fasta" -o "rRNA/${BASENAME}.rRNA.fna" >/dev/null 2>&1 || true
  ' ::: "$base"/*.fasta

  local fasta500
  fasta500=$(count_fasta "$base/500")
  if (( fasta500 == 0 )); then
    msg "No *.fasta under $base/500, skipping CheckM"
    return 0
  fi

  checkm tree -x fasta "$base/500" -t "$checkm_threads" "$base/500/checkm"
  checkm lineage_set "$base/500/checkm" "$base/500/checkm/marker_file"
  checkm analyze "$base/500/checkm/marker_file" -x fasta "$base/500" "$base/500/checkm/analyze" -t "$checkm_threads"
  checkm tree_qa "$base/500/checkm" -o 2 --tab_table > "$base/500/checkm/tree_qa.txt"
  checkm qa "$base/500/checkm/marker_file" "$base/500/checkm/analyze" -o 2 -t 30 -q -f "$base/500/checkm/qa.txt" --tab_table
}

run_split_round() {
  local dirt_base="$1"
  local output_fastq="$2"
  local output_condition="$3"
  local spades_output_dir="$4"
  local overfit_flag="$5"      # yes/no
  local checkm_qa="$6"         # may be empty when overfit_flag=no
  local jobs="${7:-$threads}"

  mkdir -p "$output_fastq" "$output_condition"
  local round_dirs
  mapfile -t round_dirs < <(find "$dirt_base" -mindepth 1 -maxdepth 1 -type d -name 'round*' | sort)
  if (( ${#round_dirs[@]} == 0 )); then
    msg "No round directories found in $dirt_base, skipping split"
    return 0
  fi

  if [[ "$overfit_flag" == "yes" ]]; then
    parallel -j "$jobs" '
      sample=$(basename {});
      python "'"$SPLIT_CLUSTER"'" "$sample" \
        --dirt_base "'"$dirt_base"'" \
        --output_fastq "'"$output_fastq"'" \
        --output_condition "'"$output_condition"'" \
        --spades_output_dir "'"$spades_output_dir"'" \
        --cutoff_of_contig "'"$cutoff_of_contig"'" \
        --minimum_cell_number "'"$minimum_cell_number"'" \
        --minimum_cell_removal "'"$minimum_cell_removal"'" \
        --threshold "'"$threshold"'" \
        --cell_kept_ratio_cutoff "'"$cell_kept_ratio_cutoff"'" \
        --overfit \
        --spades_checkm "'"$checkm_qa"'"
    ' ::: "${round_dirs[@]}"
  else
    parallel -j "$jobs" '
      sample=$(basename {});
      python "'"$SPLIT_CLUSTER"'" "$sample" \
        --dirt_base "'"$dirt_base"'" \
        --output_fastq "'"$output_fastq"'" \
        --output_condition "'"$output_condition"'" \
        --cutoff_of_contig "'"$cutoff_of_contig"'" \
        --minimum_cell_number "'"$minimum_cell_number"'" \
        --minimum_cell_removal "'"$minimum_cell_removal"'" \
        --threshold "'"$threshold"'" \
        --cell_kept_ratio_cutoff "'"$cell_kept_ratio_cutoff"'"
    ' ::: "${round_dirs[@]}"
  fi
}

process_one_overfit_level() {
  # $1 current dir that contains split_output_overfit and where overfit_results should be generated
  # $2 relative current_dir passed to process_overfit.py
  # $3 path of split_output_overfit.json output
  # $4 path of split_output_overfit.reformat.json output
  # $5 overfit_results/spades_output output dir
  local current_abs="$1"
  local current_rel="$2"
  local merged_json="$3"
  local merged_reformat="$4"
  local overfit_spades_out="$5"
  local group_n

  group_n=$(count_group_txt "$current_abs/split_output_overfit")
  if (( group_n == 0 )); then
    msg "No overfit group txt files found in $current_abs/split_output_overfit"
    return 1
  fi

  python "$PROCESS_OVERFIT" --base_dir "$COSAG_DIR/" --current_dir "$current_rel"
  merge_jsons_in_dir "$current_abs/split_output_overfit" "$merged_json" "$merged_reformat" || return 1
  python "$GENERATE_FASTQ_BY_JSON" \
    --dirt_base "$current_abs" \
    --output_dir "$current_abs/fastq_split_output_overfit" \
    -j "$merged_json" \
    --dirt_fastq "$TRIM_DIR"
  copy_overfit_fastas_from_json "$merged_json" "$overfit_spades_out" "$COSAG_DIR"
  return 0
}

generate_filter_json_and_fastq() {
  local dir="$1"
  local split_round="$2"
  local split_filter="$dir/split_output_filter"
  local json_out="$dir/split_output_filter.json"
  local fastq_out="$dir/fastq_split_output_filter"
  local n

  [[ -d "$split_filter" ]] || return 1

  # Terminal round:
  # keep the surviving results, but do NOT allow *group*_cells.txt
  # to remain as candidates for another splitting round.
  if [[ "$split_round" -ge "$MAX_SPLIT_ROUND" ]]; then
    find "$split_filter" -maxdepth 1 -type f -name '*group*_cells.txt' | while read -r f; do
      bn=$(basename "$f")
      base="${bn%_cells.txt}"
      if [[ "$base" =~ ^(.*)_group_[01]$ ]]; then
        out="${BASH_REMATCH[1]}_cells.txt"
        cat "$f" >> "$split_filter/$out"
        rm -f "$f"
      fi
    done

    # Discard unclassified cells at the terminal splitting round
    find "$split_filter" -maxdepth 1 -type f -name '*unclassified*_cells.txt' -delete

    # Do not generate JSON / FASTQ for the next round once the terminal round is reached
    rm -f "$split_filter"/*.json 2>/dev/null || true
    return 1
  fi

  n=$(find "$split_filter" -maxdepth 1 -type f -name '*group*_cells.txt' | wc -l)
  if (( n == 0 )); then
    return 1
  fi

  rm -f "$split_filter"/*.json 2>/dev/null || true
  python "$GENERATE_JSON" --directory "$split_filter"
  jq -s . "$split_filter"/*.json > "$json_out"
  python "$GENERATE_FASTQ_BY_JSON" \
    --dirt_base "$dir" \
    --output_dir "$fastq_out" \
    -j "$json_out" \
    --dirt_fastq "$TRIM_DIR"
}

run_spades_from_fastq_dir() {
  local dir="$1"
  local out="$2"
  local per_job="$3"
  mkdir -p "$out"
  shopt -s nullglob
  local r1s=("$dir"/*_R1.fastq)
  shopt -u nullglob
  if (( ${#r1s[@]} == 0 )); then
    msg "No *_R1.fastq found in $dir, skipping SPAdes"
    return 0
  fi
  parallel -j "$(( threads / per_job > 0 ? threads / per_job : 1 ))" '
    FILE_TEMP={};
    BASENAME=$(basename "$FILE_TEMP" _R1.fastq);
    OUTPUT_FILE="'"$out"'/${BASENAME}.fasta";
    if [[ ! -f "$OUTPUT_FILE" ]]; then
      echo "Processing $BASENAME";
      spades.py -t '"$per_job"' --sc --careful \
        --pe1-s "'"$dir"'/${BASENAME}_R1.fastq" \
        --pe1-s "'"$dir"'/${BASENAME}_R2.fastq" \
        -o "'"$out"'/${BASENAME}_scCareful";
      mv "'"$out"'/${BASENAME}_scCareful/contigs.fasta" "$OUTPUT_FILE";
      rm -rf "'"$out"'/${BASENAME}_scCareful";
    else
      echo "Skipping $BASENAME, output exists";
    fi
  ' ::: "$dir"/*_R1.fastq
}

# ---------------------------------
# Step 1: Initial alignment + CheckM on round_all spades_output
# ---------------------------------
if step_done "Initial Alignment and JSON Generation"; then
  msg "Skipping Initial Alignment and JSON Generation"
else
  mkdir -p "$ROUND_ALL/1st_align"
  python "$BOWTIE_CLUSTER" \
    --dirt_base "$ROUND_ALL/1st_align" \
    --dirt_fastq "$TRIM_DIR" \
    --dirt_bins "$ROUND_ALL/spades_output" \
    -j "$ROUND_ALL/round_all.json" \
    --threads "$(( threads/2 > 0 ? threads/2 : 1 ))"
  log_step "Initial Alignment and JSON Generation"
fi

if step_done "Filter Assembled Contigs and Run CheckM"; then
  msg "Skipping Filter Assembled Contigs and Run CheckM"
else
  mask_mge "$ROUND_ALL/spades_output" "$ROUND_ALL/spades_output_filter" "$threads"
  build_1000_500_and_checkm "$ROUND_ALL/spades_output" "$threads" 35
  log_step "Filter Assembled Contigs and Run CheckM"
fi

# ---------------------------------
# Step 2: First split (with overfit)
# ---------------------------------
if step_done "Split Overfitted CoSAGs (1st Round)"; then
  msg "Skipping Split Overfitted CoSAGs (1st Round)"
else
  run_split_round \
    "$ROUND_ALL/1st_align" \
    "$ROUND_ALL/1st_align/output_fastq" \
    "$ROUND_ALL/1st_align/split_output" \
    "$ROUND_ALL/spades_output" \
    yes \
    "$ROUND_ALL/spades_output/500/checkm/qa.txt" \
    "$threads"
  log_step "Split Overfitted CoSAGs (1st Round)"
fi

# ---------------------------------
# Step 3: First overfit processing + 2nd overfit alignment
# ---------------------------------
FIRST_OVERFIT_OK=0
if is_round_ge 3; then
  if step_done "Process 1st-Round Overfitted Results"; then
    msg "Skipping Process 1st-Round Overfitted Results"
    if [[ -f "$ROUND_ALL/1st_align/split_output_overfit.reformat.json" ]]; then
      FIRST_OVERFIT_OK=1
    fi
  else
    if process_one_overfit_level \
      "$ROUND_ALL/1st_align" \
      "round_ending/round_all/1st_align" \
      "$ROUND_ALL/1st_align/split_output_overfit.json" \
      "$ROUND_ALL/1st_align/split_output_overfit.reformat.json" \
      "$ROUND_ALL/1st_align/overfit_results/spades_output"; then
      FIRST_OVERFIT_OK=1
      log_step "Process 1st-Round Overfitted Results"
    else
      msg "No first-round overfit outputs; Step 4/5 will be skipped"
    fi
  fi

  if (( FIRST_OVERFIT_OK == 1 )); then
    if step_done "Align SAGs to Overfitted CoSAGs (2nd Round)"; then
      msg "Skipping Align SAGs to Overfitted CoSAGs (2nd Round)"
    else
      python "$BOWTIE_CLUSTER" \
        --dirt_base "$ROUND_ALL/1st_align/overfit_results/overfit_1st_align" \
        --dirt_fastq "$TRIM_DIR" \
        --dirt_bins "$ROUND_ALL/1st_align/overfit_results/spades_output" \
        -j "$ROUND_ALL/1st_align/split_output_overfit.reformat.json" \
        --threads "$(( threads/2 > 0 ? threads/2 : 1 ))"
      log_step "Align SAGs to Overfitted CoSAGs (2nd Round)"
    fi
  fi
fi

# ---------------------------------
# Step 4: CheckM on 2nd-round overfit CoSAGs
# ---------------------------------
if is_round_ge 4; then
  local_n=$(count_fasta "$ROUND_ALL/1st_align/overfit_results/spades_output")
  if (( local_n == 0 )); then
    msg "No 2nd-round overfit CoSAG fasta found. Skipping Step 6."
  elif step_done "CheckM and Barrnap on 2nd-Round CoSAGs"; then
    msg "Skipping CheckM and Barrnap on 2nd-Round CoSAGs"
  else
    mask_mge "$ROUND_ALL/1st_align/overfit_results/spades_output" "$ROUND_ALL/1st_align/overfit_results/spades_output_filter" "$threads"
    build_1000_500_and_checkm "$ROUND_ALL/1st_align/overfit_results/spades_output" "$threads" 35
    log_step "CheckM and Barrnap on 2nd-Round CoSAGs"
  fi
fi

# ---------------------------------
# Step 5: Second split (overfit_1st_align)
# ---------------------------------
SECOND_OVERFIT_OK=0
if is_round_ge 3; then
  second_input="$ROUND_ALL/1st_align/overfit_results/overfit_1st_align"
  if (( $(count_fasta "$ROUND_ALL/1st_align/overfit_results/spades_output") == 0 )); then
    msg "No inputs for Split Overfitted CoSAGs (2nd Round), skipping"
  elif step_done "Split Overfitted CoSAGs (2nd Round)"; then
    msg "Skipping Split Overfitted CoSAGs (2nd Round)"
  else
    run_split_round \
      "$second_input" \
      "$second_input/output_fastq" \
      "$second_input/split_output" \
      "$ROUND_ALL/1st_align/overfit_results/spades_output" \
      yes \
      "$ROUND_ALL/1st_align/overfit_results/spades_output/500/checkm/qa.txt" \
      "$threads"
    log_step "Split Overfitted CoSAGs (2nd Round)"
  fi

  if step_done "Process 2nd-Round Overfitted Results"; then
    msg "Skipping Process 2nd-Round Overfitted Results"
    [[ -f "$second_input/split_output_overfit.reformat.json" ]] && SECOND_OVERFIT_OK=1
  else
    if process_one_overfit_level \
      "$second_input" \
      "round_ending/round_all/1st_align/overfit_results/overfit_1st_align" \
      "$second_input/split_output_overfit.json" \
      "$second_input/split_output_overfit.reformat.json" \
      "$second_input/overfit_results/spades_output"; then
      SECOND_OVERFIT_OK=1
      log_step "Process 2nd-Round Overfitted Results"
    else
      msg "No second-round overfit outputs; Step 8/9 will be skipped"
    fi
  fi

  if is_round_ge 4 && (( SECOND_OVERFIT_OK == 1 )); then
    if step_done "Align SAGs to Overfitted CoSAGs (3rd Round)"; then
      msg "Skipping Align SAGs to Overfitted CoSAGs (3rd Round)"
    else
      python "$BOWTIE_CLUSTER" \
        --dirt_base "$second_input/overfit_results/overfit_2nd_align" \
        --dirt_fastq "$TRIM_DIR" \
        --dirt_bins "$second_input/overfit_results/spades_output" \
        -j "$second_input/split_output_overfit.reformat.json" \
        --threads "$(( threads/2 > 0 ? threads/2 : 1 ))"
      log_step "Align SAGs to Overfitted CoSAGs (3rd Round)"
    fi
  fi
fi

# ---------------------------------
# Step 6: CheckM on 3rd-round overfit CoSAGs
# ---------------------------------
if is_round_ge 5; then
  third_overfit_base="$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output"
  if (( $(count_fasta "$third_overfit_base") == 0 )); then
    msg "No 3rd-round overfit CoSAG fasta found. Skipping Step 10."
  elif step_done "CheckM and Barrnap on 3rd-Round CoSAGs"; then
    msg "Skipping CheckM and Barrnap on 3rd-Round CoSAGs"
  else
    mask_mge "$third_overfit_base" "$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output_filter" "$threads"
    build_1000_500_and_checkm "$third_overfit_base" "$threads" 40
    log_step "CheckM and Barrnap on 3rd-Round CoSAGs"
  fi
fi

# ---------------------------------
# Step 7: Third split (overfit_2nd_align)
# ---------------------------------
THIRD_OVERFIT_OK=0
if is_round_ge 4; then
  third_input="$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align"
  third_spades="$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/overfit_results/spades_output"
  if (( $(count_fasta "$third_spades") == 0 )); then
    msg "No inputs for Split Overfitted CoSAGs (3rd Round), skipping"
  elif step_done "Split Overfitted CoSAGs (3rd Round)"; then
    msg "Skipping Split Overfitted CoSAGs (3rd Round)"
  else
    run_split_round \
      "$third_input" \
      "$third_input/output_fastq" \
      "$third_input/split_output" \
      "$third_spades" \
      yes \
      "$third_spades/500/checkm/qa.txt" \
      "$threads"
    log_step "Split Overfitted CoSAGs (3rd Round)"
  fi

  if is_round_ge 5; then
    if step_done "Process 3rd-Round Overfitted Results"; then
      msg "Skipping Process 3rd-Round Overfitted Results"
      [[ -f "$third_input/split_output_overfit.reformat.json" ]] && THIRD_OVERFIT_OK=1
    else
      if process_one_overfit_level \
        "$third_input" \
        "round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align" \
        "$third_input/split_output_overfit.json" \
        "$third_input/split_output_overfit.reformat.json" \
        "$third_input/overfit_results/spades_output"; then
        THIRD_OVERFIT_OK=1
        log_step "Process 3rd-Round Overfitted Results"
      else
        msg "No third-round overfit outputs; Step 13/14/15 will be skipped"
      fi
    fi

    if (( THIRD_OVERFIT_OK == 1 )); then
      if step_done "Align SAGs to Overfitted CoSAGs (4th Round)"; then
        msg "Skipping Align SAGs to Overfitted CoSAGs (4th Round)"
      else
        python "$BOWTIE_CLUSTER" \
          --dirt_base "$third_input/overfit_results/overfit_3rd_align" \
          --dirt_fastq "$TRIM_DIR" \
          --dirt_bins "$third_input/overfit_results/spades_output" \
          -j "$third_input/split_output_overfit.reformat.json" \
          --threads "$(( threads/2 > 0 ? threads/2 : 1 ))"
        log_step "Align SAGs to Overfitted CoSAGs (4th Round)"
      fi
    fi
  fi
fi

# ---------------------------------
# Step 8: CheckM on 4th-round overfit CoSAGs
# ---------------------------------
if is_round_ge 5; then
  fourth_overfit_base="$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output"
  if (( $(count_fasta "$fourth_overfit_base") == 0 )); then
    msg "No 4th-round overfit CoSAG fasta found. Skipping Step 14."
  elif step_done "CheckM and Barrnap on 4th-Round CoSAGs"; then
    msg "Skipping CheckM and Barrnap on 4th-Round CoSAGs"
  else
    mask_mge "$fourth_overfit_base" "$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/spades_output_filter" "$threads"
    build_1000_500_and_checkm "$fourth_overfit_base" "$threads" 40
    log_step "CheckM and Barrnap on 4th-Round CoSAGs"
  fi
fi

# ---------------------------------
# Step 9: Fourth split (overfit_3rd_align) -- no further overfit processing
# ---------------------------------
if is_round_ge 5; then
  fourth_input="$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align"
  if step_done "Split Overfitted CoSAGs (4th Round)"; then
    msg "Skipping Split Overfitted CoSAGs (4th Round)"
  else
    run_split_round \
      "$fourth_input" \
      "$fourth_input/output_fastq" \
      "$fourth_input/split_output" \
      "" \
      no \
      "" \
      "$threads"
    log_step "Split Overfitted CoSAGs (4th Round)"
  fi
fi

# ---------------------------------
# Step 10: End overfit checking -> generate 2nd_align fastq
# ---------------------------------
if step_done "End overfit checking"; then
  msg "Skipping End overfit checking"
else
  # generate split_output_filter.json at each level if group files exist
  generate_filter_json_and_fastq "$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align" 4 || true
  generate_filter_json_and_fastq "$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align" 3 || true
  generate_filter_json_and_fastq "$ROUND_ALL/1st_align/overfit_results/overfit_1st_align" 2 || true
  generate_filter_json_and_fastq "$ROUND_ALL/1st_align" 1 || true

  mkdir -p "$ROUND_ALL/2nd_align"
  mapfile -t json_files < <(find \
    "$ROUND_ALL/1st_align/split_output_filter" \
    "$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/split_output_filter" \
    "$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output_filter" \
    "$ROUND_ALL/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align/split_output_filter" \
    -maxdepth 1 -type f -name '*.json' 2>/dev/null)

  if (( ${#json_files[@]} > 0 )); then
    jq -s . "${json_files[@]}" > "$ROUND_ALL/2nd_align/split_output_filter.json"
    python "$GENERATE_FASTQ_BY_JSON" \
      --dirt_base "$ROUND_ALL/2nd_align" \
      --output_dir "$ROUND_ALL/2nd_align/fastq_split_output_filter" \
      -j "$ROUND_ALL/2nd_align/split_output_filter.json" \
      --dirt_fastq "$TRIM_DIR"
  else
    touch "$ROUND_ALL/2nd_align/no_splitting.flag"
    msg "No split_output_filter json found across overfit levels; 2nd_align will be skipped"
  fi
  log_step "End overfit checking"
fi

# ---------------------------------
# Step 11: Normal 2nd align build + split
# ---------------------------------
if is_round_ge 3; then
  if step_done "Start re-co-assembly - from round 1 to round 2"; then
    msg "Skipping Start re-co-assembly - from round 1 to round 2"
  elif [[ -f "$ROUND_ALL/2nd_align/no_splitting.flag" ]]; then
    msg "Skipping re-co-assembly — no valid split reads for 1st re-coassembly"
  else
    run_spades_from_fastq_dir "$ROUND_ALL/2nd_align/fastq_split_output_filter" "$ROUND_ALL/2nd_align/spades_output" 5
    log_step "Start re-co-assembly - from round 1 to round 2"
  fi
fi

if is_round_ge 4; then
  if step_done "Start splitting prepare - round 2"; then
    msg "Skipping Start splitting prepare - round 2"
  else
    mask_mge "$ROUND_ALL/2nd_align/spades_output" "$ROUND_ALL/2nd_align/spades_output_filter" "$threads"
    if [[ -f "$ROUND_ALL/2nd_align/split_output_filter.json" ]]; then
      jq 'reduce .[] as $item ({}; . + $item)' "$ROUND_ALL/2nd_align/split_output_filter.json" > "$ROUND_ALL/2nd_align/split_output_filter.reformat.json"
      python "$BOWTIE_CLUSTER" \
        --dirt_base "$ROUND_ALL/2nd_align_result" \
        --dirt_fastq "$TRIM_DIR" \
        --dirt_bins "$ROUND_ALL/2nd_align/spades_output" \
        -j "$ROUND_ALL/2nd_align/split_output_filter.reformat.json" \
        --threads "$(( threads/2 > 0 ? threads/2 : 1 ))"
    fi
    log_step "Start splitting prepare - round 2"
  fi

  if step_done "Start splitting - round 2"; then
    msg "Skipping Start splitting - round 2"
  else
    run_split_round \
      "$ROUND_ALL/2nd_align_result" \
      "$ROUND_ALL/2nd_align_result/output_fastq" \
      "$ROUND_ALL/2nd_align_result/split_output" \
      "" \
      no \
      "" \
      "$threads"
    if generate_filter_json_and_fastq "$ROUND_ALL/2nd_align_result" 2; then
      :
    else
      mkdir -p "$ROUND_ALL/2nd_align_result/split_output_filter"
      touch "$ROUND_ALL/2nd_align_result/split_output_filter/no_splitting.flag"
      msg "No 2nd-round normal split group files found"
    fi
    log_step "Start splitting - round 2"
  fi
fi

# ---------------------------------
# Step 12: Normal 3rd align build + split
# ---------------------------------
if is_round_ge 4; then
  if step_done "Start re-co-assembly - from round 2 to round 3"; then
    msg "Skipping Start re-co-assembly - from round 2 to round 3"
  elif [[ -f "$ROUND_ALL/2nd_align_result/split_output_filter/no_splitting.flag" ]]; then
    msg "Skipping re-co-assembly — no valid split reads for 2nd re-coassembly"
  else
    run_spades_from_fastq_dir "$ROUND_ALL/2nd_align_result/fastq_split_output_filter" "$ROUND_ALL/2nd_align_result/spades_output" 5
    log_step "Start re-co-assembly - from round 2 to round 3"
  fi
fi

if is_round_ge 5; then
  if step_done "Start splitting prepare - round 3"; then
    msg "Skipping Start splitting prepare - round 3"
  else
    mask_mge "$ROUND_ALL/2nd_align_result/spades_output" "$ROUND_ALL/2nd_align_result/spades_output_filter" "$threads"
    if [[ -f "$ROUND_ALL/2nd_align_result/split_output_filter.json" ]]; then
      jq 'reduce .[] as $item ({}; . + $item)' "$ROUND_ALL/2nd_align_result/split_output_filter.json" > "$ROUND_ALL/2nd_align_result/split_output_filter.reformat.json"
      python "$BOWTIE_CLUSTER" \
        --dirt_base "$ROUND_ALL/3rd_align" \
        --dirt_fastq "$TRIM_DIR" \
        --dirt_bins "$ROUND_ALL/2nd_align_result/spades_output" \
        -j "$ROUND_ALL/2nd_align_result/split_output_filter.reformat.json" \
        --threads "$(( threads/2 > 0 ? threads/2 : 1 ))"
    fi
    log_step "Start splitting prepare - round 3"
  fi

  if step_done "Start splitting - round 3"; then
    msg "Skipping Start splitting - round 3"
  else
    run_split_round \
      "$ROUND_ALL/3rd_align" \
      "$ROUND_ALL/3rd_align/output_fastq" \
      "$ROUND_ALL/3rd_align/split_output" \
      "" \
      no \
      "" \
      "$threads"
    if generate_filter_json_and_fastq "$ROUND_ALL/3rd_align" 3; then
      :
    else
      mkdir -p "$ROUND_ALL/3rd_align/split_output_filter"
      touch "$ROUND_ALL/3rd_align/split_output_filter/no_splitting.flag"
      msg "No 3rd-round normal split group files found"
    fi
    log_step "Start splitting - round 3"
  fi
fi

if is_round_ge 5; then
  if step_done "Start re-co-assembly - from round 3 to end"; then
    msg "Skipping Start re-co-assembly - from round 3 to end"
  elif [[ -f "$ROUND_ALL/3rd_align/split_output_filter/no_splitting.flag" ]]; then
    msg "Skipping re-co-assembly — no valid split reads for 3rd re-coassembly"
  else
    run_spades_from_fastq_dir "$ROUND_ALL/3rd_align/fastq_split_output_filter" "$ROUND_ALL/3rd_align/spades_output" 5
    log_step "Start re-co-assembly - from round 3 to end"
  fi
fi

echo -e "\033[0;33m[Scellmate] Splitting CoSAGs in second QC have been completed.\033[0m"

