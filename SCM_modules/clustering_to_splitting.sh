#!/usr/bin/env bash
set -euo pipefail

##########################################################################################
########################     Scellmate: clustering_to_splitting.sh      ##################
##########################################################################################
# Clean version:
#   - stop round: full collection
#   - earlier rounds: collect ONLY IDs recorded in next round's overfitted.log
#   - missing overfitted.log => collect nothing (NOT full fallback)
##########################################################################################

help_message() {
  cat <<EOF
Usage: $(basename "$0") --workdir DIR

Required:
  --workdir DIR

EOF
}

workdir=""

OPTS=$(getopt -o h --long workdir:,help -- "$@")
if [[ $? -ne 0 ]]; then
  help_message
  exit 1
fi

eval set -- "$OPTS"
while true; do
  case "$1" in
    --workdir)  workdir=$2; shift 2;;
    -h|--help)  help_message; exit 0;;
    --) shift; break;;
    *) break;;
  esac
done

if [[ -z "$workdir" ]]; then
  echo "[ERROR] --workdir is required" >&2
  help_message
  exit 1
fi

workdir=$(realpath -m "$workdir")

assembly_dir="$workdir/06_CoSAG_assembly"
round_ending="$assembly_dir/round_ending"
round_all="$round_ending/round_all"
spades_out="$round_all/spades_output"
spades_out_filter="$round_all/spades_output_filter"

mkdir -p "$spades_out" "$spades_out_filter"

logfile="$assembly_dir/Record-clustering_CoSAG.log"
if [[ ! -f "$logfile" ]]; then
  echo "[ERROR] Clustering log not found: $logfile" >&2
  exit 1
fi

stop_line=$(grep -n '^stop_clustering' "$logfile" | head -n1 | cut -d: -f1 || true)
if [[ -z "$stop_line" ]]; then
  echo "[ERROR] could not find stop_clustering in $logfile" >&2
  exit 1
fi

stop_step=$(sed -n "$((stop_line-1))p" "$logfile" | cut -f1)
if [[ -z "$stop_step" ]]; then
  echo "[ERROR] failed to parse stop_step from $logfile" >&2
  exit 1
fi

echo "[INFO] stop_step = $stop_step"

# Keep ordering fixed
rounds=(round_2_inconsistent round_3 round_4 round_5)

stop_idx=-1
for i in "${!rounds[@]}"; do
  if [[ "${rounds[$i]}" == "$stop_step" ]]; then
    stop_idx=$i
    break
  fi
done
if (( stop_idx < 0 )); then
  echo "[ERROR] stop round '$stop_step' not recognized." >&2
  exit 1
fi

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

copy_fastq_and_fasta_full() {
  local r="$1"
  echo "[INFO] Full collection for $r"

  shopt -s nullglob
  for file in "$assembly_dir/$r"/*.fastq; do
    cp "$file" "$round_all/${r}_$(basename "$file")"
  done
  for file in "$assembly_dir/$r/spades_output"/*.fasta; do
    cp "$file" "$spades_out/${r}_$(basename "$file")"
  done
  shopt -u nullglob
}

copy_fastq_and_fasta_overfit_only() {
  local r="$1"
  local next_round="$2"
  local round_num="$3"

  local overfit_log="$assembly_dir/$next_round/overfitted.log"
  echo "[INFO] Overfit-only collection for $r using $overfit_log"

  if [[ ! -f "$overfit_log" ]]; then
    echo "[INFO] $overfit_log not found; collect 0 clusters for $r"
    return 0
  fi

  mapfile -t overfitted_ids < <(grep -oP "Round ${round_num} Cluster \\K\\d+" "$overfit_log" || true)

  if (( ${#overfitted_ids[@]} == 0 )); then
    echo "[INFO] No overfitted IDs found in $overfit_log for $r; collect 0 clusters"
    return 0
  fi

  printf '[INFO] %s overfitted IDs: %s\n' "$r" "${overfitted_ids[*]}"

  shopt -s nullglob
  for overfitted_id in "${overfitted_ids[@]}"; do
    for file in "$assembly_dir/$r/${overfitted_id}"_R*.fastq; do
      cp "$file" "$round_all/${r}_$(basename "$file")"
    done

    fasta_file="$assembly_dir/$r/spades_output/${overfitted_id}.fasta"
    if [[ -f "$fasta_file" ]]; then
      cp "$fasta_file" "$spades_out/${r}_$(basename "$fasta_file")"
    else
      echo "[WARN] Missing fasta for $r cluster $overfitted_id: $fasta_file" >&2
    fi
  done
  shopt -u nullglob
}

collect_json_full() {
  local r="$1"
  local infile=""

  case "$r" in
    round_2_inconsistent)
      # Your original logic used cluster_members_round_2_inconsistent.json
      infile="$assembly_dir/$r/cluster_members_${r}.json"
      ;;
    *)
      infile="$assembly_dir/$r/cluster_members_${r}.final.json"
      ;;
  esac

  if [[ ! -f "$infile" ]]; then
    echo "[ERROR] JSON source not found for $r: $infile" >&2
    exit 1
  fi

  jq 'with_entries(.key |= "'"$r"'_" + .)' \
    "$infile" > "$round_all/${r}.json"
}

collect_json_overfit_only() {
  local r="$1"
  local next_round="$2"
  local round_num="$3"
  local infile=""

  case "$r" in
    round_2_inconsistent)
      infile="$assembly_dir/$r/cluster_members_${r}.json"
      ;;
    *)
      infile="$assembly_dir/$r/cluster_members_${r}.final.json"
      ;;
  esac

  if [[ ! -f "$infile" ]]; then
    echo "[ERROR] JSON source not found for $r: $infile" >&2
    exit 1
  fi

  local overfit_log="$assembly_dir/$next_round/overfitted.log"
  if [[ ! -f "$overfit_log" ]]; then
    echo "[INFO] $overfit_log not found; writing empty JSON for $r"
    echo "{}" > "$round_all/${r}.json"
    return 0
  fi

  mapfile -t overfitted_ids < <(grep -oP "Round ${round_num} Cluster \\K\\d+" "$overfit_log" || true)
  if (( ${#overfitted_ids[@]} == 0 )); then
    echo "[INFO] No overfitted IDs found in $overfit_log for $r; writing empty JSON"
    echo "{}" > "$round_all/${r}.json"
    return 0
  fi

  ids_as_json=$(printf '%s\n' "${overfitted_ids[@]}" | jq -R . | jq -s .)

  jq --argjson ids "$ids_as_json" '
    with_entries(
      select(.key as $k | $ids | map(tostring) | index($k))
      | .key |= "'"$r"'_" + .
    )
  ' "$infile" > "$round_all/${r}.json"
}

# -------------------------------------------------------------------
# Step 1: stop round full collection
# -------------------------------------------------------------------
copy_fastq_and_fasta_full "$stop_step"
collect_json_full "$stop_step"

# -------------------------------------------------------------------
# Step 2: earlier rounds overfit-only collection
# -------------------------------------------------------------------
# Relationships:
# round_2_inconsistent <- round_3/overfitted.log  (Round 2 Cluster X)
# round_3             <- round_4/overfitted.log  (Round 3 Cluster X)
# round_4             <- round_5/overfitted.log  (Round 4 Cluster X)

if (( stop_idx >= 1 )); then
  copy_fastq_and_fasta_overfit_only "round_2_inconsistent" "round_3" "2"
  collect_json_overfit_only "round_2_inconsistent" "round_3" "2"
fi

if (( stop_idx >= 2 )); then
  copy_fastq_and_fasta_overfit_only "round_3" "round_4" "3"
  collect_json_overfit_only "round_3" "round_4" "3"
fi

if (( stop_idx >= 3 )); then
  copy_fastq_and_fasta_overfit_only "round_4" "round_5" "4"
  collect_json_overfit_only "round_4" "round_5" "4"
fi

# -------------------------------------------------------------------
# Step 3: merge round jsons to round_all.json
# Order: stop round first, then earlier rounds
# -------------------------------------------------------------------
declare -a merge_list=()

for (( idx=stop_idx; idx>=0; idx-- )); do
  r="${rounds[$idx]}"
  json_file="$round_all/${r}.json"
  if [[ -f "$json_file" ]]; then
    merge_list+=("$json_file")
  fi
done

if (( ${#merge_list[@]} == 0 )); then
  echo "[ERROR] No JSON files found to merge." >&2
  exit 1
fi

jq -s 'reduce .[] as $item ({}; . * $item)' \
  "${merge_list[@]}" > "$round_all/round_all.json"

echo "[DONE] Merged ${#merge_list[@]} round JSONs into $round_all/round_all.json"

# -------------------------------------------------------------------
# Step 4: MGE masking for spades_output -> spades_output_filter
# -------------------------------------------------------------------
shopt -s nullglob
fasta_files=("$spades_out"/*.fasta)
shopt -u nullglob

if (( ${#fasta_files[@]} == 0 )); then
  echo "[WARN] No fasta files found in $spades_out; skipping MGE masking"
else
  parallel --jobs 24 --bar '
    FILE_TEMP={};
    BASENAME=$(basename "$FILE_TEMP" .fasta);

    blastn -query "'"$spades_out"'/${BASENAME}.fasta" \
           -db "'"$assembly_dir"'/MGE" \
           -outfmt 6 \
           -out "'"$spades_out"'/${BASENAME}_results.out"

    awk -f "'"$assembly_dir"'/filter_high_identity.awk" \
        "'"$spades_out"'/${BASENAME}_results.out" \
      | sort | uniq \
      > "'"$spades_out"'/${BASENAME}_results.filter.id"

    filterbyname.sh \
      in="'"$spades_out"'/${BASENAME}.fasta" \
      out="'"$spades_out_filter"'/${BASENAME}.fasta" \
      names="'"$spades_out"'/${BASENAME}_results.filter.id" \
      include=f ow=t substring=t minlen=0 \
      > /dev/null 2>&1
  ' ::: "${fasta_files[@]}"
fi

echo -e "\033[0;33m[Scellmate] Clustered CoSAGs have been pretreated for splitting in second QC.\033[0m"
