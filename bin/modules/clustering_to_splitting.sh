#!/usr/bin/env bash
set -euo pipefail

##########################################################################################
########################        Scellmate: clustering_to_splitting.sh       #########################
##########################################################################################
# Collect all the CoSAGs from the last one clustering and also select the "overfitted" CoSAG during the intermidiate CoSAG clusters (from previous clustering rounds based on logs)
# including the assembly fasta files, as well as the fastq files~

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
  echo "[ERROR] --workdir is required"
  help_message
  exit 1
fi

workdir=$(realpath -m "$workdir")
export workdir

##########################################################################################
########################        Start to sort the clustering CoSAG ~~~~       #########################
##########################################################################################
mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output"

# ====================== Step 1: Copy all CoSAGs and reads from round_5 ====================== #
for file in "$workdir/06_CoSAG_assembly/round_5"/*.fastq; do
  cp "$file" "$workdir/06_CoSAG_assembly/round_ending/round_all/round_5_$(basename "$file")"
done
# rm "$workdir/06_CoSAG_assembly/round_5"/*.fastq

for file in "$workdir/06_CoSAG_assembly/round_5/spades_output"/*.fasta; do
  cp "$file" "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/round_5_$(basename "$file")"
done

# ====================== Step 2: Collect overfitted CoSAGs from previous rounds ====================== #
# round_2_inconsistent overfitted clusters
while IFS= read -r line; do
  overfitted_id=$(echo "$line" | grep -oP 'Round 2 Cluster \K\d+')
  for file in "$workdir/06_CoSAG_assembly/round_2_inconsistent/${overfitted_id}"_R*.fastq; do
    cp "$file" "$workdir/06_CoSAG_assembly/round_ending/round_all/round_2_inconsistent_$(basename "$file")"
  done
  
  fasta_file="$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output/${overfitted_id}.fasta"
  if [ -f "$fasta_file" ]; then
    cp "$fasta_file" "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/round_2_inconsistent_$(basename "$fasta_file")"
  fi
  echo "$overfitted_id"
done < "$workdir/06_CoSAG_assembly/round_3/overfitted.log"

# rm "$workdir/06_CoSAG_assembly/round_2_inconsistent"/*.fastq

# round_3 overfitted clusters
while IFS= read -r line; do
  overfitted_id=$(echo "$line" | grep -oP 'Round 3 Cluster \K\d+')
  for file in "$workdir/06_CoSAG_assembly/round_3/${overfitted_id}"_R*.fastq; do
    cp "$file" "$workdir/06_CoSAG_assembly/round_ending/round_all/round_3_$(basename "$file")"
  done
  fasta_file="$workdir/06_CoSAG_assembly/round_3/spades_output/${overfitted_id}.fasta"
  if [ -f "$fasta_file" ]; then
    cp "$fasta_file" "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/round_3_$(basename "$fasta_file")"
  fi
  echo "$overfitted_id"
done < "$workdir/06_CoSAG_assembly/round_4/overfitted.log"

# rm "$workdir/06_CoSAG_assembly/round_3"/*.fastq

# round_4 overfitted clusters
while IFS= read -r line; do
  overfitted_id=$(echo "$line" | grep -oP 'Round 4 Cluster \K\d+')
  for file in "$workdir/06_CoSAG_assembly/round_4/${overfitted_id}"_R*.fastq; do
    cp "$file" "$workdir/06_CoSAG_assembly/round_ending/round_all/round_4_$(basename "$file")"
  done
  fasta_file="$workdir/06_CoSAG_assembly/round_4/spades_output/${overfitted_id}.fasta"
  if [ -f "$fasta_file" ]; then
    cp "$fasta_file" "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/round_4_$(basename "$fasta_file")"
  fi
  echo "$overfitted_id"
done < "$workdir/06_CoSAG_assembly/round_5/overfitted.log"

# rm "$workdir/06_CoSAG_assembly/round_4"/*.fastq

# ====================== Step 3: Collect final JSON for each round ====================== #
jq 'with_entries(.key |= "round_5_" + .)' "$workdir/06_CoSAG_assembly/round_5/cluster_members_round_5.final.json"  > "$workdir/06_CoSAG_assembly/round_ending/round_all/round_5.json"

mapfile -t overfitted_ids < <(grep -oP 'Round 2 Cluster \K\d+' "$workdir/06_CoSAG_assembly/round_3/overfitted.log")
ids_as_json=$(printf '%s\n' "${overfitted_ids[@]}" | jq -R . | jq -s .)
jq --argjson ids "$ids_as_json" '
  with_entries(select(.key as $k | $ids | map(tostring) | index($k)) | .key |= "round_2_inconsistent_" + .)
' "$workdir/06_CoSAG_assembly/round_2_inconsistent/cluster_members_round_2_inconsistent.json" > "$workdir/06_CoSAG_assembly/round_ending/round_all/round_2_inconsistent.json"

mapfile -t overfitted_ids < <(grep -oP 'Round 3 Cluster \K\d+' "$workdir/06_CoSAG_assembly/round_4/overfitted.log")
ids_as_json=$(printf '%s\n' "${overfitted_ids[@]}" | jq -R . | jq -s .)
jq --argjson ids "$ids_as_json" '
  with_entries(select(.key as $k | $ids | map(tostring) | index($k)) | .key |= "round_3_" + .)
' "$workdir/06_CoSAG_assembly/round_3/cluster_members_round_3.final.json" > "$workdir/06_CoSAG_assembly/round_ending/round_all/round_3.json"

mapfile -t overfitted_ids < <(grep -oP 'Round 4 Cluster \K\d+' "$workdir/06_CoSAG_assembly/round_5/overfitted.log")
ids_as_json=$(printf '%s\n' "${overfitted_ids[@]}" | jq -R . | jq -s .)
jq --argjson ids "$ids_as_json" '
  with_entries(select(.key as $k | $ids | map(tostring) | index($k)) | .key |= "round_4_" + .)
' "$workdir/06_CoSAG_assembly/round_4/cluster_members_round_4.final.json" > "$workdir/06_CoSAG_assembly/round_ending/round_all/round_4.json"

# ====================== Step 4: Merge all final JSONs ====================== #
jq -s '.[0] * .[1] * .[2] * .[3]' \
  "$workdir/06_CoSAG_assembly/round_ending/round_all/round_5.json" \
  "$workdir/06_CoSAG_assembly/round_ending/round_all/round_4.json" \
  "$workdir/06_CoSAG_assembly/round_ending/round_all/round_3.json" \
  "$workdir/06_CoSAG_assembly/round_ending/round_all/round_2_inconsistent.json" \
  > "$workdir/06_CoSAG_assembly/round_ending/round_all/round_all.json"

echo "[DONE] clustering_to_splitting.sh finished successfully."

# ====================== Step 5: Filter MGE contigs ====================== #
mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output_filter"

parallel -j 20 '
    FILE_TEMP={};
    BASENAME=$(basename $FILE_TEMP .fasta);

    blastn -query "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/${BASENAME}.fasta" \
           -db   "$workdir/06_CoSAG_assembly/MGE" \
           -outfmt 6 \
           -out "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/${BASENAME}_results.out"

    awk -f "$workdir/06_CoSAG_assembly/filter_high_identity.awk" \
        "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/${BASENAME}_results.out" \
      | sort | uniq \
      > "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/${BASENAME}_results.filter.id"

    filterbyname.sh \
        in="$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/${BASENAME}.fasta" \
        out="$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output_filter/${BASENAME}.fasta" \
        names="$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output/${BASENAME}_results.filter.id" \
        include=f ow=t substring=t minlen=0
' ::: "$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output"/*.fasta