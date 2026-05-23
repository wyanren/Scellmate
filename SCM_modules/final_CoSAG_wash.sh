#!/bin/bash
set -euo pipefail

shopt -s nullglob

# Default value for threads
threads=10

# fastANI merge switch (default: on)
do_fastani=1

# maximum CheckM contamination allowed for species-level curated CoSAG annotation
max_CoSAG_cont=10

# Parameter parsing
while [[ $# -gt 0 ]]; do
  case $1 in
    -w|--workdir)
      workdir=$(realpath "$2")
      shift 2
      ;;
    -s|--script)
      script=$(realpath "$2")
      shift 2
      ;;
    -t|--threads)
      threads="$2"
      shift 2
      ;;
    --max_CoSAG_cont)
      max_CoSAG_cont="$2"
      shift 2
      ;;
  --no_fastani)
      do_fastani=0
      shift 1
      ;;
    -h|--help)
      echo "Usage: $0 -w <workdir> -s <script_dir> [-t <threads>] [--max_CoSAG_cont <value>] [--no_fastani]"
      exit 0
      ;;
    *)
      echo "Unknown parameter: $1"
      echo "Usage: $0 -w <workdir> -s <script_dir> [-t <threads>] [--max_CoSAG_cont <value>] [--no_fastani]"
      exit 1
      ;;
  esac
done

# Check required parameters
if [[ -z "$workdir" || -z "$script" ]]; then
  echo "Error: --workdir and --script are required"
  echo "Usage: $0 -w <workdir> -s <script_dir> [-t <threads>] [--max_CoSAG_cont <value>] [--no_fastani]"
  exit 1
fi

RUN_SPADES="$script/spades_traversal.py"
export workdir script threads RUN_SPADES max_CoSAG_cont

safe_copy_splitting_results() {
  local split_dir="$1"
  local base_dir="$2"
  local target_dir="$3"

  if [ -d "$split_dir" ] && compgen -G "$split_dir/"'*_cells.txt' > /dev/null; then
    echo "[INFO] Copying *_cells.txt from $split_dir to $target_dir"
    find "$split_dir" -name "*_cells.txt" ! -name "*group*" -exec cp {} "$target_dir" \;

    find "$split_dir" -name "*_cells.txt" ! -name "*group*" | while read file; do
      base=$(basename "$file" "_cells.txt")
      mapping_src="$base_dir/$base/self_mapping_df.csv"
      mapping_dst="$target_dir/${base}-self_mapping_df.csv"
      if [ -f "$mapping_src" ]; then
        cp "$mapping_src" "$mapping_dst"
      else
        echo "[INFO] Missing: $mapping_src"
      fi
    done
  else
    echo "[INFO] Skipped $split_dir — not exist or no *_cells.txt"
  fi
}


# Log file for step tracking
logfile="$workdir/06_CoSAG_assembly/round_ending/round_all/Record-final_wash_CoSAG.log"

log_step() {
  echo "$1" >> "$logfile"
}

step_done() {
  grep -Fxq "$1" "$logfile" 2>/dev/null
}

# ========== Step 1: collect the splitting results during overfitting detection ==========
if step_done "step 1: collect splitting results during overfitting"; then
  echo "Skipping collect splitting results during overfitting"
else
  echo "Start collect splitting results during overfitting"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all"
  mkdir -p combined_full/1st 2>/dev/null || true
  cd combined_full/1st
  safe_copy_splitting_results "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align/split_output_filter/" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align/" ./
  safe_copy_splitting_results "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output_filter/" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/" ./
  safe_copy_splitting_results "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/split_output_filter/" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/" ./
  safe_copy_splitting_results "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output_filter/" "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/" ./
  
  combined_dir="$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/1st"
  spades_dir="$workdir/06_CoSAG_assembly/round_ending/round_all/spades_output"

  if ! ls "$combined_dir"/*.json 1>/dev/null 2>&1; then
    echo "[INFO] No JSON found in $combined_dir, generating from *_cells.txt"

    for file in "$combined_dir"/*_cells.txt; do
      [ -f "$file" ] || continue
      base=$(basename "$file" "_cells.txt")

      {
        echo "{"
        echo "  \"$base\": ["
        awk '{print "    \"" $1 "\","}' "$file" | sed '$s/,$//'
        echo "  ]"
        echo "}"
      } > "$combined_dir/${base}.json"

      fasta_path="$spades_dir/${base}.fasta"
      if [ -f "$fasta_path" ]; then
        cp "$fasta_path" "$combined_dir/"
      else
        echo "[WARN] FASTA not found for $base at $fasta_path"
      fi
    done
  fi

  log_step "step 1: collect splitting results during overfitting"
fi

# ========== Step 2: collect the splitting results ==========
if step_done "step 2: collect splitting results"; then
  echo "Skipping collect splitting results"
else
  echo "Start collect splitting results"

  mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/2nd" 2>/dev/null || true
  safe_copy_splitting_results \
    "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/split_output_filter" \
    "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result" \
    "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/2nd"

  mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/3rd" 2>/dev/null || true
  safe_copy_splitting_results \
    "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/split_output_filter" \
    "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align" \
    "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/3rd"

  mkdir -p "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/3rd_final" 2>/dev/null || true
  safe_copy_splitting_results \
    "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/split_output_filter" \
    "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align" \
    "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/3rd_final"

  split_dir="$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/split_output_filter"
  if [ -d "$split_dir" ]; then
    find "$split_dir" -name "*_cells.txt" ! -name "*group*" | while read file; do
      base=$(basename "$file" "_cells.txt")
      mapping_src="$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/$base/self_mapping_df.csv"
      mapping_dst="$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/3rd_final/${base}-self_mapping_df.csv"
      if [ -f "$mapping_src" ]; then
        cp "$mapping_src" "$mapping_dst"
      else
        echo "[INFO] Missing mapping: $mapping_src"
      fi
    done
  else
    echo "[INFO] Skipped 3rd_align/split_output_filter — directory not found"
  fi

  log_step "step 2: collect splitting results"
fi

# ========== Step 3: generate the final json and fasta files ==========
if step_done "step 3: generate the final json and fasta files"; then
  echo "Skipping step 3: generate the final json and fasta files"
else
  echo "Start step 3: generate the final json and fasta files"

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  for dir in 1st 2nd 3rd 3rd_final; do
    if [ -d "$dir" ] && compgen -G "$dir/*_cells.txt" > /dev/null; then
      python "$script/generate_json_saturate.py" --directory "$dir"
    else
      echo "[INFO] Skipped $dir — no txt files to saturate"
    fi
  done

  jsons_exist=$(find 1st 2nd 3rd 3rd_final -maxdepth 1 -name '*.json' 2>/dev/null | wc -l)
  if [ "$jsons_exist" -gt 0 ]; then
    jq -s . 1st/*.json 2nd/*.json 3rd/*.json 3rd_final/*.json > final.json || true
    jq 'reduce .[] as $item ({}; . + $item)' final.json > final.reformat.json
  else
    echo "[WARN] Skipped final.json generation — no JSON files"
    touch final.json final.reformat.json
  fi

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  # Copy round_*.fasta for 1st
  if compgen -G "1st/*.json" > /dev/null; then
    jq -s . 1st/*.json > 1st.json
    base_dir="."
    output_dir="$base_dir/1st"
    grep 'round' 1st.json | while read line; do
      cluster_id=$(echo "$line" | grep -oP '(round_\d+_inconsistent_\d+|round_3_\d+|round_4_\d+|round_5_\d+)')
      round_dir=$(echo "$cluster_id" | grep -oP 'round_\d+_inconsistent|round_3|round_4|round_5')
      file_num=$(echo "$cluster_id" | grep -oP '\d+$')
      fasta_file="$workdir/06_CoSAG_assembly/$round_dir/spades_output/$file_num.fasta"
      if [ -f "$fasta_file" ]; then
        cp "$fasta_file" "$output_dir/${cluster_id}.fasta"
      else
        echo "[INFO] File not found: $fasta_file"
      fi
    done
  else
    echo "[INFO] Skipped 1st.json → no input JSON"
  fi

  # Copy 2nd
  if compgen -G "2nd/*.json" > /dev/null; then
    jq -s . 2nd/*.json > 2nd.json
    grep 'round' 2nd.json | while read line; do
      cluster_id=$(echo "$line" | grep -oP '(round_\w*)')
      fasta_path="$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align/spades_output/$cluster_id.fasta"
      if [ -f "$fasta_path" ]; then
        cp "$fasta_path" 2nd/
      else
        echo "[INFO] Missing: $fasta_path"
      fi
    done
  else
    echo "[INFO] Skipped 2nd.json — no JSON in 2nd/"
  fi

  # Copy 3rd
  if compgen -G "3rd/*.json" > /dev/null; then
    jq -s . 3rd/*.json > 3rd.json
    grep 'round' 3rd.json | while read line; do
      cluster_id=$(echo "$line" | grep -oP '(round_\w*)')
      fasta_path="$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/spades_output/$cluster_id.fasta"
      if [ -f "$fasta_path" ]; then
        cp "$fasta_path" 3rd/
      else
        echo "[INFO] Missing: $fasta_path"
      fi
    done
  else
    echo "[INFO] Skipped 3rd.json — no JSON in 3rd/"
  fi

  # Copy 3rd_final
  if compgen -G "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/spades_output/*.fasta" > /dev/null; then
    cp "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/spades_output/"*.fasta 3rd_final/
  else
    echo "[INFO] Skipped 3rd_final copy — no fasta files found"
  fi

  log_step "step 3: generate the final json and fasta files"
fi


# ========== Step 4: sorting CoSAG genome ==========
if step_done "step 4: sorting CoSAG genome"; then
  echo "Skipping step 4: sorting CoSAG genome"
else
  echo "Start step 4: sorting CoSAG genome"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  mkdir -p SAG_BINs 2>/dev/null
  for dir in 1st 2nd 3rd 3rd_final; do
    if compgen -G "${dir}/*.fasta" > /dev/null; then
      cp "${dir}"/*.fasta SAG_BINs/
    else
      echo "[INFO] No fasta files found in $dir, skipping"
    fi
  done
  mkdir -p SAG_BINs_500_masking 2>/dev/null || true; mkdir -p SAG_BINs_500 2>/dev/null || true
  parallel -j 16 '
      FILE_TEMP={};
      BASENAME=$(basename $FILE_TEMP .fasta);
      echo $BASENAME;
      cd ./SAG_BINs
      blastn -query ${BASENAME}.fasta -db "$workdir/06_CoSAG_assembly/MGE" -outfmt 6 -out ${BASENAME}_results.out
      awk -f "$workdir/06_CoSAG_assembly/filter_high_identity.awk" ${BASENAME}_results.out | sort | uniq > ${BASENAME}_results.filter.id
      filterbyname.sh substring=name ow=t include=f minlen=500 names=${BASENAME}_results.filter.id in=${BASENAME}.fasta out=../SAG_BINs_500_masking/${BASENAME}.fasta > /dev/null 2>&1
      filterbyname.sh substring=name ow=t minlen=500 in=${BASENAME}.fasta out=../SAG_BINs_500/${BASENAME}.fasta > /dev/null 2>&1
  ' ::: ./SAG_BINs/*.fasta

  mkdir -p SAG_BINs_1000 2>/dev/null || true
  parallel -j 16 '
      FILE_TEMP={};
      BASENAME=$(basename $FILE_TEMP .fasta);
      echo $BASENAME;
      cd ./SAG_BINs
      filterbyname.sh substring=name ow=t minlen=1000 in=${BASENAME}.fasta out=../SAG_BINs_1000/${BASENAME}.fasta > /dev/null 2>&1
  ' ::: ./SAG_BINs/*.fasta

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/SAG_BINs_1000"
  if [ ! -f checkm/qa.txt ]; then
    checkm tree -x fasta ./ -t "$threads" checkm > /dev/null 2>&1
    checkm lineage_set checkm checkm/marker_file > /dev/null 2>&1
    checkm analyze checkm/marker_file -x fasta ./ checkm/analyze -t "$threads" > /dev/null 2>&1
    checkm tree_qa checkm -o 2 --tab_table > checkm/tree_qa.txt > /dev/null 2>&1
    checkm qa checkm/marker_file checkm/analyze -o 2 -t 30 -q -f checkm/qa.txt --tab_table
  fi

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/"
  array=( "$SCELLMATE_GTDB_PATH"/GTDBtk_DB* )
  echo "[INFO] Using GTDB-Tk data path: ${array[0]}"
  export GTDBTK_DATA_PATH="${array[0]}"

  if [ ! -f SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv ]; then
    gtdbtk classify_wf --genome_dir SAG_BINs_1000/ --extension fasta --out_dir SAG_BINs_1000/gtdbtk_output --cpus "$threads" --pplacer_cpus "$threads" --skip_ani_screen > /dev/null 2>&1
  fi
  log_step "step 4: sorting CoSAG genome"
fi


# ========== Step 5: clustering CoSAG genome ==========
if step_done "step 5: clustering CoSAG genome"; then
  echo "Skipping step 5: clustering CoSAG genome"
else
  echo "Start step 5: clustering CoSAG genome"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/"
  find SAG_BINs_500_masking/ -maxdepth 1 -type f -exec realpath {} \; > rl_MGE_mask.txt
  cp rl_MGE_mask.txt ql_MGE_mask.txt
  fastANI --ql rl_MGE_mask.txt --rl ql_MGE_mask.txt -t "$threads" --fragLen 500 -o fastANI_MGE_mask.txt > /dev/null 2>&1

  # awk -F',' 'NR>1 {count[$2]++} END {for (id in count) if (count[id] > 1) print id, count[id]}' sequence_clusters.csv
  # 155, 280, 395, 527 are the cluster ids that we want to keep
  # awk -F',' 'NR==1 || $2=="155" || $2=="280" || $2=="395" || $2=="527"' sequence_clusters.csv


  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  python "$script/clustering_fastANI.py" -i ./fastANI_MGE_mask.txt -o ./sequence_clusters.csv

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  python "$script/fastANI_cluster.py" -c ./sequence_clusters.csv -j ./final.reformat.json -o ./merged_clusters.json -x ".fasta" -t clustered_sequences.txt
  python "$script/cluster.py" --dirt_base . --dirt_fastq "$workdir/01_trim_SAGs/" -j ./merged_clusters.json
  log_step "step 5: clustering CoSAG genome"
fi

if [[ "$do_fastani" -eq 1 ]] && compgen -G "cluster/*_R1.fastq" > /dev/null; then
  # ========== Step 6: re-co-assembly the clustered CoSAGs ==========
  if step_done "step 6: re-co-assembly the clustered CoSAGs"; then
    echo "Skipping step 6: re-co-assembly the clustered CoSAGs"
  else
    echo "Start step 6: re-co-assembly the clustered CoSAGs"
    cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
    mkdir -p cluster/spades_output 2>/dev/null || true
    
    $RUN_SPADES \
      -i "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster" \
      -o "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster/spades_output" \
      --threads "$threads" \
      --per-job 8 \
      --suffix_R1 _R1.fastq
    log_step "step 6: re-co-assembly the clustered CoSAGs"
  fi

  # ========== Step 7: sorting the clustered CoSAGs ==========
  if step_done "step 7: sorting the clustered CoSAGs"; then
    echo "Skipping step 7: sorting the clustered CoSAGs"
  else
    echo "Start step 7: sorting the clustered CoSAGs"
    cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster/spades_output"
    mkdir ./SAG_BINs_1000 2>/dev/null || true
    parallel -j 16 '
        FILE_TEMP={};
        BASENAME=$(basename $FILE_TEMP .fasta);
        echo $BASENAME;
        filterbyname.sh substring=name ow=t minlen=1000 in=${BASENAME}.fasta out=./SAG_BINs_1000/${BASENAME}.fasta > /dev/null 2>&1
    ' ::: ./*.fasta

    cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster/spades_output/SAG_BINs_1000"
    if [ ! -f checkm/qa.txt ]; then
      checkm tree -x fasta ./ -t "$threads" checkm > /dev/null 2>&1
      checkm lineage_set checkm checkm/marker_file > /dev/null 2>&1
      checkm analyze checkm/marker_file -x fasta ./ checkm/analyze -t "$threads" > /dev/null 2>&1
      checkm tree_qa checkm -o 2 --tab_table > checkm/tree_qa.txt > /dev/null 2>&1
      checkm qa checkm/marker_file checkm/analyze -o 2 -t 30 -q -f checkm/qa.txt --tab_table
    fi

    cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster/spades_output/"
    
    if [ ! -f SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv ]; then
      gtdbtk classify_wf --genome_dir ./SAG_BINs_1000/ -x fasta --out_dir SAG_BINs_1000/gtdbtk_output --cpus "$threads" --pplacer_cpus "$threads" --skip_ani_screen > /dev/null 2>&1
    fi
    log_step "step 7: sorting the clustered CoSAGs"
  fi

  # ========== Step 8: add the cluster information to the final json ==========
  if step_done "step 8: add the cluster information to the final json"; then
    echo "Skipping step 8: add the cluster information to the final json"
  else
    echo "Start step 8: add the cluster information to the final json"
    cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/"
    if [[ -s clustered_sequences.txt ]] && compgen -G "cluster/spades_output/SAG_BINs_1000/*.fasta" > /dev/null; then
  keys=$(awk -F': ' '{print $2}' clustered_sequences.txt | sed 's/, /\n/g' | sed '/^$/d' | sed 's/^[ \t]*//;s/[ \t]*$//' | sed 's/"/\\"/g' | jq -R -s 'split("\n") | map(select(length > 0))' | jq -c '.')
else
  keys='[]'
fi

    jq --argjson keys "$keys" 'del( .[ $keys[] ] )' final.reformat.json > temp.json
    jq -s '.[0] * .[1]' temp.json merged_clusters.json > final.reformat.add_cluster.json
    rm temp.json

    log_step "step 8: add the cluster information to the final json"
  fi
else
  echo "[INFO] No _R1.fastq found in cluster/, skipping Step 6–8 for clustered CoSAGs."
  # log_step "step 6: re-co-assembly the clustered CoSAGs (skipped: no input)"
  # log_step "step 7: sorting the clustered CoSAGs (skipped: no input)"
  cp "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/final.reformat.json" "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/final.reformat.add_cluster.json"
  # log_step "step 8: add the cluster information to the final json (skipped: no input)"
fi


# ========== Step 9: kill the SAG that exist in other SAG_bin --- Final decontamination ==========
# in this step, we will kill the SAG that show inconsistency in the genus level of SAG and speceis level of CoSAG
if step_done "step 9: kill the SAG that exist in other SAG_bin --- Final decontamination"; then
  echo "Skipping step 9: kill the SAG that exist in other SAG_bin --- Final decontamination"
else
  echo "Start step 9: kill the SAG that exist in other SAG_bin --- Final decontamination"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/"
  cat SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv > SAG_BINs_1000/gtdbtk_output/gtdbtk.summary.tsv
  [ -s SAG_BINs_1000/gtdbtk_output/gtdbtk.ar53.summary.tsv ] && tail -n +2 SAG_BINs_1000/gtdbtk_output/gtdbtk.ar53.summary.tsv >> SAG_BINs_1000/gtdbtk_output/gtdbtk.summary.tsv
  cp SAG_BINs_1000/gtdbtk_output/gtdbtk.summary.tsv ./gtdbtk_1.summary.tsv
  if [[ -s clustered_sequences.txt ]] && compgen -G "cluster/spades_output/SAG_BINs_1000/*.fasta" > /dev/null; then
  keys=$(awk -F': ' '{print $2}' clustered_sequences.txt | sed 's/, /\n/g' | sed '/^$/d' | sed 's/^[ \t]*//;s/[ \t]*$//' | sed 's/"/\\"/g' | jq -R -s 'split("\n") | map(select(length > 0))' | jq -c '.')
else
  keys='[]'
fi
  echo "$keys" | jq -r '.[]' > temp
  cp gtdbtk_1.summary.tsv gtdbtk_1.summary.tsv.bak
  if [ -s temp ]; then
    awk 'NR==FNR {exclude[$1]; next} !($1 in exclude)' temp gtdbtk_1.summary.tsv.bak > gtdbtk_1.summary.tsv
  else
    cp gtdbtk_1.summary.tsv.bak gtdbtk_1.summary.tsv
  fi
  rm temp gtdbtk_1.summary.tsv.bak
  [ -s cluster/spades_output/SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv ] && tail -n +2 cluster/spades_output/SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv >> gtdbtk_1.summary.tsv
  jq -r 'to_entries | .[] | .key as $k | .value[] | "\($k) \(.)"' final.reformat.add_cluster.json > final.reformat.add_cluster.txt
  awk 'NR==FNR {map[$2] = $1; next} {print $0, (map[$2] ? map[$2] : "NA")}' "$workdir/05_first_QC/QC_non_mock-genus.tsv" <(cat final.reformat.add_cluster.txt) > final.reformat.add_cluster.annotate.txt

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/"
  gtdb_file="gtdbtk_1.summary.tsv"
  output_file="contamination.txt"
  rm -f "$output_file"
  awk '
      BEGIN {
          FS = "[ \t]+"  # Field separator is one or more spaces or tabs
      }
      NR==FNR {
          if (NR == 1) next  # Skip header line
          bin_name = $1
          taxonomy = $2
          # Concatenate remaining fields in case taxonomy includes spaces
          for (i = 3; i <= NF; i++) {
              taxonomy = taxonomy " " $i
          }
          # Extract genus from taxonomy
          match(taxonomy, /g__([^;]*)/, arr)
          genus = arr[1]
          # Skip if genus is empty
          if (genus == "") next
          bin_genus[bin_name] = genus
          next
      }
      {
          # Now processing final_output.tsv
          bin_name = $1
          genus_final = $3
          # Proceed only if the bin has a genus from GTDB
          if (bin_genus[bin_name]) {
              genus = bin_genus[bin_name]
              if (genus_final != "NA" && genus_final != "N/A" && genus_final != "" && genus_final != genus) {
                  print $0 >> "'"$output_file"'"
              }
          }
      }
  ' "$gtdb_file" final.reformat.add_cluster.annotate.txt

  if [ -s "$output_file" ]; then
  # construct the contamination_map.json
  awk '{key_values[$1]=key_values[$1]" "$2} END {for (key in key_values) print key, key_values[key]}' contamination.txt | \
    jq -Rn '
      (inputs | split(" ")) as $line |
      ($line[0]) as $key |
      ($line[1:] | map(select(length > 0))) as $values |
      {($key): $values}
    ' | jq -s 'add' > contamination_map.json

  cp "$script/edit_cont.py" ./edit_cont.py
  python edit_cont.py
  # this step would generate magic edit_cont .json file --- final.reformat.add_cluster.edit_cont.json

  python "$script/cluster_cont.py" --dirt_base . --dirt_fastq "$workdir/01_trim_SAGs/" -j ./edit_cont.json
  else
    echo "[INFO] contamination.txt not found, skipping edit and re-clustering..."
  fi
  log_step "step 9: kill the SAG that exist in other SAG_bin --- Final decontamination"
fi

# ========== Step 10: re-co-assembly the clean CoSAGs ==========
if [ -s "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/contamination.txt" ]; then
  if step_done "step 10: re-co-assembly the clean CoSAGs"; then
    echo "Skipping step 10: re-co-assembly the clean CoSAGs"
  else
    echo "Start step 10: re-co-assembly the clean CoSAGs"
    cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
    mkdir -p cluster_cont/spades_output 2>/dev/null || true

    if compgen -G "cluster_cont/*_R1.fastq" > /dev/null; then
      $RUN_SPADES \
        -i "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster_cont" \
        -o "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster_cont/spades_output" \
        --threads "$threads" \
        --per-job 8 \
        --suffix_R1 _R1.fastq
    else
      echo "[INFO] No _R1.fastq found in cluster_cont/, skip SPAdes"
    fi

    log_step "step 10: re-co-assembly the clean CoSAGs"
  fi
fi

if [ -s $workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/contamination.txt ]; then
  # ========== Step 11: sorting the clean CoSAGs ==========
  if step_done "step 11: sorting the clean CoSAGs"; then
    echo "Skipping step 11: sorting the clean CoSAGs"
  else
    echo "Start step 11: sorting the clean CoSAGs"
    cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster_cont/spades_output"
    mkdir ./SAG_BINs_1000 2>/dev/null || true
    parallel -j 16 '
        FILE_TEMP={};
        BASENAME=$(basename $FILE_TEMP .fasta);
        echo $BASENAME;
        cd ./
        filterbyname.sh substring=name ow=t minlen=1000 in=${BASENAME}.fasta out=./SAG_BINs_1000/${BASENAME}.fasta
    ' ::: ./*.fasta

    cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster_cont/spades_output/SAG_BINs_1000"
    if [ ! -f checkm/qa.txt ]; then
      checkm tree -x fasta ./ -t "$threads" checkm > /dev/null 2>&1
      checkm lineage_set checkm checkm/marker_file > /dev/null 2>&1
      checkm analyze checkm/marker_file -x fasta ./ checkm/analyze -t "$threads" > /dev/null 2>&1
      checkm tree_qa checkm -o 2 --tab_table > checkm/tree_qa.txt > /dev/null 2>&1
      checkm qa checkm/marker_file checkm/analyze -o 2 -t 30 -q -f checkm/qa.txt --tab_table
    fi

    cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster_cont/spades_output/"
    
    if [ ! -f SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv ]; then
      gtdbtk classify_wf --genome_dir ./SAG_BINs_1000/ -x fasta --out_dir SAG_BINs_1000/gtdbtk_output --cpus "$threads" --pplacer_cpus "$threads" --skip_ani_screen > /dev/null 2>&1
    fi
    log_step "step 11: sorting the clean CoSAGs"
  fi
else
  cp "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/final.reformat.add_cluster.json" "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/final.reformat.add_cluster.edit_cont.json" 
  echo "[INFO] contamination.txt not found or empty, skipping re-coassembly"
fi


# ========== Step 12: sort the final json table ==========
if step_done "step 12: sort the final json table"; then
  echo "Skipping step 12: sort the final json table"
else
  echo "Start step 12: sort the final json table"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  jq -r 'to_entries[] | [.key, .value[]] | @tsv' final.reformat.add_cluster.edit_cont.json > temp.tsv
  sort -k1,1 -k2,2n temp.tsv > sorted_temp.tsv
  # jq -R 'split("\t") | {key: .[0], value: .[1]}' sorted_temp.tsv | jq -s 'group_by(.key) | map({key: .[0].key, value: [.[] | .value]}) | from_entries' > final_updated.json
  jq -R 'split("\t") as $fields | { ($fields[0]): ($fields[1:] | map(select(length > 0))) }' sorted_temp.tsv | jq -s 'add' > final_updated.json
  jq -r 'to_entries[] | .key as $k | .value[] | "\($k)\t\(.)"' final_updated.json > final_updated.tsv
  sed -i 's/ \+/\t/g' final_updated.tsv
  log_step "step 12: sort the final json table"
fi


# ========== Step 13: update the gtdb and checkm summary table ==========
cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
if [[ -s contamination.txt ]]; then
  awk 'NR==FNR {exclude[$1]; next} !($1 in exclude)' contamination.txt gtdbtk_1.summary.tsv > gtdbtk_2.summary.tsv
else
  cp gtdbtk_1.summary.tsv gtdbtk_2.summary.tsv
fi

if [[ -s cluster_cont/spades_output/SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv ]]; then
  tail -n +2 cluster_cont/spades_output/SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv >> gtdbtk_2.summary.tsv
fi

# update the qa.txt
cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"

cp SAG_BINs_1000/checkm/qa.txt ./qa.txt.bak

# add clustered CoSAG
if [[ -s clustered_sequences.txt ]] && compgen -G "cluster/spades_output/SAG_BINs_1000/*.fasta" > /dev/null; then
  keys=$(awk -F': ' '{print $2}' clustered_sequences.txt | sed 's/, /\n/g' | sed '/^$/d' | sed 's/^[ \t]*//;s/[ \t]*$//' | sed 's/"/\\"/g' | jq -R -s 'split("\n") | map(select(length > 0))' | jq -c '.')
else
  keys='[]'
fi
echo "$keys" | jq -r '.[]' > keys_list.txt

if [[ -s keys_list.txt ]]; then
  awk 'NR==FNR {keys[$1]; next} !($1 in keys)' keys_list.txt qa.txt.bak > qa_1.txt
else
  cp qa.txt.bak qa_1.txt
fi

if [[ -s cluster/spades_output/SAG_BINs_1000/checkm/qa.txt ]]; then
  tail -n +2 cluster/spades_output/SAG_BINs_1000/checkm/qa.txt >> qa_1.txt
fi

# add contaminated CoSAG
if [[ -s contamination.txt ]]; then
  awk 'NR==FNR {exclude[$1]; next} !($1 in exclude)' contamination.txt qa_1.txt > qa_2.txt
else
  cp qa_1.txt qa_2.txt
fi

if [[ -s cluster_cont/spades_output/SAG_BINs_1000/checkm/qa.txt ]]; then
  tail -n +2 cluster_cont/spades_output/SAG_BINs_1000/checkm/qa.txt >> qa_2.txt
fi


# ========== Step 14: annotate the final CoSAGs ==========
cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
cp "$script/curated_cosag_annotate.py" ./curated_cosag_annotate.py
python curated_cosag_annotate.py --max_CoSAG_cont "$max_CoSAG_cont"
# CoSAG with specific annotation and contamination lower than max_CoSAG_cont are retained as species-level CoSAGs;
# CoSAGs without species-level annotation are treated as putative novel species only when they pass the CheckM criteria.

# ========== Step 15: summary the final CoSAGs ==========
cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
annotated_bins="annotated_bins.tsv"
final_updated="final_updated.tsv"
output_file="annotated_bins_with_counts.tsv"
awk '{print $1}' "$final_updated" | sort | uniq -c > bin_counts_raw.txt
awk '{print $2 "\t" $1}' bin_counts_raw.txt > bin_counts.txt
head -n1 "$annotated_bins" > header.txt
tail -n +2 "$annotated_bins" | sort -k1,1 > annotated_bins_sorted.tsv
sort -k1,1 bin_counts.txt > bin_counts_sorted.txt
join -t $'\t' -1 1 -2 1 annotated_bins_sorted.tsv bin_counts_sorted.txt > annotated_bins_counts.tsv
echo -e "bin_name\ttaxonomy\tannotate\tlevel\tcount" > "$output_file"
cat annotated_bins_counts.tsv >> "$output_file"
rm bin_counts_raw.txt bin_counts.txt annotated_bins_sorted.tsv bin_counts_sorted.txt annotated_bins_counts.tsv header.txt
sed -i 's/\r//g' annotated_bins_with_counts.tsv
tail -n +2 annotated_bins_with_counts.tsv | sort -k5,5nr > annotated_bins_with_counts_sorted.tsv

# annotate the final CoSAGs
cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
jq -r 'to_entries | .[] | .key as $k | .value[] | "\($k) \(.)"' final.reformat.add_cluster.edit_cont.json > final.reformat.add_cluster.edit_cont.txt
awk 'NR==FNR {map[$2] = $1; next} {print $0, (map[$2] ? map[$2] : "NA")}' \
"$workdir/05_first_QC/QC_non_mock-genus.tsv" <(cat final.reformat.add_cluster.edit_cont.txt) > final.reformat.add_cluster.edit_cont.annotate.txt

# count the number of CoSAGs for each species
cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
awk -F'\t' '$4 == "Species" {gsub(" ", "_", $3); print $1, $3}' annotated_bins_with_counts_sorted.tsv > species.tsv
awk 'FNR==NR { species[$1]=$2; next } { print $0, (species[$1] ? species[$1] : "no_species") }' species.tsv final.reformat.add_cluster.edit_cont.annotate.txt > temp
awk '$4 != "no_species"' temp > final.reformat.add_cluster.edit_cont.annotate.species.txt

awk '{temp=$1; $1=$4; $4=temp; print}' final.reformat.add_cluster.edit_cont.annotate.species.txt > temp
awk '{print $1, $2}' temp | sort -k1,1 | jq -Rn '
  [inputs | select(length > 0) | split(" ") | {key: .[0], value: .[1]}] | 
  group_by(.key) | 
  map({ (.[0].key): map(.value) }) | 
  add
' > final.reformat.add_cluster.edit_cont.annotate.species.json
rm temp

# collect all the final CoSAGs
cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
mkdir post_final_decont_genome 2>/dev/null || true
cp SAG_BINs_1000/*.fasta post_final_decont_genome/
echo "[INFO] The number of CoSAGs without refinement: $(find post_final_decont_genome -maxdepth 1 -type f -name "*.fasta" | wc -l)"

# update the clustered CoSAGs with new version
if [[ -s clustered_sequences.txt ]] && compgen -G "cluster/spades_output/SAG_BINs_1000/*.fasta" > /dev/null; then
  if [[ -s clustered_sequences.txt ]] && compgen -G "cluster/spades_output/SAG_BINs_1000/*.fasta" > /dev/null; then
  keys=$(awk -F': ' '{print $2}' clustered_sequences.txt | sed 's/, /\n/g' | sed '/^$/d' | sed 's/^[ \t]*//;s/[ \t]*$//' | sed 's/"/\\"/g' | jq -R -s 'split("\n") | map(select(length > 0))' | jq -c '.')
else
  keys='[]'
fi
  echo "$keys" | jq -r '.[]' | xargs -I {} rm -v post_final_decont_genome/{}.fasta
  echo "The number of CoSAGs after clustering without adding it: $(find post_final_decont_genome -maxdepth 1 -type f -name "*.fasta" | wc -l)"
  cp cluster/spades_output/SAG_BINs_1000/*.fasta post_final_decont_genome/
  echo "The number of CoSAGs after clustering with adding it: $(find post_final_decont_genome -maxdepth 1 -type f -name "*.fasta" | wc -l)"
else
  echo "[INFO] no clustering CoSAG generated after splitting"
fi

if [[ -s contamination.txt ]]; then
  # update the contaminated CoSAGs with new version
  cut -d ' ' -f 1 contamination.txt | sort | uniq | xargs -I {} rm -v post_final_decont_genome/{}.fasta
  echo "The number of CoSAGs after decontaminated without adding it: $(find post_final_decont_genome -maxdepth 1 -type f -name "*.fasta" | wc -l)"
  cp cluster_cont/spades_output/SAG_BINs_1000/*.fasta post_final_decont_genome/
  echo "The number of CoSAGs after decontaminated with adding it: $(find post_final_decont_genome -maxdepth 1 -type f -name "*.fasta" | wc -l)"
else
  echo "[INFO] no contaminated CoSAG detected after splitting"
fi

# ========== Step 16: collect the well-annotated CoSAGs ==========
mkdir well_annotated_genome 2>/dev/null || true
cp post_final_decont_genome/*.fasta well_annotated_genome/
cut -d ' ' -f 1 final.reformat.add_cluster.edit_cont.annotate.species.txt | sort | uniq | sed 's/$/.fasta/' > keep_files.txt
find well_annotated_genome -type f -name "*.fasta" | grep -v -f keep_files.txt | xargs rm
rm keep_files.txt

# ========== Step 17: collect the representative CoSAGs ==========
cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"

awk '{print $1, $2}' final.reformat.add_cluster.edit_cont.annotate.species.txt | sort -k1,1 | jq -Rn '
  [inputs | select(length > 0) | split(" ") | {key: .[0], value: .[1]}] | 
  group_by(.key) | 
  map({ (.[0].key): map(.value) }) | 
  add
' > final.reformat.add_cluster.edit_cont.annotate.species.nocluster.json

awk '{print $1, $4}' final.reformat.add_cluster.edit_cont.annotate.species.txt| sort -k2,2 |awk '{print $1 "_" $2, $0}' | sort -u -k1,1 | awk '{$1=""; print substr($0,2)}'| sort -k2,2 > species_to_nocluster_deduplicated.tsv

python "$script/decide_rep.py" -s "species_to_nocluster_deduplicated.tsv" -q qa_2.txt -o qa_rep.txt
awk '{print $1}' qa_rep.txt | tail -n+2 > rep.id
mkdir -p rep_fasta 2>/dev/null || true
while read id; do
  cp well_annotated_genome/${id}.fasta ./rep_fasta/
done < rep.id

echo -e "\033[0;33m[Scellmate] Curation and annotation of CoSAGs after splitting in second QC has been completed.\033[0m"
