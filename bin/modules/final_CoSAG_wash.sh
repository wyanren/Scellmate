#!/bin/bash
set -euo pipefail

# Default value for threads
threads=10

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
    -h|--help)
      echo "Usage: $0 -w <workdir> -s <script_dir> [-t <threads>]"
      exit 0
      ;;
    *)
      echo "Unknown parameter: $1"
      echo "Usage: $0 -w <workdir> -s <script_dir> [-t <threads>]"
      exit 1
      ;;
  esac
done

# Check required parameters
if [[ -z "$workdir" || -z "$script" ]]; then
  echo "Error: --workdir and --script are required"
  echo "Usage: $0 -w <workdir> -s <script_dir> [-t <threads>]"
  exit 1
fi

RUN_SPADES="$script/spades_traversal.py"
export workdir script threads RUN_SPADES

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
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align/split_output_filter/" -name "*_cells.txt" ! -name "*group*" -exec cp {} ./ \;
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align/split_output_filter/" -name "*_cells.txt" ! -name "*group*" | while read file; do
      base=$(basename "$file" "_cells.txt")
      cp "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/overfit_results/overfit_3rd_align/$base/self_mapping_df.csv" "./${base}-self_mapping_df.csv"
  done

  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output_filter/" -name "*_cells.txt" ! -name "*group*" -exec cp {} ./ \;
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/split_output_filter" -name "*_cells.txt" ! -name "*group*" | while read file; do
      base=$(basename "$file" "_cells.txt")
      cp "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/overfit_results/overfit_2nd_align/$base/self_mapping_df.csv" "./${base}-self_mapping_df.csv"
  done

  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/split_output_filter/" -name "*_cells.txt" ! -name "*group*" -exec cp {} ./ \;
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/split_output_filter/" -name "*_cells.txt" ! -name "*group*" | while read file; do
      base=$(basename "$file" "_cells.txt")
      cp "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/overfit_results/overfit_1st_align/$base/self_mapping_df.csv" "./${base}-self_mapping_df.csv"
  done

  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output_filter/" -name "*_cells.txt" ! -name "*group*" -exec cp {} ./ \;
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/split_output_filter/" -name "*_cells.txt" ! -name "*group*" | while read file; do
      base=$(basename "$file" "_cells.txt")
      cp "$workdir/06_CoSAG_assembly/round_ending/round_all/1st_align/$base/self_mapping_df.csv" "./${base}-self_mapping_df.csv"
  done
  log_step "step 1: collect splitting results during overfitting"
fi

# ========== Step 2: collect the splitting results ==========
if step_done "step 2: collect splitting results"; then
  echo "Skipping collect splitting results"
else
  echo "Start collect splitting results"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  mkdir -p 2nd 2>/dev/null || true; cd 2nd
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/split_output_filter" -name "*_cells.txt" ! -name "*group*" -exec cp {} ./ \;
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/split_output_filter" -name "*_cells.txt" ! -name "*group*" | while read file; do
      base=$(basename "$file" "_cells.txt")
      cp "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/$base/self_mapping_df.csv" "./${base}-self_mapping_df.csv"
  done

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  mkdir -p 3rd 2>/dev/null || true; cd 3rd
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/split_output_filter" -name "*_cells.txt" ! -name "*group*" -exec cp {} ./ \;
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/split_output_filter" -name "*_cells.txt" ! -name "*group*" | while read file; do
      base=$(basename "$file" "_cells.txt")
      cp "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/$base/self_mapping_df.csv" "./${base}-self_mapping_df.csv"
  done

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  mkdir -p 3rd_final 2>/dev/null || true; cd 3rd_final
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/split_output_filter" -name "*_cells.txt" -name "*group*" -exec cp {} ./ \;
  find "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/split_output_filter" -name "*_cells.txt" ! -name "*group*" | while read file; do
      base=$(basename "$file" "_cells.txt")
      cp "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/$base/self_mapping_df.csv" "./${base}-self_mapping_df.csv"
  done
  cp "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/split_output_filter/"*.json .
  log_step "step 2: collect splitting results"
fi

# ========== Step 3: generate the final json and fasta files ==========
if step_done "step 3: generate the final json and fasta files"; then
  echo "Skipping step 3: generate the final json and fasta files"
else
  echo "Start step 3: generate the final json and fasta files"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  python "$script/generate_json_saturate.py" --directory 1st
  python "$script/generate_json_saturate.py" --directory 2nd
  python "$script/generate_json_saturate.py" --directory 3rd
  jq -s . 1st/*.json 2nd/*.json 3rd/*.json 3rd_final/*.json > final.json
  jq 'reduce .[] as $item ({}; . + $item)' final.json > final.reformat.json

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  jq -s . 1st/*json > 1st.json
  base_dir="."
  output_dir="$base_dir/1st"
  grep 'round' 1st.json | while read line; do
      # Extract the full cluster ID, e.g., round_2_inconsistent_1291
      cluster_id=$(echo $line | grep -oP '(round_\d+_inconsistent_\d+|round_3_\d+|round_4_\d+|round_5_\d+)')
      # Extract round directory and file number from the cluster ID
      round_dir=$(echo $cluster_id | grep -oP 'round_\d+_inconsistent|round_3|round_4|round_5')
      file_num=$(echo $cluster_id | grep -oP '\d+$')
      # Define the source fasta file location based on the extracted IDs
      fasta_file="$workdir/06_CoSAG_assembly/$round_dir/spades_output/$file_num.fasta"
      # Check if the fasta file exists
      if [ -f "$fasta_file" ]; then
          # Copy the fasta file to the output directory with a new name
          cp "$fasta_file" "$output_dir/${cluster_id}.fasta"
      else
          echo "File not found: $fasta_file"
      fi
  done
  echo "Copy operation completed."

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  jq -s . 2nd/*json > 2nd.json
  grep 'round' 2nd.json | while read line; do
    cluster_id=$(echo $line | grep -oP '(round_\w*)')
    echo $cluster_id
    cp "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align/spades_output/$cluster_id.fasta" 2nd/
  done

  jq -s . 3rd/*json > 3rd.json
  grep 'round' 3rd.json | while read line; do
    cluster_id=$(echo $line | grep -oP '(round_\w*)')
    echo $cluster_id
    cp "$workdir/06_CoSAG_assembly/round_ending/round_all/2nd_align_result/spades_output/$cluster_id.fasta" 3rd/
  done

  cp "$workdir/06_CoSAG_assembly/round_ending/round_all/3rd_align/spades_output/"*.fasta 3rd_final/
  log_step "step 3: generate the final json and fasta files"
fi

# ========== Step 4: sorting CoSAG genome ==========
if step_done "step 4: sorting CoSAG genome"; then
  echo "Skipping step 4: sorting CoSAG genome"
else
  echo "Start step 4: sorting CoSAG genome"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  mkdir -p SAG_BINs 2>/dev/null || true; cp 1st/*.fasta SAG_BINs; cp 2nd/*.fasta SAG_BINs; cp 3rd/*.fasta SAG_BINs; cp 3rd_final/*.fasta SAG_BINs
  mkdir -p SAG_BINs_500_masking 2>/dev/null || true; mkdir -p SAG_BINs_500 2>/dev/null || true
  parallel -j 16 '
      FILE_TEMP={};
      BASENAME=$(basename $FILE_TEMP .fasta);
      echo $BASENAME;
      cd ./SAG_BINs
      blastn -query ${BASENAME}.fasta -db "$workdir/06_CoSAG_assembly/MGE" -outfmt 6 -out ${BASENAME}_results.out
      awk -f "$workdir/06_CoSAG_assembly/filter_high_identity.awk" ${BASENAME}_results.out | sort | uniq > ${BASENAME}_results.filter.id
      filterbyname.sh substring=name ow=t include=f minlen=500 names=${BASENAME}_results.filter.id in=${BASENAME}.fasta out=../SAG_BINs_500_masking/${BASENAME}.fasta
      filterbyname.sh substring=name ow=t minlen=500 in=${BASENAME}.fasta out=../SAG_BINs_500/${BASENAME}.fasta
  ' ::: ./SAG_BINs/*.fasta

  mkdir -p SAG_BINs_1000 2>/dev/null || true
  parallel -j 16 '
      FILE_TEMP={};
      BASENAME=$(basename $FILE_TEMP .fasta);
      echo $BASENAME;
      cd ./SAG_BINs
      filterbyname.sh substring=name ow=t minlen=1000 in=${BASENAME}.fasta out=../SAG_BINs_1000/${BASENAME}.fasta
  ' ::: ./SAG_BINs/*.fasta

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/SAG_BINs_1000"
  if [ ! -f checkm/qa.txt ]; then
    checkm tree -x fasta ./ -t "$threads" checkm
    checkm lineage_set checkm checkm/marker_file
    checkm analyze checkm/marker_file -x fasta ./ checkm/analyze -t "$threads"
    checkm tree_qa checkm -o 2 --tab_table > checkm/tree_qa.txt
    checkm qa checkm/marker_file checkm/analyze -o 2 -t 30 -q -f checkm/qa.txt --tab_table
  fi

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/"
  export GTDBTK_DATA_PATH="/mnt/md0/CX/ARGO/argo-evaluation/db/gtdbtk/release220"
  if [ ! -f SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv ]; then
    gtdbtk classify_wf --genome_dir SAG_BINs_1000/ --extension fasta --out_dir SAG_BINs_1000/gtdbtk_output --cpus "$threads" --pplacer_cpus "$threads" --skip_ani_screen
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
  conda run -n scellmate_fastANI --no-capture-output fastANI --ql rl_MGE_mask.txt --rl ql_MGE_mask.txt -t "$threads" --fragLen 500 -o fastANI_MGE_mask.txt

  awk -F',' 'NR>1 {count[$2]++} END {for (id in count) if (count[id] > 1) print id, count[id]}' sequence_clusters.csv
  # 155, 280, 395, 527 are the cluster ids that we want to keep
  awk -F',' 'NR==1 || $2=="155" || $2=="280" || $2=="395" || $2=="527"' sequence_clusters.csv


  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  python "$script/clustering_fastANI.py" -i ./fastANI_MGE_mask.txt -o ./sequence_clusters.csv

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  python "$script/fastANI_cluster.py" -c ./sequence_clusters.csv -j ./final.reformat.json -o ./merged_clusters.json -x ".fasta" -t clustered_sequences.txt
  python "$script/cluster.py" --dirt_base . --dirt_fastq "$workdir/01_trim_SAGs/" -j ./merged_clusters.json
  log_step "step 5: clustering CoSAG genome"
fi

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
      filterbyname.sh substring=name ow=t minlen=1000 in=${BASENAME}.fasta out=./SAG_BINs_1000/${BASENAME}.fasta
  ' ::: ./*.fasta

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster/spades_output/SAG_BINs_1000"
  if [ ! -f checkm/qa.txt ]; then
    checkm tree -x fasta ./ -t "$threads" checkm
    checkm lineage_set checkm checkm/marker_file
    checkm analyze checkm/marker_file -x fasta ./ checkm/analyze -t "$threads"
    checkm tree_qa checkm -o 2 --tab_table > checkm/tree_qa.txt
    checkm qa checkm/marker_file checkm/analyze -o 2 -t 30 -q -f checkm/qa.txt --tab_table
  fi

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster/spades_output/"
  export GTDBTK_DATA_PATH="/mnt/md0/CX/ARGO/argo-evaluation/db/gtdbtk/release220"
  if [ ! -f SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv ]; then
    gtdbtk classify_wf --genome_dir ./SAG_BINs_1000/ -x fasta --out_dir SAG_BINs_1000/gtdbtk_output --cpus "$threads" --pplacer_cpus "$threads" --skip_ani_screen
  fi
  log_step "step 7: sorting the clustered CoSAGs"
fi

# ========== Step 8: add the cluster information to the final json ==========
if step_done "step 8: add the cluster information to the final json"; then
  echo "Skipping step 8: add the cluster information to the final json"
else
  echo "Start step 8: add the cluster information to the final json"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/"
  keys=$(awk -F': ' '{print $2}' clustered_sequences.txt | sed 's/, /\n/g' | sed '/^$/d' | sed 's/^[ \t]*//;s/[ \t]*$//' | sed 's/"/\\"/g' | jq -R -s 'split("\n") | map(select(length > 0))' | jq -c '.')

  jq --argjson keys "$keys" 'del( .[ $keys[] ] )' final.reformat.json > temp.json
  jq -s '.[0] * .[1]' temp.json merged_clusters.json > final.reformat.add_cluster.json
  rm temp.json

  log_step "step 8: add the cluster information to the final json"
fi

# ========== Step 9: kill the SAG that exist in other SAG_bin --- Final decontamination ==========
# in this step, we will kill the SAG that show inconsistency in the genus level of SAG and speceis level of CoSAG
if step_done "step 9: kill the SAG that exist in other SAG_bin --- Final decontamination"; then
  echo "Skipping step 9: kill the SAG that exist in other SAG_bin --- Final decontamination"
else
  echo "Start step 9: kill the SAG that exist in other SAG_bin --- Final decontamination"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/"
  cat SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv > SAG_BINs_1000/gtdbtk_output/gtdbtk.summary.tsv
  tail -n +2 SAG_BINs_1000/gtdbtk_output/gtdbtk.ar53.summary.tsv >> SAG_BINs_1000/gtdbtk_output/gtdbtk.summary.tsv
  cp SAG_BINs_1000/gtdbtk_output/gtdbtk.summary.tsv ./gtdbtk_1.summary.tsv
  keys=$(awk -F': ' '{print $2}' clustered_sequences.txt | sed 's/, /\n/g' | sed '/^$/d' | sed 's/^[ \t]*//;s/[ \t]*$//' | sed 's/"/\\"/g' | jq -R -s 'split("\n") | map(select(length > 0))' | jq -c '.')
  echo "$keys" | jq -r '.[]' > temp
  cp gtdbtk_1.summary.tsv gtdbtk_1.summary.tsv.bak
  awk 'NR==FNR {exclude[$1]; next} !($1 in exclude)' temp gtdbtk_1.summary.tsv.bak > gtdbtk_1.summary.tsv
  rm temp gtdbtk_1.summary.tsv.bak
  tail -n +2 cluster/spades_output/SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv >> gtdbtk_1.summary.tsv
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
              if (genus_final != "NA" && genus_final != genus) {
                  print $0 >> "'"$output_file"'"
              }
          }
      }
  ' "$gtdb_file" final.reformat.add_cluster.annotate.txt

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

  python "$script/cluster_cont.py" --dirt_base . --dirt_fastq "$workdir/01_trim_SAGs/" -j ./edit_cont.json
  log_step "step 9: kill the SAG that exist in other SAG_bin --- Final decontamination"
fi

# ========== Step 10: re-co-assembly the clean CoSAGs ==========
if step_done "step 10: re-co-assembly the clean CoSAGs"; then
  echo "Skipping step 10: re-co-assembly the clean CoSAGs"
else
  echo "Start step 10: re-co-assembly the clean CoSAGs"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  mkdir -p cluster_cont/spades_output 2>/dev/null || true
  $RUN_SPADES \
    -i "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster_cont" \
    -o "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster_cont/spades_output" \
    --threads "$threads" \
    --per-job 8 \
    --suffix_R1 _R1.fastq
  log_step "step 10: re-co-assembly the clean CoSAGs"
fi

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
    checkm tree -x fasta ./ -t "$threads" checkm
    checkm lineage_set checkm checkm/marker_file
    checkm analyze checkm/marker_file -x fasta ./ checkm/analyze -t "$threads"
    checkm tree_qa checkm -o 2 --tab_table > checkm/tree_qa.txt
    checkm qa checkm/marker_file checkm/analyze -o 2 -t 30 -q -f checkm/qa.txt --tab_table
  fi

  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/cluster_cont/spades_output/"
  export GTDBTK_DATA_PATH="/mnt/md0/CX/ARGO/argo-evaluation/db/gtdbtk/release220"
  if [ ! -f SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv ]; then
    gtdbtk classify_wf --genome_dir ./SAG_BINs_1000/ -x fasta --out_dir SAG_BINs_1000/gtdbtk_output --cpus "$threads" --pplacer_cpus "$threads" --skip_ani_screen
  fi
  log_step "step 11: sorting the clean CoSAGs"
fi

# ========== Step 12: sort the final json table ==========
if step_done "step 12: sort the final json table"; then
  echo "Skipping step 12: sort the final json table"
else
  echo "Start step 12: sort the final json table"
  cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
  jq -r 'to_entries[] | [.key, .value[]] | @tsv' final.reformat.add_cluster.edit_cont.json > temp.tsv
  sort -k1,1 -k2,2n temp.tsv > sorted_temp.tsv
  jq -R 'split("\t") | {key: .[0], value: .[1]}' sorted_temp.tsv | jq -s 'group_by(.key) | map({key: .[0].key, value: [.[] | .value]}) | from_entries' > final_updated.json
  jq -r 'to_entries[] | [.key, .value[]] | @tsv' final_updated.json > final_updated.tsv
  sed -i 's/ \+/\t/g' final_updated.tsv
  log_step "step 12: sort the final json table"
fi

# ========== Step 13: update the gtdb and checkm summary table ==========
cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
awk 'NR==FNR {exclude[$1]; next} !($1 in exclude)' contamination.txt gtdbtk_1.summary.tsv > gtdbtk_2.summary.tsv
wc -l gtdbtk_*.summary.tsv
tail -n +2 cluster_cont/spades_output/SAG_BINs_1000/gtdbtk_output/gtdbtk.bac120.summary.tsv >> gtdbtk_2.summary.tsv
wc -l gtdbtk_*.summary.tsv

# update the qa.txt
cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
cp SAG_BINs_1000/checkm/qa.txt ./qa.txt.bak
keys=$(awk -F': ' '{print $2}' clustered_sequences.txt | sed 's/, /\n/g' | sed '/^$/d' \
| sed 's/^[ \t]*//;s/[ \t]*$//' | sed 's/"/\\"/g' \
| jq -R -s 'split("\n") | map(select(length > 0))' | jq -c '.')
echo "$keys" | jq -r '.[]'
echo "$keys" | jq -r '.[]' | awk 'NR==FNR {keys[$1]; next} !($1 in keys)' - qa.txt.bak > qa_1.txt
tail -n +2 cluster/spades_output/SAG_BINs_1000/checkm/qa.txt >> qa_1.txt
awk 'NR==FNR {exclude[$1]; next} !($1 in exclude)' contamination.txt qa_1.txt > qa_2.txt
wc -l qa_2.txt
tail -n +2 cluster_cont/spades_output/SAG_BINs_1000/checkm/qa.txt >> qa_2.txt
wc -l qa_*.txt

# ========== Step 14: annotate the final CoSAGs ==========
cd "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full"
cp "$script/annotate2.py" ./annotate2.py
python annotate2.py
# CoSAG with specific annotation, and for those CoSAG without, if they are medium quality then count them as novel species;
# Also, if there are already a species, count the latter CoSAG as novel species too.
head annotated_bins.tsv

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
wc -l final.reformat.add_cluster.edit_cont.annotate.species.txt
# 4198

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
echo "The number of CoSAGs without refinement: $(ls post_final_decont_genome/*.fasta | wc -l)"

# update the clustered CoSAGs with new version
keys=$(awk -F': ' '{print $2}' clustered_sequences.txt | sed 's/, /\n/g' | sed '/^$/d' \
| sed 's/^[ \t]*//;s/[ \t]*$//' | sed 's/"/\\"/g' \
| jq -R -s 'split("\n") | map(select(length > 0))' | jq -c '.')
echo "$keys" | jq -r '.[]' | xargs -I {} rm -v post_final_decont_genome/{}.fasta
echo "The number of CoSAGs after clustering without adding it: $(ls post_final_decont_genome/*.fasta | wc -l)"
cp cluster/spades_output/SAG_BINs_1000/*.fasta post_final_decont_genome/
echo "The number of CoSAGs after clustering with adding it: $(ls post_final_decont_genome/*.fasta | wc -l)"

# update the contaminated CoSAGs with new version
cut -d ' ' -f 1 contamination.txt | sort | uniq | xargs -I {} rm -v post_final_decont_genome/{}.fasta
echo "The number of CoSAGs after decontaminated without adding it: $(ls post_final_decont_genome/*.fasta | wc -l)"
cp cluster_cont/spades_output/SAG_BINs_1000/*.fasta post_final_decont_genome/
echo "The number of CoSAGs after decontaminated with adding it: $(ls post_final_decont_genome/*.fasta | wc -l)"

# ========== Step 16: collect the well-annotated CoSAGs ==========
mkdir well_annotated_genome 2>/dev/null || true
cp post_final_decont_genome/*.fasta well_annotated_genome/
cut -d ' ' -f 1 final.reformat.add_cluster.edit_cont.annotate.species.txt | sort | uniq | sed 's/$/.fasta/' > keep_files.txt
find well_annotated_genome -type f -name "*.fasta" | grep -v -f keep_files.txt | xargs rm -v
rm keep_files.txt
ls well_annotated_genome/*.fasta | wc -l

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
  cp well_annotated_genome//${id}.fasta ./rep_fasta/
done < rep.id