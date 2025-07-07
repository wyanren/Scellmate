#!/usr/bin/env bash
set -euo pipefail
# currently different workstation acts differently, so we use the current environment still, and wait for final package to be released

##########################################################################################
########################        Scellmate: CoSAG Iterative Assembly & Clustering        #################
##########################################################################################

# ========== Default parameters ==========
threads=48
workdir=""
script=""
overwrite=false

# ========== Help message ==========
help_message() {
  cat <<EOF

Usage: $(basename "$0") --workdir DIR --script DIR [options]

Required
  --workdir DIR               Working directory
  --script DIR                Script directory (should contain parse_to_round*.py, etc.)

Optional
  -t, --threads N             Number of parallel jobs (default: $threads)
  --overwrite                 Force rerun all steps and overwrite log
  -h, --help                  Show this help message

EOF
}

# ========== Parse arguments ==========
OPTS=$(getopt -o t:h --long workdir:,script:,threads:,overwrite,help -- "$@")
if [[ $? -ne 0 ]]; then
  help_message
  exit 1
fi
eval set -- "$OPTS"
while true; do
  case "$1" in
    --workdir)  workdir=$2; shift 2 ;;
    --script)   script=$2; shift 2 ;;
    -t|--threads) threads=$2; shift 2 ;;
    --overwrite) overwrite=true; shift ;;
    -h|--help)  help_message; exit 0 ;;
    --) shift; break ;;
    *) break ;;
  esac
done

if [[ -z "$workdir" || -z "$script" ]]; then
  echo "[ERROR] --workdir and --script are required." >&2
  help_message
  exit 1
fi

# Ensure main working directory exists
mkdir -p "$workdir/06_CoSAG_assembly"

# ========== Stepwise Logging Setup ==========
logfile="$workdir/06_CoSAG_assembly/Record-clustering_CoSAG.log"
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


# ========== Constants ==========
MGE_DB="${workdir}/06_CoSAG_assembly/MGE"
FILTER_AWK="${workdir}/06_CoSAG_assembly/filter_high_identity.awk"
PARSE2="$script/parse_to_round2.py"
PARSE3="$script/parse_to_round3.py"
PARSE4="$script/parse_to_round4.py"
PARSE5="$script/parse_to_round5.py"
RUN_SPADES="$script/spades_traversal.py"
LOGFILE="${LOGFILE:-${workdir}/06_CoSAG_assembly/clustering.running.log}"

##########################################################################################
##########            DIR                      ##########################################
##########################################################################################
echo '$3 >= 99 && $4 >= 500 {print $1}' > "${workdir}/06_CoSAG_assembly/filter_high_identity.awk"

##########################################################################################
##########            round_1 → assembly               #########
##########################################################################################

if step_done "round_1"; then
  echo "Skipping round_1"
else
  mkdir -p "$workdir/06_CoSAG_assembly/round_1/" 2>/dev/null || true
  
  find "$workdir/01_trim_SAGs" -maxdepth 1 -type f -name '*_R1_paired.fastq' -printf '%f\n' | sed 's/_R1_paired\.fastq$//' |sort -u > "$workdir/all.id"
  awk -F'\t' 'NR>1 && $3!="fail-1st-QC" {print $1}' "$workdir/05_first_QC/QC_1st.tsv" > "$workdir/06_CoSAG_assembly/round_1/exclude_fail_1st_QC.id"
  comm -12 \
  <(sort -u "$workdir/all.id") \
  <(sort -u "$workdir/06_CoSAG_assembly/round_1/exclude_fail_1st_QC.id") \
  > "$workdir/06_CoSAG_assembly/round_1/intersect.id"

  mkdir -p "$workdir/06_CoSAG_assembly/round_1/spades_output_filtered_sample/" 2>/dev/null || true
  while read -r id; do
    src="$workdir/04_SAG_assembly/spades_output/${id}.fasta"
    if [ -f "$src" ]; then
        cp "$src" "$workdir/06_CoSAG_assembly/round_1/spades_output_filtered_sample/"
    else
        echo "⚠️  $src not found" >&2
    fi
  done < "$workdir/06_CoSAG_assembly/round_1/intersect.id"

  mkdir -p "$workdir/06_CoSAG_assembly/round_1/spades_output_filter/" 2>/dev/null || true

  workdir=$(realpath "$workdir")
  LOGFILE=$(realpath "$LOGFILE")
  MGE_DB=$(realpath "$MGE_DB")
  FILTER_AWK=$(realpath "$FILTER_AWK")
  export workdir LOGFILE MGE_DB FILTER_AWK
  parallel --env threads --env LOGFILE --env FILTER_AWK --env MGE_DB -j "$threads" --bar '
    infile={};
    sample=$(basename "$infile" .fasta);
     {
      echo "[round_1 MGE] $sample @ $(date)";
      blastn -query "$infile" -db "$MGE_DB" -outfmt 6 -out "$workdir/06_CoSAG_assembly/round_1/spades_output_filtered_sample/${sample}_results.out";
      awk -f "$FILTER_AWK" "$workdir/06_CoSAG_assembly/round_1/spades_output_filtered_sample/${sample}_results.out" | sort | uniq > "$workdir/06_CoSAG_assembly/round_1/spades_output_filtered_sample/${sample}_results.filter.id";
      filterbyname.sh substring=name ow=t include=f minlen=0 \
        names="$workdir/06_CoSAG_assembly/round_1/spades_output_filtered_sample/${sample}_results.filter.id" \
        in="$workdir/06_CoSAG_assembly/round_1/spades_output_filtered_sample/${sample}.fasta" \
        out="$workdir/06_CoSAG_assembly/round_1/spades_output_filter/${sample}.fasta"
     } &>> "$LOGFILE"
  ' ::: "$workdir/06_CoSAG_assembly/round_1/spades_output_filtered_sample"/*.fasta


  mkdir -p "$workdir/06_CoSAG_assembly/round_1/sourmash/" 2>/dev/null || true
  parallel --env threads --env LOGFILE --env workdir -j "$threads" --bar '
    sample=$(basename {} .fasta);
    {
      echo "[round_1 sourmash] $sample @ $(date)";
      sourmash compute --track-abundance "$workdir/06_CoSAG_assembly/round_1/spades_output_filter/${sample}.fasta" --output "$workdir/06_CoSAG_assembly/round_1/sourmash/${sample}.sig";
    } &>> "$LOGFILE"
  ' ::: "$workdir/06_CoSAG_assembly/round_1/spades_output_filter"/*.fasta

  echo "[round_1 compare] @ $(date)";
  sourmash compare -p "$threads" \
      "$workdir/06_CoSAG_assembly/round_1/sourmash"/*.sig -k 51 \
      -o "$workdir/06_CoSAG_assembly/round_1/1st_round_cmp.npy" \
      --csv "$workdir/06_CoSAG_assembly/round_1/1st_round_cmp.csv"
  log_step "round_1"
fi


##########################################################################################
##########            round_2_inconsistent → assembly               #########
##########################################################################################
if step_done "round_2_inconsistent"; then
  echo "Skipping round_2_inconsistent"
else
  mkdir -p "$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output" \
           "$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output_filter" \
           "$workdir/06_CoSAG_assembly/round_2_inconsistent/sourmash"
  
  # Check if cluster_members_round_2_inconsistent.json exists, if so skip PARSE2
  if [[ -f "$workdir/06_CoSAG_assembly/round_2_inconsistent/cluster_members_round_2_inconsistent.json" ]]; then
    echo "[INFO] cluster_members_round_2_inconsistent.json exists, skipping PARSE2 command"
  else
    python $PARSE2 --dirt-base "$workdir/06_CoSAG_assembly" --dirt-fastq "$workdir/01_trim_SAGs" --qc-file "$workdir/05_first_QC/QC_non_mock-genus.tsv"
  fi
  
  echo "Generating round_2 temporary CoSAG cluster"
  set +e
  num_input=$(ls "$workdir/06_CoSAG_assembly/round_2_inconsistent"/*_R1.fastq 2>/dev/null | wc -l)
  num_fasta=$(ls "$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output"/*.fasta 2>/dev/null | wc -l)
  set -e
  if [[ $num_input -gt 0 && $num_input -eq $num_fasta ]]; then
    echo "[INFO] round_2_inconsistent: All fasta files present, skipping SPAdes."
  else
    $RUN_SPADES \
      -i "$workdir/06_CoSAG_assembly/round_2_inconsistent" \
      -o "$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output" \
      --threads "$threads" \
      --per-job 2 \
      --suffix_R1 _R1.fastq
  fi


  echo "Generating round_2 temporary CoSAG assemblies"

  workdir=$(realpath "$workdir")
  LOGFILE=$(realpath "$LOGFILE")
  MGE_DB=$(realpath "$MGE_DB")
  FILTER_AWK=$(realpath "$FILTER_AWK")
  export workdir LOGFILE MGE_DB FILTER_AWK
  parallel --env threads --env LOGFILE --env FILTER_AWK --env MGE_DB -j "$threads" --bar '
    infile={};
    sample=$(basename "$infile" .fasta);
     {
      echo "[round_2_inconsistent MGE] $sample @ $(date)";
      blastn -query "$infile" -db "$MGE_DB" -outfmt 6 -out "$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output/${sample}_results.out";
      awk -f "$FILTER_AWK" "$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output/${sample}_results.out" | sort | uniq > "$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output/${sample}_results.filter.id";
      filterbyname.sh substring=name ow=t include=f minlen=0 \
        names="$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output/${sample}_results.filter.id" \
        in="$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output/${sample}.fasta" \
        out="$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output_filter/${sample}.fasta"
     } &>> "$LOGFILE"
  ' ::: "$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output"/*.fasta

  parallel --env threads --env LOGFILE --env workdir -j "$threads" --bar '
    sample=$(basename {} .fasta);
    {
      echo "[round_2_inconsistent sourmash] $sample @ $(date)";
      sourmash compute --track-abundance "$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output_filter/${sample}.fasta" --output "$workdir/06_CoSAG_assembly/round_2_inconsistent/sourmash/${sample}.sig";
    } &>> "$LOGFILE"
  ' ::: "$workdir/06_CoSAG_assembly/round_2_inconsistent/spades_output_filter"/*.fasta

  echo "[round_2_inconsistent compare] @ $(date)";
  sourmash compare -p "$threads" \
      "$workdir/06_CoSAG_assembly/round_2_inconsistent/sourmash"/*.sig -k 51 \
      -o "$workdir/06_CoSAG_assembly/round_2_inconsistent/2nd_round_cmp.npy" \
      --csv "$workdir/06_CoSAG_assembly/round_2_inconsistent/2nd_round_cmp.csv"
  log_step "round_2_inconsistent"
fi

##########################################################################################
##########            round_3 → assembly               #########
##########################################################################################

if step_done "round_3"; then
  echo "Skipping round_3"
else
  mkdir -p "$workdir/06_CoSAG_assembly/round_3/spades_output" \
           "$workdir/06_CoSAG_assembly/round_3/spades_output_filter" \
           "$workdir/06_CoSAG_assembly/round_3/sourmash"
  
  # Check if cluster_members_round_3.*.json exists, if so skip PARSE3
  if [[ -f "$workdir/06_CoSAG_assembly/round_3/cluster_members_round_3.final.json" ]]; then
    echo "[INFO] cluster_members_round_3.final.json exists, skipping PARSE3 command"
  else
    python $PARSE3 --dirt-base "$workdir/06_CoSAG_assembly" --dirt-fastq "$workdir/01_trim_SAGs" --qc-file "$workdir/05_first_QC/QC_non_mock-genus.tsv"
  fi
  
  echo "Generating round_3 temporary CoSAG cluster"
  set +e
  num_input3=$(ls "$workdir/06_CoSAG_assembly/round_3"/*_R1.fastq 2>/dev/null | wc -l)
  num_fasta3=$(ls "$workdir/06_CoSAG_assembly/round_3/spades_output"/*.fasta 2>/dev/null | wc -l)
  set -e
  if [[ $num_input3 -gt 0 && $num_input3 -eq $num_fasta3 ]]; then
    echo "[INFO] round_3: All fasta files present, skipping SPAdes."
  else
    $RUN_SPADES \
      -i "$workdir/06_CoSAG_assembly/round_3" \
      -o "$workdir/06_CoSAG_assembly/round_3/spades_output" \
      --threads "$threads" \
      --per-job 2 \
      --suffix_R1 _R1.fastq
  fi

  echo "Generating round_3 temporary CoSAG assemblies"

  workdir=$(realpath "$workdir")
  LOGFILE=$(realpath "$LOGFILE")
  MGE_DB=$(realpath "$MGE_DB")
  FILTER_AWK=$(realpath "$FILTER_AWK")
  export workdir LOGFILE MGE_DB FILTER_AWK
  parallel --env threads --env LOGFILE --env FILTER_AWK --env MGE_DB -j "$threads" --bar '
    infile={};
    sample=$(basename "$infile" .fasta);
    {
      echo "[round_3 MGE] $sample @ $(date)";
      blastn -query "$infile" -db "$MGE_DB" -outfmt 6 \
      -out "$workdir/06_CoSAG_assembly/round_3/spades_output/${sample}_results.out";
      awk -f "$FILTER_AWK" "$workdir/06_CoSAG_assembly/round_3/spades_output/${sample}_results.out" | sort | uniq > "$workdir/06_CoSAG_assembly/round_3/spades_output/${sample}_results.filter.id";
      filterbyname.sh substring=name ow=t include=f minlen=0 \
        names="$workdir/06_CoSAG_assembly/round_3/spades_output/${sample}_results.filter.id" \
        in="$workdir/06_CoSAG_assembly/round_3/spades_output/${sample}.fasta" \
        out="$workdir/06_CoSAG_assembly/round_3/spades_output_filter/${sample}.fasta";
    } &>> "$LOGFILE"
  ' ::: "$workdir/06_CoSAG_assembly/round_3/spades_output"/*.fasta
  parallel --env threads --env LOGFILE -j "$threads" --bar '
    sample=$(basename {} .fasta);
    {
      echo "[round_3 sourmash] $sample @ $(date)";
      sourmash compute --track-abundance "$workdir/06_CoSAG_assembly/round_3/spades_output_filter/${sample}.fasta" \
        --output "$workdir/06_CoSAG_assembly/round_3/sourmash/${sample}.sig";
    } &>> "$LOGFILE"
  ' ::: "$workdir/06_CoSAG_assembly/round_3/spades_output_filter"/*.fasta
  echo "[round_3 compare] @ $(date)";
  sourmash compare -p "$threads" \
      "$workdir/06_CoSAG_assembly/round_3/sourmash"/*.sig -k 51 \
      -o "$workdir/06_CoSAG_assembly/round_3/3rd_round_cmp.npy" \
      --csv "$workdir/06_CoSAG_assembly/round_3/3rd_round_cmp.csv"
  log_step "round_3"
fi


##########################################################################################
##########            round_4 → assembly               #########
##########################################################################################

if step_done "round_4"; then
  echo "Skipping round_4"
else
  mkdir -p "$workdir/06_CoSAG_assembly/round_4/spades_output" \
           "$workdir/06_CoSAG_assembly/round_4/spades_output_filter" \
           "$workdir/06_CoSAG_assembly/round_4/sourmash"
  
  # Check if cluster_members_round_4.*.json exists, if so skip PARSE4
  if [[ -f "$workdir/06_CoSAG_assembly/round_4/cluster_members_round_4.final.json" ]]; then
    echo "[INFO] cluster_members_round_4.final.json exists, skipping PARSE4 command"
  else
    python $PARSE4 --dirt-base "$workdir/06_CoSAG_assembly" --dirt-fastq "$workdir/01_trim_SAGs" --qc-file "$workdir/05_first_QC/QC_non_mock-genus.tsv"
  fi
  
  echo "Generating round_4 temporary CoSAG cluster"
  set +e
  num_input4=$(ls "$workdir/06_CoSAG_assembly/round_4"/*_R1.fastq 2>/dev/null | wc -l)
  num_fasta4=$(ls "$workdir/06_CoSAG_assembly/round_4/spades_output"/*.fasta 2>/dev/null | wc -l)
  set -e
  if [[ $num_input4 -gt 0 && $num_input4 -eq $num_fasta4 ]]; then
    echo "[INFO] round_4: All fasta files present, skipping SPAdes."
  else
    $RUN_SPADES \
      -i "$workdir/06_CoSAG_assembly/round_4" \
      -o "$workdir/06_CoSAG_assembly/round_4/spades_output" \
      --threads "$threads" \
      --per-job 2 \
      --suffix_R1 _R1.fastq
  fi

  echo "Generating round_4 temporary CoSAG assemblies"
  workdir=$(realpath "$workdir")
  LOGFILE=$(realpath "$LOGFILE")
  MGE_DB=$(realpath "$MGE_DB")
  FILTER_AWK=$(realpath "$FILTER_AWK")
  export workdir LOGFILE MGE_DB FILTER_AWK
  parallel --env threads --env LOGFILE --env FILTER_AWK --env MGE_DB -j "$threads" --bar '
    infile={};
    sample=$(basename "$infile" .fasta);
    {
      echo "[round_4 MGE] $sample @ $(date)";
      blastn -query "$infile" -db "$MGE_DB" -outfmt 6 \
      -out "$workdir/06_CoSAG_assembly/round_4/spades_output/${sample}_results.out";
      awk -f "$FILTER_AWK" "$workdir/06_CoSAG_assembly/round_4/spades_output/${sample}_results.out" | sort | uniq > "$workdir/06_CoSAG_assembly/round_4/spades_output/${sample}_results.filter.id";
      filterbyname.sh substring=name ow=t include=f minlen=0 \
        names="$workdir/06_CoSAG_assembly/round_4/spades_output/${sample}_results.filter.id" \
        in="$workdir/06_CoSAG_assembly/round_4/spades_output/${sample}.fasta" \
        out="$workdir/06_CoSAG_assembly/round_4/spades_output_filter/${sample}.fasta";
    } &>> "$LOGFILE"
  ' ::: "$workdir/06_CoSAG_assembly/round_4/spades_output"/*.fasta
  parallel --env threads --env LOGFILE -j "$threads" --bar '
    sample=$(basename {} .fasta);
    {
      echo "[round_4 sourmash] $sample @ $(date)";
      sourmash compute --track-abundance "$workdir/06_CoSAG_assembly/round_4/spades_output_filter/${sample}.fasta" \
        --output "$workdir/06_CoSAG_assembly/round_4/sourmash/${sample}.sig";
    } &>> "$LOGFILE"
  ' ::: "$workdir/06_CoSAG_assembly/round_4/spades_output_filter"/*.fasta
  echo "[Round4 compare] @ $(date)";
  sourmash compare -p "$threads" \
      "$workdir/06_CoSAG_assembly/round_4/sourmash"/*.sig -k 51 \
      -o "$workdir/06_CoSAG_assembly/round_4/4th_round_cmp.npy" \
      --csv "$workdir/06_CoSAG_assembly/round_4/4th_round_cmp.csv"
  log_step "round_4"
fi

##########################################################################################
##########            round_5 → assembly               #########
##########################################################################################

if step_done "round_5"; then
  echo "Skipping round_5"
else
  mkdir -p "$workdir/06_CoSAG_assembly/round_5/spades_output" \
           "$workdir/06_CoSAG_assembly/round_5/spades_output_filter" \
           "$workdir/06_CoSAG_assembly/round_5/sourmash"
  
  # Check if cluster_members_round_5.*.json exists, if so skip PARSE5
  if [[ -f "$workdir/06_CoSAG_assembly/round_5/cluster_members_round_5.final.json" ]]; then
    echo "[INFO] cluster_members_round_5.final.json exists, skipping PARSE5 command"
  else
    python $PARSE5 --dirt-base "$workdir/06_CoSAG_assembly" --dirt-fastq "$workdir/01_trim_SAGs" --qc-file "$workdir/05_first_QC/QC_non_mock-genus.tsv"
  fi
  
  echo "Generating round_5 temporary CoSAG cluster"
  set +e
  num_input5=$(ls "$workdir/06_CoSAG_assembly/round_5"/*_R1.fastq 2>/dev/null | wc -l)
  num_fasta5=$(ls "$workdir/06_CoSAG_assembly/round_5/spades_output"/*.fasta 2>/dev/null | wc -l)
  set -e
  if [[ $num_input5 -gt 0 && $num_input5 -eq $num_fasta5 ]]; then
    echo "[INFO] round_5: All fasta files present, skipping SPAdes."
  else
    $RUN_SPADES \
      -i "$workdir/06_CoSAG_assembly/round_5" \
      -o "$workdir/06_CoSAG_assembly/round_5/spades_output" \
      --threads "$threads" \
      --per-job 2 \
      --suffix_R1 _R1.fastq
  fi

  echo "Generating round_5 temporary CoSAG assemblies"
  parallel --env threads --env LOGFILE --env FILTER_AWK --env MGE_DB -j "$threads" --bar '
    infile={};
    sample=$(basename "$infile" .fasta);
    {
      echo "[round_5 MGE] $sample @ $(date)";
      blastn -query "$infile" -db "$MGE_DB" -outfmt 6 \
      -out "$workdir/06_CoSAG_assembly/round_5/spades_output/${sample}_results.out";
      awk -f "$FILTER_AWK" "$workdir/06_CoSAG_assembly/round_5/spades_output/${sample}_results.out" | sort | uniq > "$workdir/06_CoSAG_assembly/round_5/spades_output/${sample}_results.filter.id";
      filterbyname.sh substring=name ow=t include=f minlen=0 \
        names="$workdir/06_CoSAG_assembly/round_5/spades_output/${sample}_results.filter.id" \
        in="$workdir/06_CoSAG_assembly/round_5/spades_output/${sample}.fasta" \
        out="$workdir/06_CoSAG_assembly/round_5/spades_output_filter/${sample}.fasta";
    } &>> "$LOGFILE"
  ' ::: "$workdir/06_CoSAG_assembly/round_5/spades_output"/*.fasta
  parallel --env threads --env LOGFILE -j "$threads" --bar '
    sample=$(basename {} .fasta);
    {
      echo "[round_5 sourmash] $sample @ $(date)";
      sourmash compute --track-abundance "$workdir/06_CoSAG_assembly/round_5/spades_output_filter/${sample}.fasta" \
        --output "$workdir/06_CoSAG_assembly/round_5/sourmash/${sample}.sig";
    } &>> "$LOGFILE"
  ' ::: "$workdir/06_CoSAG_assembly/round_5/spades_output_filter"/*.fasta
  echo "[round_5 compare] @ $(date)";
  sourmash compare -p "$threads" \
      "$workdir/06_CoSAG_assembly/round_5/sourmash"/*.sig -k 51 \
      -o "$workdir/06_CoSAG_assembly/round_5/4th_round_cmp.npy" \
      --csv "$workdir/06_CoSAG_assembly/round_5/4th_round_cmp.csv"
  log_step "round_5"
fi

