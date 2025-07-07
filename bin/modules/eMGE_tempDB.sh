#!/usr/bin/env bash
set -euo pipefail
ulimit -n 100000
################################################################################
########################        Scellmate: eMGE_tempDB MODULE       ###########
################################################################################

help_message () {
cat <<EOF

Usage: $(basename "$0") --workdir DIR [options]

Required
  --workdir DIR          Scellmate working directory (expects 01_trim_SAGs inside)

Optional
  -t,  --threads INT     total CPU threads for SPAdes / GNU parallel      [40]
  --min_contig INT       minimum contig length kept before genoMAD        [500]
  --genomad_db DIR       geNomad database directory                       [/mnt/md0/wangyanren/database/genomad/genomad_db]
  --skip-assemble        skip SPAdes assembly if FASTA files already exist
  -h,  --help            show this help and exit

EOF
}

###############################################################################
# defaults
###############################################################################
threads=40
min_contig=500
genomad_db="/mnt/md0/wangyanren/database/genomad/genomad_db"
skip_assemble=false
###############################################################################

# CLI parsing
OPTS=$(getopt -o ht: --long workdir:,threads:,min_contig:,genomad_db:,skip-assemble,help -- "$@")
[[ $? -ne 0 ]] && { help_message; exit 1; }
eval set -- "$OPTS"
while true; do
  case "$1" in
    --workdir)      workdir=$2; shift 2 ;;
    -t|--threads)   threads=$2; shift 2 ;;
    --min_contig)   min_contig=$2; shift 2 ;;
    --genomad_db)   genomad_db=$2; shift 2 ;;
    --skip-assemble) skip_assemble=true; shift ;;
    -h|--help)      help_message; exit 0 ;;
    --) shift; break ;;
    *)  break ;;
  esac
done
[[ -z $workdir ]] && { echo "[ERROR] --workdir required"; exit 1; }

# Define log file path
log_file="$workdir/03_reference_db/eMGE_tempDB.log"

# Function to log step completion
log_step() {
    local step_name="$1"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $step_name completed" >> "$log_file"
}

# Function to check if step is done
step_done() {
    local step_name="$1"
    if [ -f "$log_file" ] && grep -q "$step_name completed" "$log_file"; then
        return 0  # Step is done
    else
        return 1  # Step is not done
    fi
}

# Create log directory
mkdir -p "$(dirname "$log_file")"

log="$workdir/04_SAG_assembly/eMGE.log"
mkdir -p "$workdir/04_SAG_assembly/spades_output"
: > "$log"


###############################################################################
# 1)  Assemble each SAG with SPAdes
if step_done "step 1: SPAdes assembly"; then
  echo "Skipping step 1: SPAdes assembly"
else
  echo "Start step 1: SPAdes assembly"

  set +e
  num_trim=$(ls "$workdir"/01_trim_SAGs/*_R1_paired.fastq 2>/dev/null | wc -l)
  num_fasta=$(ls "$workdir"/04_SAG_assembly/spades_output/*.fasta 2>/dev/null | wc -l)
  if [[ $num_trim -gt 0 && $num_trim -eq $num_fasta ]]; then
    skip_assemble=true
  fi
  set -e

  # ───────────────────── run or skip ──────────────────
  if ! $skip_assemble; then
    echo "[INFO] Launching spades_traversal.py for $num_trim samples …"

    python "bin/scripts/spades_traversal.py" \
          -i "$workdir/01_trim_SAGs" \
          -o "$workdir/04_SAG_assembly/spades_output" \
          --threads "$threads" \
          --per-job 1 \
          --suffix_R1 _R1_paired.fastq

    echo "[INFO] Assembly finished; logs are under spades_output/"
  else
    echo "[INFO] Assembly step skipped (files already present)."
  fi
  log_step "step 1: SPAdes assembly"
fi

###############################################################################
# 2) Merge all SAG assemblies with unique contig IDs
if step_done "step 2: merge SAG assemblies"; then
  echo "Skipping step 2: merge SAG assemblies"
else
  echo "Start step 2: merge SAG assemblies"
  combine="$workdir/04_SAG_assembly/SAG_assembly_combine.fasta"
  > "$combine"
  for fasta in "$workdir"/04_SAG_assembly/spades_output/*.fasta; do
    base=$(basename "$fasta" .fasta)
    awk -v pre="$base" '{
          if(/^>/){sub(/^>/,""); print ">"pre"_"$0}
          else print
        }' "$fasta" >> "$combine"
  done
  log_step "step 2: merge SAG assemblies"
fi

###############################################################################
# 3) Filter contigs by length
if step_done "step 3: filter contigs by length"; then
  echo "Skipping step 3: filter contigs by length"
else
  echo "Start step 3: filter contigs by length"
  filtered="$workdir/04_SAG_assembly/SAG_assembly_combine.${min_contig}.fasta"
  filterbyname.sh in="$combine" out="$filtered" ow=t minlen="$min_contig" &>> "$log"
  log_step "step 3: filter contigs by length"
fi

###############################################################################
# 4) Run geNomad (skip if already finished once)
if step_done "step 4: run geNomad"; then
  echo "Skipping step 4: run geNomad"
else
  echo "Start step 4: run geNomad"

  if ls "$workdir/04_SAG_assembly/genomad_output"/*_summary.log &>/dev/null; then
    echo "[INFO] geNomad output already present → skip" | tee -a "$log"
  else
    echo "[INFO] Running geNomad …" | tee -a "$log"
    conda run -n scellmate_genomad --no-capture-output genomad end-to-end --cleanup --splits 8 \
          "$filtered" \
          "$workdir/04_SAG_assembly/genomad_output" \
          "$genomad_db" &>> "$log"
  fi
  log_step "step 4: run geNomad"
fi


###############################################################################
# 5) Make blast database
if step_done "step 5: make blast database"; then
  echo "Skipping step 5: make blast database"
else
  mkdir -p "$workdir/06_CoSAG_assembly" 2>/dev/null || true
  cat "$workdir/04_SAG_assembly/genomad_output/SAG_assembly_combine.500_summary/SAG_assembly_combine.500_plasmid.fna" \
      "$workdir/04_SAG_assembly/genomad_output/SAG_assembly_combine.500_summary/SAG_assembly_combine.500_virus.fna" > "$workdir/06_CoSAG_assembly/MGE.fna"
  makeblastdb -in "$workdir/06_CoSAG_assembly/MGE.fna" -dbtype nucl -out "$workdir/06_CoSAG_assembly/MGE"
  log_step "step 5: make blast database"
fi

echo "eMGE calling pipeline finished" | tee -a "$log"
