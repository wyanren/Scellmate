#!/usr/bin/env bash
set -euo pipefail

##########################################################################################
########################        Scellmate: KRAKEN2 DB SUBSETTING MODULE        #######################
##########################################################################################

##########################################################################################
# default thresholds & general options
MIN_GENUS_READS=5000      # keep genus row in .report if clade_reads > xxx —— species-level
MIN_SPECIES_DIRECT=500    # keep species if direct_reads > xxx  —— species-level
RATIO_CUTOFF=0.50         # direct/clade ratio to define “unclassified” —— genus-level
MIN_GENUS_OCC=6           # ratio-qualified genus must appear in > this many SAGs —— genus-level
threads=30                # parallel jobs for Kraken2 / GNU parallel
##########################################################################################

help_message () {
cat <<EOF

Usage: $(basename "$0") --workdir DIR [options]

Required
  --workdir DIR              working directory (trimmed reads are expected in DIR/01_trim_SAGs)

Optional
  -t, --threads N            parallel jobs   (default: $threads)
      --min_genus_reads N    genus row cutoff (default: $MIN_GENUS_READS)
      --min_species_direct N species direct-read cutoff (default: $MIN_SPECIES_DIRECT)
      --ratio_cutoff FLOAT   direct/clade ratio (default: $RATIO_CUTOFF)
      --min_genus_occ N      min #SAGs a ratio-qualified genus must appear in (default: $MIN_GENUS_OCC)
  -h, --help                 show this help

EOF
}

##########################################################################################
########################                Load config                         ##########################
##########################################################################################
CONFIG_FILE=$(realpath "$(dirname "$0")/../../config/scellmate-config")
[[ -f $CONFIG_FILE ]] || { echo "[ERROR] config missing: $CONFIG_FILE"; exit 1; }
source "$CONFIG_FILE"

workdir=""
OPTS=$(getopt -o t:h --long workdir:,threads:,help,\
min_genus_reads:,min_species_direct:,ratio_cutoff:,min_genus_occ: -- "$@")
[[ $? -ne 0 ]] && { help_message; exit 1; }
eval set -- "$OPTS"
while true; do
  case "$1" in
    --workdir)             workdir=$2; shift 2 ;;
    -t|--threads)          threads=$2; shift 2 ;;
    --min_genus_reads)     MIN_GENUS_READS=$2; shift 2 ;;
    --min_species_direct)  MIN_SPECIES_DIRECT=$2; shift 2 ;;
    --ratio_cutoff)        RATIO_CUTOFF=$2; shift 2 ;;
    --min_genus_occ)       MIN_GENUS_OCC=$2; shift 2 ;;
    -h|--help)             help_message; exit 0 ;;
    --) shift; break ;;
    *)  break ;;
  esac
done
[[ -z $workdir ]] && { echo "[ERROR] --workdir required"; exit 1; }

##########################################################################################
########################     constants (edit if needed)     ##########################
##########################################################################################
KRAKEN2_DB_SRC="/mnt/md0/wangyanren/database/release220/kraken2.gtdb_rep_release220/" # temporary!!!!!!!!
GTDB_TAXONOMY="/mnt/md0/wangyanren/database/release220/taxonomy_r220.processed.tsv" # temporary!!!!!!!!
GTDB_GENOMES_ID="/mnt/md0/wangyanren/database/release220/all_genomes.id.txt" # temporary!!!!!!!!

########################  prepare directories  ################################
mkdir -p "$workdir"/{02_kraken2_output,00_reference_db/chromosome,03_reference_db/chromosome}

########################  copy Kraken2 DB to /dev/shm (if absent) #############
if [[ ! -d "/dev/shm/kraken2.gtdb_rep_release220" || $(ls /dev/shm/kraken2.gtdb_rep_release220/*.k2d 2>/dev/null | wc -l) -ne $(ls "$KRAKEN2_DB_SRC"/*.k2d | wc -l) ]]; then
    echo "[INFO] Copying Kraken2 DB to /dev/shm..."
    mkdir -p /dev/shm/kraken2.gtdb_rep_release220
    cp "$KRAKEN2_DB_SRC"/*.k2d /dev/shm/kraken2.gtdb_rep_release220/
else
    echo "[INFO] /dev/shm/kraken2.gtdb_rep_release220 already exist. Skipping copy."
fi

##########################################################################################
########################                Run kraken2                           ##########################
##########################################################################################
mkdir -p "$workdir/02_kraken2_output"
logfile="$workdir/02_kraken2_output/kraken2.log"; : > "$logfile"

export workdir threads logfile
parallel --env workdir --env logfile -j "$threads" '
  sample=$(basename {} _R2_paired.fastq)
  {
    kraken2 --memory-mapping \
      --db /dev/shm/kraken2.gtdb_rep_release220 \
      "${workdir}/01_trim_SAGs/${sample}_R1_paired.fastq" \
      "${workdir}/01_trim_SAGs/${sample}_R2_paired.fastq" \
      --report "${workdir}/02_kraken2_output/${sample}.report" \
      --output "${workdir}/02_kraken2_output/${sample}.output"
  } &>> "$logfile"
#  echo "[INFO] ${sample} finished @ $(date)"
' ::: "$workdir"/01_trim_SAGs/*_R2_paired.fastq

##########################################################################################
########################         Filter report & get tax ids              ##########################
##########################################################################################

# Build an AWK filter to pull genus‐level hits (col 4 == "G") with ≥5000 reads:
cat > "$workdir/03_reference_db/chromosome/temp-filter.awk" <<'AWK'
BEGIN{FS=OFS="\t"}
FNR==1{sample=FILENAME;sub(".*/","",sample)}
$4=="G" && $2>min_genus_reads{
  tax=$6;for(i=7;i<=NF;i++)tax=tax OFS $i
  print $1,$2,$3,$4,$5,tax,sample
}
AWK
export MIN_GENUS_READS
awk_script="$workdir/03_reference_db/chromosome/temp-filter.awk"
sed -i "s/min_genus_reads/$MIN_GENUS_READS/" "$awk_script"
# Run that filter in parallel over all .report files:
> "$workdir/03_reference_db/chromosome/temp-filter_ori.tsv"
parallel --env workdir -j "$threads" \
  'awk -f '"$awk_script"' {}' \
  ::: "$workdir"/02_kraken2_output/*.report \
  >> "$workdir/03_reference_db/chromosome/temp-filter_ori.tsv"


# From temp-filter_ori.tsv (cols: pct, clade_count, direct_count, rank, taxid, taxon_name, filename),
# build temp-species.tsv by pulling matching lines from each sample’s .output:
> "${workdir}/03_reference_db/chromosome/temp-species.tsv"
while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7; do
  # col6 = genus (e.g. Escherichia), col7 = XXX.report
  report_path="${workdir}/02_kraken2_output/${col7}"
  if [[ ! -f "$report_path" ]]; then
    echo "[WARN] missing report file: $report_path" >&2
    continue
  fi
  # Match “genus + space” exactly (rank==S lines) and append "genus" back to original SAG sample "species" —— to make sure that grep all species in abundant genus
  grep -F "${col6} " "$report_path" \
    | awk -v fname="$col7" 'BEGIN{OFS="\t"} {print $0, fname}' \
    >> "${workdir}/03_reference_db/chromosome/temp-species.tsv"
done < "${workdir}/03_reference_db/chromosome/temp-filter_ori.tsv"


# Species-level filtering & GTDB mapping to ID
# Replace single spaces with underscores as in the legacy script, then remove any remaining spaces (tabs are kept)
sed -E \
  's/([0-9]) ([0-9])/\1_\2/g;
   s/([A-Za-z]) ([A-Za-z])/\1_\2/g;
   s/([0-9]) ([A-Za-z])/\1_\2/g;
   s/([A-Za-z]) ([0-9])/\1_\2/g;
   s/ //g' \
  "${workdir}/03_reference_db/chromosome/temp-species.tsv" \
  > "${workdir}/03_reference_db/chromosome/temp-species_name.tsv"

# Keep rows where direct_reads > MIN_SPECIES_DIRECT, take column 6 (species name), deduplicate
awk -v c="$MIN_SPECIES_DIRECT" '$2>c' \
  "$workdir/03_reference_db/chromosome/temp-species_name.tsv" |
  cut -f6 | awk '!a[$0]++' \
  > "$workdir/03_reference_db/chromosome/temp-species-list.txt"

# 5-3) Map species name → GTDB genome ID
awk 'NR==FNR { a[$1]; next } $2 in a { print $1 }' \
    "${workdir}/03_reference_db/chromosome/temp-species-list.txt" \
    "$GTDB_TAXONOMY" \
  | sed 's/^[^_]*_//' \
  > "${workdir}/03_reference_db/chromosome/temp-grep.species.txt"

# 5-4) Extract matching genome IDs
grep -Fxf \
  "${workdir}/03_reference_db/chromosome/temp-grep.species.txt" \
  "$GTDB_GENOMES_ID" \
  > "${workdir}/03_reference_db/chromosome/temp-id.gtdb_rep-species"


# Genus-level “unclassified species” filtering & GTDB mapping
# ratio > RATIO_CUTOFF (0.50) → take column 6 (genus) → count occurrences
# occurrences filter should be noticed, default is higher than MIN_GENUS_OCC (> 6), which means that only genera that appear multiple times of unclassified reads would be considered to expand
awk -v r="$RATIO_CUTOFF" '$3/$2>r{print $6}' \
  "$workdir/03_reference_db/chromosome/temp-filter_ori.tsv" |
  sort | uniq -c |
  awk -v m="$MIN_GENUS_OCC" '$1>m{print $2}' \
  > "$workdir/03_reference_db/chromosome/genus_unclassified_species.txt"

# 6-2) Map those genera → genome IDs
awk 'NR==FNR { g[$1]; next } $3 in g { print $1 }' \
    "${workdir}/03_reference_db/chromosome/genus_unclassified_species.txt" \
    "$GTDB_TAXONOMY" |
  sed 's/^[^_]*_//' \
  > "${workdir}/03_reference_db/chromosome/temp-grep.genus.txt"

# 6-3) Extract the corresponding genome IDs
grep -Fxf \
  "${workdir}/03_reference_db/chromosome/temp-grep.genus.txt" \
  "$GTDB_GENOMES_ID" \
  > "${workdir}/03_reference_db/chromosome/temp-id.gtdb_rep-genus"

# 7. Combine species- and genus-level ID lists
cat \
  "${workdir}/03_reference_db/chromosome/temp-id.gtdb_rep-species" \
  "${workdir}/03_reference_db/chromosome/temp-id.gtdb_rep-genus" \
| sort -u \
> "${workdir}/03_reference_db/chromosome/id.gtdb_rep.txt"

echo "[✓] Filter report & tax-ID extraction complete."
