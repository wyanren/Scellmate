#!/usr/bin/env bash
set -euo pipefail

################################################################################
########################      Scellmate: MARKER_GENE_BUILD MODULE        #######
################################################################################

help_message () {
cat <<EOF

Usage: $(basename "$0") --workdir DIR [options]

Required
  --workdir DIR            workflow base directory (same as main Scellmate workdir)

#Reference-DB modes (choose **one**, default = --gtdb-only)
  #--gtdb-only              use GTDB subset only (ids in \$workdir/03_reference_db/.../id.gtdb_rep.txt)
  #--gtdb-add               GTDB subset **plus** genomes in --add DIR
  #--add-only               skip GTDB; use only genomes in --add DIR
  #--add DIR                directory containing *.fna[.gz] files to merge (if choose "--gtdb-add" or "add-only")

General options
  -t INT, --threads INT    CPU threads for copy/decompress & CheckM   (default: 40)
  -h, --help               show this help and exit

EOF
}

################################################################################
## Default parameters
################################################################################
threads=20
workdir=""

################################################################################
## Parse CLI
################################################################################
# gtdb-only,gtdb-add,add-only,add:,\
OPTS=$(getopt -o ht: \
  --long workdir:,threads:,help,min_genus_reads:,min_species_direct:,ratio_cutoff:,min_genus_occ: \
  -- "$@")

[[ $? -ne 0 ]] && { help_message; exit 1; }
eval set -- "$OPTS"
while true; do
  case "$1" in
    --workdir)              workdir=$2; shift 2 ;;
    -t|--threads)           threads=$2; shift 2 ;;
    -h|--help)              help_message; exit 0 ;;
    --) shift; break ;;
    *)  break ;;
  esac
done

[[ -z $workdir ]] && { echo "[ERROR] --workdir required"; exit 1; }

################################################################################
## Constants (edit if needed)
################################################################################
GTDB_REP_DIR="$SCELLMATE_GTDB_PATH/gtdb_rep"          # *_genomic.fna.gz
GTDB_TAXONOMY="$SCELLMATE_GTDB_PATH/taxonomy_r220.processed.tsv"
GTDB_GENOMES_ID="$SCELLMATE_GTDB_PATH/all_genomes.id.txt"

ref_root="$workdir/03_reference_db/chromosome"
mkdir -p "$ref_root"/{gtdb_rep_sub,reference_DB,checkm}

# 1) Fetch GTDB genomes (unless add-only)
# if [[ $DB_MODE != "add-only" ]]; then
  echo "[INFO] Copying GTDB genomes to reference_DB …"
  xargs -P "$threads" -I {} sh -c \
    'cp "'"$GTDB_REP_DIR"'/{}_genomic.fna.gz" "'"$ref_root/gtdb_rep_sub"'" && \
     gunzip -f "'"$ref_root/gtdb_rep_sub"'"/{}_genomic.fna.gz' \
    < "$ref_root/id.gtdb_rep.txt"
  cp "$ref_root"/gtdb_rep_sub/*.fna "$ref_root/reference_DB/"
# fi

# 2) Add user-supplied genomes if requested
#if [[ $DB_MODE != "gtdb-only" ]]; then
#  echo "[INFO] Adding user genomes from $ADD_DIR"
#  find "$ADD_DIR" -maxdepth 1 -type f -name '*.fna*' -print0 |
#    xargs -0 -I{} cp {} "$ref_root/reference_DB/"
#  gunzip -fq "$ref_root"/reference_DB/*.gz 2>/dev/null || true
#fi

# 3) Concatenate all reference genomes (optional but handy)
cat "$ref_root"/reference_DB/*.fna > "$ref_root/combine_reference_DB.fasta"

checkm data setRoot "$SCELLMATE_DB"/checkM_db/ > /dev/null 2>&1

# 4) Run CheckM (tree → lineage_set → analyze → location_info)
echo "[INFO] Running CheckM (threads = $threads)…"
(
  cd "$ref_root" || exit
  checkm tree -x fna reference_DB -t "$threads" checkm > /dev/null 2>&1
  checkm lineage_set checkm checkm/marker_file > /dev/null 2>&1
  checkm analyze checkm/marker_file -x fna reference_DB checkm/analyze -t "$threads" > /dev/null 2>&1
  checkm tree_qa checkm -o 2 --tab_table > checkm/tree_qa.txt > /dev/null 2>&1
  checkm qa checkm/marker_file checkm/analyze -o 9 -q -f checkm/qa.9-marker_location.txt --tab_table
)

# 5) Format single-copy gene metadata
awk 'BEGIN{OFS="\t";print "GeneID","Chr","Start","End","Strand"}
     NR>1{strand=($6==1?"+":($6==-1?"-":""));print $2"_"$3,$2,$4,$5,strand}' \
  "$ref_root/checkm/qa.9-marker_location.txt" \
  > "$ref_root/metadata-SCG.txt"


echo "    → subsetting prokaryotic genomes from GTDB representative:  $ref_root/reference_DB/"
echo "    → Single-Copy-Gene (SCG) metadata of genome dataset: $ref_root/metadata-SCG.txt"
echo -e "\033[0;32m[Scellmate] Construction of marker-gene metadata in first QC have been completed.\033[0m"
