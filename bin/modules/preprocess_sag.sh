#!/usr/bin/env bash
set -euo pipefail

########################################################################################################
########################                Scellmate: PREPROCESS MODULE                 ###################
########################################################################################################

help_message () {
cat <<EOF

Usage: $(basename "$0") -i <input_dir> -o <output_dir> --prefix <sample_prefix>

Required
  -i STR              input directory of raw SAGs (fastq.gz)
  -o STR              output base directory (workdir)
  --prefix STR        sample prefix (e.g. ST_Inf), will be formatted as scDNA_<prefix>_#####"

Optional
  -t, INF                    parallel jobs   (default: 30)
  -h, --help                 show this help

EOF
}

##########################################################################################
########################                parameter check                ###################
##########################################################################################
input=""
outdir=""
prefix=""
keep_raw=false # default to delete the SAG reads in 00_raw_SAGs
threads=30

OPTS=$(getopt -o i:o:t:h --long input:,output:,prefix:,keep-raw,help -- "$@")
if [ $? -ne 0 ]; then help_message; exit 1; fi
eval set -- "$OPTS"

while true; do
    case "$1" in
        -i | --input) input=$2; shift 2 ;;
        -o | --output) outdir=$2; shift 2 ;;
        --prefix) prefix=$2; shift 2 ;;
        --keep-raw) keep_raw=true; shift 1 ;;
        -t) threads=$2; shift 2 ;;
        -h | --help) help_message; exit 0 ;;
        --) shift; break ;;
        *) break ;;
    esac
done

if [[ -z "$input" || -z "$outdir" || -z "$prefix" ]]; then
    echo "[ERROR] -i, -o, and --prefix are required."
    help_message
    exit 1
fi

#######################################################################################
########################                run — rename                ###################
#######################################################################################
echo "[INFO] Starting renaming SAGs based on their sequencing depth and adding the prefix."

mkdir -p "$outdir/00_raw_SAGs"

prefix="scDNA_${prefix}"

find "$input" -name "*.fastq.gz" | parallel -j "$threads" 'gunzip -n {}'

echo "Filename,Count" > "$outdir/00_raw_SAGs/sequence_counts.csv"
for file in "$input"/*R1_001.fastq; do
    echo "${file},"$(wc -l < "$file" | awk '{print $1/4}') >> "$outdir/00_raw_SAGs/sequence_counts.csv"
done

(head -n 1 "$outdir/00_raw_SAGs/sequence_counts.csv" && tail -n +2 "$outdir/00_raw_SAGs/sequence_counts.csv" | sort -t, -k2,2nr) > "$outdir/00_raw_SAGs/temp1.csv"
awk -F, -v prefix="$prefix" 'BEGIN {OFS=","} NR==1 {print $0, "new_name"} NR>1 {print $1, $2, prefix "_" sprintf("%05d", NR-1) "_R1.fastq"}' \
    "$outdir/00_raw_SAGs/temp1.csv" > "$outdir/00_raw_SAGs/temp2.csv"
sed 's/_R1/_R2/g' "$outdir/00_raw_SAGs/temp2.csv" > "$outdir/00_raw_SAGs/temp3.csv"

sag_count=$(($(wc -l < "$outdir/00_raw_SAGs/temp2.csv") - 1))
echo "[INFO] Detecting $sag_count SAGs."

tail -n +2 "$outdir/00_raw_SAGs/temp2.csv" | while IFS=, read -r orig_filename count new_name; do
    cp "$orig_filename" "$outdir/00_raw_SAGs/$new_name"
done
tail -n +2 "$outdir/00_raw_SAGs/temp3.csv" | while IFS=, read -r orig_filename count new_name; do
    cp "$orig_filename" "$outdir/00_raw_SAGs/$new_name"
done

# detele the temp files
mv "$outdir/00_raw_SAGs/temp2.csv" "$outdir/00_raw_SAGs/rename_map_R1.csv"
mv "$outdir/00_raw_SAGs/temp3.csv" "$outdir/00_raw_SAGs/rename_map_R2.csv"
rm "$outdir/00_raw_SAGs/temp1.csv"

echo "[INFO] Done renaming, logs have been saved in $outdir/00_raw_SAGs/rename_map_R1.csv"
###########################################################################################
########################                run — trimmomatic                ###################
############################################################################################
echo "[INFO] Starting reads QC using Trimmomatic. This may take some time..."

mkdir -p "$outdir/01_trim_SAGs"
mkdir -p "$outdir/01_trim_SAGs/log_files"

# Monitor progress in background
(
    while true; do
        trimmed=$(ls "$outdir"/01_trim_SAGs/log_files/*.log 2>/dev/null | wc -l)
        echo -ne "[INFO] Trimming in progress: $trimmed / $sag_count SAGs done...\r"
        if [[ "$trimmed" -ge "$sag_count" ]]; then
            echo -e "\n[INFO] Trimming complete: $sag_count SAGs."
            break
        fi
        sleep 20
    done
) &

parallel -j "$threads" '
    FILE={};
    BASE=$(basename $FILE _R2.fastq);
    java -jar '"$TRIMMOMATIC"' PE -threads 1 -phred33 \
        -trimlog '"$outdir"'/01_trim_SAGs/log_files/${BASE}.log \
        '"$outdir"'/00_raw_SAGs/${BASE}_R1.fastq '"$outdir"'/00_raw_SAGs/${BASE}_R2.fastq \
        '"$outdir"'/01_trim_SAGs/${BASE}_R1_paired.fastq '"$outdir"'/01_trim_SAGs/${BASE}_R1_unpaired.fastq \
        '"$outdir"'/01_trim_SAGs/${BASE}_R2_paired.fastq '"$outdir"'/01_trim_SAGs/${BASE}_R2_unpaired.fastq \
        ILLUMINACLIP:'"$TRIMMOMATIC_ADAPTERS"':2:30:10:3:TRUE \
        LEADING:25 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30 \
        2> '"$outdir"'/01_trim_SAGs/log_files/${BASE}.stderr.txt
' ::: "$outdir"/00_raw_SAGs/*_R2.fastq

wait

# detele the fastq files in 00_raw_SAGs
if [[ "$keep_raw" == false ]]; then
    find "$outdir/00_raw_SAGs/" -name "*.fastq" -delete
fi

echo "Preprocess complete: logs of renaming SAGs in 00_raw_SAGs/*csv, trimmed SAGs in 01_trim_SAGs/*.fastq"
