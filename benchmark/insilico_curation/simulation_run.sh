#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash resume_simulation.sh 62
# Optional:
#   bash resume_simulation.sh 62 run
#   bash resume_simulation.sh 62 dryrun
#
# Default mode = dryrun
#   dryrun : only generate/check task tables
#   run    : actually run pending tasks

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: bash $0 <set_id> [run|dryrun]" >&2
    exit 1
fi

SET_ID="$1"
MODE="${2:-dryrun}"

if [[ "${MODE}" != "run" && "${MODE}" != "dryrun" ]]; then
    echo "Error: mode must be 'run' or 'dryrun'" >&2
    exit 1
fi

# =========================
# Basic paths and settings
# =========================
WORKDIR=/mnt/md0/wangyanren/workdir/202603_MDA
INPUT_DIR="${WORKDIR}/genome_dir_${SET_ID}"
TSV="${WORKDIR}/selected_4_genus_triplets.no_placeholder__${SET_ID}.tsv"
MAKE_SCRIPT="${WORKDIR}/make_mda_template3.py"

THREADS=20
SINGLET_CONTAMS=("0" "5" "10" "20" "30" "40")
DOUBLET_CONTAMS=("0")

# singlet: 6 seeds
SINGLET_SEEDS=$(seq 1 6)
# doublet: 3 seeds
DOUBLET_SEEDS=$(seq 1 3)

NUM_GENOMES=12
PRIMERS=../202503_MDAsimulation/preliminary_test/mdasim/examples/primerList.fasta

INTERMEDIATE_DIR="${WORKDIR}/MDA_Amplicon-intermediate-files-${SET_ID}"

# =========================
# Basic checks
# =========================
[[ -d "${INPUT_DIR}" ]] || { echo "Error: INPUT_DIR not found: ${INPUT_DIR}" >&2; exit 1; }
[[ -f "${TSV}" ]] || { echo "Error: TSV not found: ${TSV}" >&2; exit 1; }
[[ -f "${MAKE_SCRIPT}" ]] || { echo "Error: MAKE_SCRIPT not found: ${MAKE_SCRIPT}" >&2; exit 1; }

cd "${WORKDIR}"
mkdir -p "${INTERMEDIATE_DIR}"

echo "[Info] SET_ID          : ${SET_ID}"
echo "[Info] MODE            : ${MODE}"
echo "[Info] INPUT_DIR       : ${INPUT_DIR}"
echo "[Info] TSV             : ${TSV}"
echo "[Info] INTERMEDIATE_DIR: ${INTERMEDIATE_DIR}"

# =========================
# Build singlet genome list
# =========================
cut -f4 "${TSV}" | tail -n +2 | sed 's/$/_genomic/' > singlet_genomes.list

# =========================
# Build intra-genus doublet pairs:
# within each genus, all 3 choose 2 pairs
# output format: genome1<TAB>genome2
# =========================
awk -F'\t' '
NR>1 {
    key=$1 FS $2 FS $3;   # Group, Domain, Genus
    genomes[key]=genomes[key] $4 "_genomic\t";
}
END {
    for (k in genomes) {
        n=split(genomes[k], arr, "\t");
        cnt=0;
        delete g;

        for (i=1; i<=n; i++) {
            if (arr[i] != "") {
                cnt++;
                g[cnt]=arr[i];
            }
        }

        if (cnt >= 3) {
            print g[1] "\t" g[2];
            print g[1] "\t" g[3];
            print g[2] "\t" g[3];
        }
    }
}
' "${TSV}" > doublet_pairs.list

# =========================
# Main worker
# =========================
pipeline_task() {
    local mode=$1
    local genome1=$2
    local genome2=$3
    local contaminant=$4
    local seed=$5

    local prefix=""

    source ~/conda/etc/profile.d/conda.sh

    if [[ "${mode}" == "singlet" ]]; then
        prefix="${genome1}.singlet.c${contaminant}.s$(printf '%02d' "${seed}")"
        echo "[${prefix}] start"

        conda activate ~/conda/envs/microbe-seq/
        python "${MAKE_SCRIPT}" \
            --input_dir "${INPUT_DIR}" \
            --number "${NUM_GENOMES}" \
            --mode singlet \
            --genome "${genome1}" \
            --seed "${seed}" \
            --contaminated "${contaminant}" \
            --tsv \
            --output_prefix "${prefix}"

    elif [[ "${mode}" == "doublet" ]]; then
        prefix="${genome1}__${genome2}.doublet.c${contaminant}.s$(printf '%02d' "${seed}")"
        echo "[${prefix}] start"

        conda activate ~/conda/envs/microbe-seq/
        python "${MAKE_SCRIPT}" \
            --input_dir "${INPUT_DIR}" \
            --number "${NUM_GENOMES}" \
            --mode doublet \
            --genome "${genome1}" \
            --doublet_genome "${genome2}" \
            --seed "${seed}" \
            --contaminated "${contaminant}" \
            --tsv \
            --output_prefix "${prefix}"
    else
        echo "Unknown mode: ${mode}" >&2
        return 1
    fi

    conda activate ~/conda/envs/mdasim/
    mdasim \
        --input="./${prefix}_merged.fasta" \
        --coverage=20 \
        --primers="${PRIMERS}" \
        --output="${prefix}-" \
        --log="${prefix}_mdasim_errors.log" \
        > "${prefix}_mdasim_run.log"

    art_illumina \
        -ss HS20 \
        -i "${prefix}-Amplicons.fasta" \
        -o "${prefix}-" \
        -l 100 \
        -f 6 \
        -p \
        -m 500 \
        -s 10

    reformat.sh \
        in1="${prefix}-1.fq" \
        in2="${prefix}-2.fq" \
        out1="${prefix}-sub_1.fq" \
        out2="${prefix}-sub_2.fq" \
        samplebasestarget=5000000 \
        ow=t

    rm -f "${prefix}-1.fq" "${prefix}-2.fq"
    rm -f "${prefix}-1.aln" "${prefix}-2.aln"

    # move outputs into intermediate dir
    mkdir -p "${INTERMEDIATE_DIR}"
    mv -f "${prefix}_merged.fasta" "${INTERMEDIATE_DIR}/" 2>/dev/null || true
    mv -f "${prefix}_contamination_counts.tsv" "${INTERMEDIATE_DIR}/" 2>/dev/null || true
    mv -f "${prefix}_mdasim_errors.log" "${INTERMEDIATE_DIR}/" 2>/dev/null || true
    mv -f "${prefix}_mdasim_run.log" "${INTERMEDIATE_DIR}/" 2>/dev/null || true
    mv -f "${prefix}-Amplicons.fasta" "${INTERMEDIATE_DIR}/" 2>/dev/null || true
    mv -f "${prefix}-sub_1.fq" "${INTERMEDIATE_DIR}/" 2>/dev/null || true
    mv -f "${prefix}-sub_2.fq" "${INTERMEDIATE_DIR}/" 2>/dev/null || true

    echo "[${prefix}] done"
}

export -f pipeline_task
export WORKDIR INPUT_DIR TSV MAKE_SCRIPT THREADS NUM_GENOMES PRIMERS INTERMEDIATE_DIR

# =========================
# Prepare task tables
# =========================
: > singlet_tasks.tsv
while read -r genome; do
    for contaminant in "${SINGLET_CONTAMS[@]}"; do
        for seed in ${SINGLET_SEEDS}; do
            printf "singlet\t%s\tNA\t%s\t%s\n" "${genome}" "${contaminant}" "${seed}" >> singlet_tasks.tsv
        done
    done
done < singlet_genomes.list

: > doublet_tasks.tsv
while IFS=$'\t' read -r genome1 genome2; do
    for contaminant in "${DOUBLET_CONTAMS[@]}"; do
        for seed in ${DOUBLET_SEEDS}; do
            printf "doublet\t%s\t%s\t%s\t%s\n" "${genome1}" "${genome2}" "${contaminant}" "${seed}" >> doublet_tasks.tsv
        done
    done
done < doublet_pairs.list

cat singlet_tasks.tsv doublet_tasks.tsv > all_tasks.tsv

# =========================
# Check completed tasks
# Completion standard:
#   ${prefix}-sub_1.fq and ${prefix}-sub_2.fq both exist in INTERMEDIATE_DIR
# =========================
: > pending_tasks.tsv
: > completed_tasks.tsv

while IFS=$'\t' read -r mode genome1 genome2 contaminant seed; do
    if [[ "${mode}" == "singlet" ]]; then
        prefix="${genome1}.singlet.c${contaminant}.s$(printf '%02d' "${seed}")"
    elif [[ "${mode}" == "doublet" ]]; then
        prefix="${genome1}__${genome2}.doublet.c${contaminant}.s$(printf '%02d' "${seed}")"
    else
        echo "Unknown mode in task table: ${mode}" >&2
        exit 1
    fi

    fq1="${INTERMEDIATE_DIR}/${prefix}-sub_1.fq"
    fq2="${INTERMEDIATE_DIR}/${prefix}-sub_2.fq"

    if [[ -s "${fq1}" && -s "${fq2}" ]]; then
        printf "%s\t%s\t%s\t%s\t%s\n" "${mode}" "${genome1}" "${genome2}" "${contaminant}" "${seed}" >> completed_tasks.tsv
    else
        printf "%s\t%s\t%s\t%s\t%s\n" "${mode}" "${genome1}" "${genome2}" "${contaminant}" "${seed}" >> pending_tasks.tsv
    fi
done < all_tasks.tsv

echo "[Summary]"
echo "Singlet tasks   : $(wc -l < singlet_tasks.tsv)"
echo "Doublet tasks   : $(wc -l < doublet_tasks.tsv)"
echo "All tasks       : $(wc -l < all_tasks.tsv)"
echo "Completed tasks : $(wc -l < completed_tasks.tsv)"
echo "Pending tasks   : $(wc -l < pending_tasks.tsv)"

# =========================
# Run only pending tasks
# =========================
if [[ "${MODE}" == "run" ]]; then
    if [[ -s pending_tasks.tsv ]]; then
        parallel --colsep '\t' --jobs "${THREADS}" \
            pipeline_task {1} {2} {3} {4} {5} \
            :::: pending_tasks.tsv
    else
        echo "No pending tasks. Nothing to run."
    fi
else
    echo "[Dryrun] pending tasks were not executed."
fi

# =========================
# Move bookkeeping files into intermediate dir
# =========================
mv -f all_tasks.tsv "${INTERMEDIATE_DIR}/all_tasks.tsv"
mv -f completed_tasks.tsv "${INTERMEDIATE_DIR}/completed_tasks.tsv"
mv -f pending_tasks.tsv "${INTERMEDIATE_DIR}/pending_tasks.tsv"
mv -f doublet_pairs.list "${INTERMEDIATE_DIR}/doublet_pairs.list"
mv -f doublet_tasks.tsv "${INTERMEDIATE_DIR}/doublet_tasks.tsv"
mv -f singlet_genomes.list "${INTERMEDIATE_DIR}/singlet_genomes.list"
mv -f singlet_tasks.tsv "${INTERMEDIATE_DIR}/singlet_tasks.tsv"

echo "Bookkeeping files moved to ${INTERMEDIATE_DIR}"
