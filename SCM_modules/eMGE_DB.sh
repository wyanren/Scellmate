#!/usr/bin/env bash
set -euo pipefail

help_message() {
cat <<EOF
Usage: $(basename "$0") --workdir <workdir> --script <script_dir> [options]

Required:
    --workdir STR              working directory
    --script STR               script directory

Optional:
    --threads INT              number of threads (default: 20)
    --supplement_assembly STR  optional fasta used for geNomad-based length supplementation
    -h, --help                 show this help message
EOF
}

workdir=""
script=""
threads="20"
supplement_assembly=""

OPTS=$(getopt -o h --long workdir:,script:,threads:,supplement_assembly:,help -- "$@")
if [ $? -ne 0 ]; then
  help_message >&2
  exit 1
fi
eval set -- "$OPTS"

while true; do
  case "$1" in
    --workdir) workdir="$2"; shift 2 ;;
    --script)  script="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --supplement_assembly) supplement_assembly="$2"; shift 2 ;;
    -h|--help) help_message; exit 0 ;;
    --) shift; break ;;
    *) break ;;
  esac
done

if [[ -z "$workdir" || -z "$script" ]]; then
  echo "Error: --workdir and --script are required." >&2
  help_message >&2
  exit 1
fi

workdir="$(realpath -m "$workdir")"
script="$(realpath "$script")"

if [[ -n "$supplement_assembly" ]]; then
  if [[ ! -s "$supplement_assembly" ]]; then
    echo "Error: --supplement_assembly file not found or empty: $supplement_assembly" >&2
    exit 1
  fi
  supplement_assembly="$(realpath "$supplement_assembly")"
fi


pick_representatives() (
  set -euo pipefail
  shopt -s nullglob

  local input="" summary="" output=""
  if [[ $# -eq 3 ]]; then
      input=$1; summary=$2; output=$3
  else
      while [[ $# -gt 0 ]]; do
        case $1 in
          -i|--input)    input=$2; shift 2 ;;
          -s|--summary)  summary=$2; shift 2 ;;
          -o|--output)   output=$2; shift 2 ;;
          -h|--help) echo "pick_representatives -i INPUT -s SUMMARY -o OUTPUT"; return 0 ;;
          *) echo "[pick_rep] unknown arg: $1" >&2; return 1 ;;
        esac
      done
  fi

  [[ -n $input && -n $summary && -n $output ]] || { echo "[pick_rep] ERR: input/summary/output not specified" >&2; return 1; }
  [[ -s $input   ]] || { echo "[pick_rep] ERR: $input not found" >&2; return 1; }
  [[ -s $summary ]] || { echo "[pick_rep] ERR: $summary not found" >&2; return 1; }

  mkdir -p "$(dirname "$output")" 2>/dev/null || true
  : > "$output"

  awk -F'\t' -v OFS='\t' '
    function tag_rep(r) { return (r ~ /No terminal/ ? "No_terminal_repeats" : r) }

    NR==FNR { info[$1]=$0; next }

    {
      split($2, seq, /,/); n=length(seq)

      if (n == 1) {
        if (seq[1] in info) {
          split(info[seq[1]], f)
          print seq[1], $2, "Solely_" tag_rep(f[3])
        }
        next
      }

      delete id; delete len; delete rep; delete dtr; delete cov
      for (i = 1; i <= n; i++) {
        k = seq[i]
        if (!(k in info)) next
        split(info[k], f)
        id[i]=f[1]
        len[i]=f[2]+0
        rep[i]=f[3]
        dtr[i]=f[4]
        match(f[1], /_cov_([0-9.]+)/, m)
        cov[i]=m[1]+0
      }

      hasDTR = 0
      for (i in rep) if (rep[i] == "DTR") hasDTR = 1

      if (!hasDTR) {
        pick=1
        for (i in len) if (len[i] > len[pick]) pick=i
        print id[pick], $2, "No_terminal_repeats"
        next
      }

      delete d; nd=0; delete freq
      for (i in rep) if (rep[i] == "DTR") { d[++nd]=i; freq[dtr[i]]++ }

      if (nd == 1) {
        print id[d[1]], $2, "DTR_only"
        next
      }

      maxf=0
      for (l in freq) if (freq[l] > maxf) maxf=freq[l]

      delete cand; nc=0
      for (j=1; j<=nd; j++) {
        i=d[j]
        if (freq[dtr[i]] == maxf) cand[++nc]=i
      }

      if (nc == 1) {
        print id[cand[1]], $2, "DTR_Most_Frequent"
        next
      }

      mcov=-1; delete best; nb=0
      for (j=1; j<=nc; j++) {
        i=cand[j]
        if (cov[i] > mcov) mcov=cov[i]
      }
      for (j=1; j<=nc; j++) {
        i=cand[j]
        if (cov[i] == mcov) best[++nb]=i
      }

      if (nb == 1) {
        print id[best[1]], $2, "DTR_Most_Frequent"
        next
      }

      pick=best[1]
      for (j=2; j<=nb; j++) if (len[best[j]] > len[pick]) pick=best[j]
      print id[pick], $2, "DTR_Longest_Most_Frequent"
    }
  ' "$summary" "$input" >> "$output"

  echo "[pick_rep] output → $output"
)


drop_clusters_by_ids() {
  local ids="$1"
  local input="$2"
  local output="$3"

  if [[ ! -s "$input" ]]; then
    : > "$output"
    return 0
  fi

  if [[ ! -s "$ids" ]]; then
    cp "$input" "$output"
    return 0
  fi

  awk -F'\t' '
    NR==FNR { ban[$1]=1; next }
    {
      n=split($2, a, /,/)
      hit=0
      for (i=1; i<=n; i++) {
        if (a[i] in ban) { hit=1; break }
      }
      if (!hit) print
    }
  ' "$ids" "$input" > "$output"
}


keep_clusters_by_ids() {
  local ids="$1"
  local input="$2"
  local output="$3"

  if [[ ! -s "$input" || ! -s "$ids" ]]; then
    : > "$output"
    return 0
  fi

  awk -F'\t' '
    NR==FNR { keep[$1]=1; next }
    {
      n=split($2, a, /,/)
      hit=0
      for (i=1; i<=n; i++) {
        if (a[i] in keep) { hit=1; break }
      }
      if (hit) print
    }
  ' "$ids" "$input" > "$output"
}


run_self_cluster() {
  local fna="$1"
  local prefix="$2"

  makeblastdb -in "$fna" -dbtype nucl -out "$prefix"
  blastn -query "$fna" \
         -db "$prefix" \
         -outfmt '6 std qlen slen' \
         -max_target_seqs 10000 \
         -out "${prefix}.tsv" \
         -num_threads "$threads"

  python "$script/anicalc.py" \
         -i "${prefix}.tsv" \
         -o "ani-$(basename "$prefix").tsv"

  python "$script/aniclust.py" \
         --fna "$fna" \
         --ani "ani-$(basename "$prefix").tsv" \
         --out "cluster-$(basename "$prefix").tsv" \
         --min_ani 99 --min_tcov 80 --min_qcov 0 --min_length 1000
}


mkdir -p "$workdir/07_eMGE_linkage" 2>/dev/null || true
mkdir -p "$workdir/07_eMGE_linkage/MGE_db" 2>/dev/null || true
mkdir -p "$workdir/07_eMGE_linkage/chromosome" 2>/dev/null || true

##########################
# step 0: geNomad on well-annotated genomes for chromosome masking
##########################
find "$workdir/07_eMGE_linkage/chromosome" -maxdepth 1 -type f -name '*.fasta' -delete

cp "$workdir"/06_CoSAG_assembly/round_ending/round_all/combined_full/well_annotated_genome/*.fasta \
   "$workdir/07_eMGE_linkage/chromosome/" 2>/dev/null || true

combine="$workdir/07_eMGE_linkage/chromosome_combine.fasta"
: > "$combine"

for fasta in "$workdir"/07_eMGE_linkage/chromosome/*.fasta; do
  [[ -e "$fasta" ]] || continue
  base=$(basename "$fasta" .fasta)
  awk -v pre="$base" '{
        if(/^>/){sub(/^>/,""); print ">" pre "__" $0}
        else print
      }' "$fasta" >> "$combine"
done

filterbyname.sh in="$combine" \
                out="$workdir/07_eMGE_linkage/chromosome_combine.filtered.fasta" \
                ow=t minlen=500 > /dev/null 2>&1

if ls "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad"/*_summary.log &>/dev/null; then
  echo "[INFO] clean_bin geNomad output already present → skip"
else
  echo "[INFO] Running geNomad on well-annotated genomes for masking ..."
  genomad end-to-end --cleanup --splits 8 --threads "$threads" \
        "$workdir/07_eMGE_linkage/chromosome_combine.filtered.fasta" \
        "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/" \
        "$SCELLMATE_DB/genomad_db" > /dev/null 2>&1
fi

# representative_sequences_bin.id for link_mge.sh masking
{
  if [[ -s "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/chromosome_combine.filtered_summary/chromosome_combine.filtered_plasmid.fna" ]]; then
    awk '/^>/{sub(/^>/,""); print}' \
      "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/chromosome_combine.filtered_summary/chromosome_combine.filtered_plasmid.fna"
  fi

  if [[ -s "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/chromosome_combine.filtered_summary/chromosome_combine.filtered_virus.fna" ]]; then
    awk '/^>/{sub(/^>/,""); print}' \
      "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/chromosome_combine.filtered_summary/chromosome_combine.filtered_virus.fna"
  fi
} | sort -u > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_bin_all.id"

awk '!/\|provirus_/' "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_bin_all.id" \
  > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_bin.id"

# clean SAG ids
cut -f 2 -d ' ' "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/final.reformat.add_cluster.edit_cont.annotate.species.txt" \
  > "$workdir/07_eMGE_linkage/clean_SAGs.id"

cd "$workdir/07_eMGE_linkage/MGE_db"


##########################
# optional supplement geNomad
##########################
supplement_plasmid_fna=""
supplement_plasmid_summary=""
supplement_virus_fna=""
supplement_virus_summary=""

if [[ -n "$supplement_assembly" ]]; then
  echo "[INFO] supplement assembly detected: $supplement_assembly"

  mkdir -p supplement_genomad
  filterbyname.sh in="$supplement_assembly" \
                  out="supplement.filtered.fasta" \
                  ow=t minlen=500 > /dev/null 2>&1

  if ls supplement_genomad/*_summary.log &>/dev/null; then
    echo "[INFO] supplement geNomad output already present → skip"
  else
    echo "[INFO] Running geNomad on supplement assembly ..."
    genomad end-to-end --cleanup --splits 8 --threads "$threads" \
      "supplement.filtered.fasta" \
      "supplement_genomad/" \
      "$SCELLMATE_DB/genomad_db" > /dev/null 2>&1
  fi

  supplement_plasmid_fna="supplement_genomad/supplement.filtered_summary/supplement.filtered_plasmid.fna"
  supplement_plasmid_summary="supplement_genomad/supplement.filtered_summary/supplement.filtered_plasmid_summary.tsv"
  supplement_virus_fna="supplement_genomad/supplement.filtered_summary/supplement.filtered_virus.fna"
  supplement_virus_summary="supplement_genomad/supplement.filtered_summary/supplement.filtered_virus_summary.tsv"
fi


##########################
# step 1: clean SAG plasmids (round 1)
##########################
filterbyname.sh in="$workdir/04_SAG_assembly/genomad_output/SAG_assembly_combine.500_summary/SAG_assembly_combine.500_plasmid.fna" \
                out="spades_output_modified_plasmid_clean.fna" \
                names="$workdir/07_eMGE_linkage/clean_SAGs.id" include=t substring=t ow=t > /dev/null 2>&1

run_self_cluster "spades_output_modified_plasmid_clean.fna" "blastn_plasmid_clean"

pick_representatives \
  --input   "cluster-blastn_plasmid_clean.tsv" \
  --summary "$workdir/04_SAG_assembly/genomad_output/SAG_assembly_combine.500_summary/SAG_assembly_combine.500_plasmid_summary.tsv" \
  --output  "representative_sequences_plasmid_1.tsv"

awk '/DTR/'  "representative_sequences_plasmid_1.tsv" > "representative_sequences_plasmid_1_DTR.tsv" || true
awk '!/DTR/' "representative_sequences_plasmid_1.tsv" > "representative_sequences_plasmid_1_no_DTR.tsv" || true

cut -f 1 "representative_sequences_plasmid_1_DTR.tsv" > "representative_sequences_plasmid_1_DTR.id" || true
cut -f 1 "representative_sequences_plasmid_1_no_DTR.tsv" > "representative_sequences_plasmid_1_no_DTR.id" || true

if [[ -s "representative_sequences_plasmid_1_DTR.id" ]]; then
  filterbyname.sh substring=name ow=t include=t minlen=0 \
                  names="representative_sequences_plasmid_1_DTR.id" \
                  in="spades_output_modified_plasmid_clean.fna" \
                  out="rep-plasmid_1_DTR.fna" > /dev/null 2>&1
else
  : > "rep-plasmid_1_DTR.fna"
fi

if [[ -s "representative_sequences_plasmid_1_no_DTR.id" ]]; then
  filterbyname.sh substring=name ow=t include=t minlen=0 \
                  names="representative_sequences_plasmid_1_no_DTR.id" \
                  in="spades_output_modified_plasmid_clean.fna" \
                  out="rep-plasmid_1_no_DTR.fna" > /dev/null 2>&1
else
  : > "rep-plasmid_1_no_DTR.fna"
fi


##########################
# step 2: optional supplement plasmids (round 2, only optimize non-DTR)
##########################
if [[ -n "$supplement_assembly" && -s "representative_sequences_plasmid_1_no_DTR.id" && -s "$supplement_plasmid_fna" ]]; then
  cat "rep-plasmid_1_DTR.fna" "rep-plasmid_1_no_DTR.fna" "$supplement_plasmid_fna" > "add_supplement_plasmid.fna"

  run_self_cluster "add_supplement_plasmid.fna" "blastn_add_supplement_plasmid"

  drop_clusters_by_ids \
    "representative_sequences_plasmid_1_DTR.id" \
    "cluster-blastn_add_supplement_plasmid.tsv" \
    "cluster-blastn_add_supplement_plasmid.noDTR.tsv"

  keep_clusters_by_ids \
    "representative_sequences_plasmid_1_no_DTR.id" \
    "cluster-blastn_add_supplement_plasmid.noDTR.tsv" \
    "cluster-blastn_add_supplement_plasmid.filter.tsv"

  if [[ -s "cluster-blastn_add_supplement_plasmid.filter.tsv" ]]; then
    cp "$workdir/04_SAG_assembly/genomad_output/SAG_assembly_combine.500_summary/SAG_assembly_combine.500_plasmid_summary.tsv" \
       "add_supplement_plasmid_summary.tsv"
    tail -n +2 "$supplement_plasmid_summary" >> "add_supplement_plasmid_summary.tsv"

    pick_representatives \
      --input   "cluster-blastn_add_supplement_plasmid.filter.tsv" \
      --summary "add_supplement_plasmid_summary.tsv" \
      --output  "representative_sequences_plasmid_2.tsv"

    cut -f 1 "representative_sequences_plasmid_2.tsv" > "representative_sequences_plasmid_2.id"

    filterbyname.sh substring=name ow=t include=t minlen=0 \
                    names="representative_sequences_plasmid_2.id" \
                    in="add_supplement_plasmid.fna" \
                    out="rep-plasmid_2.fna" > /dev/null 2>&1

    cat "rep-plasmid_1_DTR.fna" "rep-plasmid_2.fna" > "rep-plasmid.fna"
  else
    cat "rep-plasmid_1_DTR.fna" "rep-plasmid_1_no_DTR.fna" > "rep-plasmid.fna"
  fi
else
  cat "rep-plasmid_1_DTR.fna" "rep-plasmid_1_no_DTR.fna" > "rep-plasmid.fna"
fi


##########################
# step 3: clean SAG viruses (round 1)
##########################
filterbyname.sh in="$workdir/04_SAG_assembly/genomad_output/SAG_assembly_combine.500_summary/SAG_assembly_combine.500_virus.fna" \
                out="spades_output_modified_virus_clean.fna" \
                names="$workdir/07_eMGE_linkage/clean_SAGs.id" include=t substring=t ow=t > /dev/null 2>&1

run_self_cluster "spades_output_modified_virus_clean.fna" "blastn_virus_clean"

pick_representatives \
  --input   "cluster-blastn_virus_clean.tsv" \
  --summary "$workdir/04_SAG_assembly/genomad_output/SAG_assembly_combine.500_summary/SAG_assembly_combine.500_virus_summary.tsv" \
  --output  "representative_sequences_virus_1.tsv"

awk '/DTR/'  "representative_sequences_virus_1.tsv" > "representative_sequences_virus_1_DTR.tsv" || true
awk '!/DTR/' "representative_sequences_virus_1.tsv" > "representative_sequences_virus_1_no_DTR.tsv" || true

cut -f 1 "representative_sequences_virus_1_DTR.tsv" > "representative_sequences_virus_1_DTR.id" || true
cut -f 1 "representative_sequences_virus_1_no_DTR.tsv" > "representative_sequences_virus_1_no_DTR.id" || true

if [[ -s "representative_sequences_virus_1_DTR.id" ]]; then
  filterbyname.sh substring=name ow=t include=t minlen=0 \
                  names="representative_sequences_virus_1_DTR.id" \
                  in="spades_output_modified_virus_clean.fna" \
                  out="rep-virus_1_DTR.fna" > /dev/null 2>&1
else
  : > "rep-virus_1_DTR.fna"
fi

if [[ -s "representative_sequences_virus_1_no_DTR.id" ]]; then
  filterbyname.sh substring=name ow=t include=t minlen=0 \
                  names="representative_sequences_virus_1_no_DTR.id" \
                  in="spades_output_modified_virus_clean.fna" \
                  out="rep-virus_1_no_DTR.fna" > /dev/null 2>&1
else
  : > "rep-virus_1_no_DTR.fna"
fi


##########################
# step 4: optional supplement viruses (round 2, only optimize non-DTR)
##########################
if [[ -n "$supplement_assembly" && -s "representative_sequences_virus_1_no_DTR.id" && -s "$supplement_virus_fna" ]]; then
  cat "rep-virus_1_DTR.fna" "rep-virus_1_no_DTR.fna" "$supplement_virus_fna" > "add_supplement_virus.fna"

  run_self_cluster "add_supplement_virus.fna" "blastn_add_supplement_virus"

  drop_clusters_by_ids \
    "representative_sequences_virus_1_DTR.id" \
    "cluster-blastn_add_supplement_virus.tsv" \
    "cluster-blastn_add_supplement_virus.noDTR.tsv"

  keep_clusters_by_ids \
    "representative_sequences_virus_1_no_DTR.id" \
    "cluster-blastn_add_supplement_virus.noDTR.tsv" \
    "cluster-blastn_add_supplement_virus.filter.tsv"

  if [[ -s "cluster-blastn_add_supplement_virus.filter.tsv" ]]; then
    cp "$workdir/04_SAG_assembly/genomad_output/SAG_assembly_combine.500_summary/SAG_assembly_combine.500_virus_summary.tsv" \
       "add_supplement_virus_summary.tsv"
    tail -n +2 "$supplement_virus_summary" >> "add_supplement_virus_summary.tsv"

    pick_representatives \
      --input   "cluster-blastn_add_supplement_virus.filter.tsv" \
      --summary "add_supplement_virus_summary.tsv" \
      --output  "representative_sequences_virus_2.tsv"

    cut -f 1 "representative_sequences_virus_2.tsv" > "representative_sequences_virus_2.id"

    filterbyname.sh substring=name ow=t include=t minlen=0 \
                    names="representative_sequences_virus_2.id" \
                    in="add_supplement_virus.fna" \
                    out="rep-virus_2.fna" > /dev/null 2>&1

    cat "rep-virus_1_DTR.fna" "rep-virus_2.fna" > "rep-virus.fna"
  else
    cat "rep-virus_1_DTR.fna" "rep-virus_1_no_DTR.fna" > "rep-virus.fna"
  fi
else
  cat "rep-virus_1_DTR.fna" "rep-virus_1_no_DTR.fna" > "rep-virus.fna"
fi


##########################
# step 5: separate phage / prophage
##########################
grep 'provirus' "rep-virus.fna" | sed 's/>//' > "provirus.id" || true

if [[ -s "provirus.id" ]]; then
  filterbyname.sh substring=name ow=t include=t minlen=0 \
      names="provirus.id" \
      in="rep-virus.fna" \
      out="rep-prophage.fna" > /dev/null 2>&1

  filterbyname.sh substring=name ow=t include=f minlen=0 \
      names="provirus.id" \
      in="rep-virus.fna" \
      out="rep-phage.fna" > /dev/null 2>&1
else
  cp "rep-virus.fna" "rep-phage.fna"
  : > "rep-prophage.fna"
fi

echo -e "\033[0;35m[Scellmate] Construction of the eMGE database from clean SAGs has been completed.\033[0m"
