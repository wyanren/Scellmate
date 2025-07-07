#!/usr/bin/env bash
set -euo pipefail

help_message() {
cat <<EOF
Usage: $(basename "$0") --workdir <workdir> --script <script_dir>

Required:
    --workdir STR      working directory (e.g. /path/to/workdir)
    --script STR       script directory (e.g. /path/to/scripts)
Optional:
    --genomad_db DIR   geNomad database directory                       [/mnt/md0/wangyanren/database/genomad/genomad_db]
    -h, --help         show this help message
EOF
}

# parameter parsing
default_workdir="workdir"
default_script="bin/scripts"
workdir=""
script=""
genomad_db="/mnt/md0/wangyanren/database/genomad/genomad_db"

OPTS=$(getopt -o h --long workdir:,script:,genomad_db:,help -- "$@")
if [ $? -ne 0 ]; then
  help_message >&2
  exit 1
fi
eval set -- "$OPTS"

while true; do
  case "$1" in
    --workdir) workdir="$2"; shift 2 ;;
    --script)  script="$2"; shift 2 ;;
    --genomad_db) genomad_db="$2"; shift 2 ;;
    -h|--help) help_message; exit 0 ;;
    --) shift; break ;;
    *) break ;;
  esac
done

# 检查必需参数
if [[ -z "$workdir" || -z "$script" ]]; then
  echo "Error: --workdir and --script are required." >&2
  help_message >&2
  exit 1
fi

pick_representatives() (
  set -euo pipefail
  shopt -s nullglob

  # ---------- CLI parsing ----------
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
  [[ -n $input && -n $summary && -n $output ]] ||
      { echo "[pick_rep] ERR: input/summary/output not specified" >&2; return 1; }
  [[ -s $input   ]] || { echo "[pick_rep] ERR: $input not found" >&2;   return 1; }
  [[ -s $summary ]] || { echo "[pick_rep] ERR: $summary not found" >&2; return 1; }

  mkdir -p "$(dirname "$output")" 2>/dev/null || true
  : > "$output"

  # ---------- AWK one-pass processing ----------
  awk -F'\t' -v OFS='\t' '
    # helper function must be at top-level (portable across mawk / gawk)
    function tag_rep(r) { return (r ~ /No terminal/ ? "No_terminal_repeats" : r) }

    # --- load SUMMARY into array ---
    NR==FNR { info[$1]=$0; next }

    # --- process each cluster line ---
    {
      split($2, seq, /,/); n=length(seq)

      # single-member cluster
      if (n == 1) {
        if (seq[1] in info) {
          split(info[seq[1]], f)
          print seq[1], $2, "Solely_" tag_rep(f[3])
        }
        next
      }

      # multi-member cluster: collect per-seq stats
      delete id; delete len; delete rep; delete dtr; delete cov
      for (i = 1; i <= n; i++) {
        k = seq[i]; if (!(k in info)) next
        split(info[k], f)
        id[i]=f[1]; len[i]=f[2]+0; rep[i]=f[3]; dtr[i]=f[4]
        match(f[1], /_cov_([0-9.]+)/, m); cov[i]=m[1]+0
      }

      # no DTR in cluster  -> choose longest
      hasDTR = 0
      for (i in rep) if (rep[i]=="DTR") hasDTR=1
      if (!hasDTR) {
        pick=1; for (i in len) if (len[i] > len[pick]) pick=i
        print id[pick], $2, "No_terminal_repeats"; next
      }

      # collect DTR sequences
      delete d; nd=0; delete freq
      for (i in rep) if (rep[i]=="DTR"){ d[++nd]=i; freq[dtr[i]]++ }
      if (nd==1){ print id[d[1]], $2, "DTR_only"; next }

      # most frequent DTR length
      maxf=0; for (l in freq) if (freq[l]>maxf) maxf=freq[l]
      delete cand; nc=0
      for (j=1; j<=nd; j++){ i=d[j]; if(freq[dtr[i]]==maxf) cand[++nc]=i }
      if (nc==1){ print id[cand[1]], $2, "DTR_Most_Frequent"; next }

      # highest coverage among candidates
      mcov=-1; delete best; nb=0
      for (j=1; j<=nc; j++){ i=cand[j]; if(cov[i]>mcov) mcov=cov[i] }
      for (j=1; j<=nc; j++){ i=cand[j]; if(cov[i]==mcov) best[++nb]=i }
      if (nb==1){ print id[best[1]], $2, "DTR_Most_Frequent"; next }

      # longest among remaining ties
      pick=best[1]; for (j=2; j<=nb; j++) if (len[best[j]]>len[pick]) pick=best[j]
      print id[pick], $2, "DTR_Longest_Most_Frequent"
    }
  ' "$summary" "$input" >> "$output"

  echo "[pick_rep] output → $output"
)

##########################--------------Start--------------#########################
# step 0: generate MGE from well annotated SAGs and CoSAGs
mkdir -p "$workdir/07_eMGE_linkage/MGE_db" 2>/dev/null || true
# combine all CoSAGs into one fasta file
mkdir -p "$workdir/07_eMGE_linkage/chromosome/" 2>/dev/null || true
cp "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/well_annotated_genome/"*.fasta "$workdir/07_eMGE_linkage/chromosome/"
combine="$workdir/07_eMGE_linkage/chromosome_combine.fasta"
> "$combine"
for fasta in "$workdir"/07_eMGE_linkage/chromosome/*.fasta; do
  base=$(basename "$fasta" .fasta)
  awk -v pre="$base" '{
        if(/^>/){sub(/^>/,""); print ">"pre"_"$0}
        else print
      }' "$fasta" >> "$combine"
done

filterbyname.sh in="$combine" out="$workdir/07_eMGE_linkage/chromosome_combine.filtered.fasta" ow=t minlen=500

if ls "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad"/*_summary.log &>/dev/null; then
  echo "[INFO] geNomad output already present → skip"
else
  echo "[INFO] Running geNomad …"
  conda run -n genomad --no-capture-output genomad end-to-end --cleanup --splits 8 \
        "$workdir/07_eMGE_linkage/chromosome_combine.filtered.fasta" \
        "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/" \
        "$genomad_db"
fi

ll workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/chromosome_combine.filtered_summary/chromosome_combine.filtered_plasmid.fna
ll workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/chromosome_combine.filtered_summary/chromosome_combine.filtered_virus.fna


# step1 : Cataloguing individual well-annotated SAGs
# cp "$workdir/07_eMGE_linkage/MGE_db_copy/individual_genomad/spades_output_modified_plasmid_clean.fna" "$workdir/07_eMGE_linkage/MGE_db"
cut -f 2 -d ' ' "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/final.reformat.add_cluster.edit_cont.annotate.species.txt" > "$workdir/07_eMGE_linkage/clean_SAGs.id"
filterbyname.sh in="$workdir/04_SAG_assembly/genomad_output/SAG_assembly_combine.500_summary/SAG_assembly_combine.500_plasmid.fna" \
                out="$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_plasmid_clean.fna" \
                names="$workdir/07_eMGE_linkage/clean_SAGs.id" include=t substring=t ow=t

makeblastdb -in "$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_plasmid_clean.fna" -dbtype nucl -out "$workdir/07_eMGE_linkage/MGE_db/blastn_plasmid_clean"
blastn -query "$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_plasmid_clean.fna" \
       -db "$workdir/07_eMGE_linkage/MGE_db/blastn_plasmid_clean" \
       -outfmt '6 std qlen slen' -max_target_seqs 10000 \
       -out "$workdir/07_eMGE_linkage/MGE_db/blastn_plasmid_clean.tsv" -num_threads 48
python "$script/anicalc.py" -i "$workdir/07_eMGE_linkage/MGE_db/blastn_plasmid_clean.tsv" -o "$workdir/07_eMGE_linkage/MGE_db/ani-blastn_plasmid_clean.tsv"
python "$script/aniclust.py" --fna "$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_plasmid_clean.fna" \
                             --ani "$workdir/07_eMGE_linkage/MGE_db/ani-blastn_plasmid_clean.tsv" \
                             --out "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_plasmid_clean.tsv" --min_ani 99 --min_tcov 90 --min_qcov 0 --min_length 1000


pick_representatives \
  --input   "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_plasmid_clean.tsv" \
  --summary "$workdir/07_eMGE_linkage/MGE_db_copy/individual_genomad/spades_output_modified_plasmid_summary.tsv" \
  --output  "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_1.tsv"

grep 'DTR' "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_1.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_1_DTR.tsv"
grep -v 'DTR' "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_1.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_1_no_DTR.tsv"

cut -f 1 "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_1_DTR.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_1_DTR.id"
cut -f 1 "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_1_no_DTR.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_1_no_DTR.id"

filterbyname.sh substring=name ow=t include=t minlen=0 \
                names="$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_1_DTR.id" \
                in="$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_plasmid_clean.fna" \
                out="$workdir/07_eMGE_linkage/MGE_db/rep-plasmid_1_DTR.fna"

filterbyname.sh substring=name ow=t include=t minlen=0 \
                names="$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_1_no_DTR.id" \
                in="$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_plasmid_clean.fna" \
                out="$workdir/07_eMGE_linkage/MGE_db/rep-plasmid_1_no_DTR.fna"



# step2 : Cataloguing CoSAGs
cp "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/chromosome_combine.filtered_summary/chromosome_combine.filtered_plasmid.fna" "$workdir/07_eMGE_linkage/MGE_db/temp.fna"
cp "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/chromosome_combine.filtered_summary/chromosome_combine.filtered_plasmid_summary.tsv" "$workdir/07_eMGE_linkage/MGE_db/merged-well_annotated_genome_plasmid_summary.tsv"
awk '{print $1}' "$workdir/07_eMGE_linkage/final.reformat.add_cluster.edit_cont.annotate.species.txt" | sort | uniq > "$workdir/07_eMGE_linkage/MGE_db/remain.cluster"
filterbyname.sh substring=name ow=t include=t minlen=0 names="$workdir/07_eMGE_linkage/MGE_db/remain.cluster" \
                in="$workdir/07_eMGE_linkage/MGE_db/temp.fna" \
                out="$workdir/07_eMGE_linkage/MGE_db/merged-well_annotated_genome_plasmid.fna"

# DTR sequences are added here to avoid including those already present; this is only for assigning clusters to the no_DTR round of sequences
cat "$workdir/07_eMGE_linkage/MGE_db/rep-plasmid_1_DTR.fna" "$workdir/07_eMGE_linkage/MGE_db/rep-plasmid_1_no_DTR.fna" "$workdir/07_eMGE_linkage/MGE_db/merged-well_annotated_genome_plasmid.fna" > "$workdir/07_eMGE_linkage/MGE_db/add_bin_plasmid.fna"
makeblastdb -in "$workdir/07_eMGE_linkage/MGE_db/add_bin_plasmid.fna" -dbtype nucl -out "$workdir/07_eMGE_linkage/MGE_db/add_bin_plasmid"
blastn -query "$workdir/07_eMGE_linkage/MGE_db/add_bin_plasmid.fna" -db "$workdir/07_eMGE_linkage/MGE_db/add_bin_plasmid" -outfmt '6 std qlen slen' -max_target_seqs 10000 -out "$workdir/07_eMGE_linkage/MGE_db/blastn_add_bin_plasmid.tsv" -num_threads 48
python "$script/anicalc.py" -i "$workdir/07_eMGE_linkage/MGE_db/blastn_add_bin_plasmid.tsv" -o "$workdir/07_eMGE_linkage/MGE_db/ani-blastn_add_bin_plasmid.tsv"
python "$script/aniclust.py" --fna "$workdir/07_eMGE_linkage/MGE_db/add_bin_plasmid.fna" --ani "$workdir/07_eMGE_linkage/MGE_db/ani-blastn_add_bin_plasmid.tsv" --out "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_add_bin_plasmid.tsv" --min_ani 99 --min_tcov 90 --min_qcov 0 --min_length 1000
grep '>' "$workdir/07_eMGE_linkage/MGE_db/rep-plasmid_1_DTR.fna" | sed 's/>//' > "$workdir/07_eMGE_linkage/MGE_db/temp"
grep -F -v -f "$workdir/07_eMGE_linkage/MGE_db/temp" "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_add_bin_plasmid.tsv" > "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_add_bin_plasmid.filter.tsv"

cat "$workdir/07_eMGE_linkage/MGE_db/merged-well_annotated_genome_plasmid_summary.tsv" > "$workdir/07_eMGE_linkage/MGE_db/add_bin_plasmid_summary.tsv"; tail -n+2 "$workdir/07_eMGE_linkage/MGE_db_copy/individual_genomad/spades_output_modified_plasmid_summary.tsv" >> "$workdir/07_eMGE_linkage/MGE_db/add_bin_plasmid_summary.tsv"

pick_representatives \
  --input   "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_add_bin_plasmid.filter.tsv" \
  --summary "$workdir/07_eMGE_linkage/MGE_db/add_bin_plasmid_summary.tsv" \
  --output  "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_2.tsv"

# Here, only clusters containing scDNA elements are extracted, and clusters composed entirely of bin sequences are discarded; this is because I am not sure whether these long sequences truly exist or are chromosome fragments formed by binning clusters
# Note: At the same time, solely DTR MGEs from binning are also removed
grep 'scDNA' "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_2.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_2.filter.tsv"
cut -f 1 "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_2.filter.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_2.id"
grep 'round' "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_2.id" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_2.bin.id"

filterbyname.sh substring=name ow=t include=t minlen=0 names="$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_2.id" in="$workdir/07_eMGE_linkage/MGE_db/add_bin_plasmid.fna" out="$workdir/07_eMGE_linkage/MGE_db/rep-plasmid_2.fna"

cat "$workdir/07_eMGE_linkage/MGE_db/rep-plasmid_1_DTR.fna" "$workdir/07_eMGE_linkage/MGE_db/rep-plasmid_2.fna" > "$workdir/07_eMGE_linkage/MGE_db/rep-plasmid.fna"



# step 3 virus
cut -f 2 -d ' ' "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/final.reformat.add_cluster.edit_cont.annotate.species.txt" > "$workdir/07_eMGE_linkage/clean_SAGs.id"
filterbyname.sh in="$workdir/04_SAG_assembly/genomad_output/SAG_assembly_combine.500_summary/SAG_assembly_combine.500_virus.fna" \
                out="$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_virus_clean.fna" \
                names="$workdir/07_eMGE_linkage/clean_SAGs.id" include=t substring=t ow=t

makeblastdb -in "$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_virus_clean.fna" -dbtype nucl -out "$workdir/07_eMGE_linkage/MGE_db/blastn_virus_clean"
blastn -query "$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_virus_clean.fna" \
       -db "$workdir/07_eMGE_linkage/MGE_db/blastn_virus_clean" \
       -outfmt '6 std qlen slen' -max_target_seqs 10000 \
       -out "$workdir/07_eMGE_linkage/MGE_db/blastn_virus_clean.tsv" -num_threads 48
python "$script/anicalc.py" -i "$workdir/07_eMGE_linkage/MGE_db/blastn_virus_clean.tsv" -o "$workdir/07_eMGE_linkage/MGE_db/ani-blastn_virus_clean.tsv"
python "$script/aniclust.py" --fna "$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_virus_clean.fna" \
                             --ani "$workdir/07_eMGE_linkage/MGE_db/ani-blastn_virus_clean.tsv" \
                             --out "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_virus_clean.tsv" --min_ani 99 --min_tcov 90 --min_qcov 0 --min_length 1000

pick_representatives \
  --input   "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_virus_clean.tsv" \
  --summary "$workdir/07_eMGE_linkage/MGE_db_copy/individual_genomad/spades_output_modified_virus_summary.tsv" \
  --output  "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_1.tsv"

grep 'DTR' "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_1.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_1_DTR.tsv"
grep -v 'DTR' "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_1.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_1_no_DTR.tsv"

cut -f 1 "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_1_DTR.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_1_DTR.id"
cut -f 1 "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_1_no_DTR.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_1_no_DTR.id"

filterbyname.sh substring=name ow=t include=t minlen=0 \
                names="$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_1_DTR.id" \
                in="$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_virus_clean.fna" \
                out="$workdir/07_eMGE_linkage/MGE_db/rep-virus_1_DTR.fna"

filterbyname.sh substring=name ow=t include=t minlen=0 \
                names="$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_1_no_DTR.id" \
                in="$workdir/07_eMGE_linkage/MGE_db/spades_output_modified_virus_clean.fna" \
                out="$workdir/07_eMGE_linkage/MGE_db/rep-virus_1_no_DTR.fna"



# step4 virus : Cataloguing CoSAGs
cp "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/chromosome_combine.filtered_summary/chromosome_combine.filtered_virus.fna" "$workdir/07_eMGE_linkage/MGE_db/temp.fna"
cp "$workdir/07_eMGE_linkage/MGE_db/clean_bin_genomad/chromosome_combine.filtered_summary/chromosome_combine.filtered_virus_summary.tsv" "$workdir/07_eMGE_linkage/MGE_db/merged-well_annotated_genome_virus_summary.tsv"
awk '{print $1}' "$workdir/07_eMGE_linkage/final.reformat.add_cluster.edit_cont.annotate.species.txt" | sort | uniq > "$workdir/07_eMGE_linkage/MGE_db/remain.cluster"
filterbyname.sh substring=name ow=t include=t minlen=0 names="$workdir/07_eMGE_linkage/MGE_db/remain.cluster" \
                in="$workdir/07_eMGE_linkage/MGE_db/temp.fna" \
                out="$workdir/07_eMGE_linkage/MGE_db/merged-well_annotated_genome_virus.fna"

# DTR sequences are added here to avoid including those already present; this is only for assigning clusters to the no_DTR round of sequences
cat "$workdir/07_eMGE_linkage/MGE_db/rep-virus_1_DTR.fna" "$workdir/07_eMGE_linkage/MGE_db/rep-virus_1_no_DTR.fna" "$workdir/07_eMGE_linkage/MGE_db/merged-well_annotated_genome_virus.fna" > "$workdir/07_eMGE_linkage/MGE_db/add_bin_virus.fna"
makeblastdb -in "$workdir/07_eMGE_linkage/MGE_db/add_bin_virus.fna" -dbtype nucl -out "$workdir/07_eMGE_linkage/MGE_db/add_bin_virus"
blastn -query "$workdir/07_eMGE_linkage/MGE_db/add_bin_virus.fna" -db "$workdir/07_eMGE_linkage/MGE_db/add_bin_virus" -outfmt '6 std qlen slen' -max_target_seqs 10000 -out "$workdir/07_eMGE_linkage/MGE_db/blastn_add_bin_virus.tsv" -num_threads 48
python "$script/anicalc.py" -i "$workdir/07_eMGE_linkage/MGE_db/blastn_add_bin_virus.tsv" -o "$workdir/07_eMGE_linkage/MGE_db/ani-blastn_add_bin_virus.tsv"
python "$script/aniclust.py" --fna "$workdir/07_eMGE_linkage/MGE_db/add_bin_virus.fna" --ani "$workdir/07_eMGE_linkage/MGE_db/ani-blastn_add_bin_virus.tsv" --out "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_add_bin_virus.tsv" --min_ani 99 --min_tcov 90 --min_qcov 0 --min_length 1000
grep '>' "$workdir/07_eMGE_linkage/MGE_db/rep-virus_1_DTR.fna" | sed 's/>//' > "$workdir/07_eMGE_linkage/MGE_db/temp"
grep -F -v -f "$workdir/07_eMGE_linkage/MGE_db/temp" "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_add_bin_virus.tsv" > "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_add_bin_virus.filter.tsv"

cat "$workdir/07_eMGE_linkage/MGE_db/merged-well_annotated_genome_virus_summary.tsv" > "$workdir/07_eMGE_linkage/MGE_db/add_bin_virus_summary.tsv"; tail -n+2 "$workdir/07_eMGE_linkage/MGE_db_copy/individual_genomad/spades_output_modified_virus_summary.tsv" >> "$workdir/07_eMGE_linkage/MGE_db/add_bin_virus_summary.tsv"

pick_representatives \
  --input   "$workdir/07_eMGE_linkage/MGE_db/cluster-blastn_add_bin_virus.filter.tsv" \
  --summary "$workdir/07_eMGE_linkage/MGE_db/add_bin_virus_summary.tsv" \
  --output  "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_2.tsv"

# Here, only clusters containing scDNA elements are extracted, and clusters composed entirely of bin sequences are discarded; this is because I am not sure whether these long sequences truly exist or are chromosome fragments formed by binning clusters
# Note: At the same time, solely DTR MGEs from binning are also removed
grep 'scDNA' "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_2.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_2.filter.tsv"
cut -f 1 "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_2.filter.tsv" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_2.id"
grep 'round' "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_2.id" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_2.bin.id"

filterbyname.sh substring=name ow=t include=t minlen=0 names="$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_2.id" in="$workdir/07_eMGE_linkage/MGE_db/add_bin_virus.fna" out="$workdir/07_eMGE_linkage/MGE_db/rep-virus_2.fna"

cat "$workdir/07_eMGE_linkage/MGE_db/rep-virus_1_DTR.fna" "$workdir/07_eMGE_linkage/MGE_db/rep-virus_2.fna" > "$workdir/07_eMGE_linkage/MGE_db/rep-virus.fna"

# For virus, separate phage and prophage
grep 'provirus' "$workdir/07_eMGE_linkage/MGE_db/rep-virus.fna" | sed 's/>//' > "$workdir/07_eMGE_linkage/MGE_db/provirus.id"
filterbyname.sh substring=name ow=t include=t minlen=0 names="$workdir/07_eMGE_linkage/MGE_db/provirus.id" in="$workdir/07_eMGE_linkage/MGE_db/rep-virus.fna" out="$workdir/07_eMGE_linkage/MGE_db/rep-prophage.fna"
filterbyname.sh substring=name ow=t include=f minlen=0 names="$workdir/07_eMGE_linkage/MGE_db/provirus.id" in="$workdir/07_eMGE_linkage/MGE_db/rep-virus.fna" out="$workdir/07_eMGE_linkage/MGE_db/rep-phage.fna"

cat "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_plasmid_2.bin.id" "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_virus_2.bin.id" > "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_bin.id"