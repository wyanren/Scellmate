# ====================================================================================================
# 0. Preparation
# ====================================================================================================
# Illumina dataset from isolates
fastq/
├── E_coli_113_1.fastq
├── E_coli_113_2.fastq
├── E_coli_114_1.fastq
├── E_coli_114_2.fastq
├── E_coli_115_1.fastq
├── E_coli_115_2.fastq
├── K_pneumoniae_1.fastq
├── K_pneumoniae_2.fastq
├── S_aureus_1.fastq
└── S_aureus_2.fastq
# NCBI plasmid assembly from isolates
NCBI_assembly/
├── CP015135.1__K_pneumoniae_ATCC_35657_plasmid1.fasta
├── CP034485.1__B_subtilis_ATCC_6051_plasmid1.fasta
├── E_coli_ATCC_25922_2plasmids.fasta
└── NZ_CP158285.1__S_aureus_ATCC_6538_plasmid1.fasta

# ====================================================================================================
# 1. geNomad from SPAdes-derived assemblies of publicly available Illumina sequencing
# ====================================================================================================
workdir="/mnt/md0/wangyanren/workdir/202604_mock/mock_illumina"

mkdir -p "$workdir/isolate_assembly/spades_work" \
         "$workdir/isolate_assembly/spades_output" \
         "$workdir/isolate_assembly/logs"

# Running SPAdes for contig assembly
for r1 in "$workdir"/fastq/*_1.fastq; do
    sample=$(basename "$r1" _1.fastq)
    r2="$workdir/fastq/${sample}_2.fastq"
    outdir="$workdir/isolate_assembly/spades_work/$sample"
    echo "[INFO] Running SPAdes for $sample"
    rm -rf "$outdir"
    spades.py \
        -t 60 \
        --phred-offset 33 \
        --pe1-1 "$r1" \
        --pe1-2 "$r2" \
        -o "$outdir" \
        > "$workdir/isolate_assembly/logs/${sample}.spades.log" 2>&1
    if [[ -s "$outdir/contigs.fasta" ]]; then
        cp "$outdir/contigs.fasta" "$workdir/isolate_assembly/spades_output/${sample}.fasta"
        echo "[INFO] Finished $sample"
    else
        echo "[WARNING] No contigs.fasta for $sample"
    fi
done

# Combine and renamed the contigs
combine="$workdir/isolate_assembly/isolate_assembly_combine.fasta"
> "$combine"

for fasta in "$workdir"/isolate_assembly/spades_output/*.fasta; do
    base=$(basename "$fasta" .fasta)
    awk -v pre="$base" '{
        if(/^>/){sub(/^>/,""); print ">"pre"_"$0}
        else print
    }' "$fasta" >> "$combine"
done

# Length filteration
min_contig=1000
reformat.sh \
    in="$workdir/isolate_assembly/isolate_assembly_combine.fasta" \
    out="$workdir/isolate_assembly/isolate_assembly_combine.${min_contig}.fasta" \
    minlength="$min_contig" \
    overwrite=t \
    > /dev/null 2>&1

# geNomad predict eMGE-like contigs
mkdir -p "$workdir/isolate_assembly/genomad_output"
genomad end-to-end --cleanup --splits 8 \
    "$workdir/isolate_assembly/isolate_assembly_combine.${min_contig}.fasta" \
    "$workdir/isolate_assembly/genomad_output" \
    "$SCELLMATE_DB/genomad_db" \
    > "$workdir/isolate_assembly/logs/genomad.log" 2>&1


# ====================================================================================================
# 2. plasmid catalogue construction, run separately for E. coli, K. pneumoniae, and S. aureus
# ====================================================================================================
cd /mnt/md0/wangyanren/workdir/202604_mock/mock_illumina

# define cataloging function
pick_rep_ncbi_first() {
  input="$1"
  summary="$2"
  output="$3"

  awk -F'\t' -v OFS='\t' '
    function is_ncbi(x) {
      return (x ~ /^(CP|NZ_CP)/)
    }
    function tag_rep(r) {
      return (r ~ /No terminal/ ? "No_terminal_repeats" : r)
    }

    NR==FNR { info[$1]=$0; next }

    {
      split($2, seq, /,/)
      n=length(seq)

      # 1) if any NCBI plasmid exists in this cluster, directly choose it
      ncbi_n=0
      for (i=1; i<=n; i++) {
        if (is_ncbi(seq[i])) {
          ncbi[++ncbi_n]=seq[i]
        }
      }

      if (ncbi_n >= 1) {
        # if multiple NCBI refs exist in same cluster, pick the longest
        pick=ncbi[1]
        for (i=2; i<=ncbi_n; i++) {
          a=ncbi[i]
          if ((a in info) && (pick in info)) {
            split(info[a], fa)
            split(info[pick], fp)
            if ((fa[2]+0) > (fp[2]+0)) pick=a
          }
        }
        print pick, $2, "NCBI_reference"
        delete ncbi
        next
      }
      delete ncbi

      # 2) otherwise, fall back to previous local logic
      if (n == 1) {
        if (seq[1] in info) {
          split(info[seq[1]], f)
          print seq[1], $2, "Solely_" tag_rep(f[3])
        }
        next
      }

      delete id; delete len; delete rep; delete dtr; delete cov
      valid=1
      for (i = 1; i <= n; i++) {
        k = seq[i]
        if (!(k in info)) { valid=0; break }
        split(info[k], f)
        id[i]=f[1]
        len[i]=f[2]+0
        rep[i]=f[3]
        dtr[i]=f[4]
        if (match(f[1], /_cov_([0-9.]+)/, m)) cov[i]=m[1]+0
        else cov[i]=0
      }
      if (!valid) next

      hasDTR = 0
      for (i in rep) if (rep[i]=="DTR") hasDTR=1

      if (!hasDTR) {
        pick=1
        for (i in len) if (len[i] > len[pick]) pick=i
        print id[pick], $2, "No_terminal_repeats"
        next
      }

      delete d; nd=0; delete freq
      for (i in rep) if (rep[i]=="DTR") {
        d[++nd]=i
        freq[dtr[i]]++
      }

      if (nd==1) {
        print id[d[1]], $2, "DTR_only"
        next
      }

      maxf=0
      for (l in freq) if (freq[l]>maxf) maxf=freq[l]

      delete cand; nc=0
      for (j=1; j<=nd; j++) {
        i=d[j]
        if (freq[dtr[i]]==maxf) cand[++nc]=i
      }

      if (nc==1) {
        print id[cand[1]], $2, "DTR_Most_Frequent"
        next
      }

      mcov=-1
      delete best; nb=0
      for (j=1; j<=nc; j++) {
        i=cand[j]
        if (cov[i]>mcov) mcov=cov[i]
      }
      for (j=1; j<=nc; j++) {
        i=cand[j]
        if (cov[i]==mcov) best[++nb]=i
      }

      if (nb==1) {
        print id[best[1]], $2, "DTR_Most_Frequent"
        next
      }

      pick=best[1]
      for (j=2; j<=nb; j++) if (len[best[j]]>len[pick]) pick=best[j]
      print id[pick], $2, "DTR_Longest_Most_Frequent"
    }
  ' "$summary" "$input" > "$output"
}


# 1）E. coli cataloging
mkdir -p strain_rep_pick/E_coli
cd strain_rep_pick/E_coli

conda activate microbe-seq
seqkit grep -r -p '^CP009073\.1__E_coli_ATCC_25922_plasmid_1$|^CP009074\.1__E_coli_ATCC_25922_plasmid_2$' \
  ../../NCBI_assembly/E_coli_ATCC_25922_2plasmids.fasta \
  > E_coli_ncbi.fna

seqkit grep -r -p '^E_coli_.*NODE_' \
  ../../isolate_assembly/genomad_output/isolate_assembly_combine.1000_summary/isolate_assembly_combine.1000_plasmid.fna \
  > E_coli_local.fna

cat E_coli_ncbi.fna E_coli_local.fna > E_coli_all.fna

conda activate Scellmate_env3
awk -F'\t' 'NR==1 || $1 ~ /^E_coli_.*NODE_/' \
  ../../isolate_assembly/genomad_output/isolate_assembly_combine.1000_summary/isolate_assembly_combine.1000_plasmid_summary.tsv \
  > E_coli_local_summary.tsv

{
  echo -e "seq_name\tlength\tterminal_repeat_type\tterminal_repeat_length"
  seqkit fx2tab -n -l E_coli_ncbi.fna | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,"No terminal repeats","NA"}'
} > E_coli_ncbi_summary.tsv

cat E_coli_ncbi_summary.tsv <(tail -n +2 E_coli_local_summary.tsv) > E_coli_summary.tsv

makeblastdb -in E_coli_all.fna -dbtype nucl -out E_coli_db

blastn \
  -query E_coli_all.fna \
  -db E_coli_db \
  -outfmt '6 std qlen slen' \
  -max_target_seqs 10000 \
  -out E_coli_blastn.tsv \
  -num_threads 20

python ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_scripts/anicalc.py \
  -i E_coli_blastn.tsv \
  -o E_coli_ani.tsv

python ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_scripts/aniclust.py \
  --fna E_coli_all.fna \
  --ani E_coli_ani.tsv \
  --out E_coli_cluster.tsv \
  --min_ani 99 \
  --min_tcov 80 \
  --min_qcov 0 \
  --min_length 1000

pick_rep_ncbi_first E_coli_cluster.tsv E_coli_summary.tsv E_coli_representatives.tsv
cut -f 1 E_coli_representatives.tsv > E_coli_representatives.id
filterbyname.sh substring=name ow=t include=t minlen=0 \
  names=E_coli_representatives.id \
  in=E_coli_all.fna \
  out=rep-E_coli.fna > /dev/null 2>&1


# 2）K. pneumoniae cataloging
cd ..
mkdir -p K_pneumoniae
cd K_pneumoniae

cp ../../NCBI_assembly/CP015135.1__K_pneumoniae_ATCC_35657_plasmid1.fasta K_pneumoniae_ncbi.fna
conda activate microbe-seq
seqkit grep -r -p '^K_pneumoniae_.*NODE_' \
  ../../isolate_assembly/genomad_output/isolate_assembly_combine.1000_summary/isolate_assembly_combine.1000_plasmid.fna \
  > K_pneumoniae_local.fna

cat K_pneumoniae_ncbi.fna K_pneumoniae_local.fna > K_pneumoniae_all.fna

awk -F'\t' 'NR==1 || $1 ~ /^K_pneumoniae_.*NODE_/' \
  ../../isolate_assembly/genomad_output/isolate_assembly_combine.1000_summary/isolate_assembly_combine.1000_plasmid_summary.tsv \
  > K_pneumoniae_local_summary.tsv

{
  echo -e "seq_name\tlength\tterminal_repeat_type\tterminal_repeat_length"
  seqkit fx2tab -n -l K_pneumoniae_ncbi.fna | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,"No terminal repeats","NA"}'
} > K_pneumoniae_ncbi_summary.tsv

cat K_pneumoniae_ncbi_summary.tsv <(tail -n +2 K_pneumoniae_local_summary.tsv) > K_pneumoniae_summary.tsv

conda activate Scellmate_env3

makeblastdb -in K_pneumoniae_all.fna -dbtype nucl -out K_pneumoniae_db

blastn \
  -query K_pneumoniae_all.fna \
  -db K_pneumoniae_db \
  -outfmt '6 std qlen slen' \
  -max_target_seqs 10000 \
  -out K_pneumoniae_blastn.tsv \
  -num_threads 20

python ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_scripts/anicalc.py \
  -i K_pneumoniae_blastn.tsv \
  -o K_pneumoniae_ani.tsv

python ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_scripts/aniclust.py \
  --fna K_pneumoniae_all.fna \
  --ani K_pneumoniae_ani.tsv \
  --out K_pneumoniae_cluster.tsv \
  --min_ani 99 \
  --min_tcov 80 \
  --min_qcov 0 \
  --min_length 1000

pick_rep_ncbi_first K_pneumoniae_cluster.tsv K_pneumoniae_summary.tsv K_pneumoniae_representatives.tsv
cut -f 1 K_pneumoniae_representatives.tsv > K_pneumoniae_representatives.id
filterbyname.sh substring=name ow=t include=t minlen=0 \
  names=K_pneumoniae_representatives.id \
  in=K_pneumoniae_all.fna \
  out=rep-K_pneumoniae.fna > /dev/null 2>&1

# 3）S. aureus cataloging
cd ..
mkdir -p S_aureus
cd S_aureus

cp ../../NCBI_assembly/NZ_CP158285.1__S_aureus_ATCC_6538_plasmid1.fasta S_aureus_ncbi.fna
conda activate Scellmate_env3
seqkit grep -r -p '^S_aureus_.*NODE_' \
  ../../isolate_assembly/genomad_output/isolate_assembly_combine.1000_summary/isolate_assembly_combine.1000_plasmid.fna \
  > S_aureus_local.fna

cat S_aureus_ncbi.fna S_aureus_local.fna > S_aureus_all.fna

awk -F'\t' 'NR==1 || $1 ~ /^S_aureus_.*NODE_/' \
  ../../isolate_assembly/genomad_output/isolate_assembly_combine.1000_summary/isolate_assembly_combine.1000_plasmid_summary.tsv \
  > S_aureus_local_summary.tsv

{
  echo -e "seq_name\tlength\tterminal_repeat_type\tterminal_repeat_length"
  seqkit fx2tab -n -l S_aureus_ncbi.fna | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,"No terminal repeats","NA"}'
} > S_aureus_ncbi_summary.tsv

cat S_aureus_ncbi_summary.tsv <(tail -n +2 S_aureus_local_summary.tsv) > S_aureus_summary.tsv

conda activate Scellmate_env3
makeblastdb -in S_aureus_all.fna -dbtype nucl -out S_aureus_db

blastn \
  -query S_aureus_all.fna \
  -db S_aureus_db \
  -outfmt '6 std qlen slen' \
  -max_target_seqs 10000 \
  -out S_aureus_blastn.tsv \
  -num_threads 20

python ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_scripts/anicalc.py \
  -i S_aureus_blastn.tsv \
  -o S_aureus_ani.tsv

python ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_scripts/aniclust.py \
  --fna S_aureus_all.fna \
  --ani S_aureus_ani.tsv \
  --out S_aureus_cluster.tsv \
  --min_ani 99 \
  --min_tcov 80 \
  --min_qcov 0 \
  --min_length 1000

pick_rep_ncbi_first S_aureus_cluster.tsv S_aureus_summary.tsv S_aureus_representatives.tsv
cut -f 1 S_aureus_representatives.tsv > S_aureus_representatives.id
filterbyname.sh substring=name ow=t include=t minlen=0 \
  names=S_aureus_representatives.id \
  in=S_aureus_all.fna \
  out=rep-S_aureus.fna > /dev/null 2>&1


# ====================================================================================================
# 3. summary the plasmid catalogue
# ====================================================================================================
cd ..
cat E_coli/rep-E_coli.fna K_pneumoniae/rep-K_pneumoniae.fna S_aureus/rep-S_aureus.fna > rep-plasmid.mock3strain.fna

# check
for x in E_coli K_pneumoniae S_aureus; do
  echo "====== $x ======"
  cat $x/${x}_representatives.tsv
  echo
done

# add Bacillus subtilis NCBI assembly-derived plasmid
cp ../NCBI_assembly/CP034485.1__B_subtilis_ATCC_6051_plasmid1.fasta B_subtilis_only_ref.fna
cat rep-plasmid.mock3strain.fna B_subtilis_only_ref.fna > rep-plasmid.mock4strain_with_Bsubtilis_ref.fna



# ====================================================================================================
# 2. phage catalogue construction, run separately for E. coli, K. pneumoniae, and S. aureus
# ====================================================================================================
# Only found phage-predicted sequences from E coli
cd /mnt/md0/wangyanren/workdir/202604_mock/mock_illumina/strain_rep_pick
mkdir -p virus_E_coli_local
cd virus_E_coli_local

# Remove prophage
conda activate microbe-seq
seqkit grep -r -p '^E_coli_.*NODE_' \
  ../../isolate_assembly/genomad_output/isolate_assembly_combine.1000_summary/isolate_assembly_combine.1000_virus.fna \
  | seqkit grep -v -r -p 'provirus' \
  > E_coli_nonprovirus_virus.fna

# Summary the genomad topology information
awk -F'\t' 'NR==1 || ($1 ~ /^E_coli_.*NODE_/ && $1 !~ /provirus/)' \
  ../../isolate_assembly/genomad_output/isolate_assembly_combine.1000_summary/isolate_assembly_combine.1000_virus_summary.tsv \
  > E_coli_nonprovirus_virus_summary.tsv

# Start cataloging
conda activate Scellmate_env3
makeblastdb -in E_coli_nonprovirus_virus.fna -dbtype nucl -out E_coli_nonprovirus_virus_db

# 
pick_representatives() {
  input="$1"
  summary="$2"
  output="$3"

  awk -F'\t' -v OFS='\t' '
    function tag_rep(r) { return (r ~ /No terminal/ ? "No_terminal_repeats" : r) }

    NR==FNR { info[$1]=$0; next }

    {
      split($2, seq, /,/)
      n=length(seq)

      if (n == 1) {
        if (seq[1] in info) {
          split(info[seq[1]], f)
          print seq[1], $2, "Solely_" tag_rep(f[3])
        }
        next
      }

      delete id; delete len; delete rep; delete dtr; delete cov
      valid=1
      for (i = 1; i <= n; i++) {
        k = seq[i]
        if (!(k in info)) { valid=0; break }
        split(info[k], f)
        id[i]=f[1]
        len[i]=f[2]+0
        rep[i]=f[3]
        dtr[i]=f[4]
        if (match(f[1], /_cov_([0-9.]+)/, m)) cov[i]=m[1]+0
        else cov[i]=0
      }
      if (!valid) next

      hasDTR = 0
      for (i in rep) if (rep[i]=="DTR") hasDTR=1

      if (!hasDTR) {
        pick=1
        for (i in len) if (len[i] > len[pick]) pick=i
        print id[pick], $2, "No_terminal_repeats"
        next
      }

      delete d; nd=0; delete freq
      for (i in rep) if (rep[i]=="DTR") {
        d[++nd]=i
        freq[dtr[i]]++
      }

      if (nd==1) {
        print id[d[1]], $2, "DTR_only"
        next
      }

      maxf=0
      for (l in freq) if (freq[l]>maxf) maxf=freq[l]

      delete cand; nc=0
      for (j=1; j<=nd; j++) {
        i=d[j]
        if (freq[dtr[i]]==maxf) cand[++nc]=i
      }

      if (nc==1) {
        print id[cand[1]], $2, "DTR_Most_Frequent"
        next
      }

      mcov=-1
      delete best; nb=0
      for (j=1; j<=nc; j++) {
        i=cand[j]
        if (cov[i]>mcov) mcov=cov[i]
      }
      for (j=1; j<=nc; j++) {
        i=cand[j]
        if (cov[i]==mcov) best[++nb]=i
      }

      if (nb==1) {
        print id[best[1]], $2, "DTR_Most_Frequent"
        next
      }

      pick=best[1]
      for (j=2; j<=nb; j++) if (len[best[j]]>len[pick]) pick=best[j]
      print id[pick], $2, "DTR_Longest_Most_Frequent"
    }
  ' "$summary" "$input" > "$output"
}

blastn \
  -query E_coli_nonprovirus_virus.fna \
  -db E_coli_nonprovirus_virus_db \
  -outfmt '6 std qlen slen' \
  -max_target_seqs 10000 \
  -out E_coli_nonprovirus_virus_blastn.tsv \
  -num_threads 20

python ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_scripts/anicalc.py \
  -i E_coli_nonprovirus_virus_blastn.tsv \
  -o E_coli_nonprovirus_virus_ani.tsv

python ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_scripts/aniclust.py \
  --fna E_coli_nonprovirus_virus.fna \
  --ani E_coli_nonprovirus_virus_ani.tsv \
  --out E_coli_nonprovirus_virus_cluster.tsv \
  --min_ani 99 \
  --min_tcov 90 \
  --min_qcov 0 \
  --min_length 1000

pick_representatives \
  E_coli_nonprovirus_virus_cluster.tsv \
  E_coli_nonprovirus_virus_summary.tsv \
  E_coli_nonprovirus_virus_representatives.tsv

cut -f 1 E_coli_nonprovirus_virus_representatives.tsv > E_coli_nonprovirus_virus_representatives.id

filterbyname.sh \
  substring=name \
  ow=t \
  include=t \
  minlen=0 \
  names=E_coli_nonprovirus_virus_representatives.id \
  in=E_coli_nonprovirus_virus.fna \
  out=rep-virus.E_coli.nonprovirus.fna \
  > /dev/null 2>&1

# Only found phage-predicted sequences from E coli
cp rep-virus.E_coli.nonprovirus.fna rep-phage.mock4strain.fna


# ====================================================================================================
# 3. summary the eMGE catalogue
# ====================================================================================================
cd /mnt/md0/wangyanren/workdir/202604_mock/mock_illumina/strain_rep_pick
cp rep-plasmid.mock4strain_with_Bsubtilis_ref.fna rep-plasmid.fna
cp virus_E_coli_local/rep-virus.E_coli.nonprovirus.fna ./rep-phage.fna

mv -f scellmate_mock/07_eMGE_linkage scellmate_mock/07_eMGE_linkage-test2
mkdir -p scellmate_mock/07_eMGE_linkage/MGE_db
cp /mnt/md0/wangyanren/workdir/202604_mock/mock_illumina/strain_rep_pick/rep-plasmid.fna scellmate_mock/07_eMGE_linkage/MGE_db
cp /mnt/md0/wangyanren/workdir/202604_mock/mock_illumina/strain_rep_pick/rep-phage.fna scellmate_mock/07_eMGE_linkage/MGE_db
cp scellmate_mock/07_eMGE_linkage-test2/MGE_db/representative_sequences_bin.id scellmate_mock/07_eMGE_linkage/MGE_db

for f in scellmate_mock/07_eMGE_linkage/MGE_db/rep-plasmid.fna \
         scellmate_mock/07_eMGE_linkage/MGE_db/rep-phage.fna
do
    cp "$f" "${f}.bak_$(date +%Y%m%d_%H%M%S)"
done

# Unify naming format
python3 - <<'PY'
from pathlib import Path
import re

files = [
    Path("scellmate_mock/07_eMGE_linkage/MGE_db/rep-plasmid.fna"),
    Path("scellmate_mock/07_eMGE_linkage/MGE_db/rep-phage.fna"),
]

for fp in files:
    lines = fp.read_text().splitlines()
    out = []

    header = None
    seq_chunks = []

    def flush_record(h, seqs):
        if h is None:
            return
        seq = "".join(seqs)
        h2 = h
        if not re.search(r'length_\d+', h):
            h2 = f"{h}_length_{len(seq)}"
        out.append(">" + h2)
        if seq:
            out.append(seq)

    for line in lines:
        if line.startswith(">"):
            flush_record(header, seq_chunks)
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())

    flush_record(header, seq_chunks)

    fp.write_text("\n".join(out) + "\n")
    print(f"[OK] updated {fp}")
PY


# ====================================================================================================
# 4. Running Scellmate curation for mock dataset
# ====================================================================================================
scellmate preprocess -i mock_dataset_SAG_reads -o scellmate_mock --prefix mock -t 48
scellmate first_qc --workdir scellmate_mock -t 48
scellmate second_qc --workdir scellmate_mock -t 48

# ====================================================================================================
# 5. Running Scellmate host–eMGE linking for mock dataset
# ====================================================================================================
bash ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_modules/eMGE_DB.sh --workdir scellmate_mock --script ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_scripts --threads 40

# replace the eMGE cataologue
cp mock_illumina/strain_rep_pick/rep-plasmid.fna scellmate_mock/07_eMGE_linkage/MGE_db/rep-plasmid.fna
cp mock_illumina/strain_rep_pick/rep-phage.fna scellmate_mock/07_eMGE_linkage/MGE_db/rep-phage.fna

bash ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_modules/link_mge.sh \
  -w scellmate_mock \
  -s ~/conda/envs/Scellmate_env3/share/scellmate/bin/SCM_scripts \
  -t 40 \
  --score_cutoff 0.01
