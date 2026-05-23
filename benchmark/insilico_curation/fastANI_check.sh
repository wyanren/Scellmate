for d in genome_dir_*; do
    [ -d "$d" ] || continue

    seed="${d#genome_dir_}"
    meta="selected_4_genus_triplets.no_placeholder__${seed}.tsv"
    list="./.genome_list_${seed}.txt"
    raw="./.fastani_${seed}.tsv"
    species_tsv="./.fastani_${seed}.species.tsv"
    out="./fastANI_matrix_species_${seed}.tsv"

    if [ -s "$out" ]; then
        echo "[skip] $d : $out already exists"
        continue
    fi

    if [ ! -f "$meta" ]; then
        echo "[skip] $d : missing $meta"
        continue
    fi

    echo "[run] $d"

    find "$d" -maxdepth 1 -type f -name "*.fna" | sort > "$list"

    if [ ! -s "$list" ]; then
        echo "[skip] $d : no .fna found"
        rm -f "$list"
        continue
    fi

    fastANI \
        --ql "$list" \
        --rl "$list" \
        -o "$raw" \
        -t 30

    if [ $? -ne 0 ] || [ ! -s "$raw" ]; then
        echo "[failed] fastANI $d"
        rm -f "$list" "$raw" "$species_tsv"
        continue
    fi

    awk 'BEGIN{FS=OFS="\t"}
         FNR==NR{
             id=$4
             sp=$5
             gsub(/^.*\//,"",id)
             gsub(/_genomic\.fna(\.gz)?$/,"",id)
             id2sp[id]=sp
             next
         }
         {
             q=$1; r=$2
             gsub(/^.*\//,"",q); gsub(/^.*\//,"",r)
             gsub(/_genomic\.fna$/,"",q); gsub(/_genomic\.fna$/,"",r)

             qn = (q in id2sp ? id2sp[q] : q)
             rn = (r in id2sp ? id2sp[r] : r)

             print qn, rn, $3, $4, $5
         }' "$meta" "$raw" > "$species_tsv"

    python3 - "$species_tsv" "$out" <<'PY'
import sys
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(
    infile,
    sep="\t",
    header=None,
    names=["query", "ref", "ANI", "frags_matched", "frags_total"]
)

species = sorted(set(df["query"]).union(df["ref"]))
mat = pd.DataFrame(index=species, columns=species, dtype=float)

for _, row in df.iterrows():
    mat.loc[row["query"], row["ref"]] = row["ANI"]

for s in species:
    mat.loc[s, s] = 100.0

for i in species:
    for j in species:
        if pd.isna(mat.loc[i, j]) and not pd.isna(mat.loc[j, i]):
            mat.loc[i, j] = mat.loc[j, i]

mat.to_csv(outfile, sep="\t", na_rep="NA")
PY

    if [ $? -ne 0 ] || [ ! -s "$out" ]; then
        echo "[failed] matrix build $d"
        rm -f "$out"
    fi

    rm -f "$list" "$raw" "$species_tsv"
done
