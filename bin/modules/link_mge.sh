#!/bin/bash
set -euo pipefail

# Default value for threads
threads=48 # here is the maximum threads that the system can use

# Parameter parsing
while [[ $# -gt 0 ]]; do
  case $1 in
    -w|--workdir)
      workdir=$(realpath "$2")
      shift 2
      ;;
    -s|--script)
      script=$(realpath "$2")
      shift 2
      ;;
    -t|--threads)
      threads="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 -w <workdir> -s <script_dir> [-t <threads>]"
      exit 0
      ;;
    *)
      echo "Unknown parameter: $1"
      echo "Usage: $0 -w <workdir> -s <script_dir> [-t <threads>]"
      exit 1
      ;;
  esac
done

# Check required parameters
if [[ -z "$workdir" || -z "$script" ]]; then
  echo "Error: --workdir and --script are required"
  echo "Usage: $0 -w <workdir> -s <script_dir> [-t <threads>]"
  exit 1
fi

# After parameter parsing and required parameter check
process_threads=$((threads / 9))
if [ "$process_threads" -lt 1 ]; then process_threads=1; fi

# Define log file path
log_file="$workdir/07_eMGE_linkage/link_mge.log"

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

########################################################################################################################
########################################################################################################################
########################################################################################################################
# step 1: preparation work
if step_done "step 1: preparation work"; then
  echo "Skipping step 1: preparation work"
else
  echo "Start step 1: preparation work"
  mkdir -p "$workdir/07_eMGE_linkage/" 2>/dev/null || true
  cd "$workdir/07_eMGE_linkage/"

mkdir -p "$workdir/07_eMGE_linkage/chromosome/masking_id" 2>/dev/null || true

shopt -s nullglob
for f in "$workdir/07_eMGE_linkage/chromosome/"*.fasta; do
    touch "$workdir/07_eMGE_linkage/chromosome/masking_id/$(basename "${f%.fasta}").filter.id"
done
shopt -u nullglob

# Split the representative_sequences_bin.id file and remove prefixes
awk -F'__' '{print $2 >> "chromosome/masking_id/"$1".filter.id"}' "$workdir/07_eMGE_linkage/MGE_db/representative_sequences_bin.id"

cp "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/well_annotated_genome/"*.fasta "$workdir/07_eMGE_linkage/chromosome/"
mkdir -p "$workdir/07_eMGE_linkage/chromosome/chromosome_masking/" 2>/dev/null || true
parallel -j "$threads" '
    sample=$(basename "{}" .fasta)
    echo "Handling ${sample}"
    filterbyname.sh substring=name ow=t include=f names="chromosome/masking_id/${sample}.filter.id" in="chromosome/${sample}.fasta" out="chromosome/chromosome_masking/${sample}.fasta"
' ::: chromosome/*.fasta


cp "$script/work_generate_masking.py" "$workdir/07_eMGE_linkage/"
cp "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/species_to_nocluster_deduplicated.tsv" "$workdir/07_eMGE_linkage/"
cp "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/final.reformat.add_cluster.edit_cont.annotate.species.nocluster.json" "$workdir/07_eMGE_linkage/"
cp "$workdir/06_CoSAG_assembly/round_ending/round_all/combined_full/final.reformat.add_cluster.edit_cont.annotate.species.json" "$workdir/07_eMGE_linkage/"
cd "$workdir/07_eMGE_linkage/"
python work_generate_masking.py
log_step "step 1: preparation work"
fi

# step 2: construct reference and perform mapping
if step_done "step 2: construct reference and perform mapping"; then
  echo "Skipping step 2: construct reference and perform mapping"
else
  echo "Start step 2: construct reference and perform mapping"
  # construct reference and perform mapping
  process_dir_mapping() { 
    dir="$1"
    if [ -d "$dir" ]; then
        echo "Handling directory: $dir"
        cd "$dir" || exit  # enter directory

        shopt -s nullglob
        cat chromosome/*.fasta "../..//MGE_db/rep-plasmid.fna" "../../MGE_db/rep-phage.fna" > reference.fasta
        shopt -u nullglob

        bowtie2-build --threads 8 reference.fasta reference

        mkdir -p mapping 2>/dev/null || true

        fastq_files=(fastq/*_R2_paired.fastq)
        parallel -j 3 '
            sample=$(basename "{}" _R2_paired.fastq)
            echo "Handling ${sample}"
            bowtie2 -p 3 -x ./reference --maxins 5000 \
                -1 "./fastq/${sample}_R1_paired.fastq" \
                -2 "./fastq/${sample}_R2_paired.fastq" \
                -S "mapping/${sample}.sam"
            samtools flagstat "mapping/${sample}.sam" > "mapping/${sample}.flagstat"
            samtools view -@ 1 -bS "mapping/${sample}.sam" > "mapping/${sample}.bam"
            samtools sort -@ 1 "mapping/${sample}.bam" -o "mapping/${sample}.sort.bam"
            samtools index "mapping/${sample}.sort.bam"
            samtools idxstats "mapping/${sample}.sort.bam" > "mapping/${sample}.raw.idxstats"
            samtools depth "mapping/${sample}.sort.bam" > "mapping/${sample}.raw.depth"
            rm -f "mapping/${sample}.sam" "mapping/${sample}.bam"
        ' ::: "${fastq_files[@]}"

        cd - || exit  # back to previous directory
    else
        echo "Warning: $dir is not a directory, skip"
    fi
}

  export -f process_dir_mapping
  parallel -j "$process_threads" process_dir_mapping ::: "$workdir/07_eMGE_linkage/work/"*
  log_step "step 2: construct reference and perform mapping"
fi

# step 3: calculate depth and breadth
if step_done "step 3: calculate depth and breadth"; then
  echo "Skipping step 3: calculate depth and breadth"
else
  echo "Start step 3: calculate depth and breadth"
  # calculate depth
  process_dir_depth() {
    dir="$1"
    if [ -d "$dir" ]; then
        echo "Handling directory: $dir"
        cd "$dir" || exit  # enter directory

        fastq_files=(fastq/*_R2_paired.fastq)
        parallel -j 6 '
            sample=$(basename "{}" _R2_paired.fastq)
            echo "Handling ${sample}"
            samtools view -b -f 2 -F 524 "mapping/${sample}.sort.bam" > "mapping/${sample}.unique.bam"
            samtools index "mapping/${sample}.unique.bam"
            samtools idxstats "mapping/${sample}.unique.bam" > "mapping/${sample}.idxstats"
            samtools depth "mapping/${sample}.unique.bam" > "mapping/${sample}.depth"
        ' ::: "${fastq_files[@]}"

        cd - || exit  # back to previous directory
    else
        echo "Warning: $dir is not a directory, skip"
    fi
}

  export -f process_dir_depth
  parallel -j "$process_threads" process_dir_depth ::: "$workdir/07_eMGE_linkage/work/"*


# calculate breadth
cp "$script/calculate_breadth.py" "$workdir/07_eMGE_linkage/"
cp "$script/combine_depth.py" "$workdir/07_eMGE_linkage/"

process_dir_breadth() {
    dir="$1"
    if [ -d "$dir" ]; then
        echo "Handling directory: $dir"
        cd "$dir" || exit  # enter directory

        fastq_files=(fastq/*_R2_paired.fastq)
        parallel -j 6 '
            sample=$(basename "{}" _R2_paired.fastq)
            echo "Handling ${sample}"
            python "../../calculate_breadth.py" "mapping/${sample}.depth" > "mapping/${sample}.breadth"
            python "../../calculate_breadth.py" "mapping/${sample}.raw.depth" > "mapping/${sample}.raw.breadth"
            mv "mapping/${sample}.raw.breadth" "mapping/${sample}.raw_breadth"
        ' ::: "${fastq_files[@]}"

        cd - || exit  # back to previous directory
    else
        echo "Warning: $dir is not a directory, skip"
    fi
}

  export -f process_dir_breadth
  parallel -j "$process_threads" process_dir_breadth ::: "$workdir/07_eMGE_linkage/work/"*

# combine depth and calculate breadth
process_dir_combine() {
    dir="$1"
    if [ -d "$dir" ]; then
        echo "Handling directory: $dir"
        cd "$dir" || exit  # enter directory
        python "../../combine_depth.py" -i mapping/scDNA_*_?????.depth -o temp
        awk '{print $0 "\t" 1}' temp > combined.depth
        python "../../calculate_breadth.py" combined.depth > combined.breadth
        rm temp
        cd - || exit  # back to previous directory
    else
        echo "Warning: $dir is not a directory, skip"
    fi
}

  export -f process_dir_combine
  parallel -j "$threads" process_dir_combine ::: "$workdir/07_eMGE_linkage/work/"*
  log_step "step 3: calculate depth and breadth"
fi

# step 4: generate presence absence table
if step_done "step 4: generate presence absence table"; then
  echo "Skipping step 4: generate presence absence table"
else
  echo "Start step 4: generate presence absence table"
  # generate presence absence table
  process_dir_presence() {
    dir="$1"
    if [ -d "$dir" ]; then
        echo "Handling directory: $dir"
        cd "$dir" || exit  # enter directory
        shopt -s nullglob
        local breadth_files=(mapping/scDNA_*_?????.breadth)
        shopt -u nullglob
        breadth_files=(mapping/scDNA_*_?????.breadth)
        ((${#breadth_files[@]})) || { echo "  - No breadth files, skip"; cd - >/dev/null; return; }
        awk '$1 ~ /^(scDNA|round|cluster)/ && $4 > 0.90' combined.breadth | cut -f 1 | sort | uniq > species_1st.id
        python "../../generate_presence_absence_table_new.py" \
            --coverage_rate_threshold 0.75 --num_covered_nt_threshold 850 \
            -i "${breadth_files[@]}" \
            -s species_1st.id \
            -o presence_absence_table.txt
        cd - || exit  # back to previous directory
    else
        echo "Warning: $dir is not a directory, skip"
    fi
}
  cp "$script/generate_presence_absence_table_new.py" "$workdir/07_eMGE_linkage/"
  export -f process_dir_presence
  parallel -j "$threads" process_dir_presence ::: "$workdir/07_eMGE_linkage/work/"*
  log_step "step 4: generate presence absence table"
fi

# step 5: sum presence absence table
if step_done "step 5: sum presence absence table"; then
  echo "Skipping step 5: sum presence absence table"
else
  echo "Start step 5: sum presence absence table"
  # sum presence absence table
  for dir in "$workdir/07_eMGE_linkage/work/"*; do
    if [ -d "$dir" ]; then
        input_file="$dir/presence_absence_table.txt"
        output_file="$dir/summed_table.txt"
        if [ -f "$input_file" ]; then
            echo "Processing $input_file"
            python "$script/sum_presence_absence.py" -i "$input_file" -o "$output_file"
        else
            echo "Warning: $input_file does not exist"
        fi
    fi
  done
log_step "step 5: sum presence absence table"
fi

# step 6: merge summed tables
if step_done "step 6: merge summed tables"; then
  echo "Skipping step 6: merge summed tables"
else
  echo "Start step 6: merge summed tables"
  # merge summed tables
  pattern="$workdir/07_eMGE_linkage/work/*/summed_table.txt"
  python "$script/merge_summed_tables.py" \
      -i "$pattern" \
      -o "$workdir/07_eMGE_linkage/work/merged_presence_absence_table.txt"
    
  # calculate sum
  awk 'BEGIN {OFS="\t"} 
  NR==1 {print $0, "Sum"} 
  NR>1 {sum=0; for(i=2;i<=NF;i++) sum += $i; print $0, sum}' "$workdir/07_eMGE_linkage/work/merged_presence_absence_table.txt" > "$workdir/07_eMGE_linkage/work/summed_table.txt"
log_step "step 6: merge summed tables"
fi
