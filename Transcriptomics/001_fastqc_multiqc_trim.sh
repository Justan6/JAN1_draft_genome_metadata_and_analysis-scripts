#!/bin/bash

# Activate the RNA-seq QC environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rnaseq_qc

# Go to your transcriptomics directory
cd /proj/Peat_soil/PacBio_peat/transcriptomics || exit

# 1. Run FastQC on all raw reads
mkdir -p fastqc_output
for file in Data/*.fq.gz; do
    echo "Running FastQC on $file"
    fastqc "$file" -o fastqc_output/
done

# 2. Summarize FastQC results with MultiQC
mkdir -p multiqc_report
multiqc fastqc_output/ -o multiqc_report/
echo "QC complete. MultiQC report is in multiqc_report/"

# 3. Trim reads with fastp
mkdir -p trimmed

for r1_file in Data/*_1.fq.gz; do
    sample_name=$(basename "$r1_file" | sed 's/_1.fq.gz//')
    echo "Trimming $sample_name"

    r2_file="${r1_file/_1.fq.gz/_2.fq.gz}"

    if [[ -f "$r2_file" ]]; then
        fastp -i "$r1_file" -I "$r2_file" \
              --out1 "./trimmed/${sample_name}_1_trimmed.fq.gz" \
              --out2 "./trimmed/${sample_name}_2_trimmed.fq.gz" \
              --unpaired1 "./trimmed/${sample_name}_unpaired-1.fq.gz" \
              --unpaired2 "./trimmed/${sample_name}_unpaired-2.fq.gz" \
              --trim_front1 12 --trim_poly_g --detect_adapter_for_pe \
              --average_qual 30 --length_required 30 -w 4 \
              -h "./trimmed/${sample_name}_report.html" \
              -j "./trimmed/${sample_name}_report.json"

        echo "Finished trimming $sample_name"
    else
        echo "R2 file for $sample_name not found, skipping..."
    fi
done


