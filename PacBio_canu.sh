#!/bin/bash

# ========================================================================================
# Title          : PacBio_canu.sh
# Description    : Extraction of PacBio reads, conversion to FASTQ, and assembly using Canu and Flye
# Author         : Justus Nweze
# Date           : 21.08.2024
# Version Info   : Referenced from Canu and PacBio documentation
# GitHub Ready   : YES (Formatted and documented for public sharing)
# ========================================================================================

# =======================================
# SECTION 1: Extract PacBio .tar Archives
# =======================================
tar -xvf X201SC23121384-Z01-F001_01.tar  # For PacBio report
tar -xvf X201SC23121384-Z01-F001_02.tar  # For LM_Mbac_comb sample
tar -xvf X201SC23121384-Z01-F001_03.tar  # For LM_Msar_1108 sample

# Notes:
# - Extracted files include:
#   - .subreads.bam: Raw PacBio reads
#   - .subreads.bam.pbi: BAM index files

# =========================================================
# SECTION 2: Convert BAM to FASTQ Using PacBio BAM Toolkit
# =========================================================

# Installation notes:
# conda create --name pbtk
# conda install -c bioconda pbtk

conda activate pbtk
cd ~/scratch/PacBio_peat/Data/BAM

# Copy required BAM and PBI files from project directory
cp ~/proj/Peat_soil/PacBio_peat/Data/Peat/*.bam ~/scratch/PacBio_peat/Data/BAM
cp ~/proj/Peat_soil/PacBio_peat/Data/Peat/*.bam.pbi ~/scratch/PacBio_peat/Data/BAM

# Batch conversion of BAM files to FASTQ
for bam_file in *.bam; do
    base_name=$(basename "$bam_file" .bam)
    bam2fastq -o "${base_name}" "${bam_file}"
done

# =====================================================
# SECTION 3: Assembly with Canu
# =====================================================

# Install Canu (if not yet installed)
# conda create --name canu
# conda activate canu
# conda install -c bioconda/label/cf201901 canu

mkdir -p Canu_output

# Example Canu command for PacBio reads
# Replace input file and genome size as needed
canu -p my_assembly -d ./Canu_output genomeSize=5m -pacbio-corrected BAM/Peat73.hifi_reads.fastq.gz

# Notes:
# - Canu runs in three stages: correction → trimming → assembly
# - Output directories: `correction/`, `trimming/`, `unitigging/`

# Reference:
# https://canu.readthedocs.io/en/latest/quick-start.html

# =====================================================
# SECTION 4: Assembly with Flye
# =====================================================

# Flye is suitable for both standard and metagenomic assemblies
# Install: conda install -c bioconda flye

cd ~/proj/Peat_soil/PacBio_peat/Data
conda activate flye

# Example commands:
# flye --pacbio-raw Mbac.fastq.gz --out-dir flye_meta -t 10 --meta
# flye --pacbio-raw Mbac.fastq.gz --genome-size 2.5m --out-dir flye2.5cov50 -t 20 --asm-coverage 50

# Batch Flye run for multiple FASTQ files
input_dir=~/scratch/PacBio_peat/Data/BAM
output_dir_base="flye_Output_meta"
threads=10

mkdir -p "$output_dir_base"

for fastq_file in "$input_dir"/*.fastq.gz; do
    base_name=$(basename "$fastq_file" .fastq.gz)
    output_dir="${output_dir_base}/${base_name}"
    
    mkdir -p "$output_dir"
    
    flye --pacbio-raw "$fastq_file" --out-dir "$output_dir" -t "$threads" --meta
done

# View resulting assembly graphs with Bandage:
# https://timkahlke.github.io/LongRead_tutorials/BAN.html

# ========================================================================================
# Additional References:
# - Nanopack (QC tools): https://github.com/wdecoster/nanopack
# - Canu Docs: https://canu.readthedocs.io/en/latest/quick-start.html
# - Flye GitHub Issues: https://github.com/fenderglass/Flye/issues/128
# ========================================================================================

