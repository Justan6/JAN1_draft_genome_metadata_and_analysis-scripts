#!/bin/bash
# title          : EzAAI_AAI_calculation.sh
# description    : Automated pipeline to extract protein databases using Prodigal and compute AAI using EzAAI.
# author         : Justus Nweze
# date           : 12.07.2025
#============================================================================================================

# CONDA ENVIRONMENT SETUP
# Only needs to be run once:
# conda create --name EzAAI -c bioconda ezaai

# Activate environment
conda activate EzAAI

# 1. Prepare directories
mkdir -p AAI/DB_output
cd ~/scratch/PacBio_peat/Data/Phylogenomics/Downloaded_Genomes/Methylocystis/Data

# 2. Extract protein databases using EzAAI + Prodigal
for file in *.fa; do
    base_name="${file%.fa}"
    EzAAI extract -i "$file" -o "../AAI/DB_output/${base_name}.db"
done

# 3. Compute AAI
cd ~/scratch/PacBio_peat/Data/Phylogenomics/Downloaded_Genomes/Methylocystis/AAI
mkdir -p AAI_out

EzAAI calculate -i DB_output/ -j DB_output/ -o AAI_out/aai.tsv

# Optional: deactivate environment
conda deactivate
