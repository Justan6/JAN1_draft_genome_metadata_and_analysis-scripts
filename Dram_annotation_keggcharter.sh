#!/bin/bash
#==================================================================================================================
# Title         : DRAM Annotation & KEGGCharter Visualization Pipeline
# Description   : Annotates MAGs with DRAM and visualizes metabolic potential using KEGGCharter.
# Author        : Justus Nweze
# Date          : 2022-12-08
# Usage         : Run from project root. Assumes conda environments and DRAM/KEGGCharter are installed.
#==================================================================================================================

set -euo pipefail
logFile="30_DRAM.log"

#--------------------------------------#
# STEP 1: Set up working environment
#--------------------------------------#
# Ensure you are in the correct directory
cd Peat_metagenome2022/Data/Single_assmbly2/Metabolic_Bins || exit 1

echo "Activating DRAM environment..." | tee -a "$logFile"
conda activate myDRAM

#--------------------------------------#
# STEP 2: Annotate MAGs with DRAM
#--------------------------------------#
echo "Running DRAM annotation..." | tee -a "$logFile"

DRAM.py annotate \
  -i Reformatted/MAG*/*.fasta \
  -o Annotation/DRAM_anno \
  --threads 5

# Optional alternative run with CheckM & GTDB taxonomy (update paths if needed)
# DRAM.py annotate \
#   -i MAGs/*.fasta \
#   -o annotation \
#   --min_contig_size 1000 \
#   --gtdb_taxonomy gtdbtk.bac120.summary.tsv \
#   --checkm_quality checkm_results

#--------------------------------------#
# STEP 3: Summarize DRAM outputs
#--------------------------------------#
echo "Summarizing annotations..." | tee -a "$logFile"

DRAM.py distill \
  -i Annotation/DRAM_anno/annotations.tsv \
  -o Annotation/DRAM_anno/genome_summaries \
  --trna_path Annotation/DRAM_anno/trnas.tsv \
  --rrna_path Annotation/DRAM_anno/rrnas.tsv

#--------------------------------------#
# STEP 4: Additional run (e.g., on older genomes)
#--------------------------------------#
cd /proj/Peat_soil/Peat_metagenome2022/Data/Single_assembly/Old/Methylomirabilis_Genomes/New_comparism/ANNOTATIONS/DRAM || exit 1

echo "Running DRAM on old genomes..." | tee -a "$logFile"
DRAM.py annotate -i ../Data/*.fasta -o DRAM_anno --threads 5

DRAM.py distill \
  -i DRAM_anno/annotations.tsv \
  -o DRAM_anno/genome_summaries \
  --trna_path DRAM_anno/trnas.tsv \
  --rrna_path DRAM_anno/rrnas.tsv

#--------------------------------------#
# STEP 5: KEGGCharter Visualization
#--------------------------------------#
echo "Activating KEGGCharter environment..." | tee -a "$logFile"
conda activate My_KEGGCharter

# Generate KEGG map for Methane Metabolism (map00680)
echo "Generating KEGG maps with KEGGCharter..." | tee -a "$logFile"

# Example with reCOGnizer results
keggcharter \
  -f reCOGnizer_results.xlsx \
  -rd resources_directory \
  -keggc 'KEGG' \
  -ecc 'EC number' \
  -cogc 'COG ID' \
  -iq -it "Methylomirabiales" \
  -mm 00680 \
  -o Methane_metabolism

# Example with UPIMAPI results
keggcharter \
  -f UPIMAPI_results.tsv \
  -rd resources_directory \
  -keggc 'KEGG' \
  -koc 'KO' \
  -ecc 'EC number' \
  -cogc 'COG ID' \
  -iq -it "Methylomirabiales" \
  -mm 00680 \
  -o Methane_metabolism

echo "Pipeline complete." | tee -a "$logFile"

