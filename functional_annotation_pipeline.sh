#!/bin/bash

################################################################################
# Functional and Taxonomic Annotation Pipeline
# Tools: reCOGnizer, UPIMAPI, KEGGCharter
# Author: [Your Name]
# Last updated: [Date]
# Description: Pipeline for annotating protein sequences using domain-based,
#              sequence similarity-based, and pathway-based tools.
################################################################################

# Exit immediately if a command exits with a non-zero status
set -e

##############################################
# 1. Preprocess: Add filename to FASTA headers
##############################################
echo "Adding filenames to FASTA headers..."
for file in *.faa; do
    perl -i -0777 -pe ' $x=$ARGV;$x=~s/\.faa//g; s/>/>${x}_/ ' "$file"
done

##############################################
# 2. Domain-based annotation using reCOGnizer
##############################################
echo "Running reCOGnizer annotations..."
conda activate reCOGnizer
mkdir -p ../recognizer_output
for file in *.faa; do
    recognizer -f "$file" \
        -o "../recognizer_output/$file" \
        -rd ~/scratch/PacBio_peat/Data/resources_directory \
        -t 6
done

##########################################################
# 3. Sequence similarity and taxonomy annotation (UPIMAPI)
##########################################################
echo "Running UPIMAPI annotations..."
conda activate upimapi
mkdir -p ../Annotation/UPIMAPI
for file in *.faa; do
    upimapi.py -i "$file" \
        -o "../Annotation/UPIMAPI/$file" \
        --database swissprot \
        -t 8
done

# Example for running on a concatenated file:
# upimapi.py -i Translated_GeneRFinder_named/All_Bins.faa -o upimapi_directory --database swissprot -t 10

# Example: annotation with existing BLAST results
# upimapi.py -i upimapi_directory/aligned.blast -o upimapi_directory/uniprotinfo2 --blast

##########################################################
# 4. KEGG Pathway Mapping with KEGGCharter
##########################################################
echo "Running KEGGCharter for pathway mapping..."
conda activate keggcharter
mkdir -p KEGGCharter_Output

# Adjust the KO and taxonomy column names based on your recognizer output structure
keggcharter.py \
    --file recognizer_output/reCOGnizer_results.xlsx \
    --output KEGGCharter_Output \
    --resources-directory KEGGCharter_resources_Dir \
    -mm map00680 \
    --kegg-column "Cross-reference (KEGG)"

##########################################################
# 5. Header reformatting (custom for Vojta or GeneRFinder)
##########################################################
echo "Standardizing FASTA headers..."
# For GeneRFinder-style headers
dir="Translated_GeneRFinder_named"
for file in "$dir"/*.faa; do
    if [[ -f "$file" ]]; then
        sed -i 's/^>\([0-9]\+__[^ ]\+ len=[0-9]\+\), \([^ ]\+\)/>\2, \1/' "$file"
    fi
done

# For Vojta's bins (adjust regex as needed)
dir="New"
for file in "$dir"/*.faa; do
    if [[ -f "$file" ]]; then
        sed -i 's/^>\([0-9]\+__[^ ]\+ len=[0-9]\+\), \(k[0-9]+_[0-9]+_[0-9]+\)/>\2, \1/' "$file"
    fi
done

echo "Annotation pipeline completed successfully."

