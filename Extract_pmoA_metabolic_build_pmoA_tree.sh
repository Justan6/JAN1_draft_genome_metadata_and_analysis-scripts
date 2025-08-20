#!/bin/bash

# ========================================================================
# Script: Extract_pmoA_metabolic_build_pmoA_tree.sh
# Author: Justus Nweze
# Affiliation: Microbial Ecologist
# Description: 
#   This script extracts pmoA (K10944) sequences from annotated MAGs, 
#   renames and consolidates the files, performs multiple sequence alignment 
#   using MAFFT, and constructs a maximum likelihood phylogenetic tree using IQ-TREE.
# Dependencies: mafft, iqtree, conda (mafft_env), GNU tools
# Date: 2025-07-12
# ========================================================================

set -euo pipefail

# ===========================
# Step 1: Extract pmoA hits
# ===========================

# Define source and target directories
SOURCE_DIR="Metab_output"
TARGET_DIR="Metabolic_pmoA"

# Create target directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Find and copy K10944.hmm.collection.faa files
find "$SOURCE_DIR" -type f -name "K10944.hmm.collection.faa" | while read -r SOURCE_FILE; do
    # Extract MAG name from the grandparent directory
    MAG_NAME=$(basename "$(dirname "$(dirname "$SOURCE_FILE")")")

    # Define destination path with renamed file
    DEST_FILE="$TARGET_DIR/${MAG_NAME}_K10944.hmm.collection.faa"

    cp "$SOURCE_FILE" "$DEST_FILE"
    echo "Copied: $SOURCE_FILE -> $DEST_FILE"
done

# ===========================
# Step 2: Concatenate all pmoA sequences
# ===========================

cat "$TARGET_DIR"/*.faa > Target_pmoA_concatenated.fa
echo "All pmoA sequences concatenated to Target_pmoA_concatenated.fa"

# ===========================
# Step 3: Extract unique FASTA headers
# ===========================

grep '^>' Target_pmoA_concatenated.fa | sed 's/^>//' | awk '{print $1}' | sort -u > pmoA_list.txt
echo "Extracted identifiers saved to pmoA_list.txt"

# ===========================
# Step 4: Multiple Sequence Alignment with MAFFT
# ===========================

conda activate mafft_env

mkdir -p iqtree_pmoA

mafft --auto Target_pmoA_concatenated.fa > iqtree_pmoA/Target_pmoA_concatenated_aln.fa
echo "Alignment complete: iqtree_pmoA/Target_pmoA_concatenated_aln.fa"

# ===========================
# Step 5: Phylogenetic Tree Inference with IQ-TREE
# ===========================

cd iqtree_pmoA

# Replace this with an actual outgroup sequence ID if appropriate
OUTGROUP="MAG1c_000000000001_1404_pmoC"

iqtree -s Target_pmoA_concatenated_aln.fa \
       -o "$OUTGROUP" \
       -nt AUTO \
       -m MFP \
       -bb 1000 \
       -alrt 1000

echo "IQ-TREE analysis complete. Results are saved in iqtree_pmoA/"


