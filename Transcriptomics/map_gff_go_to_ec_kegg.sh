#!/bin/bash
# ============================================================
# Script: map_gff_go_to_ec_kegg.sh
# Description: Extract GO terms from a GFF file and map them
#              to EC numbers and KEGG reactions.
# Author: Justus (2025)
# ============================================================

# --- Input GFF file ---
GFF="reference_genome/Methylocystis_sp_ref.gff"  # <-- change this to your filename
OUTDIR="GO_mapping_results"
mkdir -p "$OUTDIR"

# --- Step 1. Extract gene info and GO terms ---
echo "Extracting gene ID, product, and GO terms..."
awk -F'\t' '$3=="CDS" && $9 ~ /Ontology_term=/ {
    match($9, /ID=([^;]+)/, id);
    match($9, /product=([^;]+)/, prod);
    match($9, /Ontology_term=([^;]+)/, go);
    if (id[1] && go[1])
        print id[1]"\t"prod[1]"\t"go[1];
}' "$GFF" > "$OUTDIR/gene_go_raw.tsv"

# --- Step 2. Prepare GO term list ---
echo "Preparing unique GO term list..."
awk -F'\t' '{print $3}' "$OUTDIR/gene_go_raw.tsv" | tr ',' '\n' | sort -u > "$OUTDIR/go_terms.txt"

# --- Step 3. Download mapping files from GO ---
echo "Downloading mapping files from Gene Ontology..."
wget -q -nc https://current.geneontology.org/ontology/external2go/ec2go -O "$OUTDIR/ec2go"
wget -q -nc https://current.geneontology.org/ontology/external2go/kegg_reaction2go -O "$OUTDIR/kegg_reaction2go"

# --- Step 4. Map GO terms to EC numbers and KEGG reactions ---
echo "Mapping GO terms to EC and KEGG..."
grep -F -f "$OUTDIR/go_terms.txt" "$OUTDIR/ec2go" > "$OUTDIR/mapped_GO_to_EC.txt"
grep -F -f "$OUTDIR/go_terms.txt" "$OUTDIR/kegg_reaction2go" > "$OUTDIR/mapped_GO_to_KEGG.txt"

# --- Step 5. Format the mapping into simple lookup tables ---
awk -F'[>;]' '{gsub(/^ +| +$/, "", $1); gsub(/^ +| +$/, "", $3); print $3"\t"$1}' "$OUTDIR/mapped_GO_to_EC.txt" > "$OUTDIR/GO_to_EC.tsv"
awk -F'[>;]' '{gsub(/^ +| +$/, "", $1); gsub(/^ +| +$/, "", $3); print $3"\t"$1}' "$OUTDIR/mapped_GO_to_KEGG.txt" > "$OUTDIR/GO_to_KEGG.tsv"

# --- Step 6. Create final CSV table ---
echo "Creating final combined CSV table..."
{
    echo "Gene_ID,Product,GO_Term,EC_Number,KEGG_Reaction"
    while IFS=$'\t' read -r gene prod golist; do
        IFS=',' read -ra GO_ARRAY <<< "$golist"
        for go in "${GO_ARRAY[@]}"; do
            go=$(echo "$go" | xargs)
            ec=$(grep -w "$go" "$OUTDIR/GO_to_EC.tsv" | cut -f2 | paste -sd ";" -)
            kegg=$(grep -w "$go" "$OUTDIR/GO_to_KEGG.tsv" | cut -f2 | paste -sd ";" -)
            echo "$gene,\"$prod\",$go,$ec,$kegg"
        done
    done < "$OUTDIR/gene_go_raw.tsv"
} > "$OUTDIR/final_GO_EC_KEGG_table.csv"

echo "✅ Mapping complete!"
echo "Output files are in: $OUTDIR/"
echo "Main result: $OUTDIR/final_GO_EC_KEGG_table.csv"

