#!/bin/bash
#==============================================================================================
# title          : phylogenomic_pipeline.sh
# description    : A complete pipeline to (1) select genome files based on a folder list, 
#                  (2) copy and rename them, (3) reformat for Anvi'o, and (4) run phylogenomic 
#                  analyses using HMMs.
# author         : Justus Nweze
# date           : 20230809
# version        : v1
# usage          : bash Rename_many_file.sh
#==============================================================================================

#==============================================================================================
# Step 1: Use file list to copy all selected genome folders into another folder
#==============================================================================================

# Create directory to hold original .fna files
mkdir -p Original_Data

# Define paths
data_folder="Downloaded"
selected_folder="Original_Data"
list_file="List.txt"

# Loop through each line in List.txt
while IFS= read -r folder; do
    if [[ -d "$data_folder/$folder" && $folder == GCA_* ]]; then
        find "$data_folder/$folder" -type f -name "*.fna" -exec cp {} "$selected_folder" \;
    else
        echo "Folder $data_folder/$folder does not exist or does not start with GCA_"
    fi
done < "$list_file"

#==============================================================================================
# Step 2: Create rename table file (Rename_files.txt)
#==============================================================================================

ls Original_Data > Rename_files.txt
# Manually edit the file to add new names in the second column, tab-delimited

#==============================================================================================
# Step 3: Rename and copy selected files into 'Data' directory
#==============================================================================================

SOURCE_DIR="Original_Data"
Rename_Table="Rename_files.txt"
TARGET_DIR="Data"

# Ensure directories exist
if [ ! -d "$SOURCE_DIR" ]; then
  echo "Source directory '$SOURCE_DIR' not found."
  exit 1
fi

if [ ! -f "$Rename_Table" ]; then
  echo "Rename table file '$Rename_Table' not found."
  exit 1
fi

mkdir -p "$TARGET_DIR"

while IFS=$'\t' read -r old_name new_name; do
  old_file="$SOURCE_DIR/$old_name"
  new_file="$TARGET_DIR/$new_name.fa"
  
  if [ -e "$old_file" ]; then
    cp "$old_file" "$new_file"
    echo "Copied and renamed: $old_name -> $new_name.fa"
  else
    echo "File '$old_name' not found in the source directory."
  fi
done < "$Rename_Table"

#==============================================================================================
# Step 4: Reformat genomes for Anvi'o
#==============================================================================================

# Define working directory (update path to match your actual location)
wkdir="$HOME/scratch/PacBio_peat/Data/Phylogenomics/Downloaded_Genomes/Methylocystis"
# Uncomment below if you're working with Methylosinus
# wkdir="$HOME/scratch/PacBio_peat/Data/Phylogenomics/Downloaded_Genomes/Methylosinus"

cd "$wkdir" || exit

mkdir -p Data/Reformatted

for i in Data/*.fa; do    
  anvi-script-reformat-fasta "$i" \
    -o "Data/Reformatted/$(basename "$i" .fa).fa" \
    --simplify-names \
    --report "Data/Reformatted/name_conversions.txt" \
    --seq-type NT
done

#==============================================================================================
# Step 5: Generate Anvi’o contigs database
#==============================================================================================

mkdir -p Database

for i in Data/Reformatted/*.fa; do
  base_name=$(basename "$i" .fa)
  db_name="Database/${base_name}.db"

  anvi-gen-contigs-database -f "$i" -o "$db_name" -T 10
  anvi-run-hmms -c "$db_name"
done

# Optional: export stats (uncomment if needed)
# anvi-display-contigs-stats Database/contigs.db --report-as-text -o contigs-stats.csv

#==============================================================================================
# Step 6: Prepare 'external-genomes.txt' file (manually)
#==============================================================================================

# Format:
# name    contigs_db_path
# MAG1    Database/MAG1.db
# MAG2    Database/MAG2.db
# (use tab-delimited formatting)
#Generate external_genomes.txt from your .db files with something like:
echo -e "name\tcontigs_db_path" > external_genomes.txt
for f in Database/*.db; do
  fname=$(basename "$f" .db)
  echo -e "${fname}\t$f" >> external_genomes.txt
done


#==============================================================================================
# Step 7: Check available HMMs in Bacteria_71
#==============================================================================================

anvi-get-sequences-for-hmm-hits \
  --external-genomes external_genomes.txt \
  --hmm-source Bacteria_71 \
  --list-available-gene-names

# Create HMM presence/absence matrix
anvi-script-gen-hmm-hits-matrix-across-genomes \
  --external-genomes external_genomes.txt \
  --hmm-source Bacteria_71 \
  -o Bacteria_71_hmm_hits.txt

#==============================================================================================
# Step 8: Extract and concatenate 16S sequences
#==============================================================================================

anvi-get-sequences-for-hmm-hits \
  --external-genomes external_genomes.txt \
  -o 16S_concatenated-sequences.fa \
  --hmm-source Ribosomal_RNA_16S \
  --return-best-hit \
  --concatenate

# Know which genomes have 16S
anvi-script-gen-hmm-hits-matrix-across-genomes \
  --external-genomes external_genomes.txt \
  --hmm-source Ribosomal_RNA_16S \
  -o 16S_hmm_hits.txt

#==============================================================================================
# END OF SCRIPT
#==============================================================================================

