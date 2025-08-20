#!/bin/bash
#title          :Extracting and Merging Contigs from MAG Files
#description    :This script is designed to extract specific contigs from a set of MAG (Metagenome-Assembled Genome) files and merge them into a single output file. Each contig is identified by its corresponding file name and contig identifier, which are listed in an input text file (Contig_lists.txt). The script processes each file individually, extracts the specified contigs, and saves them in the same directory as the source file, ensuring that each extracted contig is clearly labeled with its original MAG file name. Finally, all extracted contigs are concatenated into a single file (All_MAGs.faa), with each contig's header updated to include the MAG file name for easy identification.
#author         :Justus Nweze
#date           :21.08.2024
#version        :v1
#usage          :
#==============================================================================================
# Directories
# Have genomes separately
input_dir="Proteins"       # Directory containing the MAG FASTA files
output_dir="Target_pMMO"   # Directory to save the output files
contig_list="pMMO_list.txt"     # File containing the list of contig IDs to extract

# Ensure the output directory exists
mkdir -p "$output_dir"

# Remove duplicate rows from the contig list
temp_contig_list=$(mktemp)
sort "$contig_list" | uniq > "$temp_contig_list"

# Loop through each contig in the list
while read -r contig_name; do
    # Extract the MAG name (e.g., MAG1) and the contig ID (e.g., c_000000000001_1405)
    mag_name=$(echo "$contig_name" | awk -F'_c_' '{print $1}')
    contig_id=$(echo "$contig_name" | awk -F'_c_' '{print "c_" $2}')

    # Define the input and output file paths
    input_file="$input_dir/$mag_name.faa"
    output_file="$output_dir/$mag_name.faa"

    # Check if the input file exists
    if [[ -f $input_file ]]; then
        # Extract the specific contig from the input file
        awk -v contig="$contig_id" -v mag="$mag_name" '
        BEGIN {found=0}
        /^>/ {
            if ($1 == ">" contig) {
                print ">" mag "_" substr($0, 2)
                found = 1
            } else {
                found = 0
            }
            next
        }
        found {print}
        ' "$input_file" >> "$output_file"
    else
        echo "Warning: Input file $input_file not found, skipping..."
    fi
done < "$temp_contig_list"

# Clean up the temporary file
rm "$temp_contig_list"

echo "Extraction complete. Results saved in $output_dir."

#==============================================================================================

#one single FASTA file (genes.faa) and you want to extract only the exact contig IDs (like MAG1_c_000000000001_3766 and not partial matches like ..._37667, ..._37666, etc.)

# Input and output setup
input_file="DRAM_anno/genes.faa"         # Single combined FASTA file with all MAGs
output_dir="Target_pMMO"       # Output directory
contig_list="pMMO_list.txt"         # List of contig IDs to extract (e.g., MAG1_c_000000000001_3766)

# Ensure output directory exists
mkdir -p "$output_dir"

# Remove duplicates from the contig list
temp_contig_list=$(mktemp)
sort "$contig_list" | uniq > "$temp_contig_list"

# Make an associative array of target contigs for exact match
declare -A targets
while read -r id; do
    targets["$id"]=1
done < "$temp_contig_list"

# Output file to store results
output_file="$output_dir/selected_genes.faa"

# Parse and extract exact matches
awk -v out="$output_file" '
    BEGIN {
        while ((getline line < "'"$temp_contig_list"'") > 0) {
            target[line] = 1
        }
        close("'"$temp_contig_list"'")
    }
    /^>/ {
        # Extract just the contig ID (up to first space or full line)
        header = substr($0, 2)
        split(header, parts, " ")
        id = parts[1]
        if (id in target) {
            printing = 1
            print ">" header >> out
        } else {
            printing = 0
        }
        next
    }
    printing {
        print >> out
    }
' "$input_file"

# Clean up
rm "$temp_contig_list"

echo "Extraction complete. Output saved to $output_file"

#==============================================================================================
# Iterate over all .fa files in the current directory
for file in Target_pMMO/*.faa; do cat "$file" >> Target_pMMO_concatenated.fa; done


for file in Target_mmoX/*.faa; do cat "$file" >> Target_mmoX_concatenated.fa; done
#==============================================================================================

# To remove the part of the header starting from the # symbol onward (instead of after the |), you can use the following sed command in a bash script. This script will retain the header up to the first # and remove everything after it:

# Input file
input_file="16S_concatenated-sequences"

# Output file
output_file="16S_concatenated-sequences_SH.fa"

# Use sed to process the file and shorten the headers
sed -E 's/(^>[^#]+).*/\1/' "$input_file" > "$output_file"

echo "Headers have been shortened and saved to $output_file"




# Activate the maft environment
conda activate mafft_env

mkdir iqtree_
# align the sequences
mafft --auto Target_pmoA_concatenated.fa_SH.fa > iqtree_file/Target_pmoA_concatenated_aln.fa

#Strategy:  FFT-NS-2 (cystis) and L-INS-i (sinus)


iqtree_file

iqtree -s Target_pmoA_concatenated_aln.fa -o MAG0 -nt 10 -m MFP -bb 1000


mkdir iqtree_mmoX
# align the sequences
mafft --auto Target_mmoX_concatenated.fa_SH.fa > iqtree_mmoX/Target_mmoX_concatenated_aln.fa

# L-INS-i (Probably most accurate, very slow)

iqtree_mmoX

iqtree -s Target_mmoX_concatenated_aln.fa -o MAG1_mmoY -nt 10 -m MFP -bb 1000
