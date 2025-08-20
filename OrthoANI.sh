#!/bin/bash -   
#title          :makeblastdb.sh
#description    :Genome assembly
#author         :Justus Nweze
#date           :20230809
#version        :v1
#usage          :Go to the Analysis folder and OrthoANI.sh
#==============================================================================================
# Run OrthoANI: OrthoANI is a metric proposed by Lee et al. in 2016 to improve computation of Average Nucleotide Identity. It uses BLASTn to find orthologous blocks in a pair of sequences, and then computes the average identity only considering alignments of reciprocal orthologs.
# Change this MAG1.fna to be the MAG that you want to compare with others.

cd /home/nwezejus/scratch/PacBio_peat/Data/Phylogenomics/Downloaded_Genomes/Methylocystis
cd /home/nwezejus/scratch/PacBio_peat/Data/Phylogenomics/Downloaded_Genomes/Methylosinus


# Run the script in Data/Reformatted folder
    conda activate anvio-7.1

    for i in *.fa; do
    db_name=$(basename "$i" .fa)  # Extract the database name
    output=$(orthoani -q MAG1.fa -r "$i")  # Run orthoani and capture its output
    echo -e "$i\t$output" >> orthoani_output.txt  # Append "$i" and its output to the output file
    done








