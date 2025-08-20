#!/bin/bash
# Anvi’o pangenomic workflow
# Author: Justus
# Date: [Add date]
# Purpose: Pangenomic analysis of Methylocystis and Methylosinus genomes, including functional enrichment and phylogenomics
#============================================================================================================================
# This workflow enables you to:
# - Identify gene clusters across multiple genomes (FASTA files),
# - Combine metagenome-assembled genomes (MAGs) from different sources,
# - Interactively visualize and partition pangenomes,
# - Estimate genome relationships,
# - Perform phylogenomic analyses based on gene clusters,
# - Annotate genes and inspect alignments,
# - Add contextual metadata,
# - Quantify homogeneity within gene clusters,
# - Perform functional enrichment,
# - Compute and visualize average nucleotide identity (ANI), among other functions.

#============================================================================================================================
conda activate anvio-8

# Create directories for outputs
mkdir -p Pangenome Contigs_db Genome_database Pangenome_analysis SPLIT_PANs ANI

# Step 1: Generate contigs databases for all genomes
cd ~/proj/Peat_soil/PacBio_peat/Data/Comparative_genomics/Data/Methylocystis/Reformatted

for f in *.fa; do
    anvi-gen-contigs-database -T 10 -f "$f" -o ../Contigs_db/"${f}_out.db"
done

# Step 2: Annotate contigs with NCBI COGs and HMMs
cd ../Contigs_db

anvi-setup-ncbi-cogs

for db in *.db; do
    anvi-run-ncbi-cogs --contigs-db "$db" --num-threads 12 --search-with blastp
done

for db in *.db; do
    anvi-run-hmms --contigs-db "$db" --num-threads 8 --just-do-it
done

# Step 3: Generate genomes storage file from external genome list
cd ~/proj/Peat_soil/PacBio_peat/Data/Comparative_genomics/Data/Methylocystis/Pangenome

anvi-gen-genomes-storage --external-genomes Genome_database_data.txt --output-file Genome_database/Genomes-GENOMES.db

# Step 4: Run pangenome analyses (Methylocystis and Methylosinus separately)

# Methylocystis
anvi-pan-genome --genomes-storage Genome_database/Genomes-GENOMES.db \
    --output-dir Pangenome_analysis/Methylocystis --num-threads 8 \
    --use-ncbi-blast --mcl-inflation 5 --minbit 0.5 --project-name Methylocystis_genomes

# Methylosinus
anvi-pan-genome --genomes-storage Genome_database/Genomes-GENOMES.db \
    --output-dir Pangenome_analysis/Methylosinus --num-threads 8 \
    --use-ncbi-blast --mcl-inflation 5 --minbit 0.5 --project-name Methylosinus_genomes

# Step 5: Compute functional enrichment in pangenomes

anvi-compute-functional-enrichment-in-pan -p Pangenome_analysis/Methylocystis/Methylocystis_genomes-PAN.db \
    --annotation-source COG_FUNCTION --categories L,S -o Pangenome_analysis/Methylocystis/enrichment_LS.txt

anvi-compute-functional-enrichment-in-pan -p Pangenome_analysis/Methylocystis/Methylocystis_genomes-PAN.db \
    --annotation-source COG_FUNCTION --categories N,O -o Pangenome_analysis/Methylocystis/enrichment_NO.txt

anvi-compute-functional-enrichment-in-pan -p Pangenome_analysis/Methylosinus/Methylosinus_genomes-PAN.db \
    --annotation-source COG_FUNCTION --categories L,S -o Pangenome_analysis/Methylosinus/enrichment_LS.txt

anvi-compute-functional-enrichment-in-pan -p Pangenome_analysis/Methylosinus/Methylosinus_genomes-PAN.db \
    --annotation-source COG_FUNCTION --categories N,O -o Pangenome_analysis/Methylosinus/enrichment_NO.txt

# Step 6: Compute genome similarities (ANI) with fastANI

anvi-compute-genome-similarity --external-genomes Genome_database_data.txt \
    -o ANI/Methylocystis --pan-db Pangenome_analysis/Methylocystis/Methylocystis_genomes-PAN.db \
    --program fastANI --num-threads 8

anvi-compute-genome-similarity --external-genomes Genome_database_data.txt \
    -o ANI/Methylosinus --pan-db Pangenome_analysis/Methylosinus/Methylosinus_genomes-PAN.db \
    --program fastANI --num-threads 8

# Step 7: Extract single-copy core genes (SCGs) for phylogenomics

anvi-get-sequences-for-gene-clusters --pan-db Pangenome_analysis/Methylocystis/Methylocystis_genomes-PAN.db \
    --genomes-storage Genome_database/Genomes-GENOMES.db \
    --min-num-genomes-gene-cluster-occurs 42 --max-num-genes-from-each-genome 1 \
    --concatenate-gene-clusters --output-file Pangenome_analysis/Methylocystis-SCGs.fa

anvi-get-sequences-for-gene-clusters --pan-db Pangenome_analysis/Methylosinus/Methylosinus_genomes-PAN.db \
    --genomes-storage Genome_database/Genomes-GENOMES.db \
    --min-num-genomes-gene-cluster-occurs 22 --max-num-genes-from-each-genome 1 \
    --concatenate-gene-clusters --output-file Pangenome_analysis/Methylosinus-SCGs.fa

# Step 8: Align SCG sequences with MAFFT

mafft --auto Pangenome_analysis/Methylocystis-SCGs.fa > Pangenome_analysis/Methylocystis-SCGs_aligned.fa
mafft --auto Pangenome_analysis/Methylosinus-SCGs.fa > Pangenome_analysis/Methylosinus-SCGs_aligned.fa

# Step 9: Trim alignments to remove poorly aligned columns (>50% gaps)

trimal -in Pangenome_analysis/Methylocystis-SCGs_aligned.fa -out Pangenome_analysis/Methylocystis-SCGs_aligned-clean.fa -gt 0.50
trimal -in Pangenome_analysis/Methylosinus-SCGs_aligned.fa -out Pangenome_analysis/Methylosinus-SCGs_aligned-clean.fa -gt 0.50

# Step 10: Build phylogenomic trees with IQ-TREE using WAG model and 1000 ultrafast bootstraps

iqtree -s Pangenome_analysis/Methylocystis-SCGs_aligned-clean.fa -nt 8 -m WAG -bb 1000
iqtree -s Pangenome_analysis/Methylosinus-SCGs_aligned-clean.fa -nt 8 -m WAG -bb 1000

# Step 11: Prepare layer order files for visualization (add your Newick tree paths manually if needed)

echo -e "item_name\tdata_type\tdata_value" > Pangenome_analysis/Methylocystis-phylogenomic-layer-order.txt
# Replace the path below with your actual .treefile or .contree file path from IQ-TREE
echo -e "SCGs_Bayesian_Tree\tnewick\t$(cat Pangenome_analysis/Methylocystis-SCGs_aligned-clean.fa.contree)" >> Pangenome_analysis/Methylocystis-phylogenomic-layer-order.txt

echo -e "item_name\tdata_type\tdata_value" > Pangenome_analysis/Methylosinus-phylogenomic-layer-order.txt
echo -e "SCGs_Bayesian_Tree\tnewick\t$(cat Pangenome_analysis/Methylosinus-SCGs_aligned-clean.fa.contree)" >> Pangenome_analysis/Methylosinus-phylogenomic-layer-order.txt

# Step 12: Prepare metadata for pangenome visualization

echo -e "genome_name\tTree\tSource\tAcidity\tSource\tOur_isolate" > Pangenome_analysis/Tree_Source_categories.txt
awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' Pangenome_analysis/pan_groups_orginal.txt >> Pangenome_analysis/Tree_Source_categories.txt

# Step 13: Import layer orders and metadata into pangenome databases

anvi-import-misc-data -p Pangenome_analysis/Methylocystis/Methylocystis_genomes-PAN.db -t layer_orders Pangenome_analysis/Methylocystis-phylogenomic-layer-order.txt
anvi-import-misc-data -p Pangenome_analysis/Methylocystis/Methylocystis_genomes-PAN.db -t layer_orders Pangenome_analysis/pan_groups.txt
anvi-import-misc-data -p Pangenome_analysis/Methylocystis/Methylocystis_genomes-PAN.db -t layers Pangenome_analysis/Tree_Source_categories.txt

anvi-import-misc-data -p Pangenome_analysis/Methylosinus/Methylosinus_genomes-PAN.db -t layer_orders Pangenome_analysis/Methylosinus-phylogenomic-layer-order.txt

# Step 14: Visualize pangenomes interactively

anvi-display-pan --pan-db Pangenome_analysis/Methylocystis/Methylocystis_genomes-PAN.db \
    --genomes-storage Genome_database/Genomes-GENOMES.db --title Methylocystis_Pangenome

anvi-display-pan --pan-db Pangenome_analysis/Methylosinus/Methylosinus_genomes-PAN.db \
    --genomes-storage Genome_database/Genomes-GENOMES.db --title Methylosinus_Pangenome

# Step 15: Split the pangenome into bins (clusters) and visualize

anvi-split -p Pangenome_analysis/Methylocystis/Methylocystis_genomes-PAN.db \
    --genomes-storage Genome_database/Genomes-GENOMES.db -C default -o SPLIT_PANs

anvi-display-pan -p SPLIT_PANs/SC_Core/PAN.db --genomes-storage Genome_database/Genomes-GENOMES.db



# Quantifying functional enrichment in a pangenome https://merenlab.org/2016/11/08/pangenomics-v2/#splitting-the-pangenome

# The passage describes how Anvi’o allows users to identify and quantify functionally enriched gene clusters in a pangenome—specifically functions that are overrepresented in certain groups (or clades) of genomes.


#The summary step gives you two important things: a static HTML web page you can compress and share with your colleagues or add into your publication as a supplementary data file, and a comprehensive TAB-delimited file in the output directory that describes every gene cluster.
# 1. Run Functional Enrichment Analysis  

# You can now proceed with functional enrichment analysis using the Tree and Source categories. For example: 
# For the Tree Category  
# Cite  Shaiber, Willis et al (https://doi.org/10.1186/s13059-020-02195-w)
anvi-compute-functional-enrichment-in-pan \
    --pan-db Pangenome_analysis/Methylocystis_genomes-PAN.db \
    --genomes-storage Genome_database/Genomes-GENOMES.db \
    -o enriched-functions_Tree.txt \
    --category-variable Tree \
    --annotation-source COG20_FUNCTION \
    --functional-occurrence-table-output functions_occurrence_Tree.txt \
    --include-gc-identity-as-function                


# For the Source Category

anvi-compute-functional-enrichment-in-pan \
    --pan-db Pangenome_analysis/Methylocystis_genomes-PAN.db \
    --genomes-storage Genome_database/Genomes-GENOMES.db \
    -o enriched-functions_Source.txt \
    --category-variable Source \
    --annotation-source COG20_FUNCTION \
    --functional-occurrence-table-output functions_occurrence_Source.txt \
    --include-gc-identity-as-function           


# For the Acidity Category

anvi-compute-functional-enrichment-in-pan \
    --pan-db Pangenome_analysis/Methylocystis_genomes-PAN.db \
    --genomes-storage Genome_database/Genomes-GENOMES.db \
    -o enriched-functions_Acidity.txt \
    --category-variable Acidity \
    --annotation-source COG20_FUNCTION \
    --functional-occurrence-table-output functions_occurrence_Acidity.txt \
    --include-gc-identity-as-function 

# For the Our_isolate Category
anvi-compute-functional-enrichment-in-pan \
    --pan-db Pangenome_analysis/Methylocystis_genomes-PAN.db \
    --genomes-storage Genome_database/Genomes-GENOMES.db \
    -o enriched-functions_Our_isolate.txt \
    --category-variable Our_isolate \
    --annotation-source COG20_FUNCTION \
    --functional-occurrence-table-output functions_occurrence_Our_isolate.txt \
    --include-gc-identity-as-function 

# For the nos_gene Category
anvi-compute-functional-enrichment-in-pan \
    --pan-db Pangenome_analysis/Methylocystis_genomes-PAN.db \
    --genomes-storage Genome_database/Genomes-GENOMES.db \
    -o enriched-functions_nos_gene.txt \
    --category-variable nos_gene \
    --annotation-source COG20_FUNCTION \
    --functional-occurrence-table-output functions_occurrence_nos_gene.txt \
    --include-gc-identity-as-function 


anvi-estimate-metabolism

anvi-get-metabolic-model-file
#=====================================================================================================
anvi-gen-phylogenomic-tree -f Pangenome_analysis/Methylomirabiales-SCGs-clean.fa -o phylogenomic-tree.txt

anvi-interactive -p Profile_db/phylogenomic-profile.db -t Pangenome_analysis/phylogenomic-tree.txt --title "Methylomirabiales_Phylogenomics" --manual



