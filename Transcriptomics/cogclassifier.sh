# COG(Cluster of Orthologous Genes) is a database that plays an important role in the annotation, classification, and analysis of microbial gene function. Functional annotation, classification, and analysis of each gene in newly sequenced bacterial genomes using the COG database is a common task. However, there was no COG functional classification command line software that is easy-to-use and capable of producing publication-ready figures. Therefore, I developed COGclassifier to fill this need. COGclassifier can automatically perform the processes from searching query sequences into the COG database, to annotation and classification of gene functions, to generation of publication-ready figures (See figure below).
# conda install -c conda-forge -c bioconda cogclassifier

mkdir -p ~//proj/Peat_soil/PacBio_peat/transcriptomics/COGclassifier_resources
cd ~/proj/Peat_soil/PacBio_peat/transcriptomics/COGclassifier_resources

wget https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab -O cog-24.fun.tab
wget https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.def.tab -O cog-24.def.tab
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -O cddid.tbl.gz
gunzip cddid.tbl.gz  # decompress

wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz -O Cog_LE.tar.gz
tar -xzvf Cog_LE.tar.gz


cd /proj/Peat_soil/PacBio_peat/transcriptomics/reference_genome

COGclassifier -i ~/proj/Peat_soil/PacBio_peat/transcriptomics/reference_genome/Methylocystis_sp_ref.faa -o Output_cogclassifier


