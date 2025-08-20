# JAN1 draft genome metadata and analysis-scripts
Metadata and analysis scripts associated with the draft genome of Candidatus Methylocystis sumavensis strain JAN1
 The workflow integrates several bioinformatics pipelines to achieve a comprehensive understanding of these organisms' genetics, metabolism, and evolutionary relationships.

1. Genome Acquisition and Preparation

The project begins by acquiring and preparing genome data from various sources, including publicly available databases and internal sequencing efforts. For long-read sequencing data (e.g., PacBio), raw BAM files are converted to FASTQ format using the PacBio BAM Toolkit (bam2fastq). The reads are then assembled into contigs using long-read assemblers such as Canu or Flye.

For comparative analyses, a dedicated pipeline is used to select, copy, and rename a list of target genomes. These genomes are then reformatted for compatibility with the Anvi'o platform.

2. Comparative and Phylogenetic Analysis

The core of the project involves comparative genomics and phylogenomic analysis using Anvi'o, GTDB-Tk, and OrthoANI.

    Pangenomics: The Anvi'o pangenomic workflow is applied to analyze gene clusters across multiple genomes from the genera Methylocystis and Methylosinus separately. This involves generating contigs databases, annotating genes with NCBI COGs and HMMs, and computing a pangenome database. The pipeline can then be used to compute functional enrichment, perform average nucleotide identity (ANI) calculations, and visualize the pangenomes.

    Phylogenomics: Single-copy core genes (SCGs) are extracted from the pangenomes, aligned with MAFFT, and used to construct phylogenetic trees with IQ-TREE. The results are prepared for visualization in Anvi'o's interactive interface. The project also uses GTDB-Tk for taxonomic classification and to build phylogenomic trees from MAGs.

    Genome Similarity: OrthoANI and EzAAI are used to calculate Average Nucleotide Identity (ANI) and Average Amino Acid Identity (AAI) between pairs of genomes, respectively, which are key metrics for species delineation and comparison.

3. Metabolic and Functional Annotation

The project delves into the metabolic potential of the genomes, particularly focusing on methane metabolism.

    DRAM & METABOLIC: The DRAM pipeline is used to perform a detailed functional annotation of MAGs, followed by summarization and distillation of metabolic information. The METABOLIC software is also used to predict metabolic and biogeochemical cycling pathways.

    Targeted Gene Analysis: Specific genes of interest, such as pmoA (particulate methane monooxygenase) and mmoX, are extracted from the annotated genomes. These genes are then aligned with MAFFT and used to build dedicated phylogenetic trees with IQ-TREE to investigate their evolutionary history.

    HMMER: HMMER is used to build and search for specific protein family profiles (e.g., NiFe hydrogenases), providing a powerful way to identify key functional genes across the genomes.

    KEGGCharter: This tool visualizes the metabolic pathways by mapping gene annotations onto KEGG metabolic maps, allowing for a clear representation of an organism's metabolic potential, such as methane metabolism.
