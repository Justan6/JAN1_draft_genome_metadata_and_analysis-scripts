JAN1 Draft Genome Metadata and Analysis Scripts

This repository contains metadata and analysis scripts associated with the draft genome of Candidatus Methylocystis sumavensis strain JAN1. The workflow integrates multiple bioinformatics pipelines to provide a comprehensive understanding of the organism's genetics, metabolism, and evolutionary relationships.

Workflow Overview
1. Genome Acquisition and Preparation

Data Sources: Genome data are acquired from public databases and internal sequencing efforts.

Long-Read Processing: Raw BAM files (e.g., from PacBio sequencing) are converted to FASTQ using the PacBio BAM Toolkit (bam2fastq).

Assembly: Reads are assembled into contigs using long-read assemblers such as Canu or Flye.

Comparative Genome Preparation: Selected genomes are reformatted and renamed for compatibility with Anvi’o for downstream analysis.

2. Comparative and Phylogenetic Analysis

Pangenomics: Gene clusters across Methylocystis and Methylosinus genomes are analyzed with Anvi’o. Steps include contig database creation, gene annotation with NCBI COGs and HMMs, pangenome computation, functional enrichment analysis, ANI calculations, and visualization.

Phylogenomics: Single-copy core genes (SCGs) are extracted, aligned with MAFFT, and used to build phylogenetic trees via IQ-TREE. GTDB-Tk is used for taxonomic classification and phylogenomic tree construction from MAGs.

Genome Similarity: OrthoANI and EzAAI are employed to calculate Average Nucleotide Identity (ANI) and Average Amino Acid Identity (AAI), key metrics for species delineation.

3. Metabolic and Functional Annotation

DRAM & METABOLIC: Functional annotation pipelines are applied to MAGs to summarize metabolic capabilities and predict biogeochemical cycling pathways.

Targeted Gene Analysis: Key genes (e.g., pmoA, mmoX) are extracted, aligned, and used to construct phylogenetic trees to study evolutionary patterns.

HMMER: Protein family profiles, such as NiFe hydrogenases, are identified across genomes for functional insights.

KEGGCharter: Gene annotations are mapped onto KEGG metabolic pathways to visualize metabolic potential, including methane metabolism.

This repository provides a fully reproducible framework for genomic analysis and functional annotation of methanotrophic bacteria, facilitating comparative genomics and evolutionary studies.
