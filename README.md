JAN1 Draft Genome Metadata and Analysis Scripts

This repository contains metadata and analysis scripts associated with the draft genome of Methylocystis sumavensis strain JAN1. The workflow integrates multiple bioinformatics pipelines to provide a comprehensive understanding of the organism's genetics, metabolism, and evolutionary relationships.

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

4. Transcriptomic Analysis (RNA-seq)

Scripts in the Transcriptomics/ folder cover the comparative transcriptomic analysis of strain JAN1 under acidic (pH 5.0), optimum (pH 6.8), and alkaline (pH 9.0) conditions:

Quality Control and Trimming: Raw paired-end reads are assessed with FastQC/MultiQC and trimmed with fastp (001_fastqc_multiqc_trim.sh).

Alignment: Trimmed reads are aligned to the JAN1 genome with HISAT2 and processed with SAMtools (002_hisat2_Methylocystis.sh).

Quantification: Gene-level read counts are generated with featureCounts based on the CDS annotation (003_counts_Methylocystis.sh).

Differential Expression: Counts are analyzed in DESeq2 to identify significantly up- and downregulated genes at pH 5.0 and pH 9.0 relative to pH 6.8 (004_analysis_DESeq_pH5_vs_pH6.8.R, 004_analysis_DESeq_pH9_vs_pH6.8.R).

Visualization: Heatmaps of the top differentially expressed genes and volcano plots are generated for each pH comparison (005_Top50_Heatmap_sig_up_or_down_regu_*.R, 006_Volcano_Plot_*.R).

Functional Categorization: Differentially expressed genes are classified into COG categories (cogclassifier.sh, COG_merge_all_Script.R), tested for over/under-representation per category with exact binomial and Fisher's exact tests (cog_within_category_binom_test.R), and summarized as bar/pie charts (Functional_category_bar_pie_plot.R, 007_pH5_vs_pH68_Functional_Category_plots.R, 007_pH9_vs_pH68_Functional_Category_plots.R).

Gene Annotation Support: CDS and GO term information is extracted from the genome GFF file and mapped to EC numbers and KEGG reactions to support interpretation of differentially expressed genes (Extract_CDS_information.sh, map_gff_go_to_ec_kegg.sh).

This repository provides a fully reproducible framework for genomic, functional, and transcriptomic analysis of methanotrophic bacteria, facilitating comparative genomics, evolutionary studies, and investigation of pH stress adaptation.
