#!/usr/bin/env Rscript

# ==========================
# LOAD REQUIRED PACKAGES
# ==========================
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(dplyr)

# ==========================
# 1. LOAD GENE COUNTS
# ==========================
counts_file <- "counts_Methylocystis.txt"

# Read featureCounts output
counts <- read.table(counts_file, header=TRUE, row.names=1, check.names=FALSE)

# Drop non-count annotation columns if present
annotation_cols <- c("Chr", "Start", "End", "Strand", "Length")
counts <- counts[, !(colnames(counts) %in% annotation_cols), drop=FALSE]

# Convert to integer matrix (DESeq2 expects integers)
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

# ==========================
# 2. CLEAN COLUMN NAMES
# ==========================
# Remove "aligned/" prefix and fix "_sorted" duplicates
colnames(counts) <- gsub("^aligned/", "", colnames(counts))
colnames(counts) <- gsub("_sorted_sorted$", "_sorted", colnames(counts))
colnames(counts) <- gsub("_sorted$", "_sorted", colnames(counts))  # ensure consistency

# ==========================
# 3. CREATE METADATA
# ==========================
sample_names <- colnames(counts)

# Extract pH condition from sample names
pH_condition <- sapply(sample_names, function(x) {
  if (grepl("^PH5_", x)) return("5.0")
  if (grepl("^PH68_", x)) return("6.8")
  if (grepl("^PH9_", x)) return("9.0")
  return("Unknown")
})

meta <- data.frame(
  sample = sample_names,
  pH = factor(pH_condition, levels = c("5.0", "6.8", "9.0")),
  stringsAsFactors = TRUE
)
rownames(meta) <- meta$sample

# Keep only pH 5.0 and 6.8 samples
keep_samples <- meta$pH %in% c("5.0", "6.8")
counts <- counts[, keep_samples, drop=FALSE]
meta <- meta[keep_samples, , drop=FALSE]

# Ensure factor levels match
meta$pH <- droplevels(meta$pH)
colnames(counts) <- rownames(meta)
stopifnot(all(colnames(counts) == rownames(meta)))  # sanity check

# ==========================
# 4. CREATE DESEQ2 OBJECT
# ==========================
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ pH
)

# ==========================
# 5. FILTER LOW-EXPRESSED GENES
# ==========================
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# ==========================
# 6. RUN DESEQ2
# ==========================
dds <- DESeq(dds)

# ==========================
# 7. COMPARE pH 5.0 vs 6.8
# ==========================
res_5vs68 <- results(dds, contrast = c("pH", "5.0", "6.8"))
res_5vs68 <- res_5vs68[order(res_5vs68$pvalue), ]  # sort by p-value

cat("\n--- SUMMARY ---\n")
cat("Total genes analyzed:\t", nrow(res_5vs68), "\n")
cat("Significant DEGs (padj < 0.001):\t", sum(res_5vs68$padj < 0.001, na.rm = TRUE), "\n")
cat("Upregulated (log2FC > 1, padj < 0.001):\t", sum(res_5vs68$log2FoldChange > 1 & res_5vs68$padj < 0.001, na.rm = TRUE), "\n")
cat("Downregulated (log2FC < -1, padj < 0.001):\t", sum(res_5vs68$log2FoldChange < -1 & res_5vs68$padj < 0.001, na.rm = TRUE), "\n")


# ==========================
# 8. SAVE RESULTS
# ==========================
#results_file <- "results/DEGs/pH5_vs_pH68.csv"
dir.create("results/DEGs", showWarnings = FALSE, recursive = TRUE)
#write.csv(as.data.frame(res_5vs68), file = results_file)
#cat("✅ DEG results saved to:", results_file, "\n")

# ==========================
# 9. SUMMARY STATISTICS + BASIC ANNOTATION (USING PARENT & CDS ONLY)
# ==========================
cat("\n--- SUMMARY ---\n")
cat("Total genes analyzed:\t", nrow(res_5vs68), "\n")
cat("Significant DEGs (padj < 0.05):\t", sum(res_5vs68$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated (log2FC > 1, padj < 0.05):\t", sum(res_5vs68$log2FoldChange > 1 & res_5vs68$padj < 0.05, na.rm = TRUE), "\n")
cat("Downregulated (log2FC < -1, padj < 0.05):\t", sum(res_5vs68$log2FoldChange < -1 & res_5vs68$padj < 0.05, na.rm = TRUE), "\n")

# ==========================================================
# 10. LOAD GFF FILE
# ==========================================================
gff_file <- "./reference_genome/Methylocystis_sp_ref.gff"
gff <- read.delim(gff_file, 
                  header=FALSE, 
                  stringsAsFactors=FALSE, 
                  comment.char="#", 
                  fill=TRUE, 
                  sep="\t", 
                  quote="", 
                  check.names=FALSE)

colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

gff <- gff[!is.na(gff$seqid), ]

# ==========================================================
# 11. EXTRACT GFF ATTRIBUTES USING REGEX
# ==========================================================
extract_attr <- function(attr_str, key) {
  matches <- regmatches(attr_str, regexec(paste0(key, "=([^;]+)"), attr_str))
  sapply(matches, function(x) if(length(x) > 1) x[2] else NA)
}

gff$parent <- extract_attr(gff$attributes, "Parent")
gff$locus_tag <- extract_attr(gff$attributes, "locus_tag")
gff$gene_name <- extract_attr(gff$attributes, "Name")
gff$product <- extract_attr(gff$attributes, "product")
gff$go_terms <- extract_attr(gff$attributes, "Ontology_term")
gff$go_function <- extract_attr(gff$attributes, "go_function")
gff$go_process <- extract_attr(gff$attributes, "go_process")
gff$inference <- extract_attr(gff$attributes, "inference")
gff$protein_id <- extract_attr(gff$attributes, "protein_id")
gff$dbxref <- extract_attr(gff$attributes, "Dbxref")
gff$gene_symbol <- extract_attr(gff$attributes, "gene")  # e.g., rplU

# For gene features, set parent = ID (so we can link CDS to gene later)
gff[gff$type == "gene", "parent"] <- extract_attr(gff[gff$type == "gene", "attributes"], "ID")

# ==========================================================
# 12. CONVERT DESEQ RESULTS TO DATAFRAME
# ==========================================================
res_df <- as.data.frame(res_5vs68)
res_df$gene_id <- rownames(res_df)
sig_degs <- res_df[which(res_df$padj < 0.05), ]

# ==========================================================
# 13. MERGE WITH GFF (CDS FEATURES ONLY)
# ==========================================================
# Filter to CDS only — they have the most complete annotation
gff_cds <- gff[gff$type == "CDS", ]

# Merge: match sig_degs$gene_id (e.g., "gene-ACNHKD_RS00045") to gff_cds$parent
top_degs_annotated <- merge(
  sig_degs,
  gff_cds[, c("seqid", "start", "end", "strand", "parent", "locus_tag", "gene_name", "gene_symbol", "product", "go_terms", "go_function", "go_process", "inference", "protein_id", "dbxref")],
  by.x = "gene_id",
  by.y = "parent",
  all.x = TRUE
)

# Remove any remaining duplicates (should be rare)
top_degs_annotated <- top_degs_annotated[!duplicated(top_degs_annotated$gene_id), ]

# ==========================================================
# 14. SAVE ANNOTATED DEG RESULTS
# ==========================================================
# Save basic annotated results
#write.csv(top_degs_annotated, file = "results/DEGs/pH5_vs_pH68_top_DEGs_annotated.csv", row.names = FALSE)
#cat("✅ Annotated DEG results saved to: results/DEGs/pH5_vs_pH68_top_DEGs_annotated.csv\n")

# ==========================================================
# 15. COMPUTE TPM AND LOG2TPM PER CONDITION (CD vs ME)
# ==========================================================
# Use annotated DEGs from Step 13
degs <- top_degs_annotated

# Ensure numeric columns
degs$baseMean <- as.numeric(degs$baseMean)
degs$padj <- as.numeric(degs$padj)
degs$log2FoldChange <- as.numeric(degs$log2FoldChange)

# Load counts and metadata (if not already in environment)
counts <- read.table("counts_Methylocystis.txt", header=TRUE, row.names=1, check.names=FALSE)
colnames(counts) <- gsub("^aligned/", "", colnames(counts))
colnames(counts) <- gsub("_sorted_sorted$", "_sorted", colnames(counts))

sample_names <- colnames(counts)
pH_condition <- sapply(sample_names, function(x) {
  if (grepl("^PH5_", x)) return("5.0")
  if (grepl("^PH68_", x)) return("6.8")
  return(NA)
})
meta <- data.frame(sample = sample_names, pH = factor(pH_condition, levels = c("5.0", "6.8")))
rownames(meta) <- meta$sample

# Keep only pH 5.0 and 6.8 samples
keep_samples <- !is.na(meta$pH)
counts <- counts[, keep_samples]
meta <- meta[keep_samples, ]

# Compute TPM per sample
total_counts <- colSums(counts)
tpm_matrix <- sweep(counts, 2, total_counts, FUN="/") * 1e6
tpm_df <- as.data.frame(tpm_matrix)
tpm_df$gene_id <- rownames(tpm_df)

# Pivot to long format and merge with metadata
tpm_long <- tpm_df %>%
  pivot_longer(cols = -gene_id, names_to = "sample", values_to = "TPM") %>%
  left_join(meta, by = c("sample" = "sample"))

# Average TPM per condition and compute log2TPM
tpm_wide <- tpm_long %>%
  group_by(gene_id, pH) %>%
  summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = pH, values_from = TPM, names_prefix = "TPM_") %>%
  mutate(
    log2TPM_5.0 = log2(TPM_5.0 + 1),
    log2TPM_6.8 = log2(TPM_6.8 + 1)
  )

# Merge condition-specific TPMs with annotated DEGs
degs_tpm <- degs %>%
  left_join(tpm_wide, by = "gene_id")

# Save full TPM/log2TPM table
write.csv(degs_tpm,
          file = "results/DEGs/pH5_vs_pH68_TPM_log2TPM_by_condition.csv",
          row.names = FALSE)
cat("✅ Full TPM/log2TPM by condition saved to results/DEGs/pH5_vs_pH68_TPM_log2TPM_by_condition.csv\n")
