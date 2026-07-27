#=================================================================================================
# Top 50 up- and downregulated genes (padj < 0.001)
# Heatmaps: TPM (pH9.0/pH6.8) and log2FoldChange with product annotation
#=================================================================================================

# Load libraries
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

#=================================================================================================
# Input and output paths
#=================================================================================================
input_file <- "results/DEGs/pH9_vs_pH68_TPM_log2TPM_by_condition.csv"
output_dir <- "results/figures/top50_heatmaps"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#=================================================================================================
# Load and prepare data
#=================================================================================================
degs <- read.csv(input_file, stringsAsFactors = FALSE)

# Ensure numeric
numeric_cols <- c("TPM_9.0", "TPM_6.8", "log2TPM_9.0", "log2TPM_6.8", "log2FoldChange", "padj")
degs[numeric_cols] <- lapply(degs[numeric_cols], as.numeric)

# Filter for significant genes (padj < 0.001 and |log2FC| > 1)
sig_genes <- degs %>%
  filter(!is.na(padj) & padj < 0.001 & abs(log2FoldChange) > 1)

# Select top 25 upregulated and top 25 downregulated
top25_up <- sig_genes %>%
  arrange(desc(log2FoldChange)) %>%
  slice(1:25)

top25_down <- sig_genes %>%
  arrange(log2FoldChange) %>%
  slice(1:25)

# Combine up- and downregulated genes
top50 <- bind_rows(top25_up, top25_down)

# Preserve consistent order
row_order <- order(top50$log2FoldChange, decreasing = TRUE)
top50 <- top50[row_order, ]

locus_labels <- top50$locus_tag
product_labels <- top50$product

#=================================================================================================
# Shared color scale across both conditions (log2TPM)
#=================================================================================================
max_log2TPM <- max(c(top50$log2TPM_9.0, top50$log2TPM_6.8), na.rm = TRUE)
col_fun <- colorRamp2(
  c(0, max_log2TPM / 3, 2 * max_log2TPM / 3, max_log2TPM),
  c("white", "lightpink", "hotpink", "deeppink")
)

#=================================================================================================
# TPM matrices
#=================================================================================================
mat_9.0_log <- matrix(top50$log2TPM_9.0, ncol = 1, dimnames = list(NULL, "9.0 (pH 9.0)"))
mat_6.8_log  <- matrix(top50$log2TPM_6.8, ncol = 1, dimnames = list(NULL, "6.8 (Control)"))
mat_9.0_raw <- matrix(top50$TPM_9.0, ncol = 1)
mat_6.8_raw  <- matrix(top50$TPM_6.8, ncol = 1)

#=================================================================================================
# TPM Heatmaps (values = raw TPM; color = log2TPM)
#=================================================================================================
ht_9.0 <- Heatmap(
  mat_9.0_log,
  name = "log2TPM",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_title = "pH 9.0 (Low pH)",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", mat_9.0_raw[i, j]), x, y, gp = gpar(fontsize = 11))
  },
  width = unit(3, "cm")
)

ht_6.8 <- Heatmap(
  mat_6.8_log,
  name = "log2TPM",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_title = "pH 6.8 (control)",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", mat_6.8_raw[i, j]), x, y, gp = gpar(fontsize = 11))
  },
  width = unit(3, "cm")
)

# Row annotation for locus tags
row_anno <- rowAnnotation(
  locus_tag = anno_text(locus_labels, location = 0, just = "left", gp = gpar(fontsize = 11)),
  width = unit(4.5, "cm")
)

# Combine TPM heatmaps
ht_list_tpm <- ht_6.8 + ht_9.0 + row_anno

#=================================================================================================
# log2FoldChange heatmap (colored by log2FC, product as annotation)
#=================================================================================================
fc_min <- min(top50$log2FoldChange, na.rm = TRUE)
fc_max <- max(top50$log2FoldChange, na.rm = TRUE)
fc_range <- max(abs(fc_min), abs(fc_max))

col_fun_fc <- colorRamp2(
  c(-fc_range, 0, fc_range),
  c("#87CEFA", "white", "#FFA07A")
)

mat_fc <- matrix(top50$log2FoldChange, ncol = 1, dimnames = list(NULL, "log2FC"))

ht_fc <- Heatmap(
  mat_fc,
  name = "log2FC",
  col = col_fun_fc,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_title = "log2FoldChange",
  width = unit(2, "cm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mat_fc[i, j]), x, y, gp = gpar(fontsize = 11, col = "black"))
  }
)

# Product annotation
row_anno_pro <- rowAnnotation(
  product = anno_text(product_labels, location = 0, just = "left", gp = gpar(fontsize = 11)),
  width = unit(4.5, "cm")
)

# Combine FC heatmap
ht_list_fc <- ht_fc + row_anno_pro

#=================================================================================================
# Combine all heatmaps (pH 9.0, pH 6.8, log2FC)
#=================================================================================================
ht_list_all <- ht_6.8 + ht_9.0 + ht_fc + row_anno + row_anno_pro

#=================================================================================================
# Save outputs
#=================================================================================================
pdf(file.path(output_dir, "2.pH9.0_vs_pH6.8_Top50_up_downregulated_log2FC_product_overlay.pdf"), width = 20, height = 10)
draw(ht_list_all, heatmap_legend_side = "bottom")
dev.off()

png(file.path(output_dir, "2.pH9.0_vs_pH6.8_Top50_up_downregulated_log2FC_product_overlay.png"), width = 1500, height = 800, res = 150)
draw(ht_list_all, heatmap_legend_side = "bottom")
dev.off()

svg(file.path(output_dir, "2.pH9.0_vs_pH6.8_Top50_up_downregulated_log2FC_product_overlay.svg"), width = 20, height = 10)
draw(ht_list_all, heatmap_legend_side = "bottom")
dev.off()

message("✅ Heatmaps successfully generated in: ", output_dir)

