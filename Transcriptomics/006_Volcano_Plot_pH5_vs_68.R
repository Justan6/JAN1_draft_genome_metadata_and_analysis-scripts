# --- Volcano Plot R Script ---

# 1. Install and Load Required Packages (if not already installed)
# install.packages(c("ggplot2", "dplyr", "readr"))
library(ggplot2)
library(dplyr)
library(readr)

# 2. Define Thresholds
LOG2FC_THRESHOLD <- 1.0
P_ADJ_THRESHOLD <- 0.001

# 3. Load Data
# IMPORTANT: Adjust the file path if your file is not in the working directory
data <- read_csv("results/DEGs/pH5_vs_pH68_TPM_log2TPM_by_condition.csv")

# 4. Prepare Data for Plotting
# Calculate -log10(padj) and categorize the genes
data_plot <- data %>%
  # Handle padj = 0 by replacing it with a small value for log transform
  mutate(
    padj_clean = ifelse(padj == 0, min(padj[padj > 0], na.rm = TRUE) / 10, padj),
    neg_log10_padj = -log10(padj_clean)
  ) %>%
  # Define the differential expression status
  mutate(
    DE_status = case_when(
      log2FoldChange >= LOG2FC_THRESHOLD & padj < P_ADJ_THRESHOLD ~ "Up-regulated (pH 5.0)",
      log2FoldChange <= -LOG2FC_THRESHOLD & padj < P_ADJ_THRESHOLD ~ "Down-regulated (pH 6.8)",
      TRUE ~ "Not significant"
    )
  )

# 5. Create the Volcano Plot
volcano_plot <- ggplot(data_plot, aes(x = log2FoldChange, y = neg_log10_padj)) +
  
  # A. Draw the points, colored by status
  geom_point(aes(color = DE_status), alpha = 0.8, size = 1.5) +
  
  # B. Set colors and legend order
  scale_color_manual(
    values = c(
      "Up-regulated (pH 5.0)" = "#1F78B4",     # Blue
      "Down-regulated (pH 6.8)" = "#E31A1C",   # Red
      "Not significant" = "grey70"              # Grey
    ),
    name = "Significance"
  ) +
  
  # C. Add threshold lines (dashed grey)
  geom_hline(yintercept = -log10(P_ADJ_THRESHOLD), linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), linetype = "dashed", color = "grey30") +
  
  # D. Customize labels and title
  labs(
    title = "Volcano Plot: Differential Gene Expression (pH 5.0 vs pH 6.8)",
    x = expression(log[2]~"Fold Change"~"(pH 5.0 / pH 6.8)"),
    y = expression(-log[10]~"(Adjusted P-value)")
  ) +
  
  # E. Set theme
  theme_minimal(base_size = 14) +
  
  # F. Adjust plot limits (optional, but good for symmetry)
  xlim(c(floor(min(data_plot$log2FoldChange)), ceiling(max(data_plot$log2FoldChange))))

# 6. Display the plot
print(volcano_plot)

dir.create("results/figures/Volcano_plot", showWarnings = FALSE, recursive = TRUE)

# 7. Save the plot to a file
ggsave("results/figures/Volcano_plot/1.Volcano_Plot_5.0.0vs6.8_Msumavensis_labeled.png", plot = volcano_plot, width = 10, height = 8, units = "in", dpi = 300)
ggsave("results/figures/Volcano_plot/1.Volcano_Plot_5.0.0vs6.8_Msumavensis_labeled.svg", plot = volcano_plot, width = 10, height = 8, units = "in", dpi = 300)
ggsave("results/figures/Volcano_plot/1.Volcano_Plot_5.0.0vs6.8_Msumavensis_labeled.pdf", plot = volcano_plot, width = 10, height = 8, units = "in", dpi = 300)
