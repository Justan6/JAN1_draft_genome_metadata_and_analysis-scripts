# Load libraries
library(ggplot2)
library(dplyr)
library(readODS)
library(grid) # for unit()
library(forcats)  # for factor manipulation, including fct_inorder()

# Read data from ODS
# NOTE: Using the user's file path.
deg_data <- read_ods("results/DEGs/Functional_Categories.ods", sheet = "Selected_pH_9.0_vs_pH6.8_Functional_Categories")

# Define the custom COG order
custom_cog_order <- c("M", "V", "C", "P", "E", "H", "O", "J", "N", "T", "K", "L", "F", "I", "Q", "G", "U", "W", "R", "S", "A", "B", "D", "Y", "Z", "X")

# --- 2. Data Wrangle and Reordering ---
df_plot <- deg_data %>%
  # 1. Filter out "hypothetical protein"
  filter(product != "hypothetical protein") %>%
  
  # 2. Assign Direction based on the log2FoldChange threshold (keep only significant ones)
  mutate(Direction = case_when(
    log2FoldChange > 1  ~ "Up-regulated (pH 9.0)",
    log2FoldChange < -1 ~ "Down-regulated (pH 9.0)",
    TRUE ~ NA_character_ # Assign NA for non-significant changes
  )) %>%
  
  # IMPORTANT: Filter out the non-significant genes (where Direction is NA) 
  # to avoid the color manual scale error and only plot DEGs.
  filter(!is.na(Direction)) %>% 
  
  # IMPORTANT: Apply the custom order to COG_LETTER for correct facet display
  mutate(COG_LETTER = factor(COG_LETTER, levels = custom_cog_order)) %>%
  
  # 3. Make the gene names unique within each COG group
  # Correction: Use n() > 1 to identify duplicates (2 or more entries).
  group_by(COG_LETTER, product) %>%
  mutate(product_unique = ifelse(n() > 1, paste0(product, " (", row_number(), ")"), product)) %>%
  ungroup() %>%
  
  # 4. Order the rows: first by COG (factor), then by log2FoldChange.
  arrange(COG_LETTER, log2FoldChange) %>%
  
  # 5. Convert the unique product names to a factor based on the arranged order (fct_inorder)
  mutate(product_factor = fct_inorder(product_unique))

# --- 3. Plotting with ggplot2 ---
plot_genes <- ggplot(df_plot, aes(x = log2FoldChange, y = product_factor, fill = Direction)) +
  
  # Add the bars (geom_col is shorthand for geom_bar(stat="identity"))
  geom_col(width = 0.8) +
  
  # Use a manual color scale (Red for down, Blue for up)
  scale_fill_manual(values = c("Up-regulated (pH 9.0)" = "blue", 
                               "Down-regulated (pH 9.0)" = "red")) +
  
  # Use facet_grid to separate the plot by COG_LETTER (now correctly ordered)
  # switch="y" moves the facet label to the left side
  facet_grid(COG_LETTER ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # ✅ Symmetric axis for log2FoldChange values
  scale_x_continuous(
    breaks = seq(-10, 10, by = 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  # Customize the theme to look like a standard barplot for DEGs
  theme_bw() +
  theme(
    # Move the gene names (y-axis labels) to the right side
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.right = element_text(size = 7, hjust = 0), # Increased size slightly
    axis.ticks.y.right = element_line(),
    
    # Customize the facet strips (COG labels)
    strip.text.y = element_text(angle = 0, size = 10, face = "bold"), # Increased size
    strip.background = element_rect(fill = "white", color = "black"),
    
    # General clean-up
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  
  # Set labels
  labs(
    x = expression(log[2]~"Fold Change"),
    y = NULL,
    title = "Differentially Expressed Genes by COG Functional Category (pH 9.0 vs pH 6.8)"
  ) +
  
  # Force the y-axis labels to appear on the right side of the facet
  scale_y_discrete(position = "right")

# Print the plot
print(plot_genes)

# Save (Directories will be created if they do not exist)
fig_dir <- "results/figures/Functional_Category"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

ggsave(file.path(fig_dir, "pH9.0_vs_pH6.8_DEGs_barplot_shaded.pdf"), plot_genes, width = 15, height = 30)
ggsave(file.path(fig_dir, "pH9.0_vs_pH6.8_DEGs_barplot_shaded.svg"), plot_genes, width = 15, height = 30)
ggsave(file.path(fig_dir, "pH9.0_vs_pH6.8_DEGs_barplot_shaded.png"), plot_genes, width = 15, height = 30)

