# Script: COG functional classification plots with per-category bias tests

library(ggplot2)
library(dplyr)
library(readODS)
library(RColorBrewer)
library(svglite)
library(Polychrome)
library(grid)
library(tidyr)

set.seed(42)

cog_order <- c(
  "J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W",
  "U", "O", "X", "C", "G", "E", "F", "H", "I", "P", "Q", "R", "S"
)

palette_colors <- createPalette(26, c("#ff0000", "#00ff00", "#0000ff"))
palette_colors_named <- setNames(palette_colors, cog_order)

output_dir <- "/proj/Peat_soil/PacBio_peat/transcriptomics/results/figures/COG_FUN_CAT"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
input_file <- "/data/Peat_soil/PacBio_peat/transcriptomics/results/DEGs/Functional_Categories.ods"

prepare_cog_data <- function(sheet_name) {
  data <- read_ods(input_file, sheet = sheet_name) %>%
    mutate(
      COG_Label = paste(COG_LETTER, ":", COG_DESCRIPTION),
      Regulation = case_when(
        log2FoldChange > 0 ~ "Upregulated",
        log2FoldChange < 0 ~ "Downregulated",
        TRUE ~ "Neutral"
      )
    ) %>%
    filter(Regulation != "Neutral")

  summary <- data %>%
    group_by(COG_LETTER, COG_DESCRIPTION, COG_Label, Regulation) %>%
    summarise(AbsCount = n(), .groups = "drop") %>%
    mutate(
      Count = if_else(Regulation == "Downregulated", -AbsCount, AbsCount),
      COG_LETTER = factor(COG_LETTER, levels = cog_order)
    )

  label_levels <- summary %>%
    distinct(COG_LETTER, COG_Label) %>%
    arrange(COG_LETTER) %>%
    pull(COG_Label)

  summary$COG_Label <- factor(summary$COG_Label, levels = label_levels)
  summary
}

compute_bias_stats <- function(summary_df) {
  counts_wide <- summary_df %>%
    transmute(
      COG_LETTER,
      COG_DESCRIPTION,
      Regulation,
      AbsCount = abs(Count)
    ) %>%
    group_by(COG_LETTER, COG_DESCRIPTION, Regulation) %>%
    summarise(AbsCount = sum(AbsCount), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = Regulation,
      values_from = AbsCount,
      values_fill = 0
    ) %>%
    rename(Up = Upregulated, Down = Downregulated) %>%
    mutate(COG_LETTER = factor(COG_LETTER, levels = cog_order)) %>%
    arrange(COG_LETTER)

  total_up <- sum(counts_wide$Up)
  total_down <- sum(counts_wide$Down)

  fisher_results <- lapply(seq_len(nrow(counts_wide)), function(i) {
    row <- counts_wide[i, ]
    contingency <- matrix(
      c(
        row$Up,
        row$Down,
        total_up - row$Up,
        total_down - row$Down
      ),
      nrow = 2,
      byrow = TRUE
    )

    test <- fisher.test(contingency)
    direction <- if (is.finite(unname(test$estimate)) && unname(test$estimate) < 1) {
      "Downregulated"
    } else if (!is.finite(unname(test$estimate)) || unname(test$estimate) > 1) {
      "Upregulated"
    } else {
      "Balanced"
    }

    data.frame(
      COG_LETTER = row$COG_LETTER,
      COG_DESCRIPTION = row$COG_DESCRIPTION,
      Up = row$Up,
      Down = row$Down,
      OddsRatio = unname(test$estimate),
      PValue = test$p.value,
      Direction = direction
    )
  })

  stats <- bind_rows(fisher_results) %>%
    mutate(
      FDR = p.adjust(PValue, method = "BH"),
      Significant = FDR < 0.05
    )

  max_count <- max(c(stats$Up, stats$Down), na.rm = TRUE)
  offset <- max(0.6, max_count * 0.08)

  stats %>%
    mutate(
      Asterisk = if_else(Significant, "*", ""),
      AnnotationY = case_when(
        Significant & Direction == "Upregulated" ~ Up + offset,
        Significant & Direction == "Downregulated" ~ -(Down + offset),
        TRUE ~ NA_real_
      )
    )
}

make_pie_data <- function(summary_df, regulation) {
  summary_df %>%
    filter(Regulation == regulation) %>%
    group_by(COG_LETTER, COG_DESCRIPTION) %>%
    summarise(Count = sum(abs(Count)), .groups = "drop") %>%
    mutate(Percent = Count / sum(Count) * 100)
}

make_bar_plot <- function(summary_df, stats_df, comparison_label) {
  sig_stats_df <- stats_df[stats_df$Significant, , drop = FALSE]

  ggplot(summary_df, aes(x = COG_LETTER, y = Count, fill = COG_LETTER)) +
    geom_bar(stat = "identity", color = "black", width = 0.8) +
    geom_hline(yintercept = 0, color = "black") +
    geom_text(
      data = sig_stats_df,
      aes(x = COG_LETTER, y = AnnotationY, label = Asterisk),
      inherit.aes = FALSE,
      size = 7,
      fontface = "bold"
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_fill_manual(
      values = palette_colors_named,
      labels = levels(summary_df$COG_Label),
      guide = guide_legend(ncol = 1)
    ) +
    labs(
      title = paste0("COG Functional Classification (", comparison_label, ")"),
      subtitle = "* FDR < 0.05 from Fisher's exact test for up/down bias within each COG category",
      x = "Functional Category",
      y = "Number of Sequences (Up vs Downregulated)",
      fill = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.text.x = element_text(vjust = 1, hjust = 1),
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      axis.line = element_line(color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

make_pie_plot <- function(pie_df, plot_title) {
  palette_local <- colorRampPalette(brewer.pal(12, "Paired"))(nrow(pie_df))

  ggplot(pie_df, aes(x = "", y = Percent, fill = COG_LETTER)) +
    geom_bar(stat = "identity", color = "black", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(
      values = palette_local,
      labels = paste(pie_df$COG_LETTER, pie_df$COG_DESCRIPTION, sep = ": ")
    ) +
    labs(title = plot_title, fill = NULL) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 10)
    )
}

save_plot_set <- function(prefix, plot_list) {
  for (plot_name in names(plot_list)) {
    file_prefix <- file.path(output_dir, paste0(prefix, "_", plot_name))
    ggsave(paste0(file_prefix, ".pdf"), plot_list[[plot_name]], width = 20, height = 15, units = "in", dpi = 300)
    ggsave(paste0(file_prefix, ".svg"), plot_list[[plot_name]], width = 20, height = 15, units = "in")
    ggsave(paste0(file_prefix, ".png"), plot_list[[plot_name]], width = 20, height = 15, units = "in", dpi = 300)
  }
}

run_comparison <- function(sheet_name, comparison_label, prefix) {
  summary_df <- prepare_cog_data(sheet_name)
  stats_df <- compute_bias_stats(summary_df)

  write.csv(
    stats_df,
    file.path(output_dir, paste0(prefix, "_COG_bias_statistics.csv")),
    row.names = FALSE
  )

  pie_up_df <- make_pie_data(summary_df, "Upregulated")
  pie_down_df <- make_pie_data(summary_df, "Downregulated")

  plots <- list(
    COG_bar = make_bar_plot(summary_df, stats_df, comparison_label),
    COG_pie_up = make_pie_plot(pie_up_df, paste0("Upregulated Genes Functional Category (", comparison_label, ")")),
    COG_pie_down = make_pie_plot(pie_down_df, paste0("Downregulated Genes Functional Category (", comparison_label, ")"))
  )

  save_plot_set(prefix, plots)
}

run_comparison(
  sheet_name = "All_pH_5.0_vs_pH6.8_Functional_Categories",
  comparison_label = "pH 5.0 vs pH 6.8",
  prefix = "5.0"
)

run_comparison(
  sheet_name = "All_pH_9.0_vs_pH6.8_Functional_Categories",
  comparison_label = "pH 9.0 vs pH 6.8",
  prefix = "9.0"
)

cat("All plots and COG bias statistics saved successfully to:", output_dir, "\n")
