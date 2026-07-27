library(dplyr)
library(readODS)

input_file <- "/data/Peat_soil/PacBio_peat/transcriptomics/results/DEGs/Functional_Categories.ods"
output_dir <- "/proj/Peat_soil/PacBio_peat/transcriptomics/results/figures/COG_FUN_CAT"

run_within_category_test <- function(sheet_name, prefix) {
  data <- read_ods(input_file, sheet = sheet_name) %>%
    mutate(
      Regulation = case_when(
        log2FoldChange > 0 ~ "Upregulated",
        log2FoldChange < 0 ~ "Downregulated",
        TRUE ~ "Neutral"
      )
    ) %>%
    filter(Regulation != "Neutral")

  counts <- data %>%
    group_by(COG_LETTER, COG_DESCRIPTION, Regulation) %>%
    summarise(Count = n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = Regulation,
      values_from = Count,
      values_fill = 0
    ) %>%
    rename(Up = Upregulated, Down = Downregulated) %>%
    mutate(Total = Up + Down)

  results <- lapply(seq_len(nrow(counts)), function(i) {
    row <- counts[i, ]
    test <- binom.test(row$Up, row$Total, p = 0.5)
    direction <- if (row$Up > row$Down) {
      "Upregulated"
    } else if (row$Up < row$Down) {
      "Downregulated"
    } else {
      "Balanced"
    }

    data.frame(
      COG_LETTER = row$COG_LETTER,
      COG_DESCRIPTION = row$COG_DESCRIPTION,
      Up = row$Up,
      Down = row$Down,
      Total = row$Total,
      ProportionUp = row$Up / row$Total,
      PValue = test$p.value,
      Direction = direction
    )
  })

  results_df <- bind_rows(results) %>%
    mutate(
      FDR = p.adjust(PValue, method = "BH"),
      Significant = FDR < 0.05
    ) %>%
    arrange(COG_LETTER)

  write.csv(
    results_df,
    file.path(output_dir, paste0(prefix, "_COG_within_category_binom_statistics.csv")),
    row.names = FALSE
  )
}

run_within_category_test(
  sheet_name = "All_pH_5.0_vs_pH6.8_Functional_Categories",
  prefix = "5.0"
)

run_within_category_test(
  sheet_name = "All_pH_9.0_vs_pH6.8_Functional_Categories",
  prefix = "9.0"
)
