# Load required packages
library(readxl)
library(readr)
library(tidyverse)

# Read files
cog_file <- read_tsv("cog_classify.tsv")

csv_file_pH5_vs_pH68 <-  read_csv("pH5_vs_pH68_TPM_log2TPM_by_condition.csv")

Joined_files <- left_join(csv_file_pH5_vs_pH68, cog_file, by = c("gene_name" = "QUERY_ID"))

# Filter significant genes
Data <- Joined_files %>%
  filter(padj < 0.001)

# Export as CSV file
write.csv(Data, "pH5_vs_pH68_Functiona_table_cog_classify.csv", row.names = FALSE)


# Read files
cog_file <- read_tsv("cog_classify.tsv")

csv_file_pH9_vs_pH68 <-  read_csv("pH9_vs_pH68_TPM_log2TPM_by_condition.csv")

Joined_files <- left_join(csv_file_pH9_vs_pH68, cog_file, by = c("gene_name" = "QUERY_ID"))

# Filter significant genes
Data <- Joined_files %>%
  filter(padj < 0.001)

# Export as CSV file
write.csv(Data, "pH9_vs_pH68_Functiona_table_cog_classify.csv", row.names = FALSE)

