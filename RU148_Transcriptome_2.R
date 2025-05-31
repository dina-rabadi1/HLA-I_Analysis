# Load required libraries
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Set working directory to the raw data folder
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata")

# Read the Excel file
gene_data <- read_excel("Normalized_Gene_counts_FLCdb_Panel_1.xlsx")

# Select columns of interest
selected_data <- gene_data %>%
  select(RU148_N, RU148_T8, RU148_T11, Mean.Normal, Mean.Tumor, 
         geneID, symbol, biotype, chromosome, gene_start, gene_end, gene_length, description)

# Calculate log2 fold change for each tumor sample compared to normal
selected_data <- selected_data %>%
  mutate(
    log2FC_T8 = log2((RU148_T8 + 0.1) / (RU148_N + 0.1)),
    log2FC_T11 = log2((RU148_T11 + 0.1) / (RU148_N + 0.1)),
    mean_log2FC = (log2((RU148_T8 + 0.1) / (RU148_N + 0.1)) + 
                     log2((RU148_T11 + 0.1) / (RU148_N + 0.1))) / 2
  )

# Create a dataset using criteria that we know should yield results based on your diagnostics
# We'll create three different files with different filtering approaches

# Option 1: Using individual sample criteria that we know will yield genes
upregulated_individual <- selected_data %>%
  filter(
    # Require log2FC > 1 in at least one tumor sample
    (log2FC_T8 > 1 | log2FC_T11 > 1),
    # Ensure expression level > 10 in the corresponding sample
    ((log2FC_T8 > 1 & RU148_T8 > 10) | (log2FC_T11 > 1 & RU148_T11 > 10))
  ) %>%
  arrange(desc(mean_log2FC))

# Option 2: Focus on genes up in both tumor samples compared to normal
upregulated_both_samples <- selected_data %>%
  filter(
    # Up in both samples (less stringent threshold)
    log2FC_T8 > 0.5 & log2FC_T11 > 0.5,
    # Reasonable expression in both
    RU148_T8 > 5 & RU148_T11 > 5
  ) %>%
  arrange(desc(mean_log2FC))

# Option 3: Broader approach focusing on the sample with stronger signal
upregulated_any_sample <- selected_data %>%
  filter(
    # Either sample shows upregulation with higher expression
    (log2FC_T8 > 1 & RU148_T8 > 10) | (log2FC_T11 > 1 & RU148_T11 > 10)
  ) %>%
  # Add metadata to identify which sample showed the stronger effect
  mutate(
    stronger_in_T8 = ifelse(log2FC_T8 > log2FC_T11, "Yes", "No"),
    significantly_upregulated = case_when(
      log2FC_T8 > 1 & log2FC_T11 > 1 ~ "Both samples",
      log2FC_T8 > 1 ~ "T8 only",
      log2FC_T11 > 1 ~ "T11 only",
      TRUE ~ "Neither"
    )
  ) %>%
  arrange(desc(mean_log2FC))

# Print counts of genes found
cat("Number of genes upregulated in either sample (Option 1):", nrow(upregulated_individual), "\n")
cat("Number of genes upregulated in both samples (Option 2):", nrow(upregulated_both_samples), "\n")
cat("Number of genes using broader approach (Option 3):", nrow(upregulated_any_sample), "\n")

# Write these results to Excel files
write_xlsx(upregulated_individual, "Upregulated_genes_either_tumor.xlsx")
write_xlsx(upregulated_both_samples, "Upregulated_genes_both_tumors.xlsx")
write_xlsx(upregulated_any_sample, "Upregulated_genes_comprehensive.xlsx")

# Additional analysis: overlap between the two tumor samples
sample_comparison <- selected_data %>%
  filter(
    (log2FC_T8 > 1 & RU148_T8 > 10) | (log2FC_T11 > 1 & RU148_T11 > 10)
  ) %>%
  mutate(
    up_in_T8 = log2FC_T8 > 1 & RU148_T8 > 10,
    up_in_T11 = log2FC_T11 > 1 & RU148_T11 > 10,
    up_in_both = up_in_T8 & up_in_T11,
    pattern = case_when(
      up_in_both ~ "Both tumors",
      up_in_T8 ~ "T8 only",
      up_in_T11 ~ "T11 only",
      TRUE ~ "Neither"
    )
  )

# Write this to a file with the pattern information
write_xlsx(sample_comparison, "Tumor_sample_comparison.xlsx")

# Create summary counts for each pattern
pattern_summary <- sample_comparison %>%
  group_by(pattern) %>%
  summarize(count = n()) %>%
  arrange(desc(count))

write_xlsx(pattern_summary, "Pattern_summary.xlsx")

# Create visualization showing the overlap
venn_data <- c(
  sum(sample_comparison$up_in_T8 & !sample_comparison$up_in_T11),  # T8 only
  sum(sample_comparison$up_in_both),  # Both
  sum(!sample_comparison$up_in_T8 & sample_comparison$up_in_T11)   # T11 only
)

# Plot as Venn diagram or as a bar chart
pdf("Upregulated_genes_overlap.pdf")
barplot(venn_data, 
        names.arg = c("T8 only", "Both tumors", "T11 only"),
        col = c("lightblue", "purple", "pink"),
        main = "Upregulated Genes in Tumor Samples",
        ylab = "Number of genes")
dev.off()

cat("Analysis complete. Multiple files created with different filtering approaches.")