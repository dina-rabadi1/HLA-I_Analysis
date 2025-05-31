# Immunopeptidomics Analysis Script
# This script:
# 1. Identifies peptides shared across a majority of samples (excluding 51S and 148N)
# 2. Compares those peptides to transcriptome, TMT proteome, and LFQ proteome data

setwd("20250516")

# Load required libraries
library(tidyverse)
library(readxl)
library(writexl)
library(openxlsx)

# Set parameters (easily adjustable)
sample_threshold_percent <- 50  # Percentage of samples a peptide must be in to be considered "shared"
log2fc_threshold <- 2          # Log2FC threshold for considering genes/proteins as upregulated

# Create output directory based on threshold
output_dir <- paste0("immunopeptidomics_analysis_", sample_threshold_percent, "pct_", log2fc_threshold, "FC")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to log messages with timestamps
log_message <- function(message) {
  cat(paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), message, "\n"))
}

log_message("Starting immunopeptidomics analysis")
log_message(paste0("Parameters: Sample threshold = ", sample_threshold_percent, "%, Log2FC threshold = ", log2fc_threshold))

# Step 1: Load immunopeptidome data
log_message("Loading immunopeptidome data")
immuno_file <- "data/combined_peptides.tsv"
immuno_data <- read_tsv(immuno_file, show_col_types = FALSE)

# Check the column names to ensure we have what we expect
log_message(paste0("Immunopeptidome data columns: ", paste(names(immuno_data), collapse = ", ")))

# Display sample IDs
sample_ids <- unique(immuno_data$SampleID)
log_message(paste0("Found ", length(sample_ids), " samples: ", paste(sample_ids, collapse = ", ")))

# Exclude specified samples
exclude_samples <- c("51S", "148N")
filtered_samples <- setdiff(sample_ids, exclude_samples)
log_message(paste0("Excluding samples: ", paste(exclude_samples, collapse = ", ")))
log_message(paste0("Remaining ", length(filtered_samples), " samples: ", paste(filtered_samples, collapse = ", ")))

# Filter the data to exclude specified samples
filtered_immuno_data <- immuno_data %>%
  filter(SampleID %in% filtered_samples)

# Step 2: Identify peptides shared across majority of samples
log_message("Identifying peptides shared across majority of samples")

# Count in how many samples each peptide appears
peptide_counts <- filtered_immuno_data %>%
  group_by(Peptide, Gene) %>%
  summarize(
    num_samples = n_distinct(SampleID),
    sample_list = paste(unique(SampleID), collapse = ","),
    .groups = "drop"
  )

# Calculate the threshold number of samples
threshold_samples <- ceiling(length(filtered_samples) * sample_threshold_percent / 100)
log_message(paste0("Sample threshold (", sample_threshold_percent, "%) = ", threshold_samples, " samples"))

# Get peptides that meet the threshold
shared_peptides <- peptide_counts %>%
  filter(num_samples >= threshold_samples)

log_message(paste0("Found ", nrow(shared_peptides), " peptides present in at least ", 
                   threshold_samples, " samples (", sample_threshold_percent, "%)"))

# Step A: Write shared peptides to file
write_csv(shared_peptides, file.path(output_dir, "shared_peptides.csv"))

# Step 3: Load transcriptome data
log_message("Loading transcriptome data")
transcriptome_file <- "data/SupplementaryData3_Requena_2024.xlsx"
transcriptome_data <- read_excel(transcriptome_file)

# Check the column names
log_message(paste0("Transcriptome data columns: ", paste(names(transcriptome_data), collapse = ", ")))

# Filter transcriptome data for upregulated genes
upregulated_transcriptome <- transcriptome_data %>%
  filter(Ave.log2FC > log2fc_threshold)

log_message(paste0("Found ", nrow(upregulated_transcriptome), 
                   " upregulated genes in transcriptome (Log2FC > ", log2fc_threshold, ")"))

# Step 4: Load TMT proteome data
log_message("Loading TMT proteome data")
tmt_file <- "data/adg7038_Table_S1_TMT.xlsx"
tmt_data <- read_excel(tmt_file, sheet = "Significant and 1.5x_2")

# Check the column names
log_message(paste0("TMT proteome data columns: ", paste(names(tmt_data), collapse = ", ")))

# Filter TMT data for upregulated proteins
upregulated_tmt <- tmt_data %>%
  filter(Log2_Difference > log2fc_threshold)

log_message(paste0("Found ", nrow(upregulated_tmt), 
                   " upregulated proteins in TMT proteome (Log2FC > ", log2fc_threshold, ")"))

# Step 5: Load LFQ proteome data
log_message("Loading LFQ proteome data")
lfq_file <- "data/adg7038_Table_S2_LFQ.xlsx"
lfq_data <- read_excel(lfq_file, sheet = "Significant and 1.5x_2")

# Check the column names
log_message(paste0("LFQ proteome data columns: ", paste(names(lfq_data), collapse = ", ")))

# Filter LFQ data for upregulated proteins
upregulated_lfq <- lfq_data %>%
  filter(Log2_Difference > log2fc_threshold)

log_message(paste0("Found ", nrow(upregulated_lfq), 
                   " upregulated proteins in LFQ proteome (Log2FC > ", log2fc_threshold, ")"))

# Step 6: Compare the shared peptides with upregulated genes/proteins
log_message("Comparing shared peptides with upregulated genes/proteins")

# Compare with transcriptome
log_message("Comparing with transcriptome")
transcriptome_matches <- shared_peptides %>%
  inner_join(upregulated_transcriptome, by = c("Gene" = "Gene_Symbol")) %>%
  select(Peptide, Gene, num_samples, sample_list, Ave.log2FC, Ave.FDR, Gene_Description)

log_message(paste0("Found ", nrow(transcriptome_matches), 
                   " peptides matching upregulated genes in transcriptome"))

# Compare with TMT proteome
log_message("Comparing with TMT proteome")
tmt_matches <- shared_peptides %>%
  inner_join(upregulated_tmt, by = c("Gene" = "Gene Name")) %>%
  select(Peptide, Gene, num_samples, sample_list, Log2_Difference, P.value, `Protein Name`)

log_message(paste0("Found ", nrow(tmt_matches), 
                   " peptides matching upregulated proteins in TMT proteome"))

# Compare with LFQ proteome
log_message("Comparing with LFQ proteome")
lfq_matches <- shared_peptides %>%
  inner_join(upregulated_lfq, by = c("Gene" = "Gene Name")) %>%
  select(Peptide, Gene, num_samples, sample_list, Log2_Difference, P.value, `Protein Name`)

log_message(paste0("Found ", nrow(lfq_matches), 
                   " peptides matching upregulated proteins in LFQ proteome"))

# Step 7: Find peptides present in all datasets (intersection)
log_message("Finding peptides present in all upregulated datasets")

# Get genes in common
genes_in_transcriptome <- transcriptome_matches$Gene
genes_in_tmt <- tmt_matches$Gene
genes_in_lfq <- lfq_matches$Gene

# Find genes in common across all datasets
genes_in_all_datasets <- Reduce(intersect, list(genes_in_transcriptome, genes_in_tmt, genes_in_lfq))

log_message(paste0("Found ", length(genes_in_all_datasets), 
                   " genes upregulated in all datasets (transcriptome, TMT and LFQ)"))

# Find peptides from those genes
peptides_in_all_datasets <- shared_peptides %>%
  filter(Gene %in% genes_in_all_datasets)

log_message(paste0("Found ", nrow(peptides_in_all_datasets), 
                   " peptides from genes upregulated in all datasets"))

# Step 8: Create Venn diagram data for overlaps between datasets
genes_only_in_transcriptome <- setdiff(genes_in_transcriptome, c(genes_in_tmt, genes_in_lfq))
genes_only_in_tmt <- setdiff(genes_in_tmt, c(genes_in_transcriptome, genes_in_lfq))
genes_only_in_lfq <- setdiff(genes_in_lfq, c(genes_in_transcriptome, genes_in_tmt))

genes_in_transcriptome_and_tmt <- intersect(genes_in_transcriptome, genes_in_tmt) %>%
  setdiff(genes_in_lfq)
genes_in_transcriptome_and_lfq <- intersect(genes_in_transcriptome, genes_in_lfq) %>%
  setdiff(genes_in_tmt)
genes_in_tmt_and_lfq <- intersect(genes_in_tmt, genes_in_lfq) %>%
  setdiff(genes_in_transcriptome)

venn_data <- data.frame(
  Category = c(
    "Only in Transcriptome", 
    "Only in TMT", 
    "Only in LFQ",
    "In Transcriptome and TMT",
    "In Transcriptome and LFQ",
    "In TMT and LFQ",
    "In All Datasets"
  ),
  Count = c(
    length(genes_only_in_transcriptome),
    length(genes_only_in_tmt),
    length(genes_only_in_lfq),
    length(genes_in_transcriptome_and_tmt),
    length(genes_in_transcriptome_and_lfq),
    length(genes_in_tmt_and_lfq),
    length(genes_in_all_datasets)
  ),
  Genes = c(
    paste(genes_only_in_transcriptome, collapse = ", "),
    paste(genes_only_in_tmt, collapse = ", "),
    paste(genes_only_in_lfq, collapse = ", "),
    paste(genes_in_transcriptome_and_tmt, collapse = ", "),
    paste(genes_in_transcriptome_and_lfq, collapse = ", "),
    paste(genes_in_tmt_and_lfq, collapse = ", "),
    paste(genes_in_all_datasets, collapse = ", ")
  )
)

# Step 9: Save all results to files
log_message("Saving results to files")

# Save individual result files
write_csv(transcriptome_matches, file.path(output_dir, "transcriptome_matches.csv"))
write_csv(tmt_matches, file.path(output_dir, "tmt_matches.csv"))
write_csv(lfq_matches, file.path(output_dir, "lfq_matches.csv"))
write_csv(peptides_in_all_datasets, file.path(output_dir, "peptides_in_all_datasets.csv"))
write_csv(venn_data, file.path(output_dir, "venn_diagram_data.csv"))

# Create a combined Excel file with multiple sheets
log_message("Creating combined Excel report")

wb <- createWorkbook()

# Add sheets for each result
addWorksheet(wb, "Shared Peptides")
writeData(wb, "Shared Peptides", shared_peptides)

addWorksheet(wb, "Transcriptome Matches")
writeData(wb, "Transcriptome Matches", transcriptome_matches)

addWorksheet(wb, "TMT Matches")
writeData(wb, "TMT Matches", tmt_matches)

addWorksheet(wb, "LFQ Matches")
writeData(wb, "LFQ Matches", lfq_matches)

addWorksheet(wb, "Peptides in All Datasets")
writeData(wb, "Peptides in All Datasets", peptides_in_all_datasets)

addWorksheet(wb, "Venn Diagram Data")
writeData(wb, "Venn Diagram Data", venn_data)

# Add summary sheet
summary_data <- data.frame(
  Metric = c(
    "Number of Total Samples",
    "Number of Filtered Samples",
    "Sample Threshold Percentage",
    "Minimum Number of Samples Required",
    "Log2FC Threshold",
    "Number of Shared Peptides",
    "Number of Upregulated Genes in Transcriptome",
    "Number of Upregulated Proteins in TMT",
    "Number of Upregulated Proteins in LFQ",
    "Number of Peptides Matching Transcriptome",
    "Number of Peptides Matching TMT",
    "Number of Peptides Matching LFQ",
    "Number of Genes in All Datasets",
    "Number of Peptides from Genes in All Datasets"
  ),
  Value = c(
    length(sample_ids),
    length(filtered_samples),
    sample_threshold_percent,
    threshold_samples,
    log2fc_threshold,
    nrow(shared_peptides),
    nrow(upregulated_transcriptome),
    nrow(upregulated_tmt),
    nrow(upregulated_lfq),
    nrow(transcriptome_matches),
    nrow(tmt_matches),
    nrow(lfq_matches),
    length(genes_in_all_datasets),
    nrow(peptides_in_all_datasets)
  )
)

addWorksheet(wb, "Summary", gridLines = TRUE)
writeData(wb, "Summary", summary_data)

# Save the workbook
saveWorkbook(wb, file.path(output_dir, "immunopeptidomics_analysis_results.xlsx"), overwrite = TRUE)

log_message(paste0("Analysis complete. Results saved to ", output_dir, " directory"))

# Print summary to console
cat("\n=== ANALYSIS SUMMARY ===\n")
print(summary_data)
cat("=======================\n")