# Simple script to analyze peptide sample distribution
# No complex formatting or indentation

# Setting directory
setwd("~/Documents/Github/HLA-I_Analysis/")

# Simple script to analyze peptide sample distribution
# No complex formatting or indentation

# Path to the combined peptides file
combined_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/combined_peptides.tsv"

# Load required libraries
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library(tidyverse)
}

# Read the combined peptides file
cat("Reading combined peptides file...\n")
peptides_data <- read_tsv(combined_file, show_col_types = FALSE)

# Filter to include only 8-12mers
cat("Filtering to include only 8-12mers...\n")
peptides_8_12 <- peptides_data %>%
  filter(`Peptide Length` >= 8 & `Peptide Length` <= 12)

cat("Number of 8-12mer peptides:", n_distinct(peptides_8_12$Peptide), "\n")

# Create a report showing which samples each peptide appears in
# Now including gene and protein information
peptide_report <- peptides_8_12 %>%
  # First get basic peptide and sample information
  select(Peptide, SampleID, `Peptide Length`, Gene, Protein) %>%
  distinct() %>%
  # Group by peptide to consolidate sample information
  group_by(Peptide, `Peptide Length`) %>%
  summarize(
    sample_list = paste(sort(SampleID), collapse=", "),
    sample_count = n_distinct(SampleID),
    # Consolidate gene information
    genes = paste(unique(na.omit(Gene)), collapse="; "),
    # Consolidate protein information
    proteins = paste(unique(na.omit(Protein)), collapse="; "),
    .groups = "drop"
  ) %>%
  arrange(desc(sample_count), Peptide)

# Write the report to a file
output_dir <- dirname(combined_file)
output_report <- paste0(output_dir, "/peptide_sample_distribution.tsv")
write_tsv(peptide_report, output_report)
cat("Peptide sample distribution report written to:", output_report, "\n")

# Print summary of the report
cat("\nSample count distribution:\n")
print(table(peptide_report$sample_count))

# Function to report samples for a specific peptide
get_peptide_samples <- function(peptide_seq) {
  peptide_info <- peptide_report %>% filter(Peptide == peptide_seq)
  
  if(nrow(peptide_info) == 0) {
    cat("Peptide", peptide_seq, "not found in dataset (or not an 8-12mer).\n")
    return(NULL)
  }
  
  cat("\nPeptide:", peptide_seq, "\n")
  cat("Length:", peptide_info$`Peptide Length`, "amino acids\n")
  cat("Found in", peptide_info$sample_count, "samples:", peptide_info$sample_list, "\n")
  
  # Get more details about this peptide
  peptide_details <- peptides_8_12 %>%
    filter(Peptide == peptide_seq) %>%
    select(Peptide, SampleID, Gene, Protein, SourceFile) %>%
    distinct()
  
  cat("\nAssociated genes and proteins:\n")
  genes_proteins <- peptide_details %>%
    select(Gene, Protein) %>%
    distinct()
  print(genes_proteins)
  
  return(peptide_info)
}

# Check specific peptide as example
get_peptide_samples("AAGPPISEGKY")

cat("\nTo check other peptides, use the get_peptide_samples() function.\n")
cat("Analysis complete.\n")