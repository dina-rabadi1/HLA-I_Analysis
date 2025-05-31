# File: diagnose_2cv_3cv_duplicates.R
# Purpose: Check for duplication issues in the 2CV/3CV combination process

library(dplyr)
library(tidyr)
library(readr)

# Set working directory to where your files are
setwd("~/Documents/Github/HLA-I_Analysis/20250516")

# Load the output file from the 2CV/3CV combination
output_file <- "unique_peptides_unmodified.tsv"  # Update with correct path if needed
peptide_data <- read_tsv(output_file)

# Check basics
cat("\n===== BASIC STATISTICS =====\n")
cat("Total rows in file:", nrow(peptide_data), "\n")
cat("Unique peptide sequences:", length(unique(peptide_data$Peptide)), "\n")

# Check if any peptide sequences appear multiple times
duplicate_count <- sum(duplicated(peptide_data$Peptide))
cat("Duplicate peptide sequences:", duplicate_count, "\n")

if (duplicate_count > 0) {
  # Get list of duplicated peptides
  duplicated_peptides <- peptide_data$Peptide[duplicated(peptide_data$Peptide)]
  
  # Look at the first few duplicates in detail
  cat("\n===== EXAMINING DUPLICATE EXAMPLES =====\n")
  for (i in 1:min(5, length(duplicated_peptides))) {
    peptide <- duplicated_peptides[i]
    cat("\nDuplicate peptide:", peptide, "\n")
    
    # Get all rows with this peptide
    peptide_rows <- peptide_data[peptide_data$Peptide == peptide, ]
    
    # Check what differs between the rows
    if (nrow(peptide_rows) > 1) {
      cat("Appears", nrow(peptide_rows), "times\n")
      
      # Find columns that differ between duplicates
      different_cols <- c()
      for (col in names(peptide_rows)) {
        values <- peptide_rows[[col]]
        if (length(unique(values)) > 1) {
          different_cols <- c(different_cols, col)
        }
      }
      
      cat("Columns with different values:", paste(different_cols, collapse=", "), "\n")
      
      # Display the values for differing columns
      if (length(different_cols) > 0) {
        cat("Values for differing columns:\n")
        print(peptide_rows[, c("Peptide", different_cols)])
      }
    }
  }
  
  # Check sample distribution
  if ("SampleID" %in% names(peptide_data)) {
    cat("\n===== CHECKING SAMPLE DISTRIBUTION OF DUPLICATES =====\n")
    
    # Get a larger sample of duplicates
    dup_sample <- unique(duplicated_peptides)[1:min(100, length(unique(duplicated_peptides)))]
    
    # Count how many duplicates appear in multiple samples
    multi_sample_count <- 0
    for (peptide in dup_sample) {
      rows <- peptide_data[peptide_data$Peptide == peptide, ]
      if ("SampleID" %in% names(rows)) {
        unique_samples <- unique(rows$SampleID)
        if (length(unique_samples) > 1) {
          multi_sample_count <- multi_sample_count + 1
        }
      }
    }
    
    cat("Of", length(dup_sample), "sampled duplicate peptides,", 
        multi_sample_count, "appear in multiple samples\n")
  }
  
  # Check for inconsistent gene/protein annotations
  cat("\n===== CHECKING ANNOTATION CONSISTENCY =====\n")
  if ("Gene" %in% names(peptide_data) && "Protein" %in% names(peptide_data)) {
    # Group by peptide and check for annotation consistency
    annotation_check <- peptide_data %>%
      group_by(Peptide) %>%
      summarize(
        Gene_Count = n_distinct(Gene, na.rm = TRUE),
        Protein_Count = n_distinct(Protein, na.rm = TRUE),
        .groups = 'drop'
      )
    
    multi_gene_count <- sum(annotation_check$Gene_Count > 1, na.rm = TRUE)
    multi_protein_count <- sum(annotation_check$Protein_Count > 1, na.rm = TRUE)
    
    cat("Peptides with multiple gene annotations:", multi_gene_count, "\n")
    cat("Peptides with multiple protein annotations:", multi_protein_count, "\n")
    
    if (multi_gene_count > 0 || multi_protein_count > 0) {
      cat("\nExample of peptides with multiple annotations:\n")
      example_peptides <- annotation_check %>%
        filter(Gene_Count > 1 | Protein_Count > 1) %>%
        head(3) %>%
        pull(Peptide)
      
      for (peptide in example_peptides) {
        cat("\nPeptide:", peptide, "\n")
        annotations <- peptide_data %>%
          filter(Peptide == peptide) %>%
          select(Peptide, Gene, Protein) %>%
          distinct()
        
        print(annotations)
      }
    }
  }
} else {
  cat("No duplicates found - each peptide sequence appears exactly once.\n")
}

# Check for duplicates by sample
if ("SampleID" %in% names(peptide_data)) {
  cat("\n===== CHECKING PEPTIDES PER SAMPLE =====\n")
  
  # Count peptides per sample
  peptides_per_sample <- peptide_data %>%
    group_by(SampleID) %>%
    summarize(Peptide_Count = n_distinct(Peptide), .groups = 'drop')
  
  cat("Peptide counts by sample:\n")
  print(peptides_per_sample)
  
  # Check peptide duplication within samples
  sample_duplicates <- peptide_data %>%
    group_by(SampleID) %>%
    summarize(
      Total_Rows = n(),
      Unique_Peptides = n_distinct(Peptide),
      Duplicates = Total_Rows - Unique_Peptides,
      .groups = 'drop'
    )
  
  any_sample_duplicates <- any(sample_duplicates$Duplicates > 0)
  
  if (any_sample_duplicates) {
    cat("\nSamples with duplicate peptides:\n")
    print(sample_duplicates %>% filter(Duplicates > 0))
  } else {
    cat("\nNo samples contain duplicate peptides.\n")
  }
}

cat("\nDiagnostic analysis complete!\n")