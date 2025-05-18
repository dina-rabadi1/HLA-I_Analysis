# Script to concatenate peptide data from 2cv and 3cv files
# Creates a comprehensive dataset with all peptides identified in both,
# while preserving all properties from original data and adding
# three columns: detected_2cv, detected_3cv, detected_both

# Setting directory
setwd("~/Documents/Github/HLA-I_Analysis/")

# Load required packages
library(tidyverse)
library(readr)
library(stringr)

# Function to list and categorize TSV files
list_tsv_files <- function(directory) {
  # List all files in directory
  all_files <- list.files(path = directory, pattern = ".*_peptides\\.tsv$", full.names = TRUE)
  
  if (length(all_files) == 0) {
    stop("No peptide TSV files found in ", directory)
  }
  
  # Categorize by CV type
  files_2cv <- all_files[grepl("_2CV_", all_files)]
  files_3cv <- all_files[grepl("_3CV_", all_files)]
  
  cat("Found", length(files_2cv), "2cv files and", length(files_3cv), "3cv files\n")
  
  return(list(files_2cv = files_2cv, files_3cv = files_3cv))
}

# Function to extract sample ID from filename
extract_sample_id <- function(filename) {
  # Extract sample ID using regex
  # Pattern like: DDA_2CV_117_peptides.tsv
  filename <- basename(filename)
  match <- str_extract(filename, "DDA_[23]CV_([^_]+)")
  sample_id <- gsub("DDA_[23]CV_", "", match)
  
  if (is.na(sample_id)) {
    return("unknown")
  }
  
  return(sample_id)
}

# Function to read peptide files
read_peptide_files <- function(file_list, cv_type) {
  if (length(file_list) == 0) {
    return(NULL)
  }
  
  all_data <- list()
  
  for (i in seq_along(file_list)) {
    file <- file_list[i]
    cat("Reading", cv_type, "file:", basename(file), "\n")
    
    # Read the TSV file
    data <- read_tsv(file, show_col_types = FALSE)
    
    # Add metadata columns
    data$SourceFile <- basename(file)
    data$SampleID <- extract_sample_id(file)
    data$CVType <- cv_type
    
    all_data[[i]] <- data
  }
  
  # Combine all data frames
  combined_data <- bind_rows(all_data)
  
  return(combined_data)
}

# Main function to process and concatenate peptide data
concatenate_peptide_data <- function(directory) {
  cat("Processing peptide files in:", directory, "\n")
  
  # Step 1: List and categorize TSV files
  file_lists <- list_tsv_files(directory)
  
  # Step 2: Read and process files
  data_2cv <- read_peptide_files(file_lists$files_2cv, "2cv")
  data_3cv <- read_peptide_files(file_lists$files_3cv, "3cv")
  
  cat("Read", nrow(data_2cv), "entries from 2cv files\n")
  cat("Read", nrow(data_3cv), "entries from 3cv files\n")
  
  # Step 3: Get all unique peptides from both datasets
  all_peptides <- unique(c(data_2cv$Peptide, data_3cv$Peptide))
  cat("Found", length(all_peptides), "unique peptides in total\n")
  
  # Step 4: Create a detection map for all peptides
  peptide_detection <- tibble(
    Peptide = all_peptides,
    detected_2cv = Peptide %in% data_2cv$Peptide,
    detected_3cv = Peptide %in% data_3cv$Peptide,
    detected_both = (Peptide %in% data_2cv$Peptide) & (Peptide %in% data_3cv$Peptide)
  )
  
  # Step 5: Create combined dataset
  # Start with 2cv data
  if (!is.null(data_2cv) && nrow(data_2cv) > 0) {
    # Add detection flags
    data_2cv <- data_2cv %>%
      left_join(peptide_detection, by = "Peptide")
  } else {
    data_2cv <- tibble() # Empty tibble
  }
  
  # Then 3cv data
  if (!is.null(data_3cv) && nrow(data_3cv) > 0) {
    # Add detection flags
    data_3cv <- data_3cv %>%
      left_join(peptide_detection, by = "Peptide")
  } else {
    data_3cv <- tibble() # Empty tibble
  }
  
  # Combine both datasets
  combined_data <- bind_rows(data_2cv, data_3cv)
  
  # Step 6: Calculate detection statistics
  detected_2cv_only <- sum(peptide_detection$detected_2cv & !peptide_detection$detected_3cv)
  detected_3cv_only <- sum(!peptide_detection$detected_2cv & peptide_detection$detected_3cv)
  detected_both <- sum(peptide_detection$detected_both)
  
  cat("\nDetection statistics:\n")
  cat("Peptides detected in 2cv only:", detected_2cv_only, "\n")
  cat("Peptides detected in 3cv only:", detected_3cv_only, "\n")
  cat("Peptides detected in both datasets:", detected_both, "\n")
  
  cat("\nCombined dataset has", nrow(combined_data), "entries\n")
  
  return(combined_data)
}

# Main execution logic
main <- function() {
  # Parse command line arguments or set default values
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript peptide_concatenation.R <directory> [output_file]\n")
    cat("Example: Rscript peptide_concatenation.R /path/to/tsv/files combined_peptides.tsv\n")
    return(invisible())
  }
  
  directory <- args[1]
  output_file <- if (length(args) >= 2) args[2] else "combined_peptides.tsv"
  
  tryCatch({
    # Process and concatenate peptide data
    combined_data <- concatenate_peptide_data(directory)
    
    # Save the output
    write_tsv(combined_data, output_file)
    cat("\nSuccessfully saved combined data to", output_file, "s\n")
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}

# Call the main function when the script is run directly
if (sys.nframe() == 0) {
  main()
}

# Set the path to your data directory
my_directory <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/imp_014_rawdata"  # Replace with your actual path

# Run the concatenation function
combined_data <- concatenate_peptide_data(my_directory)

# Save the result to a file
write_tsv(combined_data, "combined_peptides.tsv")