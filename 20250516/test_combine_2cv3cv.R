# Fixed Script to Combine Peptide Data from 2cv and 3cv Files
# This version ensures one row per unique peptide, with all information from both methods
# Added support for filtering modified peptides

# Setting directory (customize this for your environment)
setwd("~/Documents/Github/HLA-I_Analysis/20250516")

# Load required packages
library(tidyverse)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)

# Function to check if a peptide has modifications
has_modifications <- function(assigned_mods, observed_mods) {
  # Handle NA values properly
  if (is.na(assigned_mods) && is.na(observed_mods)) {
    return(FALSE)
  }
  
  # Handle cases where one value is NA
  if (is.na(assigned_mods)) {
    return(!is.na(observed_mods) && observed_mods != "")
  }
  
  if (is.na(observed_mods)) {
    return(!is.na(assigned_mods) && assigned_mods != "")
  }
  
  # Otherwise check if either field has content
  return(assigned_mods != "" || observed_mods != "")
}

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

# Main function to process and create a unique-peptide dataset
combine_peptide_data <- function(directory, filter_modified = FALSE) {
  cat("Processing peptide files in:", directory, "\n")
  if (filter_modified) {
    cat("Filtering out peptides with modifications\n")
  }
  
  # Step 1: List and categorize TSV files
  file_lists <- list_tsv_files(directory)
  
  # Step 2: Read and process files
  data_2cv <- read_peptide_files(file_lists$files_2cv, "2cv")
  data_3cv <- read_peptide_files(file_lists$files_3cv, "3cv")
  
  cat("Read", nrow(data_2cv), "entries from 2cv files\n")
  cat("Read", nrow(data_3cv), "entries from 3cv files\n")
  
  # Step 2.5: If requested, filter out peptides with modifications
  if (filter_modified) {
    # Count peptides before filtering
    peptide_count_before_2cv <- length(unique(data_2cv$Peptide))
    peptide_count_before_3cv <- length(unique(data_3cv$Peptide))
    
    # Filter 2CV data
    data_2cv <- data_2cv %>%
      rowwise() %>%
      filter(!has_modifications(`Assigned Modifications`, `Observed Modifications`)) %>%
      ungroup()
    
    # Filter 3CV data
    data_3cv <- data_3cv %>%
      rowwise() %>%
      filter(!has_modifications(`Assigned Modifications`, `Observed Modifications`)) %>%
      ungroup()
    
    # Count peptides after filtering
    peptide_count_after_2cv <- length(unique(data_2cv$Peptide))
    peptide_count_after_3cv <- length(unique(data_3cv$Peptide))
    
    cat("Filtered out", peptide_count_before_2cv - peptide_count_after_2cv, 
        "modified peptides from 2cv data\n")
    cat("Filtered out", peptide_count_before_3cv - peptide_count_after_3cv, 
        "modified peptides from 3cv data\n")
  }
  
  # Step 3: Create a dataset with all unique peptides and their base information
  # Start by choosing the base dataset 
  # (preferring 2CV for simplicity, but could be either)
  all_peptides <- unique(c(data_2cv$Peptide, data_3cv$Peptide))
  cat("Found", length(all_peptides), "unique peptides in total\n")
  
  # Step 4: For each unique peptide, get the base information from 2CV or 3CV
  unique_rows <- list()
  
  for (peptide in all_peptides) {
    # Check if peptide is in 2CV data
    row_2cv <- which(data_2cv$Peptide == peptide)
    row_3cv <- which(data_3cv$Peptide == peptide)
    
    # For each unique peptide, gets the base information from 2CV or 3CV
    # If available in 2CV, use that; otherwise use 3CV
    base_row <- if (length(row_2cv) > 0) {
      # Use the first occurrence if multiple
      data_2cv[row_2cv[1], ]
    } else if (length(row_3cv) > 0) {
      # Use the first occurrence if multiple
      data_3cv[row_3cv[1], ]
    } else {
      # This shouldn't happen if the unique peptide list was created correctly
      next
    }
    
    # Now add the cross-reference columns
    # Use numeric indexing for intensity values
    base_row$Intensity_2cv <- if (length(row_2cv) > 0) data_2cv$Intensity[row_2cv[1]] else NA
    base_row$Intensity_3cv <- if (length(row_3cv) > 0) data_3cv$Intensity[row_3cv[1]] else NA
    
    # Calculate final intensity (average if both methods detected, otherwise use the available value)
    base_row$final_intensity <- case_when(
      length(row_2cv) > 0 && length(row_3cv) > 0 ~ (base_row$Intensity_2cv + base_row$Intensity_3cv) / 2,
      length(row_2cv) > 0 ~ base_row$Intensity_2cv,
      length(row_3cv) > 0 ~ base_row$Intensity_3cv,
      TRUE ~ NA_real_  # This should never happen if peptide detection is correct
    )
    
    # Add source file information
    base_row$SourceFile_2cv <- if (length(row_2cv) > 0) data_2cv$SourceFile[row_2cv[1]] else NA
    base_row$SourceFile_3cv <- if (length(row_3cv) > 0) data_3cv$SourceFile[row_3cv[1]] else NA
    
    # Add sample ID information
    base_row$SampleID_2cv <- if (length(row_2cv) > 0) data_2cv$SampleID[row_2cv[1]] else NA
    base_row$SampleID_3cv <- if (length(row_3cv) > 0) data_3cv$SampleID[row_3cv[1]] else NA
    
    # Add detection flags
    base_row$detected_2cv <- length(row_2cv) > 0
    base_row$detected_3cv <- length(row_3cv) > 0
    base_row$detected_both <- (length(row_2cv) > 0) && (length(row_3cv) > 0)
    
    # Add to the list of unique rows
    unique_rows[[length(unique_rows) + 1]] <- base_row
  }
  
  # Combine all unique rows
  combined_unique_data <- bind_rows(unique_rows)
  
  # Step 5: Calculate detection statistics
  detected_2cv_only <- sum(combined_unique_data$detected_2cv & !combined_unique_data$detected_3cv)
  detected_3cv_only <- sum(!combined_unique_data$detected_2cv & combined_unique_data$detected_3cv)
  detected_both <- sum(combined_unique_data$detected_both)
  
  cat("\nDetection statistics:\n")
  cat("Peptides detected in 2cv only:", detected_2cv_only, "\n")
  cat("Peptides detected in 3cv only:", detected_3cv_only, "\n")
  cat("Peptides detected in both datasets:", detected_both, "\n")
  cat("Total unique peptides:", nrow(combined_unique_data), "\n")
  
  # Step 6: Check if any peptides have modifications
  has_mods <- vapply(
    seq_len(nrow(combined_unique_data)),
    function(i) has_modifications(
      combined_unique_data$`Assigned Modifications`[i], 
      combined_unique_data$`Observed Modifications`[i]
    ),
    logical(1)
  )
  
  modified_count <- sum(has_mods)
  
  if (modified_count > 0) {
    cat("\nModification information:\n")
    cat("Peptides with modifications:", modified_count, "\n")
    
    if (filter_modified) {
      cat("Note: Modified peptides should have been filtered out. If any remain, please check the filter logic.\n")
    } else {
      cat("Use filter_modified=TRUE to exclude these peptides.\n")
    }
  } else if (filter_modified) {
    cat("\nAll modified peptides have been successfully filtered out.\n")
  }
  
  # Step 7: Print statistics about final_intensity values
  avg_intensities <- sum(combined_unique_data$detected_both)
  solo_2cv_intensities <- sum(combined_unique_data$detected_2cv & !combined_unique_data$detected_3cv)
  solo_3cv_intensities <- sum(!combined_unique_data$detected_2cv & combined_unique_data$detected_3cv)
  
  cat("\nFinal intensity statistics:\n")
  cat("Peptides with averaged intensity (from both methods):", avg_intensities, "\n")
  cat("Peptides with intensity from 2cv only:", solo_2cv_intensities, "\n")
  cat("Peptides with intensity from 3cv only:", solo_3cv_intensities, "\n")
  cat("Total peptides with intensity values:", avg_intensities + solo_2cv_intensities + solo_3cv_intensities, "\n")
  
  cat("\nUnique peptide dataset has", nrow(combined_unique_data), "entries\n")
  
  return(combined_unique_data)
}

# Main execution logic
main <- function() {
  # Parse command line arguments or set default values
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript peptide_combination.R <directory> [output_file] [filter_modified]\n")
    cat("Example: Rscript peptide_combination.R /path/to/tsv/files unique_peptides.tsv TRUE\n")
    return(invisible())
  }
  
  directory <- args[1]
  output_file <- if (length(args) >= 2) args[2] else "unique_peptides.tsv"
  filter_modified <- if (length(args) >= 3) as.logical(args[3]) else FALSE
  
  tryCatch({
    # Process and create unique peptide dataset
    combined_data <- combine_peptide_data(directory, filter_modified)
    
    # Save the output
    write_tsv(combined_data, output_file)
    cat("\nSuccessfully saved unique peptide data to", output_file, "\n")
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}

# Call the main function when the script is run directly
if (sys.nframe() == 0) {
  main()
}

# Example usage - Run with and without modification filtering
my_directory <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/realistic_test_data"  # Replace with your actual path

# Run the combination function WITHOUT filtering modified peptides
cat("\n--- CREATING DATASET WITH ALL PEPTIDES ---\n")
all_peptides <- combine_peptide_data(my_directory, filter_modified = FALSE)
write_tsv(all_peptides, "test_unique_peptides_all.tsv")

# Run the combination function WITH filtering of modified peptides
cat("\n--- CREATING DATASET WITHOUT MODIFIED PEPTIDES ---\n")
unmodified_peptides <- combine_peptide_data(my_directory, filter_modified = TRUE)
write_tsv(unmodified_peptides, "test_unique_peptides_unmodified.tsv")

# # Modified Script to Combine Peptide Data from 2cv and 3cv Files
# # This version ensures one row per unique peptide, with all information from both methods
# # Steps 3-5 are the ones that are different from the Concatenate 2CV3CV.R
# 
# # Setting directory
# setwd("~/Documents/Github/HLA-I_Analysis/20250516")
# 
# # Load required packages
# library(tidyverse)
# library(readr)
# library(stringr)
# library(dplyr)
# library(tidyr)
# 
# # Function to list and categorize TSV files
# list_tsv_files <- function(directory) {
#   # List all files in directory
#   all_files <- list.files(path = directory, pattern = ".*_peptides\\.tsv$", full.names = TRUE)
#   
#   if (length(all_files) == 0) {
#     stop("No peptide TSV files found in ", directory)
#   }
#   
#   # Categorize by CV type
#   files_2cv <- all_files[grepl("_2CV_", all_files)]
#   files_3cv <- all_files[grepl("_3CV_", all_files)]
#   
#   cat("Found", length(files_2cv), "2cv files and", length(files_3cv), "3cv files\n")
#   
#   return(list(files_2cv = files_2cv, files_3cv = files_3cv))
# }
# 
# # Function to extract sample ID from filename
# extract_sample_id <- function(filename) {
#   # Extract sample ID using regex
#   # Pattern like: DDA_2CV_117_peptides.tsv
#   filename <- basename(filename)
#   match <- str_extract(filename, "DDA_[23]CV_([^_]+)")
#   sample_id <- gsub("DDA_[23]CV_", "", match)
#   
#   if (is.na(sample_id)) {
#     return("unknown")
#   }
#   
#   return(sample_id)
# }
# 
# # Function to read peptide files
# read_peptide_files <- function(file_list, cv_type) {
#   if (length(file_list) == 0) {
#     return(NULL)
#   }
#   
#   all_data <- list()
#   
#   for (i in seq_along(file_list)) {
#     file <- file_list[i]
#     cat("Reading", cv_type, "file:", basename(file), "\n")
#     
#     # Read the TSV file
#     data <- read_tsv(file, show_col_types = FALSE)
#     
#     # Add metadata columns
#     data$SourceFile <- basename(file)
#     data$SampleID <- extract_sample_id(file)
#     data$CVType <- cv_type
#     
#     all_data[[i]] <- data
#   }
#   
#   # Combine all data frames
#   combined_data <- bind_rows(all_data)
#   
#   return(combined_data)
# }
# 
# # Main function to process and create a unique-peptide dataset
# combine_peptide_data <- function(directory) {
#   cat("Processing peptide files in:", directory, "\n")
#   
#   # Step 1: List and categorize TSV files
#   file_lists <- list_tsv_files(directory)
#   
#   # Step 2: Read and process files
#   data_2cv <- read_peptide_files(file_lists$files_2cv, "2cv")
#   data_3cv <- read_peptide_files(file_lists$files_3cv, "3cv")
#   
#   cat("Read", nrow(data_2cv), "entries from 2cv files\n")
#   cat("Read", nrow(data_3cv), "entries from 3cv files\n")
#   
#   # Step 3: Create a dataset with all unique peptides and their base information
#   # Start by choosing the base dataset 
#   # (preferring 2CV for simplicity, but could be either)
#   all_peptides <- unique(c(data_2cv$Peptide, data_3cv$Peptide))
#   cat("Found", length(all_peptides), "unique peptides in total\n")
#   
#   # Create a unified dataset with one row per unique peptide
#   unique_peptide_data <- tibble(
#     Peptide = all_peptides
#   )
#   
#   # Step 4: For each unique peptide, get the base information from 2CV or 3CV
#   # If available in 2CV, use that; otherwise use 3CV
#   unique_rows <- list()
#   
#   for (peptide in all_peptides) {
#     # Check if peptide is in 2CV data
#     row_2cv <- which(data_2cv$Peptide == peptide)
#     row_3cv <- which(data_3cv$Peptide == peptide)
#     
#     # Use 2CV data as base if available, otherwise use 3CV
#     if (length(row_2cv) > 0) {
#       # Use the first occurrence if multiple
#       base_row <- data_2cv[row_2cv[1], ]
#     } else if (length(row_3cv) > 0) {
#       # Use the first occurrence if multiple
#       base_row <- data_3cv[row_3cv[1], ]
#     } else {
#       # This shouldn't happen if the unique peptide list was created correctly
#       next
#     }
#     
#     # Now add the cross-reference columns
#     base_row$Intensity_2cv <- if (length(row_2cv) > 0) data_2cv$Intensity[row_2cv[1]] else NA
#     base_row$Intensity_3cv <- if (length(row_3cv) > 0) data_3cv$Intensity[row_3cv[1]] else NA
#     
#     # Calculate final intensity (average if both methods detected, otherwise use the available value)
#     base_row$final_intensity <- case_when(
#       length(row_2cv) > 0 && length(row_3cv) > 0 ~ (base_row$Intensity_2cv + base_row$Intensity_3cv) / 2,
#       length(row_2cv) > 0 ~ base_row$Intensity_2cv,
#       length(row_3cv) > 0 ~ base_row$Intensity_3cv,
#       TRUE ~ NA_real_  # This should never happen if peptide detection is correct
#     )
#     
#     # Add source file information
#     base_row$SourceFile_2cv <- if (length(row_2cv) > 0) data_2cv$SourceFile[row_2cv[1]] else NA
#     base_row$SourceFile_3cv <- if (length(row_3cv) > 0) data_3cv$SourceFile[row_3cv[1]] else NA
#     
#     # Add sample ID information
#     base_row$SampleID_2cv <- if (length(row_2cv) > 0) data_2cv$SampleID[row_2cv[1]] else NA
#     base_row$SampleID_3cv <- if (length(row_3cv) > 0) data_3cv$SampleID[row_3cv[1]] else NA
#     
#     # Add detection flags
#     base_row$detected_2cv <- length(row_2cv) > 0
#     base_row$detected_3cv <- length(row_3cv) > 0
#     base_row$detected_both <- (length(row_2cv) > 0) && (length(row_3cv) > 0)
#     
#     # Add to the list of unique rows
#     unique_rows[[length(unique_rows) + 1]] <- base_row
#   }
#   
#   # Combine all unique rows
#   combined_unique_data <- bind_rows(unique_rows)
#   
#   # Step 5: Calculate detection statistics
#   detected_2cv_only <- sum(combined_unique_data$detected_2cv & !combined_unique_data$detected_3cv)
#   detected_3cv_only <- sum(!combined_unique_data$detected_2cv & combined_unique_data$detected_3cv)
#   detected_both <- sum(combined_unique_data$detected_both)
#   
#   cat("\nDetection statistics:\n")
#   cat("Peptides detected in 2cv only:", detected_2cv_only, "\n")
#   cat("Peptides detected in 3cv only:", detected_3cv_only, "\n")
#   cat("Peptides detected in both datasets:", detected_both, "\n")
#   cat("Total unique peptides:", nrow(combined_unique_data), "\n")
#   
#   # Step 6: Print statistics about final_intensity values
#   avg_intensities <- sum(combined_unique_data$detected_both)
#   solo_2cv_intensities <- sum(combined_unique_data$detected_2cv & !combined_unique_data$detected_3cv)
#   solo_3cv_intensities <- sum(!combined_unique_data$detected_2cv & combined_unique_data$detected_3cv)
#   
#   cat("\nFinal intensity statistics:\n")
#   cat("Peptides with averaged intensity (from both methods):", avg_intensities, "\n")
#   cat("Peptides with intensity from 2cv only:", solo_2cv_intensities, "\n")
#   cat("Peptides with intensity from 3cv only:", solo_3cv_intensities, "\n")
#   cat("Total peptides with intensity values:", avg_intensities + solo_2cv_intensities + solo_3cv_intensities, "\n")
#   
#   return(combined_unique_data)
# }
# 
# # Main execution logic
# main <- function() {
#   # Parse command line arguments or set default values
#   args <- commandArgs(trailingOnly = TRUE)
#   
#   if (length(args) < 1) {
#     cat("Usage: Rscript peptide_combination.R <directory> [output_file]\n")
#     cat("Example: Rscript peptide_combination.R /path/to/tsv/files unique_peptides.tsv\n")
#     return(invisible())
#   }
#   
#   directory <- args[1]
#   output_file <- if (length(args) >= 2) args[2] else "unique_peptides.tsv"
#   
#   tryCatch({
#     # Process and create unique peptide dataset
#     combined_data <- combine_peptide_data(directory)
#     
#     # Save the output
#     write_tsv(combined_data, output_file)
#     cat("\nSuccessfully saved unique peptide data to", output_file, "\n")
#     
#   }, error = function(e) {
#     cat("Error:", e$message, "\n")
#   })
# }
# 
# # Call the main function when the script is run directly
# if (sys.nframe() == 0) {
#   main()
# }
# 
# # Set the path to your data directory
# my_directory <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/realistic_test_data"  # Replace with your actual path
# 
# # Run the combination function
# unique_peptides <- combine_peptide_data(my_directory)
# 
# # Save the result to a file
# write_tsv(unique_peptides, "test_unique_peptides.tsv")