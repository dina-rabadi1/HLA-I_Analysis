#' Transcriptome Processing Functions
#' 
#' Functions for processing and fixing transcriptome data for immunopeptidomics analysis
#' Author: Generated for Dina Rabadi's HLA-I Analysis Pipeline
#' 

# Load required libraries
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(readxl)) install.packages("readxl")
if (!require(writexl)) install.packages("writexl")
if (!require(openxlsx)) install.packages("openxlsx")

library(tidyverse)
library(readxl)
library(writexl)
library(openxlsx)

#' Fix Transcriptome Data
#' 
#' Processes the normalized gene counts Excel file by:
#' 1. Filtering to only RU samples
#' 2. Recalculating Mean.Normal (RU samples ending in _N)
#' 3. Recalculating Mean.Tumor (all other RU samples)
#' 4. Calculating log2 fold change (log2(Tumor/Normal))
#' 5. Selecting only required output columns
#' 6. Saving processed data as Excel, CSV, and returning R object
#' 
#' @param input_file Path to input Excel file (default: auto-detected)
#' @param output_dir Directory to save processed files (default: auto-created)
#' @param verbose Print progress messages (default: TRUE)
#' 
#' @return Processed transcriptome data as tibble
#' 
#' @examples
#' processed_data <- fix_transcriptome()
#' processed_data <- fix_transcriptome(verbose = TRUE)
#' 
fix_transcriptome <- function(
    input_file = "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/directory/data/transcriptome/Normalized_Gene_counts_FLCdb_Panel_1.xlsx",
    output_dir = NULL,
    verbose = TRUE
) {
  
  if (verbose) cat("=== TRANSCRIPTOME PROCESSING STARTED ===\n")
  
  # Set up output directory
  if (is.null(output_dir)) {
    base_dir <- dirname(input_file)
    output_dir <- file.path(base_dir, "fix_transcriptome")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("Created output directory:", output_dir, "\n")
  }
  
  # Validate input file exists
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }
  
  if (verbose) cat("Reading input file:", basename(input_file), "\n")
  
  # Read the Excel file
  tryCatch({
    raw_data <- read_excel(input_file)
    if (verbose) cat("✓ Successfully loaded", nrow(raw_data), "genes\n")
  }, error = function(e) {
    stop("Error reading Excel file: ", e$message)
  })
  
  # Get all column names
  all_cols <- colnames(raw_data)
  if (verbose) cat("Total columns in dataset:", length(all_cols), "\n")
  
  # Identify RU sample columns
  ru_cols <- all_cols[grepl("^RU", all_cols)]
  if (verbose) cat("Found", length(ru_cols), "RU sample columns\n")
  
  if (length(ru_cols) == 0) {
    stop("No RU sample columns found in the dataset")
  }
  
  # Classify RU samples into Normal and Tumor
  normal_cols <- ru_cols[grepl("_N$", ru_cols)]
  tumor_cols <- ru_cols[!grepl("_N$", ru_cols)]
  
  if (verbose) {
    cat("Normal samples (ending in _N):", length(normal_cols), "\n")
    cat("Tumor samples (all other RU):", length(tumor_cols), "\n")
    cat("Sample breakdown:\n")
    cat("  Normal:", paste(head(normal_cols, 5), collapse = ", "))
    if (length(normal_cols) > 5) cat(", ...")
    cat("\n")
    cat("  Tumor:", paste(head(tumor_cols, 5), collapse = ", "))
    if (length(tumor_cols) > 5) cat(", ...")
    cat("\n")
  }
  
  if (length(normal_cols) == 0) {
    stop("No normal samples found (RU samples ending in _N)")
  }
  
  if (length(tumor_cols) == 0) {
    stop("No tumor samples found (RU samples not ending in _N)")
  }
  
  # Required metadata columns
  required_cols <- c("geneID", "symbol", "biotype", "chromosome", 
                     "gene_start", "gene_end", "gene_length", "description")
  
  # Check if required columns exist
  missing_cols <- required_cols[!required_cols %in% all_cols]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (verbose) cat("Processing gene expression data...\n")
  
  # Calculate new means and log2FC
  processed_data <- raw_data %>%
    # Select only the RU sample columns and required metadata columns
    select(all_of(c(ru_cols, required_cols))) %>%
    # Calculate mean expression for normal and tumor samples
    rowwise() %>%
    mutate(
      Mean.Normal = mean(c_across(all_of(normal_cols)), na.rm = TRUE),
      Mean.Tumor = mean(c_across(all_of(tumor_cols)), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # Calculate log2 fold change
    mutate(
      log2FC = log2(Mean.Tumor / Mean.Normal)
    ) %>%
    # Select only the final required columns
    select(Mean.Normal, Mean.Tumor, geneID, symbol, biotype, 
           chromosome, gene_start, gene_end, gene_length, description, log2FC) %>%
    # Arrange by geneID for consistency
    arrange(geneID)
  
  if (verbose) {
    cat("✓ Processed", nrow(processed_data), "genes successfully\n")
    cat("Expression summary:\n")
    cat("  Mean.Normal range:", round(min(processed_data$Mean.Normal, na.rm = TRUE), 2), 
        "to", round(max(processed_data$Mean.Normal, na.rm = TRUE), 2), "\n")
    cat("  Mean.Tumor range:", round(min(processed_data$Mean.Tumor, na.rm = TRUE), 2), 
        "to", round(max(processed_data$Mean.Tumor, na.rm = TRUE), 2), "\n")
    cat("  log2FC range:", round(min(processed_data$log2FC, na.rm = TRUE), 2), 
        "to", round(max(processed_data$log2FC, na.rm = TRUE), 2), "\n")
  }
  
  # Generate output filenames with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  base_filename <- "Processed_Normalized_Gene_counts_FLCdb_Panel_1"
  
  excel_file <- file.path(output_dir, paste0(timestamp, "_", base_filename, ".xlsx"))
  csv_file <- file.path(output_dir, paste0(timestamp, "_", base_filename, ".csv"))
  
  # Save as Excel file
  if (verbose) cat("Saving Excel file...\n")
  tryCatch({
    write_xlsx(processed_data, excel_file)
    if (verbose) cat("✓ Excel file saved:", basename(excel_file), "\n")
  }, error = function(e) {
    warning("Failed to save Excel file: ", e$message)
  })
  
  # Save as CSV file
  if (verbose) cat("Saving CSV file...\n")
  tryCatch({
    write_csv(processed_data, csv_file)
    if (verbose) cat("✓ CSV file saved:", basename(csv_file), "\n")
  }, error = function(e) {
    warning("Failed to save CSV file: ", e$message)
  })
  
  # Create a summary file with processing details
  summary_info <- data.frame(
    Metric = c(
      "Input File",
      "Processing Timestamp", 
      "Total Genes Processed",
      "Normal Samples Used",
      "Tumor Samples Used",
      "Output Directory",
      "Excel File",
      "CSV File"
    ),
    Value = c(
      basename(input_file),
      timestamp,
      nrow(processed_data),
      length(normal_cols),
      length(tumor_cols),
      output_dir,
      basename(excel_file),
      basename(csv_file)
    )
  )
  
  summary_file <- file.path(output_dir, paste0(timestamp, "_processing_summary.csv"))
  write_csv(summary_info, summary_file)
  
  if (verbose) {
    cat("\n=== TRANSCRIPTOME PROCESSING COMPLETED ===\n")
    cat("📁 Output directory:", output_dir, "\n")
    cat("📊 Files generated:\n")
    cat("   -", basename(excel_file), "\n")
    cat("   -", basename(csv_file), "\n")
    cat("   -", basename(summary_file), "\n")
    cat("📈 Data ready for immunopeptidomics cross-reference\n")
  }
  
  # Return the processed data for direct use in pipeline
  return(processed_data)
}

#' Get Transcriptome Summary Statistics
#' 
#' Generate summary statistics for the processed transcriptome data
#' 
#' @param transcriptome_data Processed transcriptome data (from fix_transcriptome)
#' @param verbose Print summary (default: TRUE)
#' 
#' @return Summary statistics as data frame
#' 
get_transcriptome_summary <- function(transcriptome_data, verbose = TRUE) {
  
  if (verbose) cat("Generating transcriptome summary statistics...\n")
  
  summary_stats <- data.frame(
    Metric = c(
      "Total Genes",
      "Genes with Valid Expression",
      "Mean Normal Expression",
      "Mean Tumor Expression", 
      "Median log2FC",
      "Upregulated Genes (log2FC > 1)",
      "Downregulated Genes (log2FC < -1)",
      "Unchanged Genes (|log2FC| <= 1)"
    ),
    Value = c(
      nrow(transcriptome_data),
      sum(!is.na(transcriptome_data$Mean.Normal) & !is.na(transcriptome_data$Mean.Tumor)),
      round(mean(transcriptome_data$Mean.Normal, na.rm = TRUE), 2),
      round(mean(transcriptome_data$Mean.Tumor, na.rm = TRUE), 2),
      round(median(transcriptome_data$log2FC, na.rm = TRUE), 2),
      sum(transcriptome_data$log2FC > 1, na.rm = TRUE),
      sum(transcriptome_data$log2FC < -1, na.rm = TRUE),
      sum(abs(transcriptome_data$log2FC) <= 1, na.rm = TRUE)
    )
  )
  
  if (verbose) {
    cat("Transcriptome Summary:\n")
    for(i in 1:nrow(summary_stats)) {
      cat(sprintf("  %-35s: %s\n", summary_stats$Metric[i], summary_stats$Value[i]))
    }
  }
  
  return(summary_stats)
}

#' Example usage:
#' 
#' # Process transcriptome data
#' processed_transcriptome <- fix_transcriptome()
#' 
#' # Get summary statistics  
#' summary_stats <- get_transcriptome_summary(processed_transcriptome)
#' 
#' # Use in immunopeptidomics pipeline
#' # The processed_transcriptome object is ready for cross-referencing with:
#' # - Immunopeptidome data (by geneID/symbol)
#' # - TMT proteome data
#' # - LFQ proteome data