# Core utility functions for peptide analysis
# src/core/peptide_utils.R

#' Load and preprocess peptide data
#' 
#' @param file_path Path to the peptide TSV file
#' @param min_length Minimum peptide length to include
#' @param max_length Maximum peptide length to include
#' @return A preprocessed data frame of peptide data
load_peptide_data <- function(file_path, min_length = 8, max_length = 12) {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("Peptide data file not found: ", file_path)
  }
  
  # Read the data
  cat("Reading peptide data from:", file_path, "\n")
  peptide_data <- readr::read_tsv(file_path, show_col_types = FALSE)
  
  # Add peptide length column if not present
  if (!"Peptide Length" %in% colnames(peptide_data)) {
    peptide_data <- peptide_data %>%
      dplyr::mutate(`Peptide Length` = nchar(Peptide))
  }
  
  # Filter by peptide length
  peptide_data <- peptide_data %>%
    dplyr::filter(`Peptide Length` >= min_length & `Peptide Length` <= max_length)
  
  cat("Loaded", nrow(peptide_data), "peptide entries with lengths", 
      min_length, "to", max_length, "\n")
  
  return(peptide_data)
}

#' Create analysis output directories
#' 
#' @param config Configuration list with output settings
#' @return A list of created directory paths
setup_output_dirs <- function(config) {
  # Base directory
  base_dir <- config$output$base_dir
  
  # Create main directory if it doesn't exist
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create subdirectories
  data_dir <- file.path(base_dir, "data")
  viz_dir <- file.path(base_dir, "visualizations")
  report_dir <- file.path(base_dir, "reports")
  
  for (dir in c(data_dir, viz_dir, report_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # Return paths
  return(list(
    base_dir = base_dir,
    data_dir = data_dir,
    viz_dir = viz_dir,
    report_dir = report_dir
  ))
}

#' Filter peptide data by sample ID
#' 
#' @param peptide_data Data frame of peptide data
#' @param sample_ids Vector of sample IDs to include
#' @param sample_id_column Name of the column containing sample IDs
#' @return Filtered data frame
filter_by_samples <- function(peptide_data, sample_ids, sample_id_column = "SampleID") {
  # Ensure the sample column exists
  if (!sample_id_column %in% colnames(peptide_data)) {
    stop("Sample ID column '", sample_id_column, "' not found in data")
  }
  
  # Filter to only included samples
  filtered_data <- peptide_data %>%
    dplyr::filter(!!dplyr::sym(sample_id_column) %in% sample_ids)
  
  cat("Filtered data to", length(unique(filtered_data[[sample_id_column]])), 
      "samples:", paste(sort(unique(filtered_data[[sample_id_column]])), collapse = ", "), "\n")
  
  return(filtered_data)
}