# Module for concatenating 2CV and 3CV peptide data files

#' List and categorize TSV files for 2CV and 3CV
#' 
#' @param directory Directory containing the TSV files
#' @param pattern_2cv Pattern to match 2CV files
#' @param pattern_3cv Pattern to match 3CV files
#' @return List of file paths, categorized by CV type
list_tsv_files <- function(directory, pattern_2cv, pattern_3cv) {
  # List all files in directory
  all_files <- list.files(path = directory, full.names = TRUE)
  
  if (length(all_files) == 0) {
    stop("No files found in ", directory)
  }
  
  # Categorize by CV type
  files_2cv <- all_files[grepl(pattern_2cv, all_files)]
  files_3cv <- all_files[grepl(pattern_3cv, all_files)]
  
  cat("Found", length(files_2cv), "2CV files and", length(files_3cv), "3CV files\n")
  
  return(list(files_2cv = files_2cv, files_3cv = files_3cv))
}

#' Extract sample ID from filename
#' 
#' @param filename Peptide data filename
#' @return Sample ID string
extract_sample_id <- function(filename) {
  # Extract sample ID using regex
  filename <- basename(filename)
  
  # Try different patterns to extract sample ID
  # First try pattern like: DDA_2CV_117_peptides.tsv
  sample_id <- gsub(".*DDA_[23]CV_([^_]+).*\\.tsv$", "\\1", filename)
  
  # For filenames with extra parts like "_01", clean up the sample ID
  if (grepl("_\\d+_peptides\\.tsv$", filename)) {
    sample_id <- gsub("(.*)_\\d+$", "\\1", sample_id)
  }
  
  # If the above didn't work, try another pattern
  if (sample_id == filename) {
    sample_id <- gsub(".*_([^_]+)_peptides\\.tsv$", "\\1", filename)
  }
  
  # If still not extracted, use a fallback
  if (sample_id == filename) {
    cat("Warning: Could not extract sample ID from", filename, "\n")
    sample_id <- "unknown"
  }
  
  return(sample_id)
}

#' Read peptide files of a specific CV type
#' 
#' @param file_list List of files to read
#' @param cv_type CV type string (2CV or 3CV)
#' @return Combined data frame of all files
read_peptide_files <- function(file_list, cv_type) {
  if (length(file_list) == 0) {
    return(NULL)
  }
  
  all_data <- list()
  
  for (i in seq_along(file_list)) {
    file <- file_list[i]
    cat("Reading", cv_type, "file:", basename(file), "\n")
    
    # Read the TSV file
    data <- readr::read_tsv(file, show_col_types = FALSE)
    
    # Add metadata columns
    data$SourceFile <- basename(file)
    data$SampleID <- extract_sample_id(file)
    data$CVType <- cv_type
    
    all_data[[i]] <- data
  }
  
  # Combine all data frames
  combined_data <- dplyr::bind_rows(all_data)
  
  return(combined_data)
}

#' Run the 2CV/3CV concatenation
#' 
#' @param config Configuration list
#' @return Path to the combined peptide file
run_concatenation <- function(config) {
  cat("\n--- Running 2CV/3CV Concatenation ---\n")
  
  # Extract configuration parameters
  raw_data_dir <- config$input$raw_data_dir
  pattern_2cv <- config$preprocessing$pattern_2cv
  pattern_3cv <- config$preprocessing$pattern_3cv
  output_file <- config$input$peptide_file
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Step 1: List and categorize files
  file_lists <- list_tsv_files(raw_data_dir, pattern_2cv, pattern_3cv)
  
  # Step 2: Read and process files
  cat("\nReading 2CV files...\n")
  data_2cv <- read_peptide_files(file_lists$files_2cv, "2CV")
  
  cat("\nReading 3CV files...\n")
  data_3cv <- read_peptide_files(file_lists$files_3cv, "3CV")
  
  # Make sure both are not NULL
  if (is.null(data_2cv) && is.null(data_3cv)) {
    stop("No 2CV or 3CV data found in the specified directory")
  }
  
  # Step 3: Combine data
  combined_data <- dplyr::bind_rows(data_2cv, data_3cv)
  
  cat("\nCombined dataset has", nrow(combined_data), "entries\n")
  
  # Step 4: Calculate detection statistics
  all_peptides <- unique(combined_data$Peptide)
  
  peptides_2cv <- if (!is.null(data_2cv)) unique(data_2cv$Peptide) else character(0)
  peptides_3cv <- if (!is.null(data_3cv)) unique(data_3cv$Peptide) else character(0)
  
  # Detection counts
  detected_2cv_only <- sum(all_peptides %in% peptides_2cv & !(all_peptides %in% peptides_3cv))
  detected_3cv_only <- sum(!(all_peptides %in% peptides_2cv) & all_peptides %in% peptides_3cv)
  detected_both <- sum(all_peptides %in% peptides_2cv & all_peptides %in% peptides_3cv)
  
  cat("\nDetection statistics:\n")
  cat("Total unique peptides:", length(all_peptides), "\n")
  cat("Peptides detected in 2CV only:", detected_2cv_only, "\n")
  cat("Peptides detected in 3CV only:", detected_3cv_only, "\n")
  cat("Peptides detected in both datasets:", detected_both, "\n")
  
  # Step 5: Add detection flags
  cat("\nAdding detection flags...\n")
  combined_data <- combined_data %>%
    dplyr::mutate(
      detected_2cv = Peptide %in% peptides_2cv,
      detected_3cv = Peptide %in% peptides_3cv,
      detected_both = Peptide %in% peptides_2cv & Peptide %in% peptides_3cv
    )
  
  # Step 6: Save the combined data
  cat("\nSaving combined data to", output_file, "...\n")
  readr::write_tsv(combined_data, output_file)
  
  cat("Concatenation complete!\n")
  
  return(output_file)
}