# Module for peptide length distribution analysis

#' Analyze peptide length distribution
#' 
#' @param peptide_data Data frame of peptide data
#' @param config Configuration list
#' @return List containing analysis results
run_length_distribution_analysis <- function(peptide_data, config) {
  cat("Starting peptide length distribution analysis...\n")
  
  # Filter by included samples first
  sample_ids <- config$samples$include
  cat("Filtering length analysis to samples:", paste(sample_ids, collapse = ", "), "\n")
  filtered_data <- filter_by_samples(peptide_data, sample_ids)
  
  # Extract peptide length data
  length_data <- filtered_data %>%
    dplyr::select(Peptide, `Peptide Length`) %>%
    dplyr::distinct()
  
  # Create summary statistics
  length_dist <- length_data %>%
    dplyr::count(`Peptide Length`) %>%
    dplyr::mutate(
      percentage = n / sum(n) * 100,
      cumulative = cumsum(n),
      cumulative_percent = cumsum(percentage)
    )
  
  # Create a more detailed analysis by sample
  length_by_sample <- filtered_data %>%
    dplyr::group_by(SampleID, `Peptide Length`) %>%
    dplyr::summarize(
      peptide_count = dplyr::n_distinct(Peptide),
      .groups = "drop"
    ) %>%
    dplyr::arrange(SampleID, `Peptide Length`)
  
  # Calculate statistics for each length
  length_stats <- length_data %>%
    dplyr::group_by(`Peptide Length`) %>%
    dplyr::summarize(
      count = n(),
      percentage = count / nrow(length_data) * 100,
      .groups = "drop"
    ) %>%
    dplyr::arrange(`Peptide Length`)
  
  # Calculate mode (most common length)
  mode_length <- length_stats %>%
    dplyr::arrange(dplyr::desc(count)) %>%
    dplyr::slice(1) %>%
    dplyr::pull(`Peptide Length`)
  
  # Calculate median length
  median_length <- median(length_data$`Peptide Length`)
  
  # Calculate mean length
  mean_length <- mean(length_data$`Peptide Length`)
  
  # Return results
  results <- list(
    length_distribution = length_dist,
    length_by_sample = length_by_sample,
    length_stats = length_stats,
    mode_length = mode_length,
    median_length = median_length,
    mean_length = mean_length,
    total_peptides = nrow(length_data)
  )
  
  cat("Peptide length distribution analysis complete!\n")
  cat("Most common peptide length:", mode_length, 
      sprintf("(%.1f%%)", length_stats$percentage[length_stats$`Peptide Length` == mode_length]), "\n")
  cat("Mean peptide length:", sprintf("%.2f", mean_length), "\n")
  cat("Median peptide length:", median_length, "\n")
  
  return(results)
}