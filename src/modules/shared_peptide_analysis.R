# Shared peptide analysis module
# shared_peptide_analysis.R

#' Perform shared peptide analysis across samples
#' 
#' @param peptide_data Data frame of peptide data
#' @param config Configuration list
#' @return List containing analysis results
run_shared_peptide_analysis <- function(peptide_data, config) {
  cat("Starting shared peptide analysis...\n")
  
  # Extract relevant config parameters
  sample_ids <- config$samples$include
  exclude_samples <- config$samples$exclude_from_shared
  
  # Remove excluded samples
  sample_ids <- setdiff(sample_ids, exclude_samples)
  
  cat("Analyzing shared peptides across", length(sample_ids), 
      "samples:", paste(sample_ids, collapse = ", "), "\n")
  
  # First filter data to only selected samples
  filtered_data <- filter_by_samples(peptide_data, sample_ids)
  
  # Create a report showing which samples each peptide appears in
  peptide_report <- filtered_data %>%
    # First get basic peptide and sample information
    dplyr::select(Peptide, SampleID, `Peptide Length`, Gene, Protein) %>%
    dplyr::distinct() %>%
    # Group by peptide to consolidate sample information
    dplyr::group_by(Peptide, `Peptide Length`) %>%
    dplyr::summarize(
      sample_list = paste(sort(SampleID), collapse = ", "),
      sample_count = dplyr::n_distinct(SampleID),
      # Consolidate gene information
      genes = paste(unique(na.omit(Gene)), collapse = "; "),
      # Consolidate protein information
      proteins = paste(unique(na.omit(Protein)), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(sample_count), Peptide)
  
  # Create summary of peptide sharing
  sharing_summary <- peptide_report %>%
    dplyr::count(sample_count) %>%
    dplyr::mutate(
      percentage = n / sum(n) * 100,
      cumulative = cumsum(n),
      cumulative_percent = cumsum(percentage)
    )
  
  # Calculate categories for unique vs shared
  unique_shared <- data.frame(
    category = c("Unique (1 sample)", "Shared (2+ samples)"),
    count = c(
      sum(sharing_summary$n[sharing_summary$sample_count == 1]),
      sum(sharing_summary$n[sharing_summary$sample_count > 1])
    )
  ) %>%
    dplyr::mutate(percentage = count / sum(count) * 100)
  
  # Get peptides by sharing level
  unique_peptides <- peptide_report %>%
    dplyr::filter(sample_count == 1)
  
  shared_peptides <- peptide_report %>%
    dplyr::filter(sample_count > 1) %>%
    dplyr::arrange(dplyr::desc(sample_count), Peptide)
  
  # Group shared peptides by sharing level
  sharing_level_groups <- list()
  for (i in 2:max(peptide_report$sample_count)) {
    sharing_level_groups[[paste0("in_", i, "_samples")]] <- 
      peptide_report %>% 
      dplyr::filter(sample_count == i)
  }
  
  # Return results
  results <- list(
    peptide_report = peptide_report,
    sharing_summary = sharing_summary,
    unique_shared = unique_shared,
    unique_peptides = unique_peptides,
    shared_peptides = shared_peptides,
    sharing_level_groups = sharing_level_groups,
    sample_ids = sample_ids,
    n_samples = length(sample_ids),
    n_total_peptides = nrow(peptide_report),
    n_unique_peptides = nrow(unique_peptides),
    n_shared_peptides = nrow(shared_peptides)
  )
  
  cat("Shared peptide analysis complete!\n")
  cat("Total peptides:", results$n_total_peptides, "\n")
  cat("Unique peptides (1 sample):", results$n_unique_peptides, 
      sprintf("(%.1f%%)", results$n_unique_peptides/results$n_total_peptides*100), "\n")
  cat("Shared peptides (2+ samples):", results$n_shared_peptides, 
      sprintf("(%.1f%%)", results$n_shared_peptides/results$n_total_peptides*100), "\n")
  
  return(results)
}