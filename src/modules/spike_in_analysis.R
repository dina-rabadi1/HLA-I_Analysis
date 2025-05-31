# Module for 51 vs 51S spike-in peptide analysis


# Add this to spike_in_analysis.R

#' Debug function to print peptide search results
#'
#' @param peptide_data All peptide data
#' @param spike_peptides List of spike-in peptides
#' @param original_sample Original sample ID
#' @param spiked_sample Spiked sample ID
debug_spike_peptide_search <- function(peptide_data, spike_peptides, original_sample, spiked_sample) {
  cat("\nDebugging spike peptide search:\n")
  
  for (peptide in spike_peptides) {
    # Check if this peptide is in either sample
    orig_match <- peptide_data %>% 
      filter(Peptide == peptide, SampleID == original_sample)
    spike_match <- peptide_data %>% 
      filter(Peptide == peptide, SampleID == spiked_sample)
    
    cat(peptide, ":\n")
    cat("  - Found in original sample (", original_sample, "): ", nrow(orig_match) > 0, 
        if(nrow(orig_match) > 0) paste0(" (Intensity: ", sum(orig_match$Intensity), ")") else "", "\n", sep="")
    cat("  - Found in spiked sample (", spiked_sample, "): ", nrow(spike_match) > 0,
        if(nrow(spike_match) > 0) paste0(" (Intensity: ", sum(spike_match$Intensity), ")") else "", "\n", sep="")
  }
  
  # Also check for any file pattern that might be related to spike-in
  spike_files <- peptide_data %>%
    filter(grepl("spike|51S|untargeted", SourceFile, ignore.case = TRUE)) %>%
    select(SourceFile) %>%
    distinct()
  
  if(nrow(spike_files) > 0) {
    cat("\nFound potential spike-in related files:\n")
    for(i in 1:nrow(spike_files)) {
      cat("  - ", spike_files$SourceFile[i], "\n", sep="")
    }
  }
}

#' Analyze spike-in peptides (51 vs 51S)
#' 
#' @param peptide_data Data frame of peptide data
#' @param config Configuration settings
#' @return List containing analysis results
run_spike_in_analysis <- function(peptide_data, config) {
  cat("Starting spike-in peptide analysis (51 vs 51S)...\n")
  
  # Extract configuration parameters
  original_sample <- config$analysis$sample_comparisons$original_sample
  spiked_sample <- config$analysis$sample_comparisons$spiked_sample
  spiked_peptides <- config$analysis$sample_comparisons$spiked_peptides
  
  cat("Comparing", original_sample, "(original) with", spiked_sample, "(spiked)\n")
  cat("Analyzing", length(spiked_peptides), "spiked peptides\n")
  
  # Filter data to the two samples of interest with improved handling
  sample_data <- peptide_data %>%
    # Check for untargeted files which might contain spike data
    mutate(is_spike_file = grepl("51S|untargeted", SourceFile, ignore.case = TRUE)) %>%
    dplyr::filter(
      # Either the sample ID matches what we want
      SampleID %in% c(original_sample, spiked_sample) |
        # Or it's from a special spike file (very important for 51S)
        (is_spike_file & grepl(spiked_sample, SourceFile, ignore.case = TRUE))
    )
  
  # Check if any data remains after filtering
  if(nrow(sample_data) == 0) {
    warning("No data found after filtering for samples ", 
            original_sample, " and ", spiked_sample, 
            ". Check sample IDs and file names.")
    
    # Print some sample IDs to help debugging
    cat("Available sample IDs in data:", 
        paste(unique(peptide_data$SampleID)[1:min(10, length(unique(peptide_data$SampleID)))], 
              collapse = ", "), 
        if(length(unique(peptide_data$SampleID)) > 10) "..." else "", "\n")
  }
  
  # Add extra debugging for filtered data
  cat("\nAfter filtering, found", nrow(sample_data), "rows of data\n")
  cat("Sample IDs in filtered data:", paste(unique(sample_data$SampleID), collapse = ", "), "\n")
  
  # Extract the spike-in peptides
  spike_peptide_data <- sample_data %>%
    dplyr::filter(Peptide %in% spiked_peptides)
  
  # Add debugging call here
  debug_spike_peptide_search(peptide_data, spiked_peptides, original_sample, spiked_sample)
  
  # Create summary of detection - with column name check
  # First check which column names are available
  spectral_count_col <- if("Spectral.Count" %in% colnames(spike_peptide_data)) {
    "Spectral.Count"
  } else if("spectral_count" %in% colnames(spike_peptide_data)) {
    "spectral_count"
  } else {
    # If neither column exists, create a dummy column with zeros
    spike_peptide_data$spectral_count <- 0
    "spectral_count"
  }
  
  intensity_col <- if("Intensity" %in% colnames(spike_peptide_data)) {
    "Intensity"
  } else if("intensity" %in% colnames(spike_peptide_data)) {
    "intensity"
  } else {
    # If neither column exists, create a dummy column with zeros
    spike_peptide_data$intensity <- 0
    "intensity"
  }
  
  # Now use the detected column names in the summarize function
  spike_summary <- spike_peptide_data %>%
    dplyr::group_by(SampleID, Peptide) %>%
    dplyr::summarize(
      spectral_count = sum(!!dplyr::sym(spectral_count_col)),
      total_intensity = sum(!!dplyr::sym(intensity_col)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Peptide, SampleID)
  
  # Create a complete matrix with all peptide-sample combinations
  all_combinations <- expand.grid(
    Peptide = spiked_peptides,
    SampleID = c(original_sample, spiked_sample),
    stringsAsFactors = FALSE
  )
  
  # Merge with actual data, filling in zeros for missing combinations
  spike_matrix_data <- all_combinations %>%
    dplyr::left_join(spike_summary, by = c("Peptide", "SampleID")) %>%
    dplyr::mutate(
      spectral_count = ifelse(is.na(spectral_count), 0, spectral_count),
      total_intensity = ifelse(is.na(total_intensity), 0, total_intensity),
      detected = total_intensity > 0
    )
  
  # Create a side-by-side comparison
  spike_comparison <- spike_matrix_data %>%
    tidyr::pivot_wider(
      names_from = SampleID,
      values_from = c(detected, spectral_count, total_intensity)
    )
  
  # Add comparison metrics
  comparison_metrics <- spike_comparison %>%
    dplyr::mutate(
      detection_status = dplyr::case_when(
        get(paste0("detected_", original_sample)) & get(paste0("detected_", spiked_sample)) ~ 
          "Detected in both",
        !get(paste0("detected_", original_sample)) & get(paste0("detected_", spiked_sample)) ~ 
          "Only in spiked sample",
        get(paste0("detected_", original_sample)) & !get(paste0("detected_", spiked_sample)) ~ 
          "Only in original sample",
        TRUE ~ "Not detected in either"
      ),
      intensity_fold_change = dplyr::case_when(
        get(paste0("total_intensity_", original_sample)) > 0 ~ 
          get(paste0("total_intensity_", spiked_sample)) / get(paste0("total_intensity_", original_sample)),
        get(paste0("total_intensity_", spiked_sample)) > 0 ~ Inf,
        TRUE ~ NA_real_
      ),
      intensity_increase = get(paste0("total_intensity_", spiked_sample)) - 
        get(paste0("total_intensity_", original_sample))
    )
  
  # Calculate recovery statistics
  recovery_stats <- list(
    total_peptides = length(spiked_peptides),
    detected_original = sum(spike_comparison[[paste0("detected_", original_sample)]]),
    detected_spiked = sum(spike_comparison[[paste0("detected_", spiked_sample)]]),
    detected_both = sum(spike_comparison[[paste0("detected_", original_sample)]] & 
                          spike_comparison[[paste0("detected_", spiked_sample)]]),
    only_in_original = sum(spike_comparison[[paste0("detected_", original_sample)]] & 
                             !spike_comparison[[paste0("detected_", spiked_sample)]]),
    only_in_spiked = sum(!spike_comparison[[paste0("detected_", original_sample)]] & 
                           spike_comparison[[paste0("detected_", spiked_sample)]]),
    not_detected = sum(!spike_comparison[[paste0("detected_", original_sample)]] & 
                         !spike_comparison[[paste0("detected_", spiked_sample)]])
  )
  
  # Return results
  results <- list(
    spike_matrix_data = spike_matrix_data,
    spike_comparison = spike_comparison,
    comparison_metrics = comparison_metrics,
    recovery_stats = recovery_stats,
    spiked_peptides = spiked_peptides,
    original_sample = original_sample,
    spiked_sample = spiked_sample
  )
  
  # Print summary
  cat("\nSpike-in peptide recovery statistics:\n")
  cat("Total spike-in peptides:", recovery_stats$total_peptides, "\n")
  cat("Detected in original sample (", original_sample, "): ", 
      recovery_stats$detected_original, " (", 
      sprintf("%.1f%%", recovery_stats$detected_original/recovery_stats$total_peptides*100), ")\n", sep="")
  cat("Detected in spiked sample (", spiked_sample, "): ", 
      recovery_stats$detected_spiked, " (", 
      sprintf("%.1f%%", recovery_stats$detected_spiked/recovery_stats$total_peptides*100), ")\n", sep="")
  cat("Detected in both samples: ", 
      recovery_stats$detected_both, " (", 
      sprintf("%.1f%%", recovery_stats$detected_both/recovery_stats$total_peptides*100), ")\n", sep="")
  cat("Only in spiked sample: ", 
      recovery_stats$only_in_spiked, " (", 
      sprintf("%.1f%%", recovery_stats$only_in_spiked/recovery_stats$total_peptides*100), ")\n", sep="")
  
  cat("Spike-in analysis complete!\n")
  
  # Add to your spike_in_analysis.R file
  cat("All sample IDs in dataset:", paste(unique(peptide_data$SampleID), collapse=", "), "\n")
  
  # List files that might contain 51S data
  potential_51S_files <- peptide_data %>%
    filter(grepl("51S|untargeted", SourceFile, ignore.case=TRUE)) %>%
    select(SourceFile, SampleID) %>%
    distinct()
  
  print(potential_51S_files)
  
  # Create a simple source file report with enhanced sample info
  source_report <- peptide_data %>%
    filter(Peptide %in% spiked_peptides) %>%
    # Add descriptive sample designation
    mutate(
      SampleWithSource = case_when(
        SampleID == original_sample & grepl("2CV", SourceFile) ~ paste0(original_sample, "-2CV"),
        SampleID == original_sample & grepl("3CV", SourceFile) ~ paste0(original_sample, "-3CV"),
        SampleID == spiked_sample ~ paste0(spiked_sample, 
                                           ifelse(grepl("untargeted", SourceFile), "-untargeted", "")),
        TRUE ~ SampleID
      ),
      # Extract acquisition method directly
      AcquisitionMethod = case_when(
        grepl("2CV", SourceFile) ~ "2CV",
        grepl("3CV", SourceFile) ~ "3CV",
        grepl("untargeted", SourceFile) ~ "untargeted",
        TRUE ~ "unknown"
      )
    ) %>%
    select(Peptide, SampleID, SampleWithSource, AcquisitionMethod, SourceFile, Intensity) %>%
    arrange(Peptide, SampleID, AcquisitionMethod)
  
  # Add this after creating the source_report to debug whether considering 51 2cv and 51 3cv
  cat("\nDetailed peptide source information:\n")
  for (peptide in spiked_peptides) {
    sources <- source_report %>%
      filter(Peptide == peptide)
    
    if (nrow(sources) > 0) {
      cat("\n", peptide, "found in:\n")
      for (i in 1:nrow(sources)) {
        cat("  - Sample:", sources$SampleWithSource[i],
            "  Method:", sources$AcquisitionMethod[i],
            "  Intensity:", format(sources$Intensity[i], scientific = TRUE, digits = 3),
            "\n")
      }
    } else {
      cat("\n", peptide, "not found in any file\n")
    }
  }
  
  # Save the source report as a TSV file
  output_dir <- file.path(config$output$base_dir, "data")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  write_tsv(source_report, file.path(output_dir, "peptide_sources.tsv"))
  
  # Include in the results list
  results$source_report <- source_report
  
  return(results)
}
