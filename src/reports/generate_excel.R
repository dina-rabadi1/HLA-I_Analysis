# Module for generating Excel reports

#' Create an Excel workbook with styles
#' 
#' @return Excel workbook object
create_styled_workbook <- function() {
  wb <- openxlsx::createWorkbook()
  
  # Create styles
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9D9D9",
    border = "bottom",
    fontSize = 12
  )
  
  number_style <- openxlsx::createStyle(
    numFmt = "#,##0"
  )
  
  percent_style <- openxlsx::createStyle(
    numFmt = "0.0%"
  )
  
  intensity_style <- openxlsx::createStyle(
    numFmt = "0.00E+00"
  )
  
  # Add styles to workbook
  wb$styles <- list(
    header = header_style,
    number = number_style,
    percent = percent_style,
    intensity = intensity_style
  )
  
  return(wb)
}

#' Add shared peptide results to Excel workbook
#' 
#' @param wb Excel workbook object
#' @param shared_results Shared peptide analysis results
#' @return Updated Excel workbook
add_shared_peptide_sheets <- function(wb, shared_results) {
  # Add main peptide report
  if (!is.null(shared_results$peptide_report) && nrow(shared_results$peptide_report) > 0) {
    openxlsx::addWorksheet(wb, "Peptide Sharing Report")
    openxlsx::writeData(wb, "Peptide Sharing Report", shared_results$peptide_report)
    openxlsx::addStyle(wb, "Peptide Sharing Report", wb$styles$header, 
                       rows = 1, cols = 1:ncol(shared_results$peptide_report))
    
    # Adjust column widths
    openxlsx::setColWidths(wb, "Peptide Sharing Report", 
                           cols = 1:ncol(shared_results$peptide_report), 
                           widths = "auto")
  }
  
  # Add sharing summary
  if (!is.null(shared_results$sharing_summary) && nrow(shared_results$sharing_summary) > 0) {
    openxlsx::addWorksheet(wb, "Sharing Summary")
    openxlsx::writeData(wb, "Sharing Summary", shared_results$sharing_summary)
    openxlsx::addStyle(wb, "Sharing Summary", wb$styles$header, 
                       rows = 1, cols = 1:ncol(shared_results$sharing_summary))
  }
  
  # Add unique vs shared summary
  if (!is.null(shared_results$unique_shared) && nrow(shared_results$unique_shared) > 0) {
    openxlsx::addWorksheet(wb, "Unique vs Shared")
    openxlsx::writeData(wb, "Unique vs Shared", shared_results$unique_shared)
    openxlsx::addStyle(wb, "Unique vs Shared", wb$styles$header, 
                       rows = 1, cols = 1:ncol(shared_results$unique_shared))
  }
  
  # Add shared peptides
  if (!is.null(shared_results$shared_peptides) && nrow(shared_results$shared_peptides) > 0) {
    openxlsx::addWorksheet(wb, "Shared Peptides")
    openxlsx::writeData(wb, "Shared Peptides", shared_results$shared_peptides)
    openxlsx::addStyle(wb, "Shared Peptides", wb$styles$header, 
                       rows = 1, cols = 1:ncol(shared_results$shared_peptides))
  }
  
  # Add sheets for different sharing levels
  if (!is.null(shared_results$sharing_level_groups)) {
    for (i in 2:min(10, length(shared_results$sharing_level_groups) + 1)) {
      group_name <- paste0("in_", i, "_samples")
      if (!is.null(shared_results$sharing_level_groups[[group_name]]) &&
          nrow(shared_results$sharing_level_groups[[group_name]]) > 0) {
        
        sheet_name <- paste0("In_", i, "_Samples")
        openxlsx::addWorksheet(wb, sheet_name)
        openxlsx::writeData(wb, sheet_name, 
                            shared_results$sharing_level_groups[[group_name]])
        openxlsx::addStyle(wb, sheet_name, wb$styles$header, 
                           rows = 1, 
                           cols = 1:ncol(shared_results$sharing_level_groups[[group_name]]))
        
        # Adjust column widths
        openxlsx::setColWidths(wb, sheet_name, 
                               cols = 1:ncol(shared_results$sharing_level_groups[[group_name]]), 
                               widths = "auto")
      }
    }
  }
  
  return(wb)
}

#' Add length distribution results to Excel workbook
#' 
#' @param wb Excel workbook object
#' @param length_results Length distribution analysis results
#' @return Updated Excel workbook
add_length_distribution_sheets <- function(wb, length_results) {
  # Add main length distribution summary
  if (!is.null(length_results$length_stats) && nrow(length_results$length_stats) > 0) {
    openxlsx::addWorksheet(wb, "Length Distribution")
    openxlsx::writeData(wb, "Length Distribution", length_results$length_stats)
    openxlsx::addStyle(wb, "Length Distribution", wb$styles$header, 
                       rows = 1, cols = 1:ncol(length_results$length_stats))
  }
  
  # Add sample-level length distribution
  if (!is.null(length_results$length_by_sample) && nrow(length_results$length_by_sample) > 0) {
    openxlsx::addWorksheet(wb, "Length by Sample")
    openxlsx::writeData(wb, "Length by Sample", length_results$length_by_sample)
    openxlsx::addStyle(wb, "Length by Sample", wb$styles$header, 
                       rows = 1, cols = 1:ncol(length_results$length_by_sample))
  }
  
  # Add summary statistics
  summary_stats <- data.frame(
    Statistic = c("Total Peptides", "Mean Length", "Median Length", "Mode Length"),
    Value = c(
      ifelse(is.null(length_results$total_peptides), 0, length_results$total_peptides),
      ifelse(is.null(length_results$mean_length), 0, length_results$mean_length),
      ifelse(is.null(length_results$median_length), 0, length_results$median_length),
      ifelse(is.null(length_results$mode_length), 0, length_results$mode_length)
    )
  )
  
  openxlsx::addWorksheet(wb, "Length Statistics")
  openxlsx::writeData(wb, "Length Statistics", summary_stats)
  openxlsx::addStyle(wb, "Length Statistics", wb$styles$header, 
                     rows = 1, cols = 1:ncol(summary_stats))
  
  return(wb)
}

#' Add spike-in analysis results to Excel workbook
#' 
#' @param wb Excel workbook object
#' @param spike_results Spike-in analysis results
#' @return Updated Excel workbook
add_spike_in_sheets <- function(wb, spike_results) {
  # Add spike-in comparison sheet
  if (!is.null(spike_results$comparison_metrics) && nrow(spike_results$comparison_metrics) > 0) {
    openxlsx::addWorksheet(wb, "Spike Comparison")
    openxlsx::writeData(wb, "Spike Comparison", spike_results$comparison_metrics)
    openxlsx::addStyle(wb, "Spike Comparison", wb$styles$header, 
                       rows = 1, cols = 1:ncol(spike_results$comparison_metrics))
  }
  
  # Add raw data sheet
  if (!is.null(spike_results$spike_matrix_data) && nrow(spike_results$spike_matrix_data) > 0) {
    openxlsx::addWorksheet(wb, "Spike Raw Data")
    openxlsx::writeData(wb, "Spike Raw Data", spike_results$spike_matrix_data)
    openxlsx::addStyle(wb, "Spike Raw Data", wb$styles$header, 
                       rows = 1, cols = 1:ncol(spike_results$spike_matrix_data))
  }
  
  # Add recovery statistics if available
  if (!is.null(spike_results$recovery_stats)) {
    # Check if all expected fields exist
    required_fields <- c("total_peptides", "detected_original", "detected_spiked", 
                         "detected_both", "only_in_original", "only_in_spiked", 
                         "not_detected")
    
    # Ensure all fields exist or provide defaults
    stats <- spike_results$recovery_stats
    for (field in required_fields) {
      if (is.null(stats[[field]])) {
        stats[[field]] <- 0
      }
    }
    
    # Only proceed if we have the original and spiked sample names
    if (!is.null(spike_results$original_sample) && !is.null(spike_results$spiked_sample)) {
      recovery_df <- data.frame(
        Statistic = c(
          "Total Spike-in Peptides",
          paste0("Detected in ", spike_results$original_sample),
          paste0("Detected in ", spike_results$spiked_sample),
          "Detected in Both",
          paste0("Only in ", spike_results$original_sample),
          paste0("Only in ", spike_results$spiked_sample),
          "Not Detected in Either"
        ),
        Count = c(
          stats$total_peptides,
          stats$detected_original,
          stats$detected_spiked,
          stats$detected_both,
          stats$only_in_original,
          stats$only_in_spiked,
          stats$not_detected
        )
      )
      
      # Add percentage column if total_peptides > 0
      if (stats$total_peptides > 0) {
        recovery_df$Percentage <- c(
          100,
          stats$detected_original / stats$total_peptides * 100,
          stats$detected_spiked / stats$total_peptides * 100,
          stats$detected_both / stats$total_peptides * 100,
          stats$only_in_original / stats$total_peptides * 100,
          stats$only_in_spiked / stats$total_peptides * 100,
          stats$not_detected / stats$total_peptides * 100
        )
      } else {
        recovery_df$Percentage <- rep(0, nrow(recovery_df))
      }
      
      openxlsx::addWorksheet(wb, "Spike Recovery Stats")
      openxlsx::writeData(wb, "Spike Recovery Stats", recovery_df)
      openxlsx::addStyle(wb, "Spike Recovery Stats", wb$styles$header, 
                         rows = 1, cols = 1:ncol(recovery_df))
    }
  }
  
  return(wb)
}

#' Generate Excel reports for analysis results
#' 
#' @param results List containing all analysis results
#' @param config Configuration list
#' @return Path to the generated Excel file
generate_excel_reports <- function(results, config) {
  # Get output directory
  dirs <- results$dirs
  
  # Create workbook
  wb <- create_styled_workbook()
  
  # Add analysis overview sheet
  overview <- data.frame(
    Parameter = c(
      "Analysis Date",
      "Analysis Type",
      "Samples Included",
      "Total Peptides Analyzed",
      "Peptide Length Range",
      "Output Directory"
    ),
    Value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      paste(names(results)[!(names(results) %in% c("dirs", "files"))], collapse = ", "),
      paste(config$samples$include, collapse = ", "),
      ifelse(is.null(results$shared_peptide), "Not analyzed", 
             as.character(results$shared_peptide$n_total_peptides)),
      paste(config$analysis$min_length, "-", config$analysis$max_length),
      dirs$base_dir
    )
  )
  
  openxlsx::addWorksheet(wb, "Analysis Overview")
  openxlsx::writeData(wb, "Analysis Overview", overview)
  openxlsx::addStyle(wb, "Analysis Overview", wb$styles$header, 
                     rows = 1, cols = 1:ncol(overview))
  
  # Add specific analysis sheets with error handling
  tryCatch({
    if (!is.null(results$shared_peptide)) {
      cat("Adding shared peptide results to Excel report...\n")
      wb <- add_shared_peptide_sheets(wb, results$shared_peptide)
    }
  }, error = function(e) {
    cat("Warning: Error adding shared peptide results to Excel report:", e$message, "\n")
  })
  
  tryCatch({
    if (!is.null(results$length_distribution)) {
      cat("Adding length distribution results to Excel report...\n")
      wb <- add_length_distribution_sheets(wb, results$length_distribution)
    }
  }, error = function(e) {
    cat("Warning: Error adding length distribution results to Excel report:", e$message, "\n")
  })
  
  tryCatch({
    if (!is.null(results$spike_in)) {
      cat("Adding spike-in analysis results to Excel report...\n")
      wb <- add_spike_in_sheets(wb, results$spike_in)
    }
  }, error = function(e) {
    cat("Warning: Error adding spike-in analysis results to Excel report:", e$message, "\n")
  })
  
  # Save workbook
  excel_path <- file.path(dirs$report_dir, "analysis_results.xlsx")
  openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
  
  cat("Excel report generated:", excel_path, "\n")
  
  return(excel_path)
}