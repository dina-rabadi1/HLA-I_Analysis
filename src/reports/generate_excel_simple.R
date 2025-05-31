# Simplified Excel report generation module

#' Generate simplified Excel reports for analysis results
#' 
#' @param results List containing all analysis results
#' @param config Configuration list
#' @return Path to the generated Excel file
generate_excel_reports_simple <- function(results, config) {
  # Get output directory
  dirs <- results$dirs
  
  # Create list of data frames for each sheet
  excel_sheets <- list()
  
  # Add analysis overview sheet
  excel_sheets[["Analysis_Overview"]] <- data.frame(
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
  
  # Add shared peptide results if available
  if (!is.null(results$shared_peptide)) {
    # Main peptide report
    if (!is.null(results$shared_peptide$peptide_report) && 
        nrow(results$shared_peptide$peptide_report) > 0) {
      excel_sheets[["Peptide_Sharing_Report"]] <- results$shared_peptide$peptide_report
    }
    
    # Sharing summary
    if (!is.null(results$shared_peptide$sharing_summary) && 
        nrow(results$shared_peptide$sharing_summary) > 0) {
      excel_sheets[["Sharing_Summary"]] <- results$shared_peptide$sharing_summary
    }
    
    # Unique vs Shared
    if (!is.null(results$shared_peptide$unique_shared) && 
        nrow(results$shared_peptide$unique_shared) > 0) {
      excel_sheets[["Unique_vs_Shared"]] <- results$shared_peptide$unique_shared
    }
    
    # Shared peptides
    if (!is.null(results$shared_peptide$shared_peptides) && 
        nrow(results$shared_peptide$shared_peptides) > 0) {
      excel_sheets[["Shared_Peptides"]] <- results$shared_peptide$shared_peptides
    }
  }
  
  # Add length distribution results if available
  if (!is.null(results$length_distribution)) {
    # Length stats
    if (!is.null(results$length_distribution$length_stats) && 
        nrow(results$length_distribution$length_stats) > 0) {
      excel_sheets[["Length_Distribution"]] <- results$length_distribution$length_stats
    }
    
    # Length by sample
    if (!is.null(results$length_distribution$length_by_sample) && 
        nrow(results$length_distribution$length_by_sample) > 0) {
      excel_sheets[["Length_by_Sample"]] <- results$length_distribution$length_by_sample
    }
    
    # Summary statistics
    excel_sheets[["Length_Statistics"]] <- data.frame(
      Statistic = c("Total Peptides", "Mean Length", "Median Length", "Mode Length"),
      Value = c(
        ifelse(is.null(results$length_distribution$total_peptides), 0, 
               results$length_distribution$total_peptides),
        ifelse(is.null(results$length_distribution$mean_length), 0, 
               results$length_distribution$mean_length),
        ifelse(is.null(results$length_distribution$median_length), 0, 
               results$length_distribution$median_length),
        ifelse(is.null(results$length_distribution$mode_length), 0, 
               results$length_distribution$mode_length)
      )
    )
  }
  
  # Add spike-in analysis results if available
  if (!is.null(results$spike_in)) {
    # Comparison metrics
    if (!is.null(results$spike_in$comparison_metrics) && 
        nrow(results$spike_in$comparison_metrics) > 0) {
      excel_sheets[["Spike_Comparison"]] <- results$spike_in$comparison_metrics
    }
    
    # Raw data
    if (!is.null(results$spike_in$spike_matrix_data) && 
        nrow(results$spike_in$spike_matrix_data) > 0) {
      excel_sheets[["Spike_Raw_Data"]] <- results$spike_in$spike_matrix_data
    }
    
    # Add the source report sheet (if available)
    if (!is.null(results$spike_in$source_report) && 
        nrow(results$spike_in$source_report) > 0) {
      excel_sheets[["Peptide_Sources"]] <- results$spike_in$source_report
    }
    
    # Recovery statistics
    if (!is.null(results$spike_in$recovery_stats) && 
        !is.null(results$spike_in$original_sample) && 
        !is.null(results$spike_in$spiked_sample)) {
      
      stats <- results$spike_in$recovery_stats
      required_fields <- c("total_peptides", "detected_original", "detected_spiked", 
                           "detected_both", "only_in_original", "only_in_spiked", 
                           "not_detected")
      
      # Check if required fields exist
      has_all_fields <- all(sapply(required_fields, function(field) !is.null(stats[[field]])))
      
      if (has_all_fields) {
        recovery_df <- data.frame(
          Statistic = c(
            "Total Spike-in Peptides",
            paste0("Detected in ", results$spike_in$original_sample),
            paste0("Detected in ", results$spike_in$spiked_sample),
            "Detected in Both",
            paste0("Only in ", results$spike_in$original_sample),
            paste0("Only in ", results$spike_in$spiked_sample),
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
        
        excel_sheets[["Spike_Recovery_Stats"]] <- recovery_df
      }
    }
  }
  
  # Add fusion analysis results if available
  if (!is.null(results$fusion)) {
    # Fusion peptide results
    if (!is.null(results$fusion$fusion_results) && 
        nrow(results$fusion$fusion_results) > 0) {
      # Just use the columns we know exist
      excel_sheets[["Fusion_Peptides"]] <- results$fusion$fusion_results %>%
        select(Peptide, SampleID, Intensity, Probability, 
               seq1_part, seq2_part, seq1_contribution, seq2_contribution, 
               length) %>%
        arrange(SampleID, desc(Intensity))
    }
    
    # Sample counts
    if (!is.null(results$fusion$fusion_metrics) && 
        !is.null(results$fusion$fusion_metrics$sample_counts) && 
        nrow(results$fusion$fusion_metrics$sample_counts) > 0) {
      excel_sheets[["Fusion_Sample_Counts"]] <- results$fusion$fusion_metrics$sample_counts
    }
    
    # Peptide length counts
    if (!is.null(results$fusion$fusion_metrics) && 
        !is.null(results$fusion$fusion_metrics$peptide_length_counts) && 
        nrow(results$fusion$fusion_metrics$peptide_length_counts) > 0) {
      excel_sheets[["Fusion_Length_Counts"]] <- results$fusion$fusion_metrics$peptide_length_counts
    }
    
    # Add summary statistics for fusion analysis
    fusion_stats <- data.frame(
      Statistic = c(
        "Total Fusion Peptides Detected",
        "Unique Fusion Peptides",
        "Samples with Fusion Peptides",
        "Total Samples Analyzed",
        "Detection Rate (%)",
        "Total Possible Fusion Peptides",
        "Coverage of Possible Peptides (%)"
      ),
      Value = c(
        results$fusion$fusion_metrics$total_fusion_peptides,
        results$fusion$fusion_metrics$unique_fusion_peptides,
        results$fusion$fusion_metrics$samples_with_fusion,
        results$fusion$fusion_metrics$samples_analyzed,
        ifelse(results$fusion$fusion_metrics$samples_analyzed > 0,
               results$fusion$fusion_metrics$samples_with_fusion / 
                 results$fusion$fusion_metrics$samples_analyzed * 100, 0),
        results$fusion$fusion_metrics$total_possible_peptides,
        ifelse(results$fusion$fusion_metrics$total_possible_peptides > 0,
               results$fusion$fusion_metrics$unique_fusion_peptides / 
                 results$fusion$fusion_metrics$total_possible_peptides * 100, 0)
      )
    )
    excel_sheets[["Fusion_Summary"]] <- fusion_stats
    
    # Add fusion peptide list - these are the theoretical peptides we searched for
    if (!is.null(results$fusion$fusion_peptides) && 
        nrow(results$fusion$fusion_peptides) > 0) {
      excel_sheets[["Theoretical_Fusion_Peptides"]] <- results$fusion$fusion_peptides
    }
    
    # Add fusion peptide matrix
    if (!is.null(results$fusion$fusion_matrix) && 
        !is.null(results$fusion$fusion_matrix$presence_matrix) && 
        nrow(results$fusion$fusion_matrix$presence_matrix) > 0) {
      # Convert presence matrix to dataframe
      presence_df <- as.data.frame(results$fusion$fusion_matrix$presence_matrix)
      presence_df$Peptide <- rownames(results$fusion$fusion_matrix$presence_matrix)
      presence_df <- presence_df %>%
        select(Peptide, everything())
      
      excel_sheets[["Fusion_Peptide_Matrix"]] <- presence_df
    }
  }
  
  # Create the Excel file path
  excel_path <- file.path(dirs$report_dir, "analysis_results.xlsx")
  
  # Use a simple approach to save the Excel file
  tryCatch({
    # Save each sheet to a separate CSV file first
    csv_files <- list()
    for (sheet_name in names(excel_sheets)) {
      csv_path <- file.path(dirs$report_dir, paste0(sheet_name, ".csv"))
      write.csv(excel_sheets[[sheet_name]], csv_path, row.names = FALSE)
      csv_files[[sheet_name]] <- csv_path
    }
    
    # Try to save as Excel if openxlsx is available
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      wb <- openxlsx::createWorkbook()
      
      for (sheet_name in names(excel_sheets)) {
        sheet_name_safe <- substr(gsub("[^A-Za-z0-9_]", "_", sheet_name), 1, 31)
        
        tryCatch({
          openxlsx::addWorksheet(wb, sheet_name_safe)
          openxlsx::writeData(wb, sheet_name_safe, excel_sheets[[sheet_name]])
        }, error = function(e) {
          cat("Warning: Could not add worksheet", sheet_name_safe, ":", e$message, "\n")
        })
      }
      
      tryCatch({
        openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
        cat("Excel report generated:", excel_path, "\n")
      }, error = function(e) {
        cat("Warning: Could not save Excel workbook:", e$message, "\n")
        cat("CSV files were saved in the report directory as a fallback.\n")
      })
    } else {
      cat("Warning: openxlsx package is not available. CSV files were saved instead.\n")
    }
    
    # Return the path even if only CSVs were saved
    return(excel_path)
  }, error = function(e) {
    cat("Warning: Error in Excel generation:", e$message, "\n")
    cat("Falling back to individual CSV files in the report directory.\n")
    return(dirs$report_dir) # Return the directory path instead
  })
}