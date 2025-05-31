# Main script for running the peptide analysis pipeline with transcriptome integration
# run_peptide_transcriptome_pipeline.R

# Load required packages
library(tidyverse)
library(pheatmap)
library(openxlsx)
library(readxl)
library(VennDiagram)

# Source core utilities
source("src/core/peptide_utils.R")
source("src/core/visualization_utils.R")
source("src/core/transcriptome_viz_utils.R")

# Source analysis modules
source("src/modules/concatenate_2cv_3cv.R")
source("src/modules/shared_peptide_analysis.R")
source("src/modules/length_distribution_analysis.R")
source("src/modules/spike_in_analysis.R")
source("src/modules/fusion_analysis.R")
source("src/modules/transcriptome_analysis.R")
source("src/modules/multi_omics_integration.R")

# Source reporting modules
if (file.exists("src/reports/generate_excel_simple.R")) {
  source("src/reports/generate_excel_simple.R")
} else {
  source("src/reports/generate_excel.R")
}

# Get config file from command line arg or use default
args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) > 0) args[1] else "config/transcriptome_config.R"
cat("Using configuration file:", config_file, "\n")
source(config_file)

# Setup output directories
dirs <- setup_output_dirs(config)

# Run concatenation if enabled
if (config$preprocessing$concatenate_2cv_3cv) {
  config$input$peptide_file <- run_concatenation(config)
}

# Load and preprocess data
peptide_data <- load_peptide_data(
  config$input$peptide_file,
  min_length = config$analysis$min_length,
  max_length = config$analysis$max_length
)

# Store results
results <- list(
  dirs = dirs,
  files = list()
)

# Run shared peptide analysis if enabled
if (config$analysis$do_shared_peptide) {
  cat("\n--- Running Shared Peptide Analysis ---\n")
  results$shared_peptide <- run_shared_peptide_analysis(peptide_data, config)
  
  # Create visualizations
  cat("\nCreating shared peptide visualizations...\n")
  viz_paths <- list()
  
  # Only create visualizations if there are actually shared peptides
  if (results$shared_peptide$n_shared_peptides > 0) {
    # Sample count distribution
    viz_paths$sample_count <- create_sample_count_distribution(
      results$shared_peptide$sharing_summary,
      dirs$viz_dir,
      config
    )
    
    # Unique vs shared pie chart
    viz_paths$unique_shared <- create_unique_shared_pie(
      results$shared_peptide$unique_shared,
      dirs$viz_dir,
      config
    )
  } else {
    cat("No shared peptides found - skipping shared peptide visualizations\n")
  }
  
  # Store visualization paths
  results$shared_peptide$visualizations <- viz_paths
  
  # Export result data
  cat("\nExporting shared peptide results...\n")
  write.csv(
    results$shared_peptide$peptide_report,
    file.path(dirs$data_dir, "peptide_sharing_report.csv"),
    row.names = FALSE
  )
  
  write.csv(
    results$shared_peptide$sharing_summary,
    file.path(dirs$data_dir, "sharing_summary.csv"),
    row.names = FALSE
  )
  
  if (results$shared_peptide$n_shared_peptides > 0) {
    write.csv(
      results$shared_peptide$shared_peptides,
      file.path(dirs$data_dir, "shared_peptides.csv"),
      row.names = FALSE
    )
  }
}

# Run peptide length distribution analysis if enabled
if (config$analysis$do_length_distribution) {
  cat("\n--- Running Peptide Length Distribution Analysis ---\n")
  results$length_distribution <- run_length_distribution_analysis(peptide_data, config)
  
  # Create visualizations
  cat("\nCreating length distribution visualizations...\n")
  viz_paths <- list()
  
  # Length distribution bar chart
  viz_paths$length_distribution <- create_length_distribution_plot(
    results$length_distribution$length_stats,
    dirs$viz_dir,
    config
  )
  
  # Length by sample heatmap
  viz_paths$length_by_sample <- create_length_by_sample_heatmap(
    results$length_distribution$length_by_sample,
    dirs$viz_dir,
    config
  )
  
  # Store visualization paths
  results$length_distribution$visualizations <- viz_paths
  
  # Export result data
  cat("\nExporting length distribution results...\n")
  write.csv(
    results$length_distribution$length_stats,
    file.path(dirs$data_dir, "length_distribution.csv"),
    row.names = FALSE
  )
  
  write.csv(
    results$length_distribution$length_by_sample,
    file.path(dirs$data_dir, "length_by_sample.csv"),
    row.names = FALSE
  )
}

# Run 51 vs 51S spike-in analysis if enabled
if (!is.null(config$analysis$sample_comparisons) && 
    !is.null(config$analysis$sample_comparisons$do_51_vs_51S) && 
    config$analysis$sample_comparisons$do_51_vs_51S) {
  cat("\n--- Running 51 vs 51S Spike-in Analysis ---\n")
  results$spike_in <- run_spike_in_analysis(peptide_data, config)
  
  # Create visualizations
  cat("\nCreating spike-in visualizations...\n")
  viz_paths <- list()
  
  # Comparison bar chart
  viz_paths$comparison <- create_spike_comparison_plot(
    results$spike_in$spike_matrix_data,
    dirs$viz_dir,
    config
  )
  
  # Detection heatmap
  viz_paths$heatmap <- create_spike_detection_heatmap(
    results$spike_in$spike_matrix_data,
    dirs$viz_dir,
    config
  )
  
  # Store visualization paths
  results$spike_in$visualizations <- viz_paths
  
  # Export result data
  cat("\nExporting spike-in analysis results...\n")
  write.csv(
    results$spike_in$comparison_metrics,
    file.path(dirs$data_dir, "spike_comparison.csv"),
    row.names = FALSE
  )
}

# Run fusion analysis if enabled
if (!is.null(config$analysis$do_fusion) && config$analysis$do_fusion) {
  cat("\n--- Running DNAJB1:PRKACA Fusion Analysis ---\n")
  results$fusion <- run_fusion_analysis(peptide_data, config)
  
  # Create visualizations
  cat("\nCreating fusion peptide visualizations...\n")
  viz_paths <- list()
  
  # Only create visualizations if there are actually fusion peptides
  if (results$fusion$fusion_metrics$total_fusion_peptides > 0) {
    # Detection heatmap
    viz_paths$heatmap <- create_fusion_detection_heatmap(
      results$fusion$fusion_matrix,
      dirs$viz_dir,
      config
    )
    
    # Distribution plot
    viz_paths$distribution <- create_fusion_distribution_plot(
      results$fusion$fusion_metrics,
      dirs$viz_dir,
      config
    )
    
    # Sequence visualization
    viz_paths$sequence <- create_fusion_sequence_plot(
      results$fusion$fusion_results,
      results$fusion$fusion_peptides,
      dirs$viz_dir,
      config
    )
  } else {
    cat("No fusion peptides found - skipping fusion peptide visualizations\n")
  }
  
  # Store visualization paths
  results$fusion$visualizations <- viz_paths
  
  # Export result data
  cat("\nExporting fusion peptide results...\n")
  
  # Export fusion peptides if found
  if (results$fusion$fusion_metrics$total_fusion_peptides > 0) {
    # Export detected fusion peptides - only use columns we know exist
    write.csv(
      results$fusion$fusion_results %>% 
        select(Peptide, SampleID, Intensity, seq1_part, seq2_part, 
               seq1_contribution, seq2_contribution, Probability),
      file.path(dirs$data_dir, "fusion_peptides.csv"),
      row.names = FALSE
    )
    
    # Export metrics
    write.csv(
      results$fusion$fusion_metrics$sample_counts,
      file.path(dirs$data_dir, "fusion_sample_counts.csv"),
      row.names = FALSE
    )
    
    write.csv(
      results$fusion$fusion_metrics$peptide_length_counts,
      file.path(dirs$data_dir, "fusion_peptide_lengths.csv"),
      row.names = FALSE
    )
  } else {
    # Create empty file to indicate analysis was performed
    write("No DNAJB1:PRKACA fusion peptides were detected in any sample.",
          file = file.path(dirs$data_dir, "fusion_peptides_not_found.txt"))
  }
  
  # Always export the theoretical peptides
  write.csv(
    results$fusion$fusion_peptides,
    file.path(dirs$data_dir, "theoretical_fusion_peptides.csv"),
    row.names = FALSE
  )
}

# Run transcriptome analysis if enabled
if (!is.null(config$analysis$do_transcriptome) && config$analysis$do_transcriptome) {
  cat("\n--- Running Transcriptome Analysis ---\n")
  results$transcriptome <- run_transcriptome_analysis(peptide_data, config)
  
  # Generate visualizations
  if (!is.null(results$transcriptome)) {
    cat("\nGenerating transcriptome visualizations...\n")
    results$transcriptome$visualizations <- generate_transcriptome_visualizations(
      results$transcriptome,
      dirs$viz_dir,
      config
    )
    
    # Export result data
    cat("\nExporting transcriptome analysis results...\n")
    
    # Export paired analysis results if available
    if (!is.null(results$transcriptome$paired_analysis)) {
      for (pair_name in names(results$transcriptome$paired_analysis)) {
        write.csv(
          results$transcriptome$paired_analysis[[pair_name]],
          file.path(dirs$data_dir, paste0("transcriptome_", pair_name, ".csv")),
          row.names = FALSE
        )
      }
    }
    
    # Export sample-specific results
    if (!is.null(results$transcriptome$sample_analysis) && 
        !is.null(results$transcriptome$sample_analysis$sample_specific)) {
      # Export summary of peptide-gene matches
      sample_match_summary <- data.frame(
        SampleID = character(),
        TotalPeptides = integer(),
        MatchedWithTranscriptome = integer(),
        MatchRate = numeric(),
        stringsAsFactors = FALSE
      )
      
      for (sample_id in names(results$transcriptome$sample_analysis$sample_specific)) {
        sample_data <- results$transcriptome$sample_analysis$sample_specific[[sample_id]]
        
        # Only add to summary if sample_data and required fields exist
        if (!is.null(sample_data) && !is.null(sample_data$summary)) {
          # Add to summary with safer access to fields
          sample_match_summary <- rbind(
            sample_match_summary,
            data.frame(
              SampleID = sample_id,
              TotalPeptides = if (!is.null(sample_data$summary$total_peptides)) sample_data$summary$total_peptides else NA,
              MatchedWithTranscriptome = if (!is.null(sample_data$summary$has_transcriptome)) sample_data$summary$has_transcriptome else NA,
              MatchRate = if (!is.null(sample_data$match_rate)) sample_data$match_rate else NA,
              stringsAsFactors = FALSE
            )
          )
        }
        
        # Export matched data if available
        if (!is.null(sample_data$matched_data) && nrow(sample_data$matched_data) > 0) {
          write.csv(
            sample_data$matched_data,
            file.path(dirs$data_dir, paste0("transcriptome_", sample_id, "_matched.csv")),
            row.names = FALSE
          )
        }
      }
      
      # Export match summary
      write.csv(
        sample_match_summary,
        file.path(dirs$data_dir, "transcriptome_match_summary.csv"),
        row.names = FALSE
      )
    }
    
    # Export shared peptide analysis results
    if (!is.null(results$transcriptome$sample_analysis) && 
        !is.null(results$transcriptome$sample_analysis$shared_analysis)) {
      
      # Export sharing thresholds
      threshold_results <- results$transcriptome$sample_analysis$shared_analysis$threshold_results
      
      for (threshold_name in names(threshold_results)) {
        threshold_data <- threshold_results[[threshold_name]]
        threshold_value <- as.numeric(gsub("shared_", "", threshold_name))
        
        # Export matched data
        if (!is.null(threshold_data$matched_data) && nrow(threshold_data$matched_data) > 0) {
          write.csv(
            threshold_data$matched_data,
            file.path(dirs$data_dir, paste0("transcriptome_shared_", threshold_value, ".csv")),
            row.names = FALSE
          )
        }
      }
    }
    
    # Export multi-omics results if available
    if (!is.null(results$transcriptome$multi_omics)) {
      # Export potential neoantigens
      if (!is.null(results$transcriptome$multi_omics$potential_neoantigens) && 
          nrow(results$transcriptome$multi_omics$potential_neoantigens) > 0) {
        write.csv(
          results$transcriptome$multi_omics$potential_neoantigens,
          file.path(dirs$data_dir, "potential_public_neoantigens.csv"),
          row.names = FALSE
        )
      }
    }
  }
}

# Run multi-omics integration if enabled
# REPLACE lines ~438-487 with this corrected code:
if (!is.null(config$analysis$do_multi_omics) && config$analysis$do_multi_omics) {
  cat("\n--- Running Multi-Omics Integration for Shared Peptides ---\n")
  
  # Load LFQ data if available
  lfq_data <- NULL
  if (!is.null(config$input$lfq_file) && file.exists(config$input$lfq_file)) {
    cat("Loading LFQ proteome data from:", config$input$lfq_file, "\n")
    tryCatch({
      lfq_data <- readxl::read_excel(config$input$lfq_file, sheet = "Significant and 1.5x_2")
      lfq_data <- process_lfq_data(lfq_data)
      cat("Loaded", nrow(lfq_data), "proteins from LFQ proteome data\n")
    }, error = function(e) {
      cat("Error loading LFQ data:", e$message, "\n")
    })
  }
  
  # Load TMT data if available
  tmt_data <- NULL
  if (!is.null(config$input$tmt_file) && file.exists(config$input$tmt_file)) {
    cat("Loading TMT proteome data from:", config$input$tmt_file, "\n")
    tryCatch({
      tmt_data <- readxl::read_excel(config$input$tmt_file, sheet = "Significant and 1.5x_2")
      tmt_data <- process_tmt_data(tmt_data)
      cat("Loaded", nrow(tmt_data), "proteins from TMT proteome data\n")
    }, error = function(e) {
      cat("Error loading TMT data:", e$message, "\n")
    })
  }
  
  # Get transcriptome data if available
  transcriptome_data <- NULL
  if (!is.null(results$transcriptome) && !is.null(results$transcriptome$transcriptome_data)) {
    transcriptome_data <- results$transcriptome$transcriptome_data
    cat("Using transcriptome data from previous analysis\n")
  }
  
  # Run multi-omics integration
  tryCatch({
    results$multi_omics <- run_multi_omics_integration(
      peptide_data,
      transcriptome_data,
      lfq_data,
      tmt_data,
      config
    )
    
    # Generate visualizations if results available
    if (!is.null(results$multi_omics)) {
      cat("\nGenerating multi-omics visualizations...\n")
      
      # Create visualization directory if it doesn't exist
      multi_omics_viz_dir <- file.path(dirs$viz_dir, "multi_omics")
      if (!dir.exists(multi_omics_viz_dir)) {
        dir.create(multi_omics_viz_dir, recursive = TRUE)
      }
      
      # Create data directory if it doesn't exist
      multi_omics_data_dir <- file.path(dirs$data_dir, "multi_omics")
      if (!dir.exists(multi_omics_data_dir)) {
        dir.create(multi_omics_data_dir, recursive = TRUE)
      }
      
      # Generate visualizations if function exists
      if (exists("generate_multi_omics_visualizations")) {
        results$multi_omics$visualizations <- generate_multi_omics_visualizations(
          results$multi_omics,
          multi_omics_viz_dir,
          config
        )
      } else {
        cat("Warning: generate_multi_omics_visualizations function not found. Skipping visualizations.\n")
      }
      
      # Export result data
      cat("\nExporting multi-omics analysis results...\n")
      if (exists("generate_multi_omics_report")) {
        results$multi_omics$report_files <- generate_multi_omics_report(
          results$multi_omics,
          multi_omics_data_dir,
          config
        )
      } else {
        cat("Warning: generate_multi_omics_report function not found. Skipping report generation.\n")
        
        # Fallback: export raw data if report generation isn't available
        if (!is.null(results$multi_omics$threshold_results)) {
          for (threshold_name in names(results$multi_omics$threshold_results)) {
            threshold_data <- results$multi_omics$threshold_results[[threshold_name]]
            write.csv(
              threshold_data,
              file.path(multi_omics_data_dir, paste0(threshold_name, "_integrated.csv")),
              row.names = FALSE
            )
          }
        }
      }
    }
  }, error = function(e) {
    cat("Error in multi-omics integration:", e$message, "\n")
    cat("Continuing with pipeline despite error\n")
  })
}

# Generate Excel reports if enabled
if (config$output$create_excel) {
  cat("\n--- Generating Excel Reports ---\n")
  if (exists("generate_excel_reports_simple")) {
    results$files$excel <- generate_excel_reports_simple(results, config)
  } else {
    results$files$excel <- generate_excel_reports(results, config)
  }
}



# Print summary
cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", dirs$base_dir, "\n")
cat("- Data files:", dirs$data_dir, "\n")
cat("- Visualizations:", dirs$viz_dir, "\n")
cat("- Reports:", dirs$report_dir, "\n")

if (config$output$create_excel) {
  cat("- Excel report:", results$files$excel, "\n")
}

cat("\nThe following analyses were performed:\n")
analyses_performed <- c()
if (config$preprocessing$concatenate_2cv_3cv) analyses_performed <- c(analyses_performed, "2CV/3CV Concatenation")
if (config$analysis$do_shared_peptide) analyses_performed <- c(analyses_performed, "Shared Peptide Analysis")
if (config$analysis$do_length_distribution) analyses_performed <- c(analyses_performed, "Peptide Length Distribution")
if (!is.null(config$analysis$sample_comparisons) && 
    !is.null(config$analysis$sample_comparisons$do_51_vs_51S) && 
    config$analysis$sample_comparisons$do_51_vs_51S) {
  analyses_performed <- c(analyses_performed, "51 vs 51S Spike-in Analysis")
}
if (!is.null(config$analysis$do_fusion) && config$analysis$do_fusion) {
  analyses_performed <- c(analyses_performed, "DNAJB1:PRKACA Fusion Analysis")
}
if (!is.null(config$analysis$do_transcriptome) && config$analysis$do_transcriptome) {
  analyses_performed <- c(analyses_performed, "Transcriptome Analysis")
  
  # Add details about transcriptome analysis
  if (!is.null(results$transcriptome)) {
    if (!is.null(results$transcriptome$paired_analysis)) {
      analyses_performed <- c(analyses_performed, "  - Paired Tumor-Normal Transcriptome")
    }
    if (!is.null(results$transcriptome$sample_analysis)) {
      analyses_performed <- c(analyses_performed, "  - Sample-Specific Transcriptome Integration")
    }
    if (!is.null(results$transcriptome$multi_omics)) {
      analyses_performed <- c(analyses_performed, "  - Multi-Omics Integration")
    }
  }
}

if (!is.null(config$analysis$do_multi_omics) && config$analysis$do_multi_omics) {
  analyses_performed <- c(analyses_performed, "Multi-Omics Integration for Shared Peptides")
  
  # Add details if available
  if (!is.null(results$multi_omics)) {
    for (threshold_name in names(results$multi_omics$threshold_results)) {
      threshold_value <- as.numeric(gsub("shared_", "", threshold_name))
      threshold_data <- results$multi_omics$threshold_results[[threshold_name]]
      
      # Check if we have neoantigen classifications
      if ("public_neoantigen_classification" %in% colnames(threshold_data)) {
        neoantigens <- threshold_data %>%
          dplyr::filter(public_neoantigen_classification != "Not a public neoantigen")
        
        if (nrow(neoantigens) > 0) {
          analyses_performed <- c(
            analyses_performed, 
            paste0("  - Identified ", nrow(neoantigens), " potential public neoantigens in peptides shared in ", 
                   threshold_value, "+ samples")
          )
        }
      }
    }
  }
}

for (analysis in analyses_performed) {
  cat("- ", analysis, "\n", sep = "")
}
