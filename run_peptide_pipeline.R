# Main script for running the peptide analysis pipeline

# Load required packages
library(tidyverse)
library(pheatmap)
library(openxlsx)

# Source core utilities
source("src/core/peptide_utils.R")
source("src/core/visualization_utils.R")

# Source analysis modules
source("src/modules/concatenate_2cv_3cv.R")
source("src/modules/shared_peptide_analysis.R")
source("src/modules/length_distribution_analysis.R")
source("src/modules/spike_in_analysis.R")
source("src/modules/fusion_analysis.R")


# Source reporting modules
# source("src/reports/generate_excel.R")
if (file.exists("src/reports/generate_excel_simple.R")) {
  source("src/reports/generate_excel_simple.R")
} else {
  source("src/reports/generate_excel.R")
}

# Get config file from command line arg or use default
args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) > 0) args[1] else "config/default_config.R"
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
if (config$analysis$do_fusion) {
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

# Generate Excel reports if enabled
if (config$output$create_excel) {
  cat("\n--- Generating Excel Reports ---\n")
  # results$files$excel <- generate_excel_reports(results, config) 
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
if (config$analysis$sample_comparisons$do_51_vs_51S) analyses_performed <- c(analyses_performed, "51 vs 51S Spike-in Analysis")
if (config$analysis$do_fusion) analyses_performed <- c(analyses_performed, "DNAJB1:PRKACA Fusion Analysis")

for (analysis in analyses_performed) {
  cat("- ", analysis, "\n", sep = "")
}

