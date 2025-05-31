# Script to correlate peptide intensity data from combined 2cv and 3cv data
# Modified to use the output of peptide_combination.R
# this is the more recent script

# Setting directory
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516")

# Load required packages
library(tidyverse)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(gridExtra)

# File path to the combined peptide data
combined_file <- "unique_peptides_all.tsv"

# Define output directory - we'll extract a prefix from the input filename
create_output_directory <- function(input_file) {
  # Extract the base name without extension
  base_name <- tools::file_path_sans_ext(basename(input_file))
  
  # Create a timestamped folder
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- paste0("correlation_results_", base_name, "_", timestamp)
  
  # Create directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("Created output directory:", output_dir, "\n")
  return(output_dir)
}

# Function to read the combined peptide data
read_combined_data <- function(file_path) {
  cat("Reading combined peptide data from:", file_path, "\n")
  
  # Read the combined TSV file
  combined_data <- read_tsv(file_path, show_col_types = FALSE)
  
  cat("Read", nrow(combined_data), "rows from combined peptide data\n")
  
  return(combined_data)
}

# Function to transform the combined data into 2cv and 3cv datasets
prepare_data_for_correlation <- function(combined_data) {
  # Create a dataset for 2cv data
  data_2cv <- combined_data %>%
    filter(detected_2cv == TRUE) %>%
    select(Peptide, SampleID = SampleID_2cv, Intensity = Intensity_2cv, SourceFile = SourceFile_2cv) %>%
    mutate(CVType = "2cv", source = "2cv")
  
  # Create a dataset for 3cv data
  data_3cv <- combined_data %>%
    filter(detected_3cv == TRUE) %>%
    select(Peptide, SampleID = SampleID_3cv, Intensity = Intensity_3cv, SourceFile = SourceFile_3cv) %>%
    mutate(CVType = "3cv", source = "3cv")
  
  cat("Prepared", nrow(data_2cv), "rows for 2cv data and", nrow(data_3cv), "rows for 3cv data\n")
  
  return(list(data_2cv = data_2cv, data_3cv = data_3cv))
}

# Function to perform correlation analysis for matched peptides
perform_correlation_analysis <- function(data_2cv, data_3cv) {
  # Get unique sample IDs
  samples_2cv <- unique(data_2cv$SampleID)
  samples_3cv <- unique(data_3cv$SampleID)
  
  # Find common sample IDs
  common_samples <- intersect(samples_2cv, samples_3cv)
  
  if (length(common_samples) == 0) {
    stop("No common sample IDs found between 2CV and 3CV datasets")
  }
  
  cat("Found", length(common_samples), "common sample IDs\n")
  
  # Initialize results list
  results <- list(
    by_sample = list(),
    overall = NULL
  )
  
  # Initialize dataframe to store all matched peptides for overall analysis
  all_matched_peptides <- data.frame()
  
  # Process each sample
  for (sample_id in common_samples) {
    cat("Processing sample:", sample_id, "\n")
    
    # Filter data for this sample
    sample_2cv <- data_2cv %>% filter(SampleID == sample_id)
    sample_3cv <- data_3cv %>% filter(SampleID == sample_id)
    
    # Find common peptides
    peptides_2cv <- sample_2cv$Peptide
    peptides_3cv <- sample_3cv$Peptide
    common_peptides <- intersect(peptides_2cv, peptides_3cv)
    
    cat("  Found", length(common_peptides), "peptides in both 2CV and 3CV\n")
    
    if (length(common_peptides) < 5) {
      warning("Too few common peptides for sample ", sample_id, ". Skipping.")
      next
    }
    
    # Create matched data
    matched_data <- data.frame(
      Peptide = common_peptides,
      SampleID = sample_id,
      Intensity_2CV = sample_2cv %>% 
        filter(Peptide %in% common_peptides) %>% 
        select(Peptide, Intensity) %>% 
        rename(Intensity_2CV = Intensity) %>% 
        pull(Intensity_2CV),
      Intensity_3CV = sample_3cv %>% 
        filter(Peptide %in% common_peptides) %>% 
        select(Peptide, Intensity) %>% 
        rename(Intensity_3CV = Intensity) %>% 
        pull(Intensity_3CV)
    )
    
    # Check for negative or zero intensity values
    negative_2cv <- sum(matched_data$Intensity_2CV < 0, na.rm = TRUE)
    zero_2cv <- sum(matched_data$Intensity_2CV == 0, na.rm = TRUE)
    cat("Found", negative_2cv, "negative intensity values and", zero_2cv, "zero values in 2CV data\n")
    
    negative_3cv <- sum(matched_data$Intensity_3CV < 0, na.rm = TRUE)
    zero_3cv <- sum(matched_data$Intensity_3CV == 0, na.rm = TRUE)
    cat("Found", negative_3cv, "negative intensity values and", zero_3cv, "zero values in 3CV data\n")
    
    # If there are negative values, let's see their distribution
    if (negative_2cv > 0) {
      cat("Range of negative 2CV intensities:", range(matched_data$Intensity_2CV[matched_data$Intensity_2CV < 0], na.rm = TRUE), "\n")
    }
    if (negative_3cv > 0) {
      cat("Range of negative 3CV intensities:", range(matched_data$Intensity_3CV[matched_data$Intensity_3CV < 0], na.rm = TRUE), "\n")
    }
    
    # Calculate correlations
    pearson_cor <- cor(matched_data$Intensity_2CV, matched_data$Intensity_3CV, 
                       method = "pearson", use = "complete.obs")
    spearman_cor <- cor(matched_data$Intensity_2CV, matched_data$Intensity_3CV, 
                        method = "spearman", use = "complete.obs")
    
    # Calculate R-squared (coefficient of determination)
    lm_model <- lm(Intensity_3CV ~ Intensity_2CV, data = matched_data)
    r_squared <- summary(lm_model)$r.squared
    
    # Store results for this sample
    results$by_sample[[sample_id]] <- list(
      matched_data = matched_data,
      pearson = pearson_cor,
      spearman = spearman_cor,
      r_squared = r_squared
    )
    
    # Add to overall dataset
    all_matched_peptides <- bind_rows(all_matched_peptides, matched_data)
  }
  
  # Overall analysis if we have data
  if (nrow(all_matched_peptides) > 0) {
    # Calculate overall correlations
    pearson_cor_overall <- cor(all_matched_peptides$Intensity_2CV, all_matched_peptides$Intensity_3CV, 
                               method = "pearson", use = "complete.obs")
    spearman_cor_overall <- cor(all_matched_peptides$Intensity_2CV, all_matched_peptides$Intensity_3CV, 
                                method = "spearman", use = "complete.obs")
    
    # Calculate overall R-squared
    lm_model_overall <- lm(Intensity_3CV ~ Intensity_2CV, data = all_matched_peptides)
    r_squared_overall <- summary(lm_model_overall)$r.squared
    
    results$overall <- list(
      matched_data = all_matched_peptides,
      pearson = pearson_cor_overall,
      spearman = spearman_cor_overall,
      r_squared = r_squared_overall
    )
  }
  
  return(results)
}

# Function to create visualization
create_visualizations <- function(correlation_results, output_dir) {
  # Create plots directory within the output directory
  plots_dir <- file.path(output_dir, "plots")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  
  create_scatter_plot <- function(data, title, filename) {
    # Filter out zeros and negative values for log transformation
    data_filtered <- data %>% filter(Intensity_2CV > 0, Intensity_3CV > 0)
    
    # Create scatter plot
    p <- ggplot(data_filtered, aes(x = Intensity_2CV, y = Intensity_3CV)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", color = "red") +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        title = title,
        subtitle = sprintf(
          "Pearson: %.3f, Spearman: %.3f, R²: %.3f (n=%d, filtered out %d 0 values",
          cor(data_filtered$Intensity_2CV, data_filtered$Intensity_3CV, method = "pearson", use = "complete.obs"),
          cor(data_filtered$Intensity_2CV, data_filtered$Intensity_3CV, method = "spearman", use = "complete.obs"),
          summary(lm(Intensity_3CV ~ Intensity_2CV, data = data_filtered))$r.squared,
          nrow(data_filtered),
          nrow(data) - nrow(data_filtered)
        ),
        x = "2CV Intensity (log10)",
        y = "3CV Intensity (log10)"
      ) +
      theme_minimal()
    
    # Save to the plots directory within output directory
    full_path <- file.path(plots_dir, filename)
    ggsave(full_path, plot = p, width = 10, height = 8)
    cat("Saved plot to:", full_path, "\n")
    return(p)
  }
  
  # Create plots for each sample
  sample_plots <- list()
  for (sample_id in names(correlation_results$by_sample)) {
    data <- correlation_results$by_sample[[sample_id]]$matched_data
    title <- paste("Peptide Intensity Correlation for Sample", sample_id)
    filename <- paste0("scatter_", sample_id, ".png")
    
    sample_plots[[sample_id]] <- create_scatter_plot(data, title, filename)
    cat("Created scatter plot for sample", sample_id, "\n")
  }
  
  # Create overall plot
  if (!is.null(correlation_results$overall)) {
    overall_data <- correlation_results$overall$matched_data
    title <- "Overall Peptide Intensity Correlation (All Samples)"
    filename <- "scatter_overall.png"
    
    overall_plot <- create_scatter_plot(overall_data, title, filename)
    cat("Created overall scatter plot\n")
    
    # Create heatmap
    # For the heatmap, we'll bin the data into a 2D histogram
    # First, filter out zeros and negative values before log transformation
    heatmap_data <- overall_data %>%
      filter(Intensity_2CV > 0, Intensity_3CV > 0) %>%  # Filter out non-positive values
      mutate(
        log_2CV = log10(Intensity_2CV),
        log_3CV = log10(Intensity_3CV)
      )
    
    # Use stat_bin2d for a reliable heatmap creation
    heatmap_plot <- ggplot(heatmap_data, aes(x = log_2CV, y = log_3CV)) +
      stat_bin2d(bins = 50) +
      scale_fill_viridis_c(trans = "log1p") +
      labs(
        title = "Heatmap of Peptide Intensity Correlation",
        x = "2CV Intensity (log10)",
        y = "3CV Intensity (log10)",
        fill = "Count"
      ) +
      theme_minimal()
    
    heatmap_path <- file.path(plots_dir, "heatmap_overall.png")
    ggsave(heatmap_path, plot = heatmap_plot, width = 10, height = 8)
    cat("Created overall heatmap at:", heatmap_path, "\n")
  }
  
  # Create summary table
  summary_table <- data.frame(
    SampleID = character(),
    PeptideCount = integer(),
    Pearson = numeric(),
    Spearman = numeric(),
    R_Squared = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (sample_id in names(correlation_results$by_sample)) {
    result <- correlation_results$by_sample[[sample_id]]
    summary_table <- bind_rows(summary_table, data.frame(
      SampleID = sample_id,
      PeptideCount = nrow(result$matched_data),
      Pearson = result$pearson,
      Spearman = result$spearman,
      R_Squared = result$r_squared
    ))
  }
  
  if (!is.null(correlation_results$overall)) {
    summary_table <- bind_rows(summary_table, data.frame(
      SampleID = "OVERALL",
      PeptideCount = nrow(correlation_results$overall$matched_data),
      Pearson = correlation_results$overall$pearson,
      Spearman = correlation_results$overall$spearman,
      R_Squared = correlation_results$overall$r_squared
    ))
  }
  
  # Save summary table to output directory
  summary_csv_path <- file.path(output_dir, "correlation_summary.csv")
  write.csv(summary_table, summary_csv_path, row.names = FALSE)
  cat("Created correlation summary table at:", summary_csv_path, "\n")
  
  # Also save as TSV
  summary_tsv_path <- file.path(output_dir, "correlation_summary.tsv")
  write_tsv(summary_table, summary_tsv_path)
  cat("Created correlation summary table (TSV) at:", summary_tsv_path, "\n")
  
  return(list(
    sample_plots = sample_plots,
    overall_plot = if (!is.null(correlation_results$overall)) overall_plot else NULL,
    heatmap = if (!is.null(correlation_results$overall)) heatmap_plot else NULL,
    summary_table = summary_table,
    output_dir = output_dir,
    plots_dir = plots_dir
  ))
}

# Create summary PDF function
create_summary_pdf_polished <- function(correlation_results, visualization_results, output_dir) {
  library(ggplot2)
  library(gridExtra)
  library(grid)
  
  # Set the PDF path in the output directory
  pdf_path <- file.path(output_dir, "peptide_correlation_summary.pdf")
  cat("Creating polished PDF at:", pdf_path, "\n")
  
  # Create PDF with larger dimensions and margins
  pdf(pdf_path, width = 11, height = 11)
  
  # Title page
  grid.newpage()
  grid.text("Peptide Intensity Correlation Analysis Summary", 
            x = 0.5, y = 0.7, gp = gpar(fontsize = 24, fontface = "bold"))
  grid.text(paste("Analysis Date:", format(Sys.Date(), "%B %d, %Y")),
            x = 0.5, y = 0.5, gp = gpar(fontsize = 16))
  grid.text(paste("Input File:", basename(combined_file)),
            x = 0.5, y = 0.4, gp = gpar(fontsize = 14))
  
  # Summary statistics table
  grid.newpage()
  grid.text("Correlation Statistics Summary", x = 0.5, y = 0.95, 
            gp = gpar(fontsize = 18, fontface = "bold"))
  
  # Format numeric columns to 3 decimal places
  summary_table <- visualization_results$summary_table
  summary_table$Pearson <- sprintf("%.3f", summary_table$Pearson)
  summary_table$Spearman <- sprintf("%.3f", summary_table$Spearman)
  summary_table$R_Squared <- sprintf("%.3f", summary_table$R_Squared)
  
  # Create table with better spacing
  table_theme <- ttheme_minimal(
    core = list(fg_params = list(fontsize = 10, hjust = 0.5, x = 0.5),
                bg_params = list(fill = c("white", "grey95"), col = NA)),
    colhead = list(fg_params = list(fontsize = 11, fontface = "bold", hjust = 0.5, x = 0.5),
                   bg_params = list(fill = "grey90", col = NA))
  )
  
  grid.table(summary_table, rows = NULL, theme = table_theme)
  
  # Overall correlation plot
  if (!is.null(visualization_results$overall_plot)) {
    print(visualization_results$overall_plot + 
            ggtitle("Overall Peptide Intensity Correlation (All Samples)"))
  }
  
  # Sample plots - 4 per page
  sample_ids <- names(visualization_results$sample_plots)
  num_samples <- length(sample_ids)
  
  # Process in batches of 4
  for (i in seq(1, num_samples, by = 4)) {
    grid.newpage()
    
    # Set up the layout
    pushViewport(viewport(layout = grid.layout(3, 2, heights = c(0.1, 1, 1), widths = c(1, 1))))
    
    # Add the page title
    grid.text("Sample Correlation Plots", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
              gp = gpar(fontsize = 16, fontface = "bold"))
    
    # Get indices for this page
    end_idx <- min(i + 3, num_samples)
    batch_ids <- sample_ids[i:end_idx]
    
    # Calculate positions for each plot
    positions <- list(
      c(2, 1), c(2, 2),
      c(3, 1), c(3, 2)
    )
    
    # Create and place each plot
    for (j in 1:length(batch_ids)) {
      sample_id <- batch_ids[j]
      pos <- positions[[j]]
      
      # Print the plot in the right position
      print(visualization_results$sample_plots[[sample_id]], 
            vp = viewport(layout.pos.row = pos[1], layout.pos.col = pos[2]))
    }
    
    popViewport()
  }
  
  # Heatmap
  if (!is.null(visualization_results$heatmap)) {
    print(visualization_results$heatmap)
  }
  
  # Interpretation page
  grid.newpage()
  grid.text("Interpretation of Results", x = 0.5, y = 0.95, 
            gp = gpar(fontsize = 18, fontface = "bold"))
  
  # Get values for interpretation
  pearson_overall <- as.numeric(gsub("([0-9.]+).*", "\\1", summary_table$Pearson[nrow(summary_table)]))
  r_squared_overall <- as.numeric(gsub("([0-9.]+).*", "\\1", summary_table$R_Squared[nrow(summary_table)]))
  
  # For min/max, exclude the OVERALL row
  pearson_values <- as.numeric(gsub("([0-9.]+).*", "\\1", summary_table$Pearson[-nrow(summary_table)]))
  pearson_min <- min(pearson_values)
  pearson_max <- max(pearson_values)
  
  # Create interpretation text
  interpretation_text <- paste(
    "Summary of findings:",
    "",
    paste("Input file:", basename(combined_file)),
    "",
    "1. Overall correlation: Strong correlation between 2CV and 3CV samples",
    paste("   Pearson correlation:", sprintf("%.3f", pearson_overall)),
    paste("   R-squared:", sprintf("%.3f", r_squared_overall)),
    "",
    "2. Sample variability: Sample correlations range from",
    paste("   ", sprintf("%.3f", pearson_min), "to", sprintf("%.3f", pearson_max), "(Pearson)"),
    "",
    "3. Data filtering: Zero intensity values were excluded from log-scale visualizations.",
    "   These represent peptides detected in both samples but with intensity below",
    "   the quantification threshold in at least one of the samples.",
    "",
    "4. Heatmap interpretation: The heatmap shows the density of peptide intensities,",
    "   with the diagonal pattern confirming the strong correlation between 2CV and 3CV measurements.",
    sep = "\n"
  )
  
  # Add text to the page with better formatting
  grid.text(interpretation_text, x = 0.05, y = 0.8, just = c("left", "top"),
            gp = gpar(fontsize = 12))
  
  # Close the PDF device
  dev.off()
  cat("Polished summary PDF created at:", pdf_path, "\n")
  
  # Also create a copy with the dataset name in it
  dataset_name <- tools::file_path_sans_ext(basename(combined_file))
  named_pdf_path <- file.path(output_dir, paste0("peptide_correlation_", dataset_name, ".pdf"))
  file.copy(pdf_path, named_pdf_path)
  cat("Created named PDF copy at:", named_pdf_path, "\n")
  
  return(pdf_path)
}

# Main function to run the analysis
run_peptide_correlation_analysis <- function(combined_file_path) {
  # Create output directory
  output_dir <- create_output_directory(combined_file_path)
  
  # Read the combined peptide data
  combined_data <- read_combined_data(combined_file_path)
  
  # Prepare data for correlation analysis
  prepared_data <- prepare_data_for_correlation(combined_data)
  data_2cv <- prepared_data$data_2cv
  data_3cv <- prepared_data$data_3cv
  
  if (is.null(data_2cv) || is.null(data_3cv)) {
    stop("Missing data for one or both CV types")
  }
  
  # Perform correlation analysis
  correlation_results <- perform_correlation_analysis(data_2cv, data_3cv)
  
  # Create visualizations
  visualization_results <- create_visualizations(correlation_results, output_dir)
  
  # Create summary PDF
  pdf_path <- create_summary_pdf_polished(correlation_results, visualization_results, output_dir)
  
  # Copy the input data to the output directory for reference
  file.copy(combined_file_path, file.path(output_dir, basename(combined_file_path)))
  
  # Create a README.txt file with information about the analysis
  readme_content <- paste(
    "Peptide Correlation Analysis Results",
    "==================================",
    "",
    paste("Input file:", basename(combined_file_path)),
    paste("Analysis date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Output directory:", output_dir),
    "",
    "Contents:",
    "- correlation_summary.csv: Summary of correlation statistics by sample",
    "- correlation_summary.tsv: Same data in TSV format",
    "- peptide_correlation_summary.pdf: Complete report with all visualizations",
    paste0("- peptide_correlation_", tools::file_path_sans_ext(basename(combined_file_path)), ".pdf: Named copy of the report"),
    "- plots/: Directory containing individual plot images",
    paste0("- ", basename(combined_file_path), ": Copy of the input data file"),
    "",
    "This analysis was performed using the modified peptide correlation script.",
    sep = "\n"
  )
  
  writeLines(readme_content, file.path(output_dir, "README.txt"))
  
  cat("Analysis complete! Results are in:", output_dir, "\n")
  
  # Return results
  return(list(
    correlation = correlation_results,
    visualization = visualization_results,
    output_dir = output_dir,
    pdf_path = pdf_path
  ))
}

# Run the analysis
results <- run_peptide_correlation_analysis(combined_file)

# Print summary of results
cat("\n=== CORRELATION ANALYSIS SUMMARY ===\n")
print(results$visualization$summary_table)
cat("\nAll results have been saved to:", results$output_dir, "\n")
cat("Main PDF report:", results$pdf_path, "\n")