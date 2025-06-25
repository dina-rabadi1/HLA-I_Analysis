# ===============================
# INTENSITY ANALYSIS FUNCTIONS
# ===============================
# Functions for comprehensive intensity analysis of immunopeptidomics data
# # immunopeptidomics_intensity_analysis_functions.R
# Author: Generated for immunopeptidomics analysis

# Function to calculate comprehensive statistics
calculate_statistics <- function(values, name) {
  # Remove NA values for calculations
  clean_values <- values[!is.na(values)]
  
  if (length(clean_values) == 0) {
    return(data.frame(
      Dataset = name,
      N = 0,
      Mean = NA,
      Median = NA,
      Mode = NA,
      Min = NA,
      Max = NA,
      Range = NA,
      Q1 = NA,
      Q3 = NA,
      IQR = NA,
      SD = NA,
      CV = NA,
      Skewness = NA
    ))
  }
  
  # Calculate mode (most frequent value)
  mode_val <- as.numeric(names(sort(table(clean_values), decreasing = TRUE))[1])
  
  # Calculate skewness manually
  skewness <- mean((clean_values - mean(clean_values))^3) / (sd(clean_values)^3)
  
  stats <- data.frame(
    Dataset = name,
    N = length(clean_values),
    Mean = mean(clean_values),
    Median = median(clean_values),
    Mode = mode_val,
    Min = min(clean_values),
    Max = max(clean_values),
    Range = max(clean_values) - min(clean_values),
    Q1 = quantile(clean_values, 0.25),
    Q3 = quantile(clean_values, 0.75),
    IQR = IQR(clean_values),
    SD = sd(clean_values),
    CV = sd(clean_values) / mean(clean_values) * 100,
    Skewness = skewness
  )
  
  return(stats)
}

# Function to prepare intensity data based on data source
prepare_intensity_data <- function(data_source, combined_data = NULL, file_path = NULL, folder_path = NULL, exclude_samples = NULL, include_samples = NULL) {
  
  if (data_source == "combined_data") {
    if (is.null(combined_data)) {
      stop("ERROR: combined_data is NULL but data_source is set to 'combined_data'")
    }
    data <- combined_data
    cat("Using combined_data from pipeline (", nrow(data), "entries)\n")
    
  } else if (data_source == "file") {
    if (is.null(file_path) || !file.exists(file_path)) {
      stop("ERROR: File path not provided or file not found: ", file_path)
    }
    cat("Loading data from file:", file_path, "\n")
    data <- read_tsv(file_path, show_col_types = FALSE)
    
  } else if (data_source == "folder") {
    if (is.null(folder_path) || !dir.exists(folder_path)) {
      stop("ERROR: Folder path not provided or folder not found: ", folder_path)
    }
    cat("Loading data from folder:", folder_path, "\n")
    
    # Get all TSV files in the folder
    tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)
    if (length(tsv_files) == 0) {
      stop("ERROR: No TSV files found in folder: ", folder_path)
    }
    
    cat("Found", length(tsv_files), "TSV files\n")
    
    # Read and combine all files
    data_list <- list()
    for (i in seq_along(tsv_files)) {
      cat("Reading file", i, "of", length(tsv_files), ":", basename(tsv_files[i]), "\n")
      temp_data <- read_tsv(tsv_files[i], show_col_types = FALSE)
      temp_data$SourceFile <- basename(tsv_files[i])
      data_list[[i]] <- temp_data
    }
    
    data <- bind_rows(data_list)
    cat("Combined data:", nrow(data), "total entries\n")
    
  } else {
    stop("ERROR: Invalid data_source. Must be 'combined_data', 'file', or 'folder'")
  }
  
  # Check required columns
  required_cols <- c("final_intensity", "Spectral Count", "SampleID", "Peptide")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("ERROR: Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Apply sample filtering
  original_count <- nrow(data)
  
  # Apply inclusion filter first
  if (!is.null(include_samples)) {
    missing_samples <- setdiff(include_samples, unique(data$SampleID))
    if (length(missing_samples) > 0) {
      warning("WARNING: Requested samples not found in data: ", paste(missing_samples, collapse = ", "))
    }
    data <- data %>% filter(SampleID %in% include_samples)
    cat("After including samples", paste(include_samples, collapse = ", "), ":", nrow(data), "peptides\n")
  }
  
  # Apply exclusion filter
  if (!is.null(exclude_samples)) {
    found_samples <- intersect(exclude_samples, unique(data$SampleID))
    if (length(found_samples) > 0) {
      data <- data %>% filter(!SampleID %in% exclude_samples)
      cat("After excluding samples", paste(found_samples, collapse = ", "), ":", nrow(data), "peptides\n")
    } else {
      cat("Note: Excluded samples not found in dataset\n")
    }
  }
  
  # Final check
  if (nrow(data) == 0) {
    stop("ERROR: No data remaining after filtering!")
  }
  
  cat("Final dataset:", nrow(data), "peptides from", length(unique(data$SampleID)), "samples\n")
  
  return(data)
}

# Function to perform comprehensive intensity analysis
perform_intensity_analysis <- function(data, output_dir, dataset_name, timestamp) {
  
  cat("\n=== INTENSITY ANALYSIS ===\n")
  
  # Create analysis subdirectory
  intensity_dir <- file.path(output_dir, "intensity_analysis")
  plots_dir <- file.path(intensity_dir, "plots")
  tables_dir <- file.path(intensity_dir, "tables")
  
  dir.create(intensity_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Sample breakdown analysis
  cat("Generating sample breakdown...\n")
  sample_summary <- data %>%
    group_by(SampleID) %>%
    summarise(
      Total_Peptides = n(),
      Zero_Intensity = sum(final_intensity == 0 | is.na(final_intensity), na.rm = TRUE),
      NonZero_Intensity = sum(final_intensity > 0, na.rm = TRUE),
      Pct_NonZero = round(NonZero_Intensity / Total_Peptides * 100, 1),
      Mean_Intensity = round(mean(final_intensity, na.rm = TRUE), 0),
      Median_Intensity = round(median(final_intensity, na.rm = TRUE), 0),
      Mean_SpectralCount = round(mean(`Spectral Count`, na.rm = TRUE), 2),
      Median_SpectralCount = median(`Spectral Count`, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(Total_Peptides))
  
  # Missing/zero value analysis
  total_peptides <- nrow(data)
  missing_intensity <- sum(is.na(data$final_intensity))
  zero_intensity <- sum(data$final_intensity == 0, na.rm = TRUE)
  non_zero_intensity <- sum(data$final_intensity > 0, na.rm = TRUE)
  
  missing_summary <- data.frame(
    Category = c("Total Peptides", "Missing Intensity", "Zero Intensity", "Non-zero Intensity"),
    Count = c(total_peptides, missing_intensity, zero_intensity, non_zero_intensity),
    Percentage = c(100, 
                   round(missing_intensity/total_peptides*100, 2),
                   round(zero_intensity/total_peptides*100, 2),
                   round(non_zero_intensity/total_peptides*100, 2))
  )
  
  # Prepare datasets
  all_data <- data
  nonzero_data <- data %>% filter(final_intensity > 0 & !is.na(final_intensity))
  
  if (nrow(nonzero_data) == 0) {
    stop("ERROR: No peptides with non-zero intensity found!")
  }
  
  # Statistical analysis
  cat("Calculating comprehensive statistics...\n")
  stats_all_raw <- calculate_statistics(all_data$final_intensity, "All_Peptides_Raw")
  stats_nonzero_raw <- calculate_statistics(nonzero_data$final_intensity, "NonZero_Peptides_Raw")
  
  # Log transform non-zero values
  nonzero_data$log_final_intensity <- log10(nonzero_data$final_intensity)
  stats_nonzero_log <- calculate_statistics(nonzero_data$log_final_intensity, "NonZero_Peptides_Log10")
  
  # Spectral count analysis
  spectral_stats_all <- calculate_statistics(all_data$`Spectral Count`, "All_Peptides_SpectralCount")
  spectral_stats_nonzero <- calculate_statistics(nonzero_data$`Spectral Count`, "NonZero_Peptides_SpectralCount")
  
  # Combine all statistics
  all_stats <- rbind(stats_all_raw, stats_nonzero_raw, stats_nonzero_log, spectral_stats_all, spectral_stats_nonzero)
  
  # Calculate percentiles
  percentiles <- quantile(nonzero_data$final_intensity, 
                          probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99), 
                          na.rm = TRUE)
  
  log_percentiles <- quantile(nonzero_data$log_final_intensity, 
                              probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99), 
                              na.rm = TRUE)
  
  # Threshold definitions
  thresholds <- data.frame(
    Threshold_Type = c("Q1 (Low)", "Q3 (High)", 
                       "10th Percentile (Low)", "90th Percentile (High)",
                       "5th Percentile (Very Low)", "95th Percentile (Very High)",
                       "Mean - 1SD (Low)", "Mean + 1SD (High)",
                       "Mean - 2SD (Very Low)", "Mean + 2SD (Very High)"),
    Raw_Intensity = c(quantile(nonzero_data$final_intensity, 0.25, na.rm = TRUE),
                      quantile(nonzero_data$final_intensity, 0.75, na.rm = TRUE),
                      quantile(nonzero_data$final_intensity, 0.10, na.rm = TRUE),
                      quantile(nonzero_data$final_intensity, 0.90, na.rm = TRUE),
                      quantile(nonzero_data$final_intensity, 0.05, na.rm = TRUE),
                      quantile(nonzero_data$final_intensity, 0.95, na.rm = TRUE),
                      mean(nonzero_data$final_intensity, na.rm = TRUE) - sd(nonzero_data$final_intensity, na.rm = TRUE),
                      mean(nonzero_data$final_intensity, na.rm = TRUE) + sd(nonzero_data$final_intensity, na.rm = TRUE),
                      mean(nonzero_data$final_intensity, na.rm = TRUE) - 2*sd(nonzero_data$final_intensity, na.rm = TRUE),
                      mean(nonzero_data$final_intensity, na.rm = TRUE) + 2*sd(nonzero_data$final_intensity, na.rm = TRUE)),
    Log10_Intensity = c(quantile(nonzero_data$log_final_intensity, 0.25, na.rm = TRUE),
                        quantile(nonzero_data$log_final_intensity, 0.75, na.rm = TRUE),
                        quantile(nonzero_data$log_final_intensity, 0.10, na.rm = TRUE),
                        quantile(nonzero_data$log_final_intensity, 0.90, na.rm = TRUE),
                        quantile(nonzero_data$log_final_intensity, 0.05, na.rm = TRUE),
                        quantile(nonzero_data$log_final_intensity, 0.95, na.rm = TRUE),
                        mean(nonzero_data$log_final_intensity, na.rm = TRUE) - sd(nonzero_data$log_final_intensity, na.rm = TRUE),
                        mean(nonzero_data$log_final_intensity, na.rm = TRUE) + sd(nonzero_data$log_final_intensity, na.rm = TRUE),
                        mean(nonzero_data$log_final_intensity, na.rm = TRUE) - 2*sd(nonzero_data$log_final_intensity, na.rm = TRUE),
                        mean(nonzero_data$log_final_intensity, na.rm = TRUE) + 2*sd(nonzero_data$log_final_intensity, na.rm = TRUE))
  )
  
  # Spectral count analysis
  correlation <- cor(nonzero_data$final_intensity, nonzero_data$`Spectral Count`, use = "complete.obs")
  
  spectral_table <- table(nonzero_data$`Spectral Count`)
  spectral_df <- data.frame(
    Spectral_Count = as.numeric(names(spectral_table)),
    Count = as.numeric(spectral_table),
    Percentage = round(as.numeric(spectral_table) / sum(spectral_table) * 100, 2)
  )
  
  # Cross-tabulation
  intensity_quartiles <- quantile(nonzero_data$final_intensity, c(0.25, 0.5, 0.75), na.rm = TRUE)
  nonzero_data$intensity_category <- cut(nonzero_data$final_intensity, 
                                         breaks = c(0, intensity_quartiles, Inf),
                                         labels = c("Low", "Medium", "High", "Very High"),
                                         include.lowest = TRUE)
  
  nonzero_data$spectral_category <- cut(nonzero_data$`Spectral Count`,
                                        breaks = c(0, 1, 2, 5, Inf),
                                        labels = c("Single", "Double", "Multiple", "High"),
                                        include.lowest = TRUE)
  
  crosstab <- table(nonzero_data$intensity_category, nonzero_data$spectral_category)
  
  # Save all tables
  cat("Saving intensity analysis tables...\n")
  write_csv(sample_summary, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_01_sample_breakdown_summary.csv")))
  write_csv(missing_summary, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_02_missing_value_summary.csv")))
  write_csv(all_stats, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_03_comprehensive_statistics.csv")))
  write_csv(data.frame(Percentile = names(percentiles), Raw_Value = as.numeric(percentiles)), 
            file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_04_raw_intensity_percentiles.csv")))
  write_csv(data.frame(Percentile = names(log_percentiles), Log10_Value = as.numeric(log_percentiles)), 
            file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_05_log_intensity_percentiles.csv")))
  write_csv(thresholds, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_06_intensity_thresholds.csv")))
  write_csv(spectral_df, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_07_spectral_count_distribution.csv")))
  write_csv(as.data.frame.matrix(crosstab), file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_08_intensity_spectral_crosstab.csv")))
  
  # Return results for plotting
  results <- list(
    data = data,
    nonzero_data = nonzero_data,
    sample_summary = sample_summary,
    missing_summary = missing_summary,
    all_stats = all_stats,
    correlation = correlation,
    total_peptides = total_peptides,
    non_zero_intensity = non_zero_intensity,
    zero_intensity = zero_intensity,
    missing_intensity = missing_intensity,
    plots_dir = plots_dir,
    tables_dir = tables_dir
  )
  
  return(results)
}

# Function to create intensity analysis visualizations
create_intensity_analysis_plots <- function(analysis_results, dataset_name, timestamp) {
  
  cat("Creating intensity analysis visualizations...\n")
  
  # Extract data from results
  data <- analysis_results$data
  nonzero_data <- analysis_results$nonzero_data
  plots_dir <- analysis_results$plots_dir
  correlation <- analysis_results$correlation
  total_peptides <- analysis_results$total_peptides
  non_zero_intensity <- analysis_results$non_zero_intensity
  zero_intensity <- analysis_results$zero_intensity
  missing_intensity <- analysis_results$missing_intensity
  
  # Create plots list
  plots <- list()
  
  # 1. Pie chart for zero vs non-zero intensity
  pie_data <- data.frame(
    Category = c("Non-zero Intensity", "Zero/Missing Intensity"),
    Count = c(non_zero_intensity, zero_intensity + missing_intensity),
    Percentage = c(round(non_zero_intensity/total_peptides*100, 1),
                   round((zero_intensity + missing_intensity)/total_peptides*100, 1))
  )
  
  plots$pie_intensity <- ggplot(pie_data, aes(x = "", y = Count, fill = Category)) +
    geom_col(width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("Non-zero Intensity" = "#4CAF50", "Zero/Missing Intensity" = "#F44336")) +
    labs(title = "Peptide Distribution by Intensity Status",
         subtitle = paste("Total peptides:", total_peptides)) +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom") +
    geom_text(aes(label = paste0(Percentage, "%\n(", Count, ")")), 
              position = position_stack(vjust = 0.5), size = 4)
  
  # 2. Histogram of raw intensity (non-zero)
  plots$hist_raw <- ggplot(nonzero_data, aes(x = final_intensity)) +
    geom_histogram(bins = 50, fill = "skyblue", alpha = 0.7, color = "black") +
    scale_x_log10(labels = scientific) +
    labs(title = "Distribution of Final Intensity (Non-zero values)",
         subtitle = paste("n =", nrow(nonzero_data), "peptides"),
         x = "Final Intensity (log10 scale)",
         y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # 3. Histogram of log-transformed intensity
  plots$hist_log <- ggplot(nonzero_data, aes(x = log_final_intensity)) +
    geom_histogram(bins = 50, fill = "lightcoral", alpha = 0.7, color = "black") +
    labs(title = "Distribution of Log10(Final Intensity)",
         subtitle = paste("n =", nrow(nonzero_data), "peptides"),
         x = "Log10(Final Intensity)",
         y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # 4. Box plot of non-zero intensities
  plots$boxplot_intensity <- ggplot(nonzero_data, aes(x = "Non-zero Intensities", y = final_intensity)) +
    geom_boxplot(fill = "lightgreen", alpha = 0.7) +
    scale_y_log10(labels = scientific) +
    labs(title = "Box Plot of Final Intensity",
         subtitle = "Non-zero values only",
         x = "",
         y = "Final Intensity (log10 scale)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # 5. Density plot of log intensity
  plots$density_log <- ggplot(nonzero_data, aes(x = log_final_intensity)) +
    geom_density(fill = "yellow", alpha = 0.5, color = "orange") +
    geom_vline(xintercept = mean(nonzero_data$log_final_intensity, na.rm = TRUE), 
               color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = median(nonzero_data$log_final_intensity, na.rm = TRUE), 
               color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Density Plot of Log10(Final Intensity)",
         subtitle = "Red line: Mean, Blue line: Median",
         x = "Log10(Final Intensity)",
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # 6. Spectral count histogram
  plots$spectral_hist <- ggplot(nonzero_data, aes(x = `Spectral Count`)) +
    geom_histogram(bins = min(30, max(nonzero_data$`Spectral Count`, na.rm = TRUE)), 
                   fill = "mediumpurple", alpha = 0.7, color = "black") +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(title = "Distribution of Spectral Counts",
         subtitle = paste("n =", nrow(nonzero_data), "peptides with non-zero intensity"),
         x = "Spectral Count",
         y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # 7. Intensity vs Spectral Count scatter plot
  plots$scatter_intensity_spectral <- ggplot(nonzero_data, aes(x = `Spectral Count`, y = final_intensity)) +
    geom_point(alpha = 0.5, color = "darkgreen") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    scale_y_log10(labels = scientific) +
    labs(title = "Final Intensity vs Spectral Count",
         subtitle = paste("Correlation =", round(correlation, 3)),
         x = "Spectral Count",
         y = "Final Intensity (log10 scale)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # 8. Sample comparison plots (if multiple samples)
  if (length(unique(data$SampleID)) > 1) {
    # Sample plot data with sampling per group
    sample_plot_data <- data %>%
      filter(!is.na(final_intensity) & final_intensity > 0)
    
    # If we have too many points, sample them
    if (nrow(sample_plot_data) > 10000) {
      sample_plot_data <- sample_plot_data %>%
        slice_sample(n = 10000)
    }
    
    plots$sample_comparison <- ggplot(sample_plot_data, aes(x = reorder(SampleID, final_intensity, median), y = final_intensity)) +
      geom_boxplot(aes(fill = SampleID), alpha = 0.7, outlier.size = 0.5) +
      scale_y_log10(labels = scientific) +
      labs(title = "Intensity Distribution by Sample (Non-zero values)",
           subtitle = paste("n =", nrow(sample_plot_data), "peptides"),
           x = "Sample ID (ordered by median intensity)",
           y = "Final Intensity (log10 scale)") +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            legend.position = "none") +
      coord_flip()
    
    # Violin plot
    plots$violin_comparison <- ggplot(sample_plot_data, aes(x = final_intensity, y = reorder(SampleID, final_intensity, median))) +
      geom_violin(aes(fill = SampleID), alpha = 0.7, trim = FALSE) +
      geom_boxplot(height = 0.1, alpha = 0.8, outlier.size = 0.3, fill = "white") +
      scale_x_log10(labels = scientific) +
      labs(title = "Intensity Distribution by Sample - Violin Plot",
           subtitle = paste("n =", nrow(sample_plot_data), "peptides"),
           y = "Sample ID (ordered by median intensity)",
           x = "Final Intensity (log10 scale)") +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 10),
            legend.position = "none")
  }
  
  # Save individual plots
  cat("Saving intensity analysis plots...\n")
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_01_intensity_status_pie_chart.png")), 
         plots$pie_intensity, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_02_raw_intensity_histogram.png")), 
         plots$hist_raw, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_03_log_intensity_histogram.png")), 
         plots$hist_log, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_04_intensity_boxplot.png")), 
         plots$boxplot_intensity, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_05_log_intensity_density.png")), 
         plots$density_log, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_06_spectral_count_histogram.png")), 
         plots$spectral_hist, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_07_intensity_vs_spectral_scatter.png")), 
         plots$scatter_intensity_spectral, width = 10, height = 8, dpi = 300)
  
  if ("sample_comparison" %in% names(plots)) {
    ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_08_sample_comparison_boxplot.png")), 
           plots$sample_comparison, width = 12, height = 8, dpi = 300)
    
    ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_09_sample_comparison_violin.png")), 
           plots$violin_comparison, width = 12, height = 8, dpi = 300)
  }
  
  # Create and save combined plots
  combined_plot1 <- grid.arrange(
    plots$pie_intensity, plots$hist_raw,
    plots$hist_log, plots$density_log,
    ncol = 2, nrow = 2
  )
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_10_combined_intensity_analysis.png")), 
         combined_plot1, width = 16, height = 12, dpi = 300)
  
  # Spectral count analysis plots
  if ("sample_comparison" %in% names(plots)) {
    spectral_plots <- grid.arrange(plots$spectral_hist, plots$scatter_intensity_spectral, plots$sample_comparison, ncol = 1)
  } else {
    spectral_plots <- grid.arrange(plots$spectral_hist, plots$scatter_intensity_spectral, ncol = 1)
  }
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_11_spectral_count_analysis.png")), 
         spectral_plots, width = 12, height = 12, dpi = 300)
  
  cat("✓ Intensity analysis plots saved!\n")
  
  return(plots)
}

# Function to create enhanced output directory structure
create_output_directory_enhanced <- function(base_dir, dataset_name, analysis_type) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- file.path(base_dir, "results", paste0(timestamp, "_", dataset_name, "_", analysis_type))
  
  # Create main directories
  dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  
  return(output_dir)
}