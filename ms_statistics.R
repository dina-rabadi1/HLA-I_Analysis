# Immunopeptidomics Data Analysis Script
# Analysis of final_intensity values, spectral counts, and peptide characteristics
# Author: Generated for immunopeptidomics analysis
# Date: May 2025

# Load required libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(scales)
library(knitr)
library(kableExtra)
library(scales)

# ===============================
# CONFIGURATION PARAMETERS
# ===============================

setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis")

# File paths
data_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis"
file_name <- "unique_peptides_all.tsv"
full_path <- file.path(data_path, file_name)

# Sample filtering options
exclude_51S <- FALSE  # Set to FALSE to include sample 51S
include_samples <- NULL  # Set to c("117", "123", "456") to analyze only specific samples, or NULL for all

# Generate timestamp and folder name
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Determine folder name based on filtering
if (!is.null(include_samples)) {
  folder_suffix <- paste0("samples_", paste(include_samples, collapse="_"))
} else if (exclude_51S) {
  folder_suffix <- "all_exclude51S"
} else {
  folder_suffix <- "all"
}

# Create output directory with timestamp and sample info
output_dir <- file.path(data_path, "intensity_analysis_results", paste0(timestamp, "_", folder_suffix))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Output directory:", output_dir, "\n")

# ===============================
# SANITY CHECKS AND DATA LOADING
# ===============================

cat("Performing sanity checks...\n")

# Check if file exists
if (!file.exists(full_path)) {
  stop("ERROR: Data file not found at: ", full_path)
}

# Load data
cat("Loading data...\n")
tryCatch({
  data <- read_tsv(full_path, show_col_types = FALSE)
}, error = function(e) {
  stop("ERROR: Failed to load data file. ", e$message)
})

# Check required columns
required_cols <- c("final_intensity", "Spectral Count", "SampleID", "Peptide")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("ERROR: Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("✓ Data loaded successfully\n")
cat("Dataset dimensions:", nrow(data), "rows x", ncol(data), "columns\n")

# ===============================
# SAMPLE FILTERING
# ===============================

original_count <- nrow(data)
cat("\nApplying sample filters...\n")
cat("Original dataset:", original_count, "peptides\n")

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
if (exclude_51S) {
  if ("51S" %in% data$SampleID) {
    data <- data %>% filter(SampleID != "51S")
    cat("After excluding 51S:", nrow(data), "peptides\n")
  } else {
    cat("Note: Sample 51S not found in dataset\n")
  }
}

# Final check
if (nrow(data) == 0) {
  stop("ERROR: No data remaining after filtering!")
}

cat("✓ Sample filtering completed\n")
cat("Final dataset:", nrow(data), "peptides from", length(unique(data$SampleID)), "samples\n")

# ===============================
# SAMPLE BREAKDOWN ANALYSIS
# ===============================

cat("\nGenerating sample breakdown...\n")

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

print(sample_summary)

# Save sample breakdown
write_csv(sample_summary, file.path(output_dir, "sample_breakdown_summary.csv"))

# ===============================
# MISSING/ZERO VALUE ANALYSIS
# ===============================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("MISSING AND ZERO VALUE ANALYSIS\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# Analyze final_intensity column
total_peptides <- nrow(data)
missing_intensity <- sum(is.na(data$final_intensity))
zero_intensity <- sum(data$final_intensity == 0, na.rm = TRUE)
non_zero_intensity <- sum(data$final_intensity > 0, na.rm = TRUE)

cat("Total peptides:", total_peptides, "\n")
cat("Missing final_intensity values:", missing_intensity, "(", round(missing_intensity/total_peptides*100, 2), "%)\n")
cat("Zero final_intensity values:", zero_intensity, "(", round(zero_intensity/total_peptides*100, 2), "%)\n")
cat("Non-zero final_intensity values:", non_zero_intensity, "(", round(non_zero_intensity/total_peptides*100, 2), "%)\n")

# Sanity check: reasonable intensity ranges
if (non_zero_intensity > 0) {
  intensity_range <- range(data$final_intensity, na.rm = TRUE)
  cat("Intensity range:", intensity_range[1], "to", intensity_range[2], "\n")
  if (intensity_range[2] > 1e10) {
    warning("WARNING: Very high intensity values detected (>10^10)")
  }
}

# Create summary table for missing/zero analysis
missing_summary <- data.frame(
  Category = c("Total Peptides", "Missing Intensity", "Zero Intensity", "Non-zero Intensity"),
  Count = c(total_peptides, missing_intensity, zero_intensity, non_zero_intensity),
  Percentage = c(100, 
                 round(missing_intensity/total_peptides*100, 2),
                 round(zero_intensity/total_peptides*100, 2),
                 round(non_zero_intensity/total_peptides*100, 2))
)

print(missing_summary)

# Save missing value summary
write_csv(missing_summary, file.path(output_dir, "missing_value_summary.csv"))

# ===============================
# DATASET PREPARATION
# ===============================

# Create datasets for analysis
all_data <- data  # All peptides
nonzero_data <- data %>% filter(final_intensity > 0 & !is.na(final_intensity))  # Only non-zero intensity

cat("\nDatasets prepared:\n")
cat("All data:", nrow(all_data), "peptides\n")
cat("Non-zero intensity data:", nrow(nonzero_data), "peptides\n")

# Sanity check
if (nrow(nonzero_data) == 0) {
  stop("ERROR: No peptides with non-zero intensity found!")
}

# ===============================
# STATISTICAL ANALYSIS FUNCTION
# ===============================

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

# ===============================
# RAW INTENSITY ANALYSIS
# ===============================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("RAW INTENSITY ANALYSIS\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# Calculate statistics for raw intensity
stats_all_raw <- calculate_statistics(all_data$final_intensity, "All_Peptides_Raw")
stats_nonzero_raw <- calculate_statistics(nonzero_data$final_intensity, "NonZero_Peptides_Raw")

# Combine statistics
raw_stats_combined <- rbind(stats_all_raw, stats_nonzero_raw)
print(raw_stats_combined)

# Calculate percentiles for non-zero data
percentiles <- quantile(nonzero_data$final_intensity, 
                        probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99), 
                        na.rm = TRUE)
cat("\nPercentiles (Non-zero intensity):\n")
print(percentiles)

# ===============================
# LOG-TRANSFORMED INTENSITY ANALYSIS
# ===============================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("LOG-TRANSFORMED INTENSITY ANALYSIS\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# Log transform non-zero values (log10)
nonzero_data$log_final_intensity <- log10(nonzero_data$final_intensity)

stats_nonzero_log <- calculate_statistics(nonzero_data$log_final_intensity, "NonZero_Peptides_Log10")
print(stats_nonzero_log)

# Calculate percentiles for log-transformed data
log_percentiles <- quantile(nonzero_data$log_final_intensity, 
                            probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99), 
                            na.rm = TRUE)
cat("\nPercentiles (Log10-transformed intensity):\n")
print(log_percentiles)

# ===============================
# THRESHOLD DEFINITIONS
# ===============================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("INTENSITY THRESHOLD DEFINITIONS\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# Define thresholds based on non-zero data
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

print(thresholds)

# ===============================
# ENHANCED SPECTRAL COUNT ANALYSIS
# ===============================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("ENHANCED SPECTRAL COUNT ANALYSIS\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# Basic spectral count statistics
spectral_stats_all <- calculate_statistics(all_data$`Spectral Count`, "All_Peptides_SpectralCount")
spectral_stats_nonzero <- calculate_statistics(nonzero_data$`Spectral Count`, "NonZero_Peptides_SpectralCount")

spectral_stats_combined <- rbind(spectral_stats_all, spectral_stats_nonzero)
print(spectral_stats_combined)

# Detailed spectral count distribution
spectral_table <- table(nonzero_data$`Spectral Count`)
cat("\nSpectral Count Distribution (Non-zero intensity peptides):\n")
spectral_df <- data.frame(
  Spectral_Count = as.numeric(names(spectral_table)),
  Count = as.numeric(spectral_table),
  Percentage = round(as.numeric(spectral_table) / sum(spectral_table) * 100, 2)
)
print(head(spectral_df, 15))

# Spectral count vs intensity correlation
correlation <- cor(nonzero_data$final_intensity, nonzero_data$`Spectral Count`, use = "complete.obs")
cat("\nCorrelation between intensity and spectral count:", round(correlation, 3), "\n")

# Cross-tabulation: Intensity quartiles vs Spectral count ranges
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
cat("\nCross-tabulation: Intensity vs Spectral Count Categories\n")
print(crosstab)
print(round(prop.table(crosstab, 1) * 100, 1))  # Row percentages

# ===============================
# SAVE STATISTICAL SUMMARIES
# ===============================

# Combine all statistics
all_stats <- rbind(raw_stats_combined, stats_nonzero_log, spectral_stats_combined)

# Save comprehensive statistics
write_csv(all_stats, file.path(output_dir, "comprehensive_statistics.csv"))
write_csv(data.frame(Percentile = names(percentiles), Raw_Value = as.numeric(percentiles)), 
          file.path(output_dir, "raw_intensity_percentiles.csv"))
write_csv(data.frame(Percentile = names(log_percentiles), Log10_Value = as.numeric(log_percentiles)), 
          file.path(output_dir, "log_intensity_percentiles.csv"))
write_csv(thresholds, file.path(output_dir, "intensity_thresholds.csv"))
write_csv(spectral_df, file.path(output_dir, "spectral_count_distribution.csv"))
write_csv(as.data.frame.matrix(crosstab), file.path(output_dir, "intensity_spectral_crosstab.csv"))

# ===============================
# VISUALIZATION
# ===============================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("GENERATING VISUALIZATIONS\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# Create plots
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

# 8. Sample comparison plot
if (length(unique(data$SampleID)) > 1) {
  # Create sample plot data with sampling per group
  sample_plot_data <- data %>%
    filter(!is.na(final_intensity) & final_intensity > 0)
  
  # If we have too many points, sample them
  if (nrow(sample_plot_data) > 10000) {
    sample_plot_data <- sample_plot_data %>%
      slice_sample(n = 10000)
  }
  
  plots$sample_comparison <- ggplot(sample_plot_data, aes(x = SampleID, y = final_intensity)) +
    geom_boxplot(aes(fill = SampleID), alpha = 0.7) +
    scale_y_log10(labels = scientific) +
    labs(title = "Intensity Distribution by Sample",
         x = "Sample ID",
         y = "Final Intensity (log10 scale)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
}

# Save individual plots
for (plot_name in names(plots)) {
  ggsave(file.path(output_dir, paste0(plot_name, ".png")), 
         plots[[plot_name]], width = 10, height = 8, dpi = 300)
}

# Create and save combined plots
combined_plot1 <- grid.arrange(
  plots$pie_intensity, plots$hist_raw,
  plots$hist_log, plots$density_log,
  ncol = 2, nrow = 2
)

ggsave(file.path(output_dir, "combined_intensity_analysis.png"), 
       combined_plot1, width = 16, height = 12, dpi = 300)

# Save spectral count analysis plots
spectral_plots <- if("sample_comparison" %in% names(plots)) {
  grid.arrange(plots$spectral_hist, plots$scatter_intensity_spectral, plots$sample_comparison, ncol = 1)
} else {
  grid.arrange(plots$spectral_hist, plots$scatter_intensity_spectral, ncol = 1)
}

ggsave(file.path(output_dir, "spectral_count_analysis.png"), 
       spectral_plots, width = 12, height = 12, dpi = 300)

# ===============================
# IMPROVED SAMPLE COMPARISON PLOTS
# ===============================

cat("Generating improved sample comparison plots...\n")

# 1. Sample comparison with non-zero values only (better formatted)
sample_plot_data_nonzero <- data %>%
  filter(!is.na(final_intensity) & final_intensity > 0)

plot_sample_nonzero <- ggplot(sample_plot_data_nonzero, aes(x = reorder(SampleID, final_intensity, median), y = final_intensity)) +
  geom_boxplot(aes(fill = SampleID), alpha = 0.7, outlier.size = 0.5) +
  scale_y_log10(labels = scientific) +
  labs(title = "Intensity Distribution by Sample (Non-zero values only)",
       subtitle = paste("n =", nrow(sample_plot_data_nonzero), "peptides"),
       x = "Sample ID (ordered by median intensity)",
       y = "Final Intensity (log10 scale)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none") +
  coord_flip()  # Flip to horizontal for better label visibility

# 2. Sample comparison including zero values
sample_plot_data_all <- data %>%
  mutate(intensity_display = ifelse(final_intensity == 0 | is.na(final_intensity), 
                                    1, final_intensity))  # Set zeros to 1 for log scale display

plot_sample_all <- ggplot(sample_plot_data_all, aes(x = reorder(SampleID, intensity_display, median), y = intensity_display)) +
  geom_boxplot(aes(fill = SampleID), alpha = 0.7, outlier.size = 0.5) +
  scale_y_log10(labels = scientific) +
  labs(title = "Intensity Distribution by Sample (All values)",
       subtitle = paste("n =", nrow(sample_plot_data_all), "peptides (zeros shown as 1 for visualization)"),
       x = "Sample ID (ordered by median intensity)",
       y = "Final Intensity (log10 scale, zeros = 1)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none") +
  coord_flip()

# Save the improved plots
ggsave(file.path(output_dir, "sample_comparison_nonzero_improved.png"), 
       plot_sample_nonzero, width = 12, height = 8, dpi = 300)

ggsave(file.path(output_dir, "sample_comparison_all_values.png"), 
       plot_sample_all, width = 12, height = 8, dpi = 300)

cat("✓ Improved sample comparison plots saved!\n")

# ===============================
# VIOLIN PLOTS FOR SAMPLE COMPARISON
# ===============================

cat("Generating violin plots for sample comparison...\n")

# 1. Violin plot with non-zero values only
violin_nonzero <- ggplot(sample_plot_data_nonzero, aes(x = reorder(SampleID, final_intensity, median), y = final_intensity)) +
  geom_violin(aes(fill = SampleID), alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.size = 0.3, fill = "white") +
  scale_y_log10(labels = scientific) +
  labs(title = "Intensity Distribution by Sample - Violin Plot (Non-zero values)",
       subtitle = paste("n =", nrow(sample_plot_data_nonzero), "peptides - Shows distribution density + boxplot"),
       x = "Sample ID (ordered by median intensity)",
       y = "Final Intensity (log10 scale)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none") +
  coord_flip()

# 2. Horizontal violin plot (better for many samples)
violin_horizontal <- ggplot(sample_plot_data_nonzero, aes(x = final_intensity, y = reorder(SampleID, final_intensity, median))) +
  geom_violin(aes(fill = SampleID), alpha = 0.7, trim = FALSE) +
  geom_boxplot(height = 0.1, alpha = 0.8, outlier.size = 0.3, fill = "white") +
  scale_x_log10(labels = scientific) +
  labs(title = "Intensity Distribution by Sample - Horizontal Violin Plot",
       subtitle = paste("n =", nrow(sample_plot_data_nonzero), "peptides - Horizontal layout for better readability"),
       y = "Sample ID (ordered by median intensity)",
       x = "Final Intensity (log10 scale)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        legend.position = "none")

# 3. Violin plot with density comparison (all samples overlaid)
violin_overlay <- ggplot(sample_plot_data_nonzero, aes(x = final_intensity, fill = SampleID)) +
  geom_density(alpha = 0.4) +
  scale_x_log10(labels = scientific) +
  labs(title = "Intensity Distribution Density Comparison by Sample",
       subtitle = paste("n =", nrow(sample_plot_data_nonzero), "peptides - Overlaid density curves"),
       x = "Final Intensity (log10 scale)",
       y = "Density",
       fill = "Sample ID") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "right")

# 4. Faceted violin plots (each sample separate)
violin_faceted <- ggplot(sample_plot_data_nonzero, aes(x = "", y = final_intensity)) +
  geom_violin(aes(fill = SampleID), alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.size = 0.2, fill = "white") +
  scale_y_log10(labels = scientific) +
  facet_wrap(~ SampleID, scales = "free_y", ncol = 4) +
  labs(title = "Intensity Distribution by Sample - Faceted Violin Plots",
       subtitle = paste("n =", nrow(sample_plot_data_nonzero), "peptides - Individual plots per sample"),
       x = "Sample",
       y = "Final Intensity (log10 scale)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(face = "bold"))

# Save violin plots
ggsave(file.path(output_dir, "violin_plot_nonzero_vertical.png"), 
       violin_nonzero, width = 12, height = 8, dpi = 300)

ggsave(file.path(output_dir, "violin_plot_nonzero_horizontal.png"), 
       violin_horizontal, width = 12, height = 8, dpi = 300)

ggsave(file.path(output_dir, "violin_plot_density_overlay.png"), 
       violin_overlay, width = 12, height = 8, dpi = 300)

ggsave(file.path(output_dir, "violin_plot_faceted.png"), 
       violin_faceted, width = 16, height = 10, dpi = 300)

# Create a combined violin plot comparison
combined_violin <- grid.arrange(
  violin_horizontal, 
  violin_overlay,
  ncol = 1, nrow = 2
)

ggsave(file.path(output_dir, "violin_plots_combined.png"), 
       combined_violin, width = 14, height = 12, dpi = 300)

cat("✓ Violin plots for sample comparison saved!\n")
cat("  - Vertical violin plot with boxplots\n")
cat("  - Horizontal violin plot (better readability)\n") 
cat("  - Density overlay comparison\n")
cat("  - Faceted individual violin plots\n")
cat("  - Combined violin plot summary\n")

# ===============================
# COMPREHENSIVE PDF REPORT
# ===============================

cat("Generating comprehensive PDF report...\n")

# Create a simple markdown report and convert to PDF
report_content <- paste0("
# Immunopeptidomics Data Analysis Report

**Generated:** ", Sys.time(), "  
**Dataset:** ", file_name, "  
**Analysis:** ", folder_suffix, "  

## Executive Summary

This analysis examined ", total_peptides, " peptides from ", length(unique(data$SampleID)), " samples. 
Key findings:

- **", round(non_zero_intensity/total_peptides*100, 1), "%** of peptides have detectable intensity values
- **Median intensity:** ", round(median(nonzero_data$final_intensity, na.rm = TRUE), 0), "
- **Intensity range:** ", round(min(nonzero_data$final_intensity, na.rm = TRUE), 0), " to ", 
                         round(max(nonzero_data$final_intensity, na.rm = TRUE), 0), "
- **Intensity-Spectral Count Correlation:** ", round(correlation, 3), "

## Recommended Intensity Thresholds

Based on quartile analysis:
- **Low intensity:** ≤ ", round(quantile(nonzero_data$final_intensity, 0.25, na.rm = TRUE), 0), " (bottom 25%)
- **Medium intensity:** ", round(quantile(nonzero_data$final_intensity, 0.25, na.rm = TRUE), 0), " to ", 
                         round(quantile(nonzero_data$final_intensity, 0.75, na.rm = TRUE), 0), " (middle 50%)
- **High intensity:** ≥ ", round(quantile(nonzero_data$final_intensity, 0.75, na.rm = TRUE), 0), " (top 25%)

## Files Generated

1. **CSV Files:** Comprehensive statistics, percentiles, thresholds, sample breakdown
2. **PNG Plots:** Individual and combined visualizations
3. **This PDF Report:** Complete analysis summary

For detailed analysis, refer to the CSV files and individual plots in the output directory.
")

# Write the report
writeLines(report_content, file.path(output_dir, "analysis_report.md"))

# Try to create PDF (requires pandoc)
tryCatch({
  rmarkdown::render(file.path(output_dir, "analysis_report.md"), 
                    output_format = "pdf_document",
                    output_file = "comprehensive_analysis_report.pdf",
                    quiet = TRUE)
}, error = function(e) {
  cat("Note: Could not generate PDF report. Markdown report available.\n")
  cat("Error:", e$message, "\n")
})

# Generate HTML report (works without LaTeX)
tryCatch({
  rmarkdown::render(file.path(output_dir, "analysis_report.md"), 
                    output_format = rmarkdown::html_document(
                      theme = "flatly",
                      highlight = "tango",
                      toc = TRUE,
                      toc_float = TRUE
                    ),
                    output_file = "comprehensive_analysis_report.html",
                    quiet = TRUE)
  cat("✓ HTML report generated successfully!\n")
}, error = function(e) {
  cat("Note: Could not generate HTML report.\n")
  cat("Error:", e$message, "\n")
})

# Also generate Word document (alternative to PDF)
tryCatch({
  rmarkdown::render(file.path(output_dir, "analysis_report.md"), 
                    output_format = "word_document",
                    output_file = "comprehensive_analysis_report.docx",
                    quiet = TRUE)
  cat("✓ Word document generated successfully!\n")
}, error = function(e) {
  cat("Note: Could not generate Word document.\n")
  cat("Error:", e$message, "\n")
})

# ===============================
# FINAL SUMMARY
# ===============================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("ANALYSIS COMPLETE - SUMMARY AND RECOMMENDATIONS\n")
cat(paste(rep("=", 60), collapse=""), "\n")

cat("Files saved to:", output_dir, "\n")
cat("✓ Sample breakdown CSV\n")
cat("✓ Comprehensive statistics CSV\n")
cat("✓ Percentiles and thresholds CSV\n")
cat("✓ Spectral count analysis CSV\n")
cat("✓ Individual and combined visualization plots\n")
cat("✓ Analysis report (Markdown and PDF if available)\n")

cat("\nKey Findings:\n")
cat("1. Total peptides analyzed:", total_peptides, "\n")
cat("2. Peptides with non-zero intensity:", non_zero_intensity, "(", round(non_zero_intensity/total_peptides*100, 1), "%)\n")
cat("3. Median intensity (non-zero):", round(median(nonzero_data$final_intensity, na.rm = TRUE), 0), "\n)")
cat("4. Mean intensity (non-zero):", round(mean(nonzero_data$final_intensity, na.rm = TRUE), 0), "\n")
cat("5. Intensity range:", round(min(nonzero_data$final_intensity, na.rm = TRUE), 0), "to", 
    round(max(nonzero_data$final_intensity, na.rm = TRUE), 0), "\n")
cat("6. Intensity-Spectral Count correlation:", round(correlation, 3), "\n")
cat("7. Samples analyzed:", paste(sort(unique(data$SampleID)), collapse = ", "), "\n")

# Suggested thresholds
q1_threshold <- quantile(nonzero_data$final_intensity, 0.25, na.rm = TRUE)
q3_threshold <- quantile(nonzero_data$final_intensity, 0.75, na.rm = TRUE)

cat("\nSuggested Intensity Classifications:\n")
cat("- Low intensity: ≤", round(q1_threshold, 0), "(bottom 25%)\n")
cat("- Medium intensity:", round(q1_threshold, 0), "to", round(q3_threshold, 0), "(middle 50%)\n")
cat("- High intensity: ≥", round(q3_threshold, 0), "(top 25%)\n")

cat("\n✓ Analysis completed successfully!\n")
