# CHOP FL57 Liver vs Lung Peptide Intensity Correlation Analysis
# This script analyzes peptide intensity correlation between liver and lung samples from patient FL57

# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)

# Set working directory (adjust as needed)
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation")

cat("=====================================================================\n")
cat("CHOP FL57 LIVER VS LUNG PEPTIDE CORRELATION ANALYSIS\n")
cat("=====================================================================\n\n")

# Define file paths
input_file <- "final_CHOP_combined_peptide.tsv"
filtered_output_file <- "FL57_liver_lung_comparison.tsv"
plot_output_file <- "FL57_liver_lung_correlation_plot.png"
stats_output_file <- "FL57_liver_lung_summary_stats.txt"

# Read the integrated CHOP data
cat("Reading integrated CHOP data...\n")
chop_data <- read.delim(input_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Check if the required columns exist
liver_col <- "H5 FL57Liver Intensity"
lung_col <- "H6 FL57Lung Intensity"

if (!liver_col %in% colnames(chop_data)) {
  # Try alternative column naming
  liver_col <- grep("FL57Liver.*Intensity", colnames(chop_data), ignore.case=TRUE, value=TRUE)[1]
}

if (!lung_col %in% colnames(chop_data)) {
  # Try alternative column naming  
  lung_col <- grep("FL57Lung.*Intensity", colnames(chop_data), ignore.case=TRUE, value=TRUE)[1]
}

cat("Using columns:\n")
cat("Liver:", liver_col, "\n")
cat("Lung:", lung_col, "\n")

# Extract relevant data for FL57 liver and lung comparison
cat("Extracting FL57 liver and lung data...\n")

# Get peptide sequence column
peptide_col <- "Peptide Sequence"
if (!peptide_col %in% colnames(chop_data)) {
  peptide_col <- grep("Peptide", colnames(chop_data), ignore.case=TRUE, value=TRUE)[1]
}

# Create filtered dataset
fl57_data <- chop_data %>%
  select(all_of(c(peptide_col, liver_col, lung_col, "Gene", "Protein"))) %>%
  rename(
    Peptide = all_of(peptide_col),
    FL57_Liver_Intensity = all_of(liver_col),
    FL57_Lung_Intensity = all_of(lung_col)
  ) %>%
  mutate(
    FL57_Liver_Intensity = as.numeric(FL57_Liver_Intensity),
    FL57_Lung_Intensity = as.numeric(FL57_Lung_Intensity)
  ) %>%
  # Replace NA with 0
  mutate(
    FL57_Liver_Intensity = ifelse(is.na(FL57_Liver_Intensity), 0, FL57_Liver_Intensity),
    FL57_Lung_Intensity = ifelse(is.na(FL57_Lung_Intensity), 0, FL57_Lung_Intensity)
  )

# Save filtered dataset
cat("Saving filtered dataset...\n")
write.table(fl57_data, file=filtered_output_file, sep="\t", quote=FALSE, row.names=FALSE)

# Calculate summary statistics
cat("Calculating summary statistics...\n")

# Peptides detected in at least one sample (intensity > 0)
peptides_at_least_one <- fl57_data %>%
  filter(FL57_Liver_Intensity > 0 | FL57_Lung_Intensity > 0)

# Peptides detected in both samples (intensity > 0 in both)
peptides_both <- fl57_data %>%
  filter(FL57_Liver_Intensity > 0 & FL57_Lung_Intensity > 0)

# Peptides detected only in liver
peptides_liver_only <- fl57_data %>%
  filter(FL57_Liver_Intensity > 0 & FL57_Lung_Intensity == 0)

# Peptides detected only in lung
peptides_lung_only <- fl57_data %>%
  filter(FL57_Liver_Intensity == 0 & FL57_Lung_Intensity > 0)

# Calculate correlations for peptides detected in both samples
if (nrow(peptides_both) > 2) {
  # Add small constant to avoid log(0)
  liver_intensities <- peptides_both$FL57_Liver_Intensity + 1
  lung_intensities <- peptides_both$FL57_Lung_Intensity + 1
  
  # Pearson correlation
  pearson_cor <- cor(liver_intensities, lung_intensities, method="pearson")
  pearson_r2 <- pearson_cor^2
  
  # Spearman correlation
  spearman_cor <- cor(liver_intensities, lung_intensities, method="spearman")
  spearman_r2 <- spearman_cor^2
  
  # Linear model for R²
  lm_model <- lm(log10(lung_intensities) ~ log10(liver_intensities))
  lm_r2 <- summary(lm_model)$r.squared
  
} else {
  pearson_cor <- pearson_r2 <- spearman_cor <- spearman_r2 <- lm_r2 <- NA
}

# Calculate descriptive statistics
liver_stats <- list(
  total_detected = sum(fl57_data$FL57_Liver_Intensity > 0),
  mean_intensity = mean(fl57_data$FL57_Liver_Intensity[fl57_data$FL57_Liver_Intensity > 0], na.rm=TRUE),
  median_intensity = median(fl57_data$FL57_Liver_Intensity[fl57_data$FL57_Liver_Intensity > 0], na.rm=TRUE),
  max_intensity = max(fl57_data$FL57_Liver_Intensity, na.rm=TRUE)
)

lung_stats <- list(
  total_detected = sum(fl57_data$FL57_Lung_Intensity > 0),
  mean_intensity = mean(fl57_data$FL57_Lung_Intensity[fl57_data$FL57_Lung_Intensity > 0], na.rm=TRUE),
  median_intensity = median(fl57_data$FL57_Lung_Intensity[fl57_data$FL57_Lung_Intensity > 0], na.rm=TRUE),
  max_intensity = max(fl57_data$FL57_Lung_Intensity, na.rm=TRUE)
)

# Create summary statistics output
summary_text <- paste0(
  "FL57 LIVER VS LUNG PEPTIDE CORRELATION ANALYSIS SUMMARY\n",
  "=====================================================\n\n",
  "DATASET OVERVIEW:\n",
  "Total peptides in dataset: ", nrow(fl57_data), "\n",
  "Peptides detected in liver only: ", nrow(peptides_liver_only), "\n",
  "Peptides detected in lung only: ", nrow(peptides_lung_only), "\n",
  "Peptides detected in both samples: ", nrow(peptides_both), "\n",
  "Peptides detected in at least one sample: ", nrow(peptides_at_least_one), "\n\n",
  
  "LIVER SAMPLE STATISTICS:\n",
  "Total peptides detected: ", liver_stats$total_detected, "\n",
  "Mean intensity (detected peptides): ", sprintf("%.2e", liver_stats$mean_intensity), "\n",
  "Median intensity (detected peptides): ", sprintf("%.2e", liver_stats$median_intensity), "\n",
  "Maximum intensity: ", sprintf("%.2e", liver_stats$max_intensity), "\n\n",
  
  "LUNG SAMPLE STATISTICS:\n",
  "Total peptides detected: ", lung_stats$total_detected, "\n",
  "Mean intensity (detected peptides): ", sprintf("%.2e", lung_stats$mean_intensity), "\n",
  "Median intensity (detected peptides): ", sprintf("%.2e", lung_stats$median_intensity), "\n",
  "Maximum intensity: ", sprintf("%.2e", lung_stats$max_intensity), "\n\n",
  
  "CORRELATION ANALYSIS (peptides detected in both samples):\n",
  "Number of peptides: ", nrow(peptides_both), "\n",
  "Pearson correlation: ", sprintf("%.3f", pearson_cor), "\n",
  "Pearson R²: ", sprintf("%.3f", pearson_r2), "\n",
  "Spearman correlation: ", sprintf("%.3f", spearman_cor), "\n",
  "Spearman R²: ", sprintf("%.3f", spearman_r2), "\n",
  "Linear model R² (log-transformed): ", sprintf("%.3f", lm_r2), "\n\n",
  
  "OVERLAP STATISTICS:\n",
  "Percentage of liver peptides also found in lung: ", 
  sprintf("%.1f%%", 100 * nrow(peptides_both) / liver_stats$total_detected), "\n",
  "Percentage of lung peptides also found in liver: ", 
  sprintf("%.1f%%", 100 * nrow(peptides_both) / lung_stats$total_detected), "\n",
  "Overall overlap (Jaccard index): ", 
  sprintf("%.3f", nrow(peptides_both) / nrow(peptides_at_least_one)), "\n"
)

# Save summary statistics
cat("Saving summary statistics...\n")
writeLines(summary_text, stats_output_file)

# Create correlation plot matching the MSKCC vs CHOP style
cat("Creating correlation plot...\n")

# Prepare data for plotting - only peptides detected in both samples
plot_data <- peptides_both %>%
  mutate(
    Liver_Log = log10(FL57_Liver_Intensity),
    Lung_Log = log10(FL57_Lung_Intensity)
  ) %>%
  filter(is.finite(Liver_Log) & is.finite(Lung_Log) & 
           FL57_Liver_Intensity > 0 & FL57_Lung_Intensity > 0)

if (nrow(plot_data) > 0) {
  
  # Print some diagnostics
  cat("Plot data range - Liver:", range(plot_data$Liver_Log), "\n")
  cat("Plot data range - Lung:", range(plot_data$Lung_Log), "\n")
  cat("Number of points:", nrow(plot_data), "\n")
  
  # Determine appropriate axis limits based on data
  x_min <- max(0, min(plot_data$Liver_Log) - 0.5)
  x_max <- max(plot_data$Liver_Log) + 0.5
  y_min <- max(0, min(plot_data$Lung_Log) - 0.5)
  y_max <- max(plot_data$Lung_Log) + 0.5
  
  # Create the plot matching your MSKCC vs CHOP style exactly
  p1 <- ggplot(plot_data, aes(x = Liver_Log, y = Lung_Log)) +
    geom_point(color = "#2E8B57", alpha = 0.2, size = 0.5) +  # Much smaller points, very transparent
    geom_smooth(method = "lm", color = "red", se = FALSE, linewidth = 1.5) +  # Red trend line
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    labs(
      title = "FL57 Liver vs Lung Intensity Correlation",
      subtitle = paste0("R² = ", sprintf("%.3f", lm_r2), ", n = ", nrow(plot_data)),
      x = "FL57 Liver Intensity (log10)",
      y = "FL57 Lung Intensity (log10)"
    ) +
    # Set axis limits based on actual data range
    scale_x_continuous(limits = c(x_min, x_max), 
                       breaks = scales::pretty_breaks(n = 6)) +
    scale_y_continuous(limits = c(y_min, y_max),
                       breaks = scales::pretty_breaks(n = 6))
  
  # Save the correlation plot
  ggsave(plot_output_file, plot = p1, width = 8, height = 6, dpi = 300, bg = "white")
  
  # Also create a subsampled version for clearer visualization
  subsample_plot_file <- "FL57_liver_lung_correlation_subsampled.png"
  
  # Subsample data for clearer visualization (every nth point)
  n_points_desired <- 1000  # Reduce to 1000 points for clarity
  if (nrow(plot_data) > n_points_desired) {
    subsample_indices <- seq(1, nrow(plot_data), length.out = n_points_desired)
    plot_data_sub <- plot_data[round(subsample_indices), ]
  } else {
    plot_data_sub <- plot_data
  }
  
  p1_sub <- ggplot(plot_data_sub, aes(x = Liver_Log, y = Lung_Log)) +
    geom_point(color = "#2E8B57", alpha = 0.6, size = 1.2) +  # Larger, less transparent points
    geom_smooth(method = "lm", color = "red", se = FALSE, linewidth = 1.5) +  # Red trend line
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    labs(
      title = "FL57 Liver vs Lung Intensity Correlation (Subsampled)",
      subtitle = paste0("R² = ", sprintf("%.3f", lm_r2), ", n = ", nrow(plot_data), " (showing ", nrow(plot_data_sub), ")"),
      x = "FL57 Liver Intensity (log10)",
      y = "FL57 Lung Intensity (log10)"
    ) +
    scale_x_continuous(limits = c(x_min, x_max), 
                       breaks = scales::pretty_breaks(n = 6)) +
    scale_y_continuous(limits = c(y_min, y_max),
                       breaks = scales::pretty_breaks(n = 6))
  
  # Save the subsampled plot
  ggsave(subsample_plot_file, plot = p1_sub, width = 8, height = 6, dpi = 300, bg = "white")
  
} else {
  cat("Warning: No peptides detected in both samples for plotting.\n")
}

# Create unique peptides analysis and plot
cat("Creating unique peptides correlation analysis...\n")

# Create data for unique peptides plot
unique_plot_file <- "FL57_unique_peptides_plot.png"

# For unique peptides, we'll show the counts in a different way
# Since we can't correlate single numbers, let's create a visualization showing the distribution

# Create a summary plot for unique peptides
unique_summary <- data.frame(
  Sample = c("Liver Only", "Lung Only", "Both Samples"),
  Count = c(nrow(peptides_liver_only), nrow(peptides_lung_only), nrow(peptides_both))
)

# Create bar plot for unique peptides
p2 <- ggplot(unique_summary, aes(x = Sample, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("Liver Only" = "#FF6B6B", "Lung Only" = "#4ECDC4", "Both Samples" = "#2E8B57")) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
    legend.position = "none"
  ) +
  labs(
    title = "FL57 Peptide Detection Distribution",
    subtitle = "Tissue-specific vs Shared Peptides",
    x = "Detection Pattern",
    y = "Number of Peptides"
  ) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4)

# Save the unique peptides plot
ggsave(unique_plot_file, plot = p2, width = 8, height = 6, dpi = 300, bg = "white")

# Additionally, create a scatter plot showing intensity vs unique detection
# This will show liver-only vs lung-only peptides in a more meaningful way
unique_scatter_file <- "FL57_unique_peptides_scatter.png"

# Prepare data for unique peptides scatter
liver_only_data <- peptides_liver_only %>%
  mutate(
    Sample_Type = "Liver Only",
    Intensity = FL57_Liver_Intensity,
    Log_Intensity = log10(FL57_Liver_Intensity)
  )

lung_only_data <- peptides_lung_only %>%
  mutate(
    Sample_Type = "Lung Only", 
    Intensity = FL57_Lung_Intensity,
    Log_Intensity = log10(FL57_Lung_Intensity)
  )

unique_combined <- bind_rows(liver_only_data, lung_only_data)

if (nrow(unique_combined) > 0) {
  # Create intensity distribution plot for unique peptides
  p3 <- ggplot(unique_combined, aes(x = Log_Intensity, fill = Sample_Type)) +
    geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
    scale_fill_manual(values = c("Liver Only" = "#FF6B6B", "Lung Only" = "#4ECDC4")) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ) +
    labs(
      title = "FL57 Tissue-Specific Peptide Intensity Distribution",
      subtitle = "Comparing Liver-Only vs Lung-Only Peptides",
      x = "Peptide Intensity (log10)",
      y = "Number of Peptides",
      fill = "Sample Type"
    )
  
  # Save the unique peptides scatter plot
  ggsave(unique_scatter_file, plot = p3, width = 10, height = 6, dpi = 300, bg = "white")
}

# Print summary to console
cat("\n=====================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=====================================================================\n\n")

cat("FILES CREATED:\n")
cat("1.", filtered_output_file, "- Filtered dataset with FL57 liver and lung data\n")
cat("2.", plot_output_file, "- Intensity correlation plot (all points)\n")
cat("3.", "FL57_liver_lung_correlation_subsampled.png", "- Correlation plot (subsampled for clarity)\n")
cat("4.", "FL57_unique_peptides_plot.png", "- Unique peptides bar chart\n")
cat("5.", "FL57_unique_peptides_scatter.png", "- Tissue-specific peptide intensity distributions\n")
cat("6.", stats_output_file, "- Summary statistics\n\n")

cat("SUMMARY STATISTICS:\n")
cat(summary_text)

cat("\nAnalysis completed successfully!\n")

# CHOP FL57 Liver vs Lung Peptide Intensity Correlation Analysis
# This script analyzes peptide intensity correlation between liver and lung samples from patient FL57

# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)

# Set working directory (adjust as needed)
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation")

cat("=====================================================================\n")
cat("CHOP FL57 LIVER VS LUNG PEPTIDE CORRELATION ANALYSIS\n")
cat("=====================================================================\n\n")

# Define file paths
input_file <- "final_CHOP_combined_peptide.tsv"
filtered_output_file <- "FL57_liver_lung_comparison.tsv"
plot_output_file <- "FL57_liver_lung_correlation_plot.png"
stats_output_file <- "FL57_liver_lung_summary_stats.txt"

# Read the integrated CHOP data
cat("Reading integrated CHOP data...\n")
chop_data <- read.delim(input_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Check if the required columns exist
liver_col <- "H5 FL57Liver Intensity"
lung_col <- "H6 FL57Lung Intensity"

if (!liver_col %in% colnames(chop_data)) {
  # Try alternative column naming
  liver_col <- grep("FL57Liver.*Intensity", colnames(chop_data), ignore.case=TRUE, value=TRUE)[1]
}

if (!lung_col %in% colnames(chop_data)) {
  # Try alternative column naming  
  lung_col <- grep("FL57Lung.*Intensity", colnames(chop_data), ignore.case=TRUE, value=TRUE)[1]
}

cat("Using columns:\n")
cat("Liver:", liver_col, "\n")
cat("Lung:", lung_col, "\n")

# Extract relevant data for FL57 liver and lung comparison
cat("Extracting FL57 liver and lung data...\n")

# Get peptide sequence column
peptide_col <- "Peptide Sequence"
if (!peptide_col %in% colnames(chop_data)) {
  peptide_col <- grep("Peptide", colnames(chop_data), ignore.case=TRUE, value=TRUE)[1]
}

# Create filtered dataset
fl57_data <- chop_data %>%
  select(all_of(c(peptide_col, liver_col, lung_col, "Gene", "Protein"))) %>%
  rename(
    Peptide = all_of(peptide_col),
    FL57_Liver_Intensity = all_of(liver_col),
    FL57_Lung_Intensity = all_of(lung_col)
  ) %>%
  mutate(
    FL57_Liver_Intensity = as.numeric(FL57_Liver_Intensity),
    FL57_Lung_Intensity = as.numeric(FL57_Lung_Intensity)
  ) %>%
  # Replace NA with 0
  mutate(
    FL57_Liver_Intensity = ifelse(is.na(FL57_Liver_Intensity), 0, FL57_Liver_Intensity),
    FL57_Lung_Intensity = ifelse(is.na(FL57_Lung_Intensity), 0, FL57_Lung_Intensity)
  )

# Save filtered dataset
cat("Saving filtered dataset...\n")
write.table(fl57_data, file=filtered_output_file, sep="\t", quote=FALSE, row.names=FALSE)

# Calculate summary statistics
cat("Calculating summary statistics...\n")

# Peptides detected in at least one sample (intensity > 0)
peptides_at_least_one <- fl57_data %>%
  filter(FL57_Liver_Intensity > 0 | FL57_Lung_Intensity > 0)

# Peptides detected in both samples (intensity > 0 in both)
peptides_both <- fl57_data %>%
  filter(FL57_Liver_Intensity > 0 & FL57_Lung_Intensity > 0)

# Peptides detected only in liver
peptides_liver_only <- fl57_data %>%
  filter(FL57_Liver_Intensity > 0 & FL57_Lung_Intensity == 0)

# Peptides detected only in lung
peptides_lung_only <- fl57_data %>%
  filter(FL57_Liver_Intensity == 0 & FL57_Lung_Intensity > 0)

# Calculate correlations for peptides detected in both samples
if (nrow(peptides_both) > 2) {
  # Add small constant to avoid log(0)
  liver_intensities <- peptides_both$FL57_Liver_Intensity + 1
  lung_intensities <- peptides_both$FL57_Lung_Intensity + 1
  
  # Pearson correlation
  pearson_cor <- cor(liver_intensities, lung_intensities, method="pearson")
  pearson_r2 <- pearson_cor^2
  
  # Spearman correlation
  spearman_cor <- cor(liver_intensities, lung_intensities, method="spearman")
  spearman_r2 <- spearman_cor^2
  
  # Linear model for R²
  lm_model <- lm(log10(lung_intensities) ~ log10(liver_intensities))
  lm_r2 <- summary(lm_model)$r.squared
  
} else {
  pearson_cor <- pearson_r2 <- spearman_cor <- spearman_r2 <- lm_r2 <- NA
}

# Calculate descriptive statistics
liver_stats <- list(
  total_detected = sum(fl57_data$FL57_Liver_Intensity > 0),
  mean_intensity = mean(fl57_data$FL57_Liver_Intensity[fl57_data$FL57_Liver_Intensity > 0], na.rm=TRUE),
  median_intensity = median(fl57_data$FL57_Liver_Intensity[fl57_data$FL57_Liver_Intensity > 0], na.rm=TRUE),
  max_intensity = max(fl57_data$FL57_Liver_Intensity, na.rm=TRUE)
)

lung_stats <- list(
  total_detected = sum(fl57_data$FL57_Lung_Intensity > 0),
  mean_intensity = mean(fl57_data$FL57_Lung_Intensity[fl57_data$FL57_Lung_Intensity > 0], na.rm=TRUE),
  median_intensity = median(fl57_data$FL57_Lung_Intensity[fl57_data$FL57_Lung_Intensity > 0], na.rm=TRUE),
  max_intensity = max(fl57_data$FL57_Lung_Intensity, na.rm=TRUE)
)

# Create summary statistics output
summary_text <- paste0(
  "FL57 LIVER VS LUNG PEPTIDE CORRELATION ANALYSIS SUMMARY\n",
  "=====================================================\n\n",
  "DATASET OVERVIEW:\n",
  "Total peptides in dataset: ", nrow(fl57_data), "\n",
  "Peptides detected in liver only: ", nrow(peptides_liver_only), "\n",
  "Peptides detected in lung only: ", nrow(peptides_lung_only), "\n",
  "Peptides detected in both samples: ", nrow(peptides_both), "\n",
  "Peptides detected in at least one sample: ", nrow(peptides_at_least_one), "\n\n",
  
  "LIVER SAMPLE STATISTICS:\n",
  "Total peptides detected: ", liver_stats$total_detected, "\n",
  "Mean intensity (detected peptides): ", sprintf("%.2e", liver_stats$mean_intensity), "\n",
  "Median intensity (detected peptides): ", sprintf("%.2e", liver_stats$median_intensity), "\n",
  "Maximum intensity: ", sprintf("%.2e", liver_stats$max_intensity), "\n\n",
  
  "LUNG SAMPLE STATISTICS:\n",
  "Total peptides detected: ", lung_stats$total_detected, "\n",
  "Mean intensity (detected peptides): ", sprintf("%.2e", lung_stats$mean_intensity), "\n",
  "Median intensity (detected peptides): ", sprintf("%.2e", lung_stats$median_intensity), "\n",
  "Maximum intensity: ", sprintf("%.2e", lung_stats$max_intensity), "\n\n",
  
  "CORRELATION ANALYSIS (peptides detected in both samples):\n",
  "Number of peptides: ", nrow(peptides_both), "\n",
  "Pearson correlation: ", sprintf("%.3f", pearson_cor), "\n",
  "Pearson R²: ", sprintf("%.3f", pearson_r2), "\n",
  "Spearman correlation: ", sprintf("%.3f", spearman_cor), "\n",
  "Spearman R²: ", sprintf("%.3f", spearman_r2), "\n",
  "Linear model R² (log-transformed): ", sprintf("%.3f", lm_r2), "\n\n",
  
  "OVERLAP STATISTICS:\n",
  "Percentage of liver peptides also found in lung: ", 
  sprintf("%.1f%%", 100 * nrow(peptides_both) / liver_stats$total_detected), "\n",
  "Percentage of lung peptides also found in liver: ", 
  sprintf("%.1f%%", 100 * nrow(peptides_both) / lung_stats$total_detected), "\n",
  "Overall overlap (Jaccard index): ", 
  sprintf("%.3f", nrow(peptides_both) / nrow(peptides_at_least_one)), "\n"
)

# Save summary statistics
cat("Saving summary statistics...\n")
writeLines(summary_text, stats_output_file)

# Create correlation plot matching the MSKCC vs CHOP style
cat("Creating correlation plot...\n")

# Prepare data for plotting - only peptides detected in both samples
plot_data <- peptides_both %>%
  mutate(
    Liver_Log = log10(FL57_Liver_Intensity),
    Lung_Log = log10(FL57_Lung_Intensity)
  ) %>%
  filter(is.finite(Liver_Log) & is.finite(Lung_Log))

if (nrow(plot_data) > 0) {
  # Create the plot matching your MSKCC vs CHOP style exactly
  p1 <- ggplot(plot_data, aes(x = Liver_Log, y = Lung_Log)) +
    geom_point(color = "#2E8B57", alpha = 0.4, size = 0.8) +  # Smaller points, more transparent
    geom_smooth(method = "lm", color = "red", se = FALSE, linewidth = 1.2) +  # Red trend line
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    labs(
      title = "FL57 Liver vs Lung Intensity Correlation",
      subtitle = paste0("R² = ", sprintf("%.3f", lm_r2), ", n = ", nrow(plot_data)),
      x = "FL57 Liver Intensity (log10)",
      y = "FL57 Lung Intensity (log10)"
    ) +
    # Set axis limits to match reference plots better
    scale_x_continuous(limits = c(0, max(plot_data$Liver_Log) * 1.1)) +
    scale_y_continuous(limits = c(0, max(plot_data$Lung_Log) * 1.1))
  
  # Save the correlation plot
  ggsave(plot_output_file, plot = p1, width = 8, height = 6, dpi = 300, bg = "white")
  
} else {
  cat("Warning: No peptides detected in both samples for plotting.\n")
}

# Create unique peptides analysis and plot
cat("Creating unique peptides correlation analysis...\n")

# Create data for unique peptides plot
unique_plot_file <- "FL57_unique_peptides_plot.png"

# For unique peptides, we'll show the counts in a different way
# Since we can't correlate single numbers, let's create a visualization showing the distribution

# Create a summary plot for unique peptides
unique_summary <- data.frame(
  Sample = c("Liver Only", "Lung Only", "Both Samples"),
  Count = c(nrow(peptides_liver_only), nrow(peptides_lung_only), nrow(peptides_both))
)

# Create bar plot for unique peptides
p2 <- ggplot(unique_summary, aes(x = Sample, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("Liver Only" = "#FF6B6B", "Lung Only" = "#4ECDC4", "Both Samples" = "#2E8B57")) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
    legend.position = "none"
  ) +
  labs(
    title = "FL57 Peptide Detection Distribution",
    subtitle = "Tissue-specific vs Shared Peptides",
    x = "Detection Pattern",
    y = "Number of Peptides"
  ) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4)

# Save the unique peptides plot
ggsave(unique_plot_file, plot = p2, width = 8, height = 6, dpi = 300, bg = "white")

# Additionally, create a scatter plot showing intensity vs unique detection
# This will show liver-only vs lung-only peptides in a more meaningful way
unique_scatter_file <- "FL57_unique_peptides_scatter.png"

# Prepare data for unique peptides scatter
liver_only_data <- peptides_liver_only %>%
  mutate(
    Sample_Type = "Liver Only",
    Intensity = FL57_Liver_Intensity,
    Log_Intensity = log10(FL57_Liver_Intensity)
  )

lung_only_data <- peptides_lung_only %>%
  mutate(
    Sample_Type = "Lung Only", 
    Intensity = FL57_Lung_Intensity,
    Log_Intensity = log10(FL57_Lung_Intensity)
  )

unique_combined <- bind_rows(liver_only_data, lung_only_data)

if (nrow(unique_combined) > 0) {
  # Create intensity distribution plot for unique peptides
  p3 <- ggplot(unique_combined, aes(x = Log_Intensity, fill = Sample_Type)) +
    geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
    scale_fill_manual(values = c("Liver Only" = "#FF6B6B", "Lung Only" = "#4ECDC4")) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ) +
    labs(
      title = "FL57 Tissue-Specific Peptide Intensity Distribution",
      subtitle = "Comparing Liver-Only vs Lung-Only Peptides",
      x = "Peptide Intensity (log10)",
      y = "Number of Peptides",
      fill = "Sample Type"
    )
  
  # Save the unique peptides scatter plot
  ggsave(unique_scatter_file, plot = p3, width = 10, height = 6, dpi = 300, bg = "white")
}

# Print summary to console
cat("\n=====================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=====================================================================\n\n")

cat("FILES CREATED:\n")
cat("1.", filtered_output_file, "- Filtered dataset with FL57 liver and lung data\n")
cat("2.", plot_output_file, "- Intensity correlation plot\n")
cat("3.", "FL57_unique_peptides_plot.png", "- Unique peptides bar chart\n")
cat("4.", "FL57_unique_peptides_scatter.png", "- Tissue-specific peptide intensity distributions\n")
cat("5.", stats_output_file, "- Summary statistics\n\n")

cat("SUMMARY STATISTICS:\n")
cat(summary_text)

cat("\nAnalysis completed successfully!\n")



# CHOP FL57 Liver vs Lung Peptide Intensity Correlation Analysis
# This script analyzes peptide intensity correlation between liver and lung samples from patient FL57

# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)

# Set working directory (adjust as needed)
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation")

cat("=====================================================================\n")
cat("CHOP FL57 LIVER VS LUNG PEPTIDE CORRELATION ANALYSIS\n")
cat("=====================================================================\n\n")

# Define file paths
input_file <- "final_CHOP_combined_peptide.tsv"
filtered_output_file <- "FL57_liver_lung_comparison.tsv"
plot_output_file <- "FL57_liver_lung_correlation_plot.png"
stats_output_file <- "FL57_liver_lung_summary_stats.txt"

# Read the integrated CHOP data
cat("Reading integrated CHOP data...\n")
chop_data <- read.delim(input_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Check if the required columns exist
liver_col <- "H5 FL57Liver Intensity"
lung_col <- "H6 FL57Lung Intensity"

if (!liver_col %in% colnames(chop_data)) {
  # Try alternative column naming
  liver_col <- grep("FL57Liver.*Intensity", colnames(chop_data), ignore.case=TRUE, value=TRUE)[1]
}

if (!lung_col %in% colnames(chop_data)) {
  # Try alternative column naming  
  lung_col <- grep("FL57Lung.*Intensity", colnames(chop_data), ignore.case=TRUE, value=TRUE)[1]
}

cat("Using columns:\n")
cat("Liver:", liver_col, "\n")
cat("Lung:", lung_col, "\n")

# Extract relevant data for FL57 liver and lung comparison
cat("Extracting FL57 liver and lung data...\n")

# Get peptide sequence column
peptide_col <- "Peptide Sequence"
if (!peptide_col %in% colnames(chop_data)) {
  peptide_col <- grep("Peptide", colnames(chop_data), ignore.case=TRUE, value=TRUE)[1]
}

# Create filtered dataset
fl57_data <- chop_data %>%
  select(all_of(c(peptide_col, liver_col, lung_col, "Gene", "Protein"))) %>%
  rename(
    Peptide = all_of(peptide_col),
    FL57_Liver_Intensity = all_of(liver_col),
    FL57_Lung_Intensity = all_of(lung_col)
  ) %>%
  mutate(
    FL57_Liver_Intensity = as.numeric(FL57_Liver_Intensity),
    FL57_Lung_Intensity = as.numeric(FL57_Lung_Intensity)
  ) %>%
  # Replace NA with 0
  mutate(
    FL57_Liver_Intensity = ifelse(is.na(FL57_Liver_Intensity), 0, FL57_Liver_Intensity),
    FL57_Lung_Intensity = ifelse(is.na(FL57_Lung_Intensity), 0, FL57_Lung_Intensity)
  )

# Save filtered dataset
cat("Saving filtered dataset...\n")
write.table(fl57_data, file=filtered_output_file, sep="\t", quote=FALSE, row.names=FALSE)

# Calculate summary statistics
cat("Calculating summary statistics...\n")

# Peptides detected in at least one sample (intensity > 0)
peptides_at_least_one <- fl57_data %>%
  filter(FL57_Liver_Intensity > 0 | FL57_Lung_Intensity > 0)

# Peptides detected in both samples (intensity > 0 in both)
peptides_both <- fl57_data %>%
  filter(FL57_Liver_Intensity > 0 & FL57_Lung_Intensity > 0)

# Peptides detected only in liver
peptides_liver_only <- fl57_data %>%
  filter(FL57_Liver_Intensity > 0 & FL57_Lung_Intensity == 0)

# Peptides detected only in lung
peptides_lung_only <- fl57_data %>%
  filter(FL57_Liver_Intensity == 0 & FL57_Lung_Intensity > 0)

# Calculate correlations for peptides detected in both samples
if (nrow(peptides_both) > 2) {
  # Add small constant to avoid log(0)
  liver_intensities <- peptides_both$FL57_Liver_Intensity + 1
  lung_intensities <- peptides_both$FL57_Lung_Intensity + 1
  
  # Pearson correlation
  pearson_cor <- cor(liver_intensities, lung_intensities, method="pearson")
  pearson_r2 <- pearson_cor^2
  
  # Spearman correlation
  spearman_cor <- cor(liver_intensities, lung_intensities, method="spearman")
  spearman_r2 <- spearman_cor^2
  
  # Linear model for R²
  lm_model <- lm(log10(lung_intensities) ~ log10(liver_intensities))
  lm_r2 <- summary(lm_model)$r.squared
  
} else {
  pearson_cor <- pearson_r2 <- spearman_cor <- spearman_r2 <- lm_r2 <- NA
}

# Calculate descriptive statistics
liver_stats <- list(
  total_detected = sum(fl57_data$FL57_Liver_Intensity > 0),
  mean_intensity = mean(fl57_data$FL57_Liver_Intensity[fl57_data$FL57_Liver_Intensity > 0], na.rm=TRUE),
  median_intensity = median(fl57_data$FL57_Liver_Intensity[fl57_data$FL57_Liver_Intensity > 0], na.rm=TRUE),
  max_intensity = max(fl57_data$FL57_Liver_Intensity, na.rm=TRUE)
)

lung_stats <- list(
  total_detected = sum(fl57_data$FL57_Lung_Intensity > 0),
  mean_intensity = mean(fl57_data$FL57_Lung_Intensity[fl57_data$FL57_Lung_Intensity > 0], na.rm=TRUE),
  median_intensity = median(fl57_data$FL57_Lung_Intensity[fl57_data$FL57_Lung_Intensity > 0], na.rm=TRUE),
  max_intensity = max(fl57_data$FL57_Lung_Intensity, na.rm=TRUE)
)

# Create summary statistics output
summary_text <- paste0(
  "FL57 LIVER VS LUNG PEPTIDE CORRELATION ANALYSIS SUMMARY\n",
  "=====================================================\n\n",
  "DATASET OVERVIEW:\n",
  "Total peptides in dataset: ", nrow(fl57_data), "\n",
  "Peptides detected in liver only: ", nrow(peptides_liver_only), "\n",
  "Peptides detected in lung only: ", nrow(peptides_lung_only), "\n",
  "Peptides detected in both samples: ", nrow(peptides_both), "\n",
  "Peptides detected in at least one sample: ", nrow(peptides_at_least_one), "\n\n",
  
  "LIVER SAMPLE STATISTICS:\n",
  "Total peptides detected: ", liver_stats$total_detected, "\n",
  "Mean intensity (detected peptides): ", sprintf("%.2e", liver_stats$mean_intensity), "\n",
  "Median intensity (detected peptides): ", sprintf("%.2e", liver_stats$median_intensity), "\n",
  "Maximum intensity: ", sprintf("%.2e", liver_stats$max_intensity), "\n\n",
  
  "LUNG SAMPLE STATISTICS:\n",
  "Total peptides detected: ", lung_stats$total_detected, "\n",
  "Mean intensity (detected peptides): ", sprintf("%.2e", lung_stats$mean_intensity), "\n",
  "Median intensity (detected peptides): ", sprintf("%.2e", lung_stats$median_intensity), "\n",
  "Maximum intensity: ", sprintf("%.2e", lung_stats$max_intensity), "\n\n",
  
  "CORRELATION ANALYSIS (peptides detected in both samples):\n",
  "Number of peptides: ", nrow(peptides_both), "\n",
  "Pearson correlation: ", sprintf("%.3f", pearson_cor), "\n",
  "Pearson R²: ", sprintf("%.3f", pearson_r2), "\n",
  "Spearman correlation: ", sprintf("%.3f", spearman_cor), "\n",
  "Spearman R²: ", sprintf("%.3f", spearman_r2), "\n",
  "Linear model R² (log-transformed): ", sprintf("%.3f", lm_r2), "\n\n",
  
  "OVERLAP STATISTICS:\n",
  "Percentage of liver peptides also found in lung: ", 
  sprintf("%.1f%%", 100 * nrow(peptides_both) / liver_stats$total_detected), "\n",
  "Percentage of lung peptides also found in liver: ", 
  sprintf("%.1f%%", 100 * nrow(peptides_both) / lung_stats$total_detected), "\n",
  "Overall overlap (Jaccard index): ", 
  sprintf("%.3f", nrow(peptides_both) / nrow(peptides_at_least_one)), "\n"
)

# Save summary statistics
cat("Saving summary statistics...\n")
writeLines(summary_text, stats_output_file)

# Create correlation plot matching the MSKCC vs CHOP style
cat("Creating correlation plot...\n")

# Prepare data for plotting - only peptides detected in both samples
plot_data <- peptides_both %>%
  mutate(
    Liver_Log = log10(FL57_Liver_Intensity + 1),
    Lung_Log = log10(FL57_Lung_Intensity + 1)
  ) %>%
  filter(is.finite(Liver_Log) & is.finite(Lung_Log))

if (nrow(plot_data) > 0) {
  # Create the plot matching your MSKCC vs CHOP style
  p <- ggplot(plot_data, aes(x = Liver_Log, y = Lung_Log)) +
    geom_point(color = "#2E8B57", alpha = 0.6) +  # Same green color as your plots
    geom_smooth(method = "lm", color = "red", se = FALSE, linewidth = 1) +  # Red trend line
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_line(color = "grey95", size = 0.25),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    ) +
    labs(
      title = "FL57 Liver vs Lung Intensity Correlation",
      subtitle = paste0("R² = ", sprintf("%.3f", lm_r2), ", n = ", nrow(plot_data)),
      x = "FL57 Liver Intensity (log10)",
      y = "FL57 Lung Intensity (log10)"
    ) +
    scale_x_continuous(trans = "identity") +
    scale_y_continuous(trans = "identity")
  
  # Save the plot
  ggsave(plot_output_file, plot = p, width = 8, height = 6, dpi = 300, bg = "white")
  
} else {
  cat("Warning: No peptides detected in both samples for plotting.\n")
}

# Print summary to console
cat("\n=====================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=====================================================================\n\n")

cat("FILES CREATED:\n")
cat("1.", filtered_output_file, "- Filtered dataset with FL57 liver and lung data\n")
cat("2.", plot_output_file, "- Correlation plot\n")
cat("3.", stats_output_file, "- Summary statistics\n\n")

cat("SUMMARY STATISTICS:\n")
cat(summary_text)

cat("\nAnalysis completed successfully!\n")