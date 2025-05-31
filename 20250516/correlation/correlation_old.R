# Immunopeptidomics Datasets Correlation Analysis
# Script for comparing MSKCC and CHOP immunopeptidomics datasets
# Date: May 16, 2025

# Load required libraries
library(tidyverse)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(corrplot)
library(scales)
library(gridExtra)

# Set working directory if needed
# setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation")

# Define file paths
mskcc_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation/MSKCC_combined_peptides.tsv"
chop_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation/CHOP_combined_peptide.tsv"


# We need to generate a new combined_CHOP file. 51 and 88 are not included in this one above.
# 

# Read data
message("Reading MSKCC data...")
mskcc_data <- read.delim(mskcc_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

message("Reading CHOP data...")
chop_data <- read.delim(chop_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Display basic information about datasets
message("Basic dataset information:")
mskcc_summary <- data.frame(
  Dataset = "MSKCC",
  Rows = nrow(mskcc_data),
  Columns = ncol(mskcc_data),
  Unique_Peptides = length(unique(mskcc_data$Peptide)),
  Sample_Count = length(unique(mskcc_data$SampleID))
)

# Extract sample IDs from CHOP columns (those with "Spectral Count" in name)
chop_sample_cols <- grep("Spectral Count", colnames(chop_data), value=TRUE)
chop_samples <- gsub(" Spectral Count", "", chop_sample_cols)

chop_summary <- data.frame(
  Dataset = "CHOP",
  Rows = nrow(chop_data),
  Columns = ncol(chop_data),
  Unique_Peptides = length(unique(chop_data$`Peptide Sequence`)),
  Sample_Count = length(chop_samples)
)

summary_df <- rbind(mskcc_summary, chop_summary)
print(summary_df)

# ============= DATA CLEANING AND PREPARATION =============

# Standardize column names for consistency
# Find the peptide sequence column in CHOP data
message("CHOP column names before standardization:")
print(head(colnames(chop_data)))

peptide_col <- grep("Peptide|peptide", colnames(chop_data), value=TRUE)[1]
message(paste("Using peptide column:", peptide_col))

# Reshape MSKCC data to have consistent format with sample IDs
message("Reshaping MSKCC data...")
# Create a list of MSKCC samples
mskcc_samples <- unique(mskcc_data$SampleID)
message(paste("MSKCC samples:", paste(mskcc_samples, collapse=", ")))

# Reshape MSKCC data using a more basic approach
message("Reshaping MSKCC data (alternative approach)...")

# First, handle non-unique peptide-sample combinations
mskcc_data_agg <- mskcc_data %>%
  group_by(Peptide, SampleID) %>%
  summarize(
    Spectral_Count = sum(as.numeric(!!sym(spectral_count_col))),
    Intensity_Value = sum(as.numeric(!!sym(intensity_col))),
    .groups = 'drop'
  )

# Now reshape for spectral counts
mskcc_sc_wide <- reshape2::dcast(
  mskcc_data_agg, 
  Peptide ~ SampleID, 
  value.var = "Spectral_Count", 
  fun.aggregate = sum,
  fill = 0
)

# Rename the columns to add prefix
colnames(mskcc_sc_wide)[-1] <- paste0("MSKCC_", colnames(mskcc_sc_wide)[-1])

# Reshape for intensity
mskcc_int_wide <- reshape2::dcast(
  mskcc_data_agg, 
  Peptide ~ SampleID, 
  value.var = "Intensity_Value", 
  fun.aggregate = sum,
  fill = 0
)

# Rename the columns to add prefix and suffix
colnames(mskcc_int_wide)[-1] <- paste0("MSKCC_", colnames(mskcc_int_wide)[-1], "_Intensity")

# Merge spectral count and intensity data
mskcc_wide <- mskcc_sc_wide %>%
  left_join(mskcc_int_wide, by = "Peptide")

# Prepare CHOP data
message("Preparing CHOP data...")
message("CHOP column names:")
print(colnames(chop_data))

# Find the peptide sequence column in CHOP data
peptide_col <- grep("Peptide", colnames(chop_data), value=TRUE)[1]
message(paste("Using peptide column:", peptide_col))

# Extract spectral count columns and intensity columns
chop_sc_cols <- grep("Spectral Count", colnames(chop_data), value=TRUE)
chop_int_cols <- grep("Intensity$", colnames(chop_data), value=TRUE)

# Rename columns for consistency
new_sc_cols <- paste0("CHOP_", gsub(" Spectral Count", "", chop_sc_cols))
new_int_cols <- paste0("CHOP_", gsub(" Intensity", "", chop_int_cols), "_Intensity")

# Create a clean CHOP dataset with renamed columns
chop_clean <- chop_data %>%
  select(!!sym(peptide_col), all_of(chop_sc_cols), all_of(chop_int_cols))

# Rename the peptide column to standardize
colnames(chop_clean)[1] <- "Peptide"

# Rename the spectral count and intensity columns
colnames(chop_clean)[colnames(chop_clean) %in% chop_sc_cols] <- new_sc_cols
colnames(chop_clean)[colnames(chop_clean) %in% chop_int_cols] <- new_int_cols

# ============= SANITY CHECKS =============
message("Performing sanity checks...")

# 1. Check peptide length distributions
# Find the peptide length column in each dataset
mskcc_length_col <- grep("Length|length", colnames(mskcc_data), value=TRUE)[1]
chop_length_col <- grep("Length|length", colnames(chop_data), value=TRUE)[1]

message(paste("MSKCC peptide length column:", mskcc_length_col))
message(paste("CHOP peptide length column:", chop_length_col))

# Create peptide length distributions
mskcc_length_dist <- mskcc_data %>%
  group_by(!!sym(mskcc_length_col)) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100,
         Dataset = "MSKCC",
         `Peptide Length` = !!sym(mskcc_length_col))

chop_length_dist <- chop_data %>%
  group_by(!!sym(chop_length_col)) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100,
         Dataset = "CHOP",
         `Peptide Length` = !!sym(chop_length_col))

# Combine the distributions for plotting
length_dist <- bind_rows(mskcc_length_dist, chop_length_dist)

# Plot peptide length distribution
p_length <- ggplot(length_dist, aes(x=factor(`Peptide Length`), y=Percentage, fill=Dataset)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  labs(title="Peptide Length Distribution", 
       x="Peptide Length", 
       y="Percentage of Peptides") +
  scale_fill_brewer(palette="Set1")

# 2. Checking common peptides between datasets
message("Checking peptide overlap between datasets...")
mskcc_peptides <- unique(mskcc_data$Peptide)
chop_peptides <- unique(chop_data$Peptide)
common_peptides <- intersect(mskcc_peptides, chop_peptides)

peptide_overlap <- data.frame(
  Category = c("MSKCC only", "CHOP only", "Common"),
  Count = c(
    length(setdiff(mskcc_peptides, chop_peptides)),
    length(setdiff(chop_peptides, mskcc_peptides)),
    length(common_peptides)
  )
)

peptide_overlap$Percentage <- peptide_overlap$Count / sum(peptide_overlap$Count) * 100

# Plot peptide overlap
p_overlap <- ggplot(peptide_overlap, aes(x="", y=Count, fill=Category)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_minimal() +
  labs(title="Peptide Overlap Between Datasets", 
       fill="Category") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5))

# 3. Data quality metrics
# Check for missing values
mskcc_missing <- colSums(is.na(mskcc_data))
chop_missing <- colSums(is.na(chop_data))

# Check for zero values in intensity columns
mskcc_zeros <- sum(mskcc_data$Intensity == 0) / nrow(mskcc_data) * 100
chop_zeros <- sapply(chop_int_cols, function(col) sum(chop_data[[col]] == 0) / nrow(chop_data) * 100)

# Create data quality summary
quality_summary <- data.frame(
  Metric = c(
    "MSKCC - Peptides with zero intensity (%)",
    paste0("CHOP - ", gsub("CHOP_|_Intensity", "", new_int_cols), " - Peptides with zero intensity (%)")
  ),
  Value = c(mskcc_zeros, chop_zeros)
)

# ============= MERGE DATASETS FOR CORRELATION ANALYSIS =============
message("Merging datasets for correlation analysis...")
# Merge by peptide sequence
merged_data <- mskcc_wide %>%
  right_join(chop_clean, by = "Peptide")

# Count peptides present in each dataset
dataset_presence <- data.frame(
  MSKCC_only = sum(!is.na(mskcc_wide$Peptide) & is.na(chop_clean$Peptide[match(mskcc_wide$Peptide, chop_clean$Peptide)])),
  CHOP_only = sum(is.na(mskcc_wide$Peptide[match(chop_clean$Peptide, mskcc_wide$Peptide)]) & !is.na(chop_clean$Peptide)),
  Common = sum(!is.na(mskcc_wide$Peptide[match(chop_clean$Peptide, mskcc_wide$Peptide)]) & !is.na(chop_clean$Peptide))
)

# ============= CORRELATION ANALYSIS =============
message("Performing correlation analysis...")

# Function to calculate correlation and generate scatter plot
create_correlation_plot <- function(x, y, x_name, y_name, log_scale = TRUE) {
  # Check if vectors are valid
  if(is.null(x) || is.null(y) || all(is.na(x)) || all(is.na(y))) {
    message(paste("Invalid data for correlation between", x_name, "and", y_name))
    return(list(
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "Insufficient data for correlation") +
        theme_void(),
      r_squared = NA,
      n = 0
    ))
  }
  
  # Convert to numeric if not already
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  # Remove zeros and NA values if using log scale
  if(log_scale) {
    valid_idx <- which(x > 0 & y > 0 & !is.na(x) & !is.na(y))
  } else {
    valid_idx <- which(!is.na(x) & !is.na(y))
  }
  
  # Check if we have enough valid data points
  if(length(valid_idx) < 3) {
    message(paste("Insufficient data points for correlation between", x_name, "and", y_name))
    return(list(
      plot = ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "Insufficient data for correlation") +
        theme_void(),
      r_squared = NA,
      n = length(valid_idx)
    ))
  }
  
  x_valid <- x[valid_idx]
  y_valid <- y[valid_idx]
  
  # Calculate correlation
  cor_value <- cor(x_valid, y_valid, method = "pearson")
  r_squared <- cor_value^2
  
  # Create scatter plot
  if(log_scale) {
    p <- ggplot(data.frame(x = x_valid, y = y_valid), aes(x = x, y = y)) +
      geom_point(alpha = 0.5, color = "blue") +
      scale_x_log10() +
      scale_y_log10() +
      geom_smooth(method = "lm", color = "red") +
      labs(
        title = paste("Correlation of", x_name, "vs", y_name),
        subtitle = paste("R² =", round(r_squared, 3), "| n =", length(x_valid)),
        x = paste(x_name, "(log10)"),
        y = paste(y_name, "(log10)")
      ) +
      theme_minimal() +
      annotation_logticks()
  } else {
    p <- ggplot(data.frame(x = x_valid, y = y_valid), aes(x = x, y = y)) +
      geom_point(alpha = 0.5, color = "blue") +
      geom_smooth(method = "lm", color = "red") +
      labs(
        title = paste("Correlation of", x_name, "vs", y_name),
        subtitle = paste("R² =", round(r_squared, 3), "| n =", length(x_valid)),
        x = x_name,
        y = y_name
      ) +
      theme_minimal()
  }
  
  return(list(plot = p, r_squared = r_squared, n = length(x_valid)))
}

# Calculate correlation matrix for spectral counts
mskcc_samples_clean <- gsub("MSKCC_", "", grep("^MSKCC_\\d+$", colnames(merged_data), value = TRUE))
chop_samples_clean <- gsub("CHOP_", "", grep("^CHOP_H\\d+_", colnames(merged_data), value = TRUE))
chop_samples_clean <- unique(gsub("_Intensity", "", chop_samples_clean))

# Create a list to store all correlation results
correlation_results <- list()
correlation_matrix <- matrix(NA, nrow = length(mskcc_samples_clean), ncol = length(chop_samples_clean))
rownames(correlation_matrix) <- mskcc_samples_clean
colnames(correlation_matrix) <- chop_samples_clean

# Create a list of plots
corr_plots <- list()

# Compare spectral counts
for (i in seq_along(mskcc_samples_clean)) {
  mskcc_col <- paste0("MSKCC_", mskcc_samples_clean[i])
  
  for (j in seq_along(chop_samples_clean)) {
    chop_col <- paste0("CHOP_", chop_samples_clean[j])
    
    # Spectral count correlation
    sc_result <- create_correlation_plot(
      merged_data[[mskcc_col]], 
      merged_data[[chop_col]], 
      paste0("MSKCC ", mskcc_samples_clean[i], " Spectral Count"), 
      paste0("CHOP ", chop_samples_clean[j], " Spectral Count"),
      log_scale = FALSE
    )
    
    corr_plots[[paste(mskcc_col, chop_col, "SC", sep = "_")]] <- sc_result$plot
    correlation_matrix[i, j] <- sc_result$r_squared
    
    # Intensity correlation
    int_result <- create_correlation_plot(
      merged_data[[paste0(mskcc_col, "_Intensity")]], 
      merged_data[[paste0(chop_col, "_Intensity")]], 
      paste0("MSKCC ", mskcc_samples_clean[i], " Intensity"), 
      paste0("CHOP ", chop_samples_clean[j], " Intensity"),
      log_scale = TRUE
    )
    
    corr_plots[[paste(mskcc_col, chop_col, "INT", sep = "_")]] <- int_result$plot
  }
}

# Create heatmap of correlation values
correlation_heatmap <- pheatmap(
  correlation_matrix,
  display_numbers = TRUE,
  number_format = "%.3f",
  main = "R² Values Between MSKCC and CHOP Samples (Spectral Counts)",
  fontsize = 10,
  fontsize_number = 8
)

# ============= ADDITIONAL ANALYSIS =============
# Analyze peptide prevalence (detection frequency) across samples
message("Analyzing peptide prevalence...")

# Function to analyze peptide prevalence
analyze_prevalence <- function(data, dataset_name) {
  # Identify samples and their spectral count columns
  if(dataset_name == "MSKCC") {
    sample_cols <- grep("^MSKCC_\\d+$", colnames(data), value = TRUE)
    # For each peptide, count in how many samples it appears
    presence_matrix <- data[, sample_cols] > 0
  } else { # CHOP
    sample_cols <- grep("^CHOP_H\\d+_[^I]", colnames(data), value = TRUE)
    # For each peptide, count in how many samples it appears
    presence_matrix <- data[, sample_cols] > 0
  }
  
  presence_count <- rowSums(presence_matrix)
  prevalence_summary <- table(presence_count) / length(presence_count) * 100
  
  return(list(
    prevalence_counts = presence_count,
    summary = prevalence_summary
  ))
}

mskcc_prevalence <- analyze_prevalence(merged_data, "MSKCC")
chop_prevalence <- analyze_prevalence(merged_data, "CHOP")

# Create a data frame for visualization
prevalence_data <- data.frame(
  Peptide = merged_data$Peptide,
  MSKCC_Samples_Count = mskcc_prevalence$prevalence_counts,
  CHOP_Samples_Count = chop_prevalence$prevalence_counts
)

# Analyze prevalence correlation
prevalence_corr <- create_correlation_plot(
  prevalence_data$MSKCC_Samples_Count,
  prevalence_data$CHOP_Samples_Count,
  "MSKCC Sample Detection Count",
  "CHOP Sample Detection Count",
  log_scale = FALSE
)

# ============= OUTPUT RESULTS =============
message("Generating result plots...")

# Create a multi-panel figure of the main findings
grid.arrange(
  p_length, p_overlap, prevalence_corr$plot,
  corr_plots[[1]], corr_plots[[2]], # Example correlation plots
  ncol = 2
)

# Print summary statistics
message("\nData Summary:")
print(summary_df)

message("\nPeptide Overlap:")
print(peptide_overlap)

message("\nData Quality Metrics:")
print(quality_summary)

message("\nCorrelation Analysis:")
print(correlation_matrix)

message("\nAnalysis completed successfully!")