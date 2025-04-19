# Simplified script to analyze DNAJB1_PRKACA fusion junction peptides across samples
# with focus on the 6 key peptides that were spiked into sample 51S

# Load required packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(writexl)  # For Excel output

# Define the path to your data files
data_path <- "/Users/dinarabadi/Desktop/2704DR"

# Define the 6 key peptides (used for spiking into 51S)
key_peptides <- c(
  "EEVKEFLAK",
  "YGEEVKEFL",
  "RYGEEVKEF",
  "RYGEEVKEFL", 
  "EIFDRYGEEV",
  "IFDRYGEEV"
)

# Read all TSV files in the directory and extract sample IDs from filenames
files <- list.files(path = data_path, pattern = ".*_peptides\\.tsv$|.*_untargeted_peptide\\.tsv$", full.names = TRUE)

if (length(files) == 0) {
  stop("No peptide TSV files found in ", data_path)
}

cat("Found", length(files), "peptide files for analysis\n")

# Extract sample IDs from filenames
all_data <- list()
sample_ids <- c()

for (file in files) {
  filename <- basename(file)
  
  # Extract sample ID from filename
  # Pattern: comes after "2704DR_IMP_014_20250403_DDA_2CV_" and before "_peptides.tsv" or "_untargeted_peptide.tsv"
  sample_id <- gsub(".*2704DR_IMP_\\d+_\\d+_DDA_2CV_([^_]+)_.*\\.tsv$", "\\1", filename)
  
  # For filenames like 2704DR_IMP_014_20250403_DDA_2CV_51_01_peptides.tsv, handle the extra "_01" part
  if (grepl("_\\d+_peptides\\.tsv$", filename)) {
    sample_id <- gsub("(.*)_\\d+$", "\\1", sample_id)
  }
  
  cat("Reading file:", filename, "- Sample ID:", sample_id, "\n")
  
  # Read the file
  data <- read.delim(file, stringsAsFactors = FALSE)
  
  # Add a column for sample ID
  data$Sample_ID <- sample_id
  data$Filename <- filename
  
  # Add to our list
  all_data[[sample_id]] <- data
  sample_ids <- c(sample_ids, sample_id)
}

# Get unique sample IDs (in case multiple files had the same sample ID)
sample_ids <- unique(sample_ids)
cat("Found data for", length(sample_ids), "unique sample IDs:", paste(sample_ids, collapse = ", "), "\n")

# Combine all data frames
combined_data <- bind_rows(all_data)

# Extract key peptides from all samples
key_peptide_data <- combined_data %>%
  filter(Peptide %in% key_peptides) %>%
  mutate(
    # Mark whether this is from the spiked sample
    is_spiked = ifelse(Sample_ID == "51S", "Spiked (10ng)", "Natural")
  )

# Create a summary of key peptide detection across samples
peptide_summary <- key_peptide_data %>%
  group_by(Sample_ID, Peptide, is_spiked) %>%
  summarize(
    detected = TRUE,
    spectral_count = sum(Spectral.Count),
    total_intensity = sum(Intensity),
    source_filename = first(Filename),
    .groups = "drop"
  ) %>%
  arrange(Peptide, Sample_ID)

# Create a wide format data frame for heatmap and Excel output
# First, get a complete matrix with all peptide-sample combinations
all_combinations <- expand.grid(
  Peptide = key_peptides,
  Sample_ID = sample_ids,
  stringsAsFactors = FALSE
)

# Then merge with actual data, filling in zeros for missing combinations
peptide_matrix_data <- all_combinations %>%
  left_join(peptide_summary, by = c("Peptide", "Sample_ID")) %>%
  mutate(
    detected = ifelse(is.na(detected), FALSE, detected),
    spectral_count = ifelse(is.na(spectral_count), 0, spectral_count),
    total_intensity = ifelse(is.na(total_intensity), 0, total_intensity),
    is_spiked = ifelse(is.na(is_spiked), "Not Detected", is_spiked)
  )

# Create a matrix for the heatmap
intensity_matrix <- peptide_matrix_data %>%
  select(Peptide, Sample_ID, total_intensity) %>%
  pivot_wider(
    names_from = Sample_ID,
    values_from = total_intensity,
    values_fill = 0
  ) %>%
  column_to_rownames("Peptide")

# Log transform intensity for better visualization (adding 1 to avoid log(0))
log_intensity_matrix <- log10(intensity_matrix + 1)

# Create a presence/absence matrix (1 if detected, 0 if not)
presence_matrix <- (intensity_matrix > 0) * 1

# Sort sample columns to keep 51 and 51S first for easy comparison
sample_order <- c(
  # Put 51 and 51S first
  sample_ids[grep("^51$", sample_ids)],
  sample_ids[grep("51S", sample_ids)],
  # Then add all other samples, sorted
  setdiff(sample_ids, c(sample_ids[grep("^51$", sample_ids)], sample_ids[grep("51S", sample_ids)]))
)

# Reorder matrices for visualization
intensity_matrix <- intensity_matrix[, sample_order]
log_intensity_matrix <- log_intensity_matrix[, sample_order]
presence_matrix <- presence_matrix[, sample_order]

# Create formatted Excel output with multiple sheets
# 1. Create a sheet with peptide detection summary
excel_summary <- peptide_matrix_data %>%
  mutate(
    Detection = ifelse(detected, "Yes", "No"),
    `MS1 Intensity` = total_intensity,
    `Spectral Count` = spectral_count,
    `Peptide Source` = is_spiked
  ) %>%
  select(
    Peptide, 
    Sample_ID, 
    Detection, 
    `MS1 Intensity`, 
    `Spectral Count`,
    `Peptide Source`
  ) %>%
  arrange(Peptide, Sample_ID)

# 2. Create a sheet with presence/absence in a wide format for easy overview
excel_presence <- peptide_matrix_data %>%
  mutate(present = ifelse(detected, "Yes", "No")) %>%
  select(Peptide, Sample_ID, present) %>%
  pivot_wider(
    names_from = Sample_ID,
    values_from = present,
    values_fill = "No"
  )

# 3. Create a sheet with intensity values in a wide format
excel_intensity <- peptide_matrix_data %>%
  select(Peptide, Sample_ID, total_intensity) %>%
  pivot_wider(
    names_from = Sample_ID,
    values_from = total_intensity,
    values_fill = 0
  )

# 4. Special comparison for samples 51 vs 51S (unspiked vs spiked)
if (all(c("51", "51S") %in% sample_ids)) {
  excel_comparison <- peptide_matrix_data %>%
    filter(Sample_ID %in% c("51", "51S")) %>%
    select(Peptide, Sample_ID, detected, total_intensity, spectral_count, is_spiked) %>%
    pivot_wider(
      names_from = Sample_ID,
      values_from = c(detected, total_intensity, spectral_count),
      values_fill = list(detected = FALSE, total_intensity = 0, spectral_count = 0)
    ) %>%
    mutate(
      detection_improvement = case_when(
        detected_51 == TRUE & detected_51S == TRUE ~ "Detected in both",
        detected_51 == FALSE & detected_51S == TRUE ~ "Only detected in spiked sample",
        detected_51 == TRUE & detected_51S == FALSE ~ "Only detected in unspiked sample",
        detected_51 == FALSE & detected_51S == FALSE ~ "Not detected in either"
      ),
      intensity_fold_change = ifelse(total_intensity_51 > 0, 
                                     total_intensity_51S / total_intensity_51,
                                     ifelse(total_intensity_51S > 0, Inf, NA)),
      intensity_increase = total_intensity_51S - total_intensity_51
    )
} else {
  excel_comparison <- data.frame(
    Peptide = key_peptides,
    Note = rep("Samples 51 and/or 51S not found for comparison", length(key_peptides))
  )
}

# Create a list of sheets for the Excel file
excel_sheets <- list(
  "Detection_Summary" = excel_summary,
  "Presence_Absence" = excel_presence,
  "Intensity_Values" = excel_intensity,
  "51_vs_51S_Comparison" = excel_comparison
)

# Write Excel file with multiple sheets
write_xlsx(excel_sheets, path = file.path(data_path, "DNAJB1_PRKACA_Junction_Peptides.xlsx"))

# Create visualizations
# 1. Heatmap of peptide presence/absence across samples
pdf(file.path(data_path, "junction_peptide_presence_heatmap.pdf"), width = 10, height = 8)
pheatmap(
  presence_matrix,
  main = "Presence of DNAJB1_PRKACA Junction Peptides Across Samples",
  color = c("white", "darkblue"),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  display_numbers = TRUE,
  number_format = "%d"
)
dev.off()

# 2. Heatmap of peptide intensities across samples
pdf(file.path(data_path, "junction_peptide_intensity_heatmap.pdf"), width = 10, height = 8)
pheatmap(
  log_intensity_matrix,
  main = "Intensity of DNAJB1_PRKACA Junction Peptides Across Samples (log10)",
  color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  display_numbers = TRUE,
  number_format = "%.1f"
)
dev.off()

# 3. Bar plot comparing 51 vs 51S for each peptide
if (all(c("51", "51S") %in% sample_ids)) {
  comparison_plot_data <- peptide_matrix_data %>%
    filter(Sample_ID %in% c("51", "51S")) %>%
    mutate(Sample_Type = ifelse(Sample_ID == "51S", "51S (Spiked 10ng)", "51 (Natural)"))
  
  pdf(file.path(data_path, "51_vs_51S_comparison.pdf"), width = 10, height = 8)
  comparison_plot <- ggplot(comparison_plot_data, aes(x = Peptide, y = total_intensity, fill = Sample_Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Comparison of Junction Peptide Detection: 51 vs 51S",
      subtitle = "51S sample has 10ng of each peptide spiked in",
      x = "Peptide Sequence",
      y = "Total MS1 Intensity",
      fill = "Sample"
    )
  
  print(comparison_plot)  # Make sure to print the plot
  dev.off()
  
  # As a backup, also save as PNG format
  png(file.path(data_path, "51_vs_51S_comparison.png"), width = 800, height = 600)
  print(comparison_plot)
  dev.off()
}

# Print summary information
cat("\nSummary of key junction peptides detected across samples:\n")
for (peptide in key_peptides) {
  detected_in <- peptide_matrix_data %>% 
    filter(Peptide == peptide, detected == TRUE) %>% 
    pull(Sample_ID)
  
  if (length(detected_in) > 0) {
    cat(paste0(peptide, ": Detected in ", length(detected_in), " samples (", 
               paste(detected_in, collapse = ", "), ")\n"))
  } else {
    cat(paste0(peptide, ": Not detected in any sample\n"))
  }
}

cat("\nAnalysis complete! Results saved to:", data_path, "\n")
cat("\nThe following files were generated:\n")
cat("1. DNAJB1_PRKACA_Junction_Peptides.xlsx - Excel file with detection summary, presence/absence matrix, intensity values, and 51 vs 51S comparison\n")
cat("2. junction_peptide_presence_heatmap.pdf - Heatmap showing presence/absence of key peptides\n")
cat("3. junction_peptide_intensity_heatmap.pdf - Heatmap showing intensities of key peptides\n")
if (all(c("51", "51S") %in% sample_ids)) {
  cat("4. 51_vs_51S_comparison.pdf - Bar plot comparing detection between 51 and 51S samples\n")
}