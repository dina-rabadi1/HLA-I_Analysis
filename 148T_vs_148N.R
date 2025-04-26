# Script to analyze peptide differences between 148T (tumor) and 148N (normal) samples
# focusing on fold change analysis

# Setting directory
setwd("~/Documents/Github/HLA-I_Analysis/")
#make data file in HLA-I_Analysis then save it

# Load required packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(writexl)  # For Excel output
library(openxlsx) # For better Excel formatting

# Define the path to your data files
data_path <- "/Users/dinarabadi/Desktop/2704DR"

# Read all TSV files in the directory and extract sample IDs from filenames
files <- list.files(path = data_path, pattern = ".*_peptides\\.tsv$|.*_untargeted_peptide\\.tsv$", full.names = TRUE)

# Filter for only 148T and 148N files
tumor_normal_files <- files[grepl("148[TN]", files)]

if (length(tumor_normal_files) == 0) {
  stop("No 148T or 148N peptide files found in ", data_path)
}

cat("Found", length(tumor_normal_files), "tumor/normal files for analysis\n")

# Read the 148T and 148N files
all_data <- list()
sample_ids <- c()

for (file in tumor_normal_files) {
  filename <- basename(file)
  
  # Extract sample ID from filename (148T or 148N)
  sample_id <- ifelse(grepl("148T", filename), "148T", "148N")
  
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

# Get unique sample IDs (should be 148T and 148N)
sample_ids <- unique(sample_ids)
cat("Found data for", length(sample_ids), "unique sample IDs:", paste(sample_ids, collapse = ", "), "\n")

# Ensure we have both tumor and normal samples
if (!all(c("148T", "148N") %in% sample_ids)) {
  stop("Missing either tumor (148T) or normal (148N) sample data")
}

# Combine all data frames
combined_data <- bind_rows(all_data)

# Create a summary of peptide detection for tumor and normal
peptide_summary <- combined_data %>%
  group_by(Sample_ID, Peptide) %>%
  summarize(
    peptide_length = first(nchar(Peptide)),
    spectral_count = sum(Spectral.Count),
    total_intensity = sum(Intensity),
    protein_ids = paste(unique(Protein.ID), collapse = "; "),
    genes = paste(unique(Gene), collapse = "; "),
    source_filename = first(Filename),
    .groups = "drop"
  ) %>%
  arrange(Peptide, Sample_ID)

# Create a wide format table with tumor and normal side by side
peptide_comparison <- peptide_summary %>%
  select(Sample_ID, Peptide, peptide_length, spectral_count, total_intensity, protein_ids, genes) %>%
  pivot_wider(
    names_from = Sample_ID,
    values_from = c(spectral_count, total_intensity, protein_ids, genes),
    values_fill = list(spectral_count = 0, total_intensity = 0, protein_ids = NA, genes = NA)
  )

# Calculate fold changes and identify tumor-specific and normal-specific peptides
peptide_analysis <- peptide_comparison %>%
  mutate(
    # Replace zero with small value to prevent division by zero or Inf
    total_intensity_148N_adj = ifelse(total_intensity_148N == 0, 0.1, total_intensity_148N),
    total_intensity_148T_adj = ifelse(total_intensity_148T == 0, 0.1, total_intensity_148T),
    
    # Calculate fold changes (log2)
    log2_fold_change = log2(total_intensity_148T_adj / total_intensity_148N_adj),
    
    # Determine if peptide is specific to tumor or normal
    detection_status = case_when(
      total_intensity_148T > 0 & total_intensity_148N == 0 ~ "Tumor-specific",
      total_intensity_148N > 0 & total_intensity_148T == 0 ~ "Normal-specific",
      total_intensity_148T > 0 & total_intensity_148N > 0 ~ "Detected in both",
      TRUE ~ "Not detected"
    ),
    
    # Add peptide length
    peptide_length = nchar(Peptide),
    
    # Simplified category for plotting
    peptide_category = case_when(
      log2_fold_change > 1 ~ "Up in Tumor (FC > 2)",
      log2_fold_change < -1 ~ "Down in Tumor (FC < 0.5)",
      TRUE ~ "Similar (-1 < log2FC < 1)"
    )
  ) %>%
  # Clean up protein and gene info
  mutate(
    genes_combined = coalesce(genes_148T, genes_148N),
    proteins_combined = coalesce(protein_ids_148T, protein_ids_148N)
  ) %>%
  # Sort by fold change for easier viewing
  arrange(desc(log2_fold_change))

# Identify DNAJB1-PRKACA fusion peptides
# Define the sequences of DNAJB1 and PRKACA parts of the fusion protein
dnajb1_seq <- "GKDYYQTLGLARGASDEEIKRAYRRQALRYHPDKNKEPGAEEKFKEIAEAYDVLSDPRKREIFDRYGEE"
prkaca_seq <- "VKEFLAKAKEDFLKKWESPAQNTAHLDQFERIKTLGTGSFGRVMLVKHKETGNHYAMKILDKQKVVKLKQIEHTLNEKRILQAVNFPFLVKLEFSFKDNSNLYMVMEYVPGGEMFSHLRRIGRFSEPHARFYAAQIVLTFEYLHSLDLIYRDLKPENLLIDQQGYIQVTDFGFAKRVKGRTWTLCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFATTDWIAIYQRKVEAPFIPKFKGPGDTSNFDDYEEEEIRVSINEKCGKEFSEF"
fusion_protein <- paste0(dnajb1_seq, prkaca_seq)

# Function to check if a peptide spans the fusion junction
is_fusion_junction_peptide <- function(peptide_seq) {
  # Define the end of DNAJB1 and start of PRKACA for the fusion
  dnajb1_end <- "IFDRYGEE"
  prkaca_start <- "VKEFLAK"
  
  # Check if the peptide spans the fusion junction
  spans_junction <- FALSE
  
  if (nchar(peptide_seq) >= 6) { # Require at least 6 AA to be meaningful
    for (i in 3:(nchar(peptide_seq) - 3)) { # At least 3 AA from each side
      left_part <- substr(peptide_seq, 1, i)
      right_part <- substr(peptide_seq, i + 1, nchar(peptide_seq))
      
      # Check if left part is in DNAJB1 (at the end) and right part in PRKACA (at the beginning)
      if (grepl(left_part, dnajb1_seq, fixed = TRUE) && 
          grepl(right_part, prkaca_seq, fixed = TRUE)) {
        
        # Additional check to ensure left part aligns with end of DNAJB1
        left_pos <- gregexpr(left_part, dnajb1_seq, fixed = TRUE)[[1]]
        if (length(left_pos) > 0 && any(left_pos + nchar(left_part) - 1 >= nchar(dnajb1_seq) - 10)) {
          
          # Additional check to ensure right part aligns with start of PRKACA
          right_pos <- gregexpr(right_part, prkaca_seq, fixed = TRUE)[[1]]
          if (length(right_pos) > 0 && any(right_pos <= 10)) {
            spans_junction <- TRUE
            break
          }
        }
      }
    }
  }
  
  # Check if peptide is from either part of the fusion protein
  from_dnajb1 <- grepl(peptide_seq, dnajb1_seq, fixed = TRUE)
  from_prkaca <- grepl(peptide_seq, prkaca_seq, fixed = TRUE)
  
  return(list(
    spans_junction = spans_junction,
    from_dnajb1 = from_dnajb1,
    from_prkaca = from_prkaca,
    from_fusion = from_dnajb1 | from_prkaca | spans_junction
  ))
}

# Add fusion protein information to the peptide analysis
peptide_analysis <- peptide_analysis %>%
  rowwise() %>%
  mutate(
    fusion_info = list(is_fusion_junction_peptide(Peptide)),
    from_fusion = fusion_info$from_fusion,
    spans_junction = fusion_info$spans_junction,
    from_dnajb1 = fusion_info$from_dnajb1,
    from_prkaca = fusion_info$from_prkaca,
    fusion_peptide_type = case_when(
      spans_junction ~ "Junction-spanning",
      from_dnajb1 ~ "DNAJB1 part",
      from_prkaca ~ "PRKACA part",
      TRUE ~ "Not from fusion"
    )
  ) %>%
  select(-fusion_info)

# Create a separate analysis specifically for fusion peptides
fusion_peptides_analysis <- peptide_analysis %>%
  filter(from_fusion) %>%
  arrange(desc(spans_junction), desc(log2_fold_change))

# Create Excel output with multiple sheets
excel_sheets <- list(
  "All_Peptides_Analysis" = peptide_analysis %>% 
    select(-total_intensity_148N_adj, -total_intensity_148T_adj),
  "Fusion_Peptides" = fusion_peptides_analysis %>% 
    select(-total_intensity_148N_adj, -total_intensity_148T_adj),
  "Tumor_Specific_Peptides" = peptide_analysis %>% 
    filter(detection_status == "Tumor-specific") %>% 
    select(-total_intensity_148N_adj, -total_intensity_148T_adj),
  "Normal_Specific_Peptides" = peptide_analysis %>% 
    filter(detection_status == "Normal-specific") %>% 
    select(-total_intensity_148N_adj, -total_intensity_148T_adj),
  "Up_In_Tumor" = peptide_analysis %>% 
    filter(log2_fold_change > 1) %>% 
    arrange(desc(log2_fold_change)) %>% 
    select(-total_intensity_148N_adj, -total_intensity_148T_adj),
  "Down_In_Tumor" = peptide_analysis %>% 
    filter(log2_fold_change < -1) %>% 
    arrange(log2_fold_change) %>% 
    select(-total_intensity_148N_adj, -total_intensity_148T_adj)
)

# Create a more formatted Excel workbook
wb <- createWorkbook()

# Add sheets with formatting
for (sheet_name in names(excel_sheets)) {
  # Add a worksheet
  addWorksheet(wb, sheet_name)
  
  # Write data
  writeData(wb, sheet_name, excel_sheets[[sheet_name]], headerStyle = createStyle(textDecoration = "bold"))
  
  # Auto-adjust column widths
  setColWidths(wb, sheet_name, cols = 1:ncol(excel_sheets[[sheet_name]]), widths = "auto")
  
  # Freeze the header row
  freezePane(wb, sheet_name, firstRow = TRUE)
}

# Save the Excel file
saveWorkbook(wb, file.path(data_path, "148T_vs_148N_Peptide_Analysis.xlsx"), overwrite = TRUE)

# Write simple Excel file as backup
write_xlsx(excel_sheets, path = file.path(data_path, "148T_vs_148N_Peptide_Analysis_simple.xlsx"))

# Create visualizations

# 1. Volcano plot of all peptides
volcano_data <- peptide_analysis %>%
  filter(!is.na(log2_fold_change))

pdf(file.path(data_path, "148T_vs_148N_volcano_plot.pdf"), width = 10, height = 8)
volcano_plot <- ggplot(volcano_data, aes(x = log2_fold_change, y = -log10(0.05), color = peptide_category)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +
  scale_color_manual(values = c("Up in Tumor (FC > 2)" = "red", 
                                "Down in Tumor (FC < 0.5)" = "blue", 
                                "Similar (-1 < log2FC < 1)" = "gray")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Peptides in 148T vs 148N",
    subtitle = "Red: Upregulated in Tumor, Blue: Downregulated in Tumor",
    x = "Log2 Fold Change (Tumor/Normal)",
    y = "-log10(p-value) [placeholder]",
    color = "Peptide Category"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )
print(volcano_plot)
dev.off()

png(file.path(data_path, "148T_vs_148N_volcano_plot.png"), width = 800, height = 600, res = 100)
print(volcano_plot)
dev.off()

# 2. Histogram of fold changes
pdf(file.path(data_path, "148T_vs_148N_fold_change_histogram.pdf"), width = 10, height = 6)
hist_plot <- ggplot(volcano_data, aes(x = log2_fold_change, fill = peptide_category)) +
  geom_histogram(bins = 50, color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("Up in Tumor (FC > 2)" = "red", 
                               "Down in Tumor (FC < 0.5)" = "blue", 
                               "Similar (-1 < log2FC < 1)" = "gray")) +
  theme_minimal() +
  labs(
    title = "Distribution of Peptide Fold Changes in 148T vs 148N",
    x = "Log2 Fold Change (Tumor/Normal)",
    y = "Count",
    fill = "Peptide Category"
  )
print(hist_plot)
dev.off()

png(file.path(data_path, "148T_vs_148N_fold_change_histogram.png"), width = 800, height = 600, res = 100)
print(hist_plot)
dev.off()

# 3. Barplot of detection status
detection_summary <- peptide_analysis %>%
  group_by(detection_status) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(count))

pdf(file.path(data_path, "148T_vs_148N_detection_status.pdf"), width = 8, height = 6)
detection_plot <- ggplot(detection_summary, aes(x = detection_status, y = count, fill = detection_status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = -0.5) +
  theme_minimal() +
  labs(
    title = "Peptide Detection Status in 148T vs 148N",
    x = "Detection Status",
    y = "Count",
    fill = "Status"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(detection_plot)
dev.off()

png(file.path(data_path, "148T_vs_148N_detection_status.png"), width = 800, height = 600, res = 100)
print(detection_plot)
dev.off()

# 4. Heatmap of fusion peptides (if any found)
if (nrow(fusion_peptides_analysis) > 0) {
  # Create a matrix for the heatmap
  fusion_intensity_data <- fusion_peptides_analysis %>%
    select(Peptide, total_intensity_148T, total_intensity_148N, fusion_peptide_type, spans_junction) %>%
    pivot_longer(
      cols = c(total_intensity_148T, total_intensity_148N),
      names_to = "Sample",
      values_to = "Intensity"
    ) %>%
    mutate(Sample = gsub("total_intensity_", "", Sample)) %>%
    pivot_wider(
      names_from = Sample,
      values_from = Intensity
    ) %>%
    arrange(desc(spans_junction), Peptide)
  
  # Log transform the values
  fusion_matrix <- as.matrix(fusion_intensity_data[, c("148T", "148N")])
  rownames(fusion_matrix) <- fusion_intensity_data$Peptide
  log_fusion_matrix <- log10(fusion_matrix + 1)
  
  # Create annotation for the rows
  row_annotation <- data.frame(
    Peptide_Type = fusion_intensity_data$fusion_peptide_type,
    row.names = fusion_intensity_data$Peptide
  )
  
  # Create a PDF of the heatmap
  pdf(file.path(data_path, "148T_vs_148N_fusion_peptides_heatmap.pdf"), width = 10, height = max(8, nrow(fusion_matrix)/3))
  pheatmap(
    log_fusion_matrix,
    main = "Intensity of DNAJB1_PRKACA Fusion Peptides in 148T vs 148N (log10)",
    color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    display_numbers = TRUE,
    number_format = "%.1f",
    fontsize_row = 10,
    fontsize_col = 10
  )
  dev.off()
  
  # Create a PNG of the heatmap
  png(file.path(data_path, "148T_vs_148N_fusion_peptides_heatmap.png"), width = 800, height = max(600, nrow(fusion_matrix)*40), res = 100)
  pheatmap(
    log_fusion_matrix,
    main = "Intensity of DNAJB1_PRKACA Fusion Peptides in 148T vs 148N (log10)",
    color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    display_numbers = TRUE,
    number_format = "%.1f",
    fontsize_row = 10,
    fontsize_col = 10
  )
  dev.off()
}

# Print summary information
cat("\nSummary of 148T vs 148N peptide analysis:\n")
cat("Total peptides analyzed:", nrow(peptide_analysis), "\n")
cat("\nDetection status:\n")
print(detection_summary)

cat("\nFold change categories:\n")
peptide_analysis %>%
  group_by(peptide_category) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count)) %>%
  print(n = Inf)

if (nrow(fusion_peptides_analysis) > 0) {
  cat("\nFusion peptides found:", nrow(fusion_peptides_analysis), "\n")
  cat("Junction-spanning peptides:", sum(fusion_peptides_analysis$spans_junction), "\n")
  
  cat("\nFusion peptide detection summary:\n")
  fusion_peptides_analysis %>%
    group_by(fusion_peptide_type, detection_status) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(fusion_peptide_type, detection_status) %>%
    print(n = Inf)
} else {
  cat("\nNo fusion peptides found in the analysis.\n")
}

cat("\nAnalysis complete! Results saved to:", data_path, "\n")
cat("\nThe following files were generated:\n")
cat("1. 148T_vs_148N_Peptide_Analysis.xlsx - Excel file with comprehensive analysis results\n")
cat("2. 148T_vs_148N_volcano_plot.pdf/png - Volcano plot showing peptide fold changes\n")
cat("3. 148T_vs_148N_fold_change_histogram.pdf/png - Histogram of fold changes\n")
cat("4. 148T_vs_148N_detection_status.pdf/png - Barplot of peptide detection status\n")
if (nrow(fusion_peptides_analysis) > 0) {
  cat("5. 148T_vs_148N_fusion_peptides_heatmap.pdf/png - Heatmap of fusion peptides\n")
}