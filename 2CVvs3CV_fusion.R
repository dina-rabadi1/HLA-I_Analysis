# Script to analyze DNAJB1-PRKACA fusion protein peptides in 2CV vs 3CV samples
# Looking for 8-12mers that span the fusion junction

# Setting directory
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/")
#make data file in HLA-I_Analysis then save it

# Load required packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(writexl)
library(stringr)
# Try to load ggrepel - install first if not available
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggrepel)

# Define the path to your data files
data_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/unique_peptides_unmodified.tsv"

# Define the fusion protein sequence - UPDATED based on provided sequences
fusion_protein <- "RKREIFDRYGEEVKEFLAKAKEDF"

# Identify the fusion junction - between E (DNAJB1) and V (PRKACA)
# FIXED: Updated to match the provided sequences
dnajb1_part <- "RKREIFDRYGEE"  # Ends with YGEE (full DNAJB1 part)
prkaca_part <- "VKEFLAKAKEDF"  # Starts with VKEF (full PRKACA part)
junction_position <- nchar(dnajb1_part)  # Position after the last E in YGEE

cat("DNAJB1 part:", dnajb1_part, "\n")
cat("PRKACA part:", prkaca_part, "\n")
cat("Junction position:", junction_position, "\n")
cat("Fusion protein:", fusion_protein, "\n")
cat("Checking junction: character at junction position is", substr(fusion_protein, junction_position, junction_position), "\n")
cat("Checking junction: next character is", substr(fusion_protein, junction_position+1, junction_position+1), "\n\n")

# Generate theoretical junction-spanning peptides (8-12mers)
theoretical_peptides <- list()
peptide_lengths <- 8:12  # Looking for 8-12mers (FIXED: explicitly excluding 7mers)

for (length in peptide_lengths) {
  for (start_pos in 1:(nchar(fusion_protein) - length + 1)) {
    peptide <- substr(fusion_protein, start_pos, start_pos + length - 1)
    
    # Check if this peptide spans the junction
    # It spans if it includes at least one residue from both proteins
    peptide_end_pos <- start_pos + length - 1
    spans_junction <- start_pos <= junction_position && peptide_end_pos > junction_position
    
    if (spans_junction) {
      theoretical_peptides[[length(theoretical_peptides) + 1]] <- list(
        Peptide = peptide,
        Length = length,
        Start_Position = start_pos,
        End_Position = peptide_end_pos,
        DNAJB1_Part = paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
        PRKACA_Part = paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = ""),
        Visualization = paste0(
          paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
          paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = "")
        )
      )
    }
  }
}

# Convert to dataframe
theoretical_junction_peptides <- do.call(rbind, lapply(theoretical_peptides, function(x) {
  data.frame(
    Peptide = x$Peptide,
    Length = x$Length,
    Start_Position = x$Start_Position,
    End_Position = x$End_Position,
    DNAJB1_Part = x$DNAJB1_Part,
    PRKACA_Part = x$PRKACA_Part,
    Visualization = x$Visualization,
    stringsAsFactors = FALSE
  )
}))

cat("Generated", nrow(theoretical_junction_peptides), "theoretical junction-spanning peptides\n")

# Count by length
for (length in peptide_lengths) {
  count <- sum(theoretical_junction_peptides$Length == length)
  cat("  Length", length, ":", count, "peptides\n")
}

# Print all theoretical peptides in a neat table
cat("\nAll theoretical junction-spanning peptides:\n")
for (i in 1:nrow(theoretical_junction_peptides)) {
  peptide <- theoretical_junction_peptides[i,]
  cat(sprintf("%2d. %s (Length: %d, Start: %d, End: %d, Vis: %s)\n", 
              i, peptide$Peptide, peptide$Length, peptide$Start_Position, 
              peptide$End_Position, peptide$Visualization))
}

# Read all TSV files in the directory and extract sample IDs and CV type
files <- list.files(path = data_path, pattern = ".*_peptides\\.tsv$|.*_untargeted_peptide\\.tsv$", full.names = TRUE)

if (length(files) == 0) {
  stop("No peptide TSV files found in ", data_path)
}

cat("Found", length(files), "peptide files for analysis\n")

# Extract sample IDs and CV type from filenames
all_data <- list()
sample_info <- data.frame(
  file = character(),
  sample_id = character(),
  cv_type = character(),
  stringsAsFactors = FALSE
)

for (file in files) {
  filename <- basename(file)
  
  # Extract CV type from filename
  cv_type <- str_extract(filename, "2CV|3CV")
  
  # Skip files that don't match either 2CV or 3CV
  if (is.na(cv_type)) {
    cat("Skipping file", filename, "- CV type not found in filename\n")
    next
  }
  
  # Extract sample ID from filename
  # Pattern matches the sample ID after DDA_2CV_ or DDA_3CV_
  sample_id <- gsub(".*DDA_[23]CV_([^_]+)_.*\\.tsv$", "\\1", filename)
  
  # For filenames with extra parts like "_01", clean up the sample ID
  if (grepl("_\\d+_peptides\\.tsv$", filename)) {
    sample_id <- gsub("(.*)_\\d+$", "\\1", sample_id)
  }
  
  cat("Reading file:", filename, "- Sample ID:", sample_id, "- CV Type:", cv_type, "\n")
  
  # Read the file
  data <- read.delim(file, stringsAsFactors = FALSE)
  
  # Add columns for sample ID and CV type
  data$Sample_ID <- sample_id
  data$CV_Type <- cv_type
  data$Filename <- filename
  
  # Add to our list
  all_data[[length(all_data) + 1]] <- data
  
  # Add to our sample info dataframe
  sample_info <- rbind(sample_info, data.frame(
    file = file,
    sample_id = sample_id,
    cv_type = cv_type,
    stringsAsFactors = FALSE
  ))
}

# Combine all data frames
combined_data <- bind_rows(all_data)

# Identify samples that have both 2CV and 3CV data for direct comparison
paired_samples <- sample_info %>%
  group_by(sample_id) %>%
  summarize(cv_count = n_distinct(cv_type)) %>%
  filter(cv_count == 2) %>%
  pull(sample_id)

cat("Found", length(paired_samples), "samples with both 2CV and 3CV data for direct comparison:", 
    paste(paired_samples, collapse = ", "), "\n")

# Create output directory for results
results_dir <- file.path(data_path, "Fusion_Junction_Analysis")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

#===============================#
# Identify all detected peptides that span the fusion junction
#===============================#

# Function to check if a peptide spans the fusion junction
is_junction_spanning <- function(peptide) {
  # First check if it's entirely within the fusion protein
  if (!grepl(peptide, fusion_protein, fixed = TRUE)) {
    return(FALSE)
  }
  
  # Find position in fusion protein
  start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
  if (start_pos == -1) {
    return(FALSE)
  }
  
  end_pos <- start_pos + nchar(peptide) - 1
  
  # Check if it spans the junction
  spans_junction <- start_pos <= junction_position && end_pos > junction_position
  
  return(spans_junction)
}

# Apply to all peptides in the dataset 
# FIXED: Added explicit filter to ensure we only include 8-12mers
combined_data <- combined_data %>%
  mutate(
    Spans_Junction = sapply(Peptide, is_junction_spanning),
    Peptide_Length = nchar(Peptide)
  ) %>%
  filter(Peptide_Length >= 8 & Peptide_Length <= 12)  # Only include 8-12mers

# Extract detected fusion-spanning peptides
detected_junction_peptides <- combined_data %>%
  filter(Spans_Junction == TRUE) %>%
  group_by(Peptide, Peptide_Length) %>%
  summarize(
    Total_Detected = n(),
    Detected_2CV = sum(CV_Type == "2CV"),
    Detected_3CV = sum(CV_Type == "3CV"),
    Avg_Intensity_2CV = mean(Intensity[CV_Type == "2CV"]),
    Avg_Intensity_3CV = mean(Intensity[CV_Type == "3CV"]),
    .groups = "drop"
  ) %>%
  arrange(Peptide_Length, Peptide)

cat("Found", nrow(detected_junction_peptides), "detected peptides that span the fusion junction\n")

if (nrow(detected_junction_peptides) > 0) {
  # Create a detailed analysis of each detected junction peptide
  detected_details <- data.frame()
  
  for (i in 1:nrow(detected_junction_peptides)) {
    peptide <- detected_junction_peptides$Peptide[i]
    start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
    end_pos <- start_pos + nchar(peptide) - 1
    
    # Calculate how many residues come from each protein
    dnajb1_residues <- max(0, min(junction_position - start_pos + 1, nchar(peptide)))
    prkaca_residues <- max(0, min(end_pos - junction_position, nchar(peptide)))
    
    # Visualization string (D for DNAJB1, P for PRKACA)
    vis_string <- paste0(
      paste(rep("D", dnajb1_residues), collapse = ""),
      paste(rep("P", prkaca_residues), collapse = "")
    )
    
    # Add to details dataframe
    detected_details <- rbind(detected_details, data.frame(
      Peptide = peptide,
      Length = nchar(peptide),
      Start_Position = start_pos,
      End_Position = end_pos,
      DNAJB1_Residues = dnajb1_residues,
      PRKACA_Residues = prkaca_residues,
      Visualization = vis_string,
      Total_Detected = detected_junction_peptides$Total_Detected[i],
      Detected_2CV = detected_junction_peptides$Detected_2CV[i],
      Detected_3CV = detected_junction_peptides$Detected_3CV[i],
      Avg_Intensity_2CV = detected_junction_peptides$Avg_Intensity_2CV[i],
      Avg_Intensity_3CV = detected_junction_peptides$Avg_Intensity_3CV[i],
      stringsAsFactors = FALSE
    ))
  }
  
  # Create a detailed heatmap for each detected junction peptide
  # showing detection across samples and CV types
  for (peptide in detected_details$Peptide) {
    # First create detected column with mutate, then select
    peptide_data <- combined_data %>%
      filter(Peptide == peptide) %>%
      mutate(detected = TRUE) %>%
      select(Sample_ID, CV_Type, detected, Intensity)
    
    # Create a complete matrix with all sample-CV combinations
    all_combinations <- expand.grid(
      Sample_ID = unique(combined_data$Sample_ID),
      CV_Type = c("2CV", "3CV"),
      stringsAsFactors = FALSE
    )
    
    peptide_matrix_data <- all_combinations %>%
      left_join(peptide_data, by = c("Sample_ID", "CV_Type")) %>%
      mutate(
        detected = ifelse(is.na(detected), FALSE, detected),
        Intensity = ifelse(is.na(Intensity), 0, Intensity)
      )
    
    # Create unique row identifiers to avoid duplicate row names
    # Create a matrix for presence/absence heatmap
    presence_matrix <- peptide_matrix_data %>%
      mutate(Row_ID = paste0(CV_Type, "_row")) %>%
      select(Row_ID, Sample_ID, present = detected) %>%
      pivot_wider(
        names_from = Sample_ID,
        values_from = present,
        values_fill = FALSE
      ) %>%
      column_to_rownames("Row_ID")
    
    # Create a matrix for intensity heatmap
    intensity_matrix <- peptide_matrix_data %>%
      mutate(Row_ID = paste0(CV_Type, "_row")) %>%
      select(Row_ID, Sample_ID, Intensity) %>%
      pivot_wider(
        names_from = Sample_ID,
        values_from = Intensity,
        values_fill = 0
      ) %>%
      column_to_rownames("Row_ID")
    
    # Convert logical presence matrix to numeric for pheatmap
    # FIX: Properly convert data frame to matrix first, then convert logical to numeric
    presence_matrix_mat <- as.matrix(presence_matrix)
    presence_matrix_numeric <- matrix(as.numeric(presence_matrix_mat), 
                                      nrow = nrow(presence_matrix_mat),
                                      dimnames = dimnames(presence_matrix_mat))
    
    # Log transform intensities
    log_intensity_matrix <- log10(intensity_matrix + 1)
    
    # Create presence/absence heatmap
    pdf(file.path(results_dir, paste0("peptide_", peptide, "_presence_heatmap.pdf")), width = 12, height = 4)
    pheatmap(
      presence_matrix_numeric,
      main = paste0("Presence of Fusion Peptide ", peptide, " Across Samples"),
      color = c("white", "steelblue"),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      # FIX: Use "%.0f" instead of "%d" for numeric format
      number_format = "%.0f"
    )
    dev.off()
    
    # Create intensity heatmap
    pdf(file.path(results_dir, paste0("peptide_", peptide, "_intensity_heatmap.pdf")), width = 12, height = 4)
    pheatmap(
      log_intensity_matrix,
      main = paste0("Intensity of Fusion Peptide ", peptide, " Across Samples (log10)"),
      color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      number_format = "%.1f"
    )
    dev.off()
  }
  
  # Compare 2CV vs 3CV detection efficiency for junction peptides
  cv_comparison_plot <- ggplot(detected_details, aes(x = Peptide, y = Length)) +
    geom_point(aes(size = Total_Detected, color = Detected_3CV / (Detected_2CV + 0.001))) +
    scale_color_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red", 
      midpoint = 1,
      name = "3CV/2CV\nDetection Ratio"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Fusion Junction Peptide Detection: 2CV vs 3CV",
      x = "Peptide Sequence",
      y = "Peptide Length",
      size = "Total\nDetections"
    )
  
  ggsave(file.path(results_dir, "junction_peptide_cv_comparison.pdf"), cv_comparison_plot, width = 12, height = 8)
  ggsave(file.path(results_dir, "junction_peptide_cv_comparison.png"), cv_comparison_plot, width = 12, height = 8)
  
  # Create a visualization showing which detected peptides span which parts of the fusion
  peptide_coverage_plot <- ggplot(detected_details, 
                                  aes(x = Start_Position, xend = End_Position, 
                                      y = reorder(Peptide, Start_Position), yend = reorder(Peptide, Start_Position))) +
    geom_segment(aes(color = Total_Detected), size = 5) +
    geom_vline(xintercept = junction_position, linetype = "dashed", color = "red") +
    annotate("text", x = junction_position - 2, y = nrow(detected_details) + 1, 
             label = "DNAJB1", color = "darkgreen", fontface = "bold") +
    annotate("text", x = junction_position + 2, y = nrow(detected_details) + 1, 
             label = "PRKACA", color = "purple", fontface = "bold") +
    scale_color_gradient(low = "lightsteelblue", high = "steelblue") +
    theme_minimal() +
    labs(
      title = "Coverage of Fusion Junction by Detected Peptides",
      x = "Position in Fusion Protein",
      y = "Peptide",
      color = "Total\nDetections"
    )
  
  ggsave(file.path(results_dir, "junction_peptide_coverage.pdf"), peptide_coverage_plot, width = 12, height = 8)
  ggsave(file.path(results_dir, "junction_peptide_coverage.png"), peptide_coverage_plot, width = 12, height = 8)
  
  # Create a visualization comparing intensities between 2CV and 3CV
  # FIX: Handle missing or NA values in the intensity data
  intensity_data <- detected_details %>%
    filter(!is.na(Avg_Intensity_2CV) & !is.na(Avg_Intensity_3CV) & 
             Avg_Intensity_2CV > 0 & Avg_Intensity_3CV > 0)
  
  if(nrow(intensity_data) > 0) {
    intensity_comparison_plot <- ggplot(intensity_data, 
                                        aes(x = Avg_Intensity_2CV + 1, y = Avg_Intensity_3CV + 1)) +
      geom_point(aes(size = Total_Detected, color = Length)) +
      geom_text_repel(aes(label = Peptide), size = 3) +
      scale_x_log10() +
      scale_y_log10() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      theme_minimal() +
      labs(
        title = "Average Intensity of Junction Peptides: 2CV vs 3CV",
        x = "Average Intensity in 2CV (log scale)",
        y = "Average Intensity in 3CV (log scale)",
        color = "Peptide\nLength",
        size = "Total\nDetections"
      )
    
    ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.pdf"), intensity_comparison_plot, width = 10, height = 8)
    ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.png"), intensity_comparison_plot, width = 10, height = 8)
  }
  
  #===============================#
  # Create combined heatmap of all detected junction peptides
  #===============================#
  
  # Create matrices for all detected junction peptides
  all_presence_2cv <- matrix(0, nrow = nrow(detected_details), 
                             ncol = length(unique(combined_data$Sample_ID[combined_data$CV_Type == "2CV"])))
  rownames(all_presence_2cv) <- detected_details$Peptide
  colnames(all_presence_2cv) <- unique(combined_data$Sample_ID[combined_data$CV_Type == "2CV"])
  
  all_presence_3cv <- matrix(0, nrow = nrow(detected_details), 
                             ncol = length(unique(combined_data$Sample_ID[combined_data$CV_Type == "3CV"])))
  rownames(all_presence_3cv) <- detected_details$Peptide
  colnames(all_presence_3cv) <- unique(combined_data$Sample_ID[combined_data$CV_Type == "3CV"])
  
  all_intensity_2cv <- matrix(0, nrow = nrow(detected_details), 
                              ncol = length(unique(combined_data$Sample_ID[combined_data$CV_Type == "2CV"])))
  rownames(all_intensity_2cv) <- detected_details$Peptide
  colnames(all_intensity_2cv) <- unique(combined_data$Sample_ID[combined_data$CV_Type == "2CV"])
  
  all_intensity_3cv <- matrix(0, nrow = nrow(detected_details), 
                              ncol = length(unique(combined_data$Sample_ID[combined_data$CV_Type == "3CV"])))
  rownames(all_intensity_3cv) <- detected_details$Peptide
  colnames(all_intensity_3cv) <- unique(combined_data$Sample_ID[combined_data$CV_Type == "3CV"])
  
  # Fill matrices with peptide detection and intensity data
  for (i in 1:nrow(detected_details)) {
    peptide <- detected_details$Peptide[i]
    
    peptide_data_2cv <- combined_data %>%
      filter(Peptide == peptide, CV_Type == "2CV")
    
    peptide_data_3cv <- combined_data %>%
      filter(Peptide == peptide, CV_Type == "3CV")
    
    # Fill 2CV matrices
    for (sample in peptide_data_2cv$Sample_ID) {
      if (sample %in% colnames(all_presence_2cv)) {
        all_presence_2cv[peptide, sample] <- 1
        all_intensity_2cv[peptide, sample] <- peptide_data_2cv$Intensity[peptide_data_2cv$Sample_ID == sample]
      }
    }
    
    # Fill 3CV matrices
    for (sample in peptide_data_3cv$Sample_ID) {
      if (sample %in% colnames(all_presence_3cv)) {
        all_presence_3cv[peptide, sample] <- 1
        all_intensity_3cv[peptide, sample] <- peptide_data_3cv$Intensity[peptide_data_3cv$Sample_ID == sample]
      }
    }
  }
  
  # Log transform intensity matrices
  log_all_intensity_2cv <- log10(all_intensity_2cv + 1)
  log_all_intensity_3cv <- log10(all_intensity_3cv + 1)
  
  # Create heatmaps for all peptides
  pdf(file.path(results_dir, "all_junction_peptides_presence_2CV.pdf"), width = 12, height = 8)
  pheatmap(
    all_presence_2cv,
    main = "Presence of All Junction Peptides in 2CV Samples",
    color = c("white", "steelblue"),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    fontsize_row = 10,
    fontsize_col = 8,
    display_numbers = TRUE,
    # FIX: Use "%.0f" instead of "%d" for numeric format
    number_format = "%.0f"
  )
  dev.off()
  
  png(file.path(results_dir, "all_junction_peptides_presence_2CV.png"), width = 1000, height = 600, res = 100)
  pheatmap(
    all_presence_2cv,
    main = "Presence of All Junction Peptides in 2CV Samples",
    color = c("white", "steelblue"),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    fontsize_row = 10,
    fontsize_col = 8,
    display_numbers = TRUE,
    # FIX: Use "%.0f" instead of "%d" for numeric format
    number_format = "%.0f"
  )
  dev.off()
  
  pdf(file.path(results_dir, "all_junction_peptides_presence_3CV.pdf"), width = 12, height = 8)
  pheatmap(
    all_presence_3cv,
    main = "Presence of All Junction Peptides in 3CV Samples",
    color = c("white", "steelblue"),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    fontsize_row = 10,
    fontsize_col = 8,
    display_numbers = TRUE,
    # FIX: Use "%.0f" instead of "%d" for numeric format
    number_format = "%.0f"
  )
  dev.off()
  
  png(file.path(results_dir, "all_junction_peptides_presence_3CV.png"), width = 1000, height = 600, res = 100)
  pheatmap(
    all_presence_3cv,
    main = "Presence of All Junction Peptides in 3CV Samples",
    color = c("white", "steelblue"),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    fontsize_row = 10,
    fontsize_col = 8,
    display_numbers = TRUE,
    # FIX: Use "%.0f" instead of "%d" for numeric format
    number_format = "%.0f"
  )
  dev.off()
  
  pdf(file.path(results_dir, "all_junction_peptides_intensity_2CV.pdf"), width = 12, height = 8)
  pheatmap(
    log_all_intensity_2cv,
    main = "Intensity of All Junction Peptides in 2CV Samples (log10)",
    color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    fontsize_row = 10,
    fontsize_col = 8,
    display_numbers = TRUE,
    number_format = "%.1f"
  )
  dev.off()
  
  png(file.path(results_dir, "all_junction_peptides_intensity_2CV.png"), width = 1000, height = 600, res = 100)
  pheatmap(
    log_all_intensity_2cv,
    main = "Intensity of All Junction Peptides in 2CV Samples (log10)",
    color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    fontsize_row = 10,
    fontsize_col = 8,
    display_numbers = TRUE,
    number_format = "%.1f"
  )
  dev.off()
  
  pdf(file.path(results_dir, "all_junction_peptides_intensity_3CV.pdf"), width = 12, height = 8)
  pheatmap(
    log_all_intensity_3cv,
    main = "Intensity of All Junction Peptides in 3CV Samples (log10)",
    color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    fontsize_row = 10,
    fontsize_col = 8,
    display_numbers = TRUE,
    number_format = "%.1f"
  )
  dev.off()
  
  png(file.path(results_dir, "all_junction_peptides_intensity_3CV.png"), width = 1000, height = 600, res = 100)
  pheatmap(
    log_all_intensity_3cv,
    main = "Intensity of All Junction Peptides in 3CV Samples (log10)",
    color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    fontsize_row = 10,
    fontsize_col = 8,
    display_numbers = TRUE,
    number_format = "%.1f"
  )
  dev.off()
  
  # If there are paired samples, create a combined 2CV-3CV heatmap for direct comparison
  if (length(paired_samples) > 0) {
    # Create data for paired samples
    # First create detected column with mutate, then select
    paired_data <- combined_data %>%
      filter(Sample_ID %in% paired_samples) %>%
      filter(Peptide %in% detected_details$Peptide) %>%
      mutate(
        Sample_CV = paste0(Sample_ID, "_", CV_Type),
        detected = TRUE
      ) %>%
      select(Peptide, Sample_CV, detected, Intensity)
    
    # Create a complete matrix with all peptide-sample-CV combinations
    all_paired_combinations <- expand.grid(
      Peptide = detected_details$Peptide,
      Sample_ID = paired_samples,
      CV_Type = c("2CV", "3CV"),
      stringsAsFactors = FALSE
    ) %>%
      mutate(Sample_CV = paste0(Sample_ID, "_", CV_Type))
    
    # Merge with actual data
    paired_matrix_data <- all_paired_combinations %>%
      left_join(
        paired_data,
        by = c("Peptide", "Sample_CV")
      ) %>%
      mutate(
        detected = ifelse(is.na(detected), FALSE, detected),
        Intensity = ifelse(is.na(Intensity), 0, Intensity)
      )
    
    # Create matrices for presence and intensity
    paired_presence <- paired_matrix_data %>%
      select(Peptide, Sample_CV, present = detected) %>%
      pivot_wider(
        names_from = Sample_CV,
        values_from = present,
        values_fill = FALSE
      ) %>%
      column_to_rownames("Peptide")
    
    paired_intensity <- paired_matrix_data %>%
      select(Peptide, Sample_CV, Intensity) %>%
      pivot_wider(
        names_from = Sample_CV,
        values_from = Intensity,
        values_fill = 0
      ) %>%
      column_to_rownames("Peptide")
    
    # Convert logical presence matrix to numeric for pheatmap
    # FIX: Properly convert data frame to matrix first, then convert logical to numeric
    paired_presence_mat <- as.matrix(paired_presence)
    paired_presence_numeric <- matrix(as.numeric(paired_presence_mat), 
                                      nrow = nrow(paired_presence_mat),
                                      dimnames = dimnames(paired_presence_mat))
    
    # Log transform intensities
    log_paired_intensity <- log10(paired_intensity + 1)
    
    # Create column annotations for CV type
    column_ann <- data.frame(
      CV_Type = str_extract(colnames(paired_presence), "2CV|3CV"),
      Sample = gsub("_[23]CV$", "", colnames(paired_presence))
    )
    rownames(column_ann) <- colnames(paired_presence)
    
    # Define colors for annotation
    ann_colors <- list(
      CV_Type = c("2CV" = "steelblue", "3CV" = "tomato"),
      Sample = setNames(
        rainbow(length(paired_samples)),
        paired_samples
      )
    )
    
    # Create paired heatmaps
    pdf(file.path(results_dir, "paired_junction_peptides_presence.pdf"), width = 14, height = 8)
    pheatmap(
      paired_presence_numeric,
      main = "Presence of Junction Peptides in Paired 2CV-3CV Samples",
      color = c("white", "steelblue"),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      # FIX: Use "%.0f" instead of "%d" for numeric format
      number_format = "%.0f",
      annotation_col = column_ann,
      annotation_colors = ann_colors
    )
    dev.off()
    
    png(file.path(results_dir, "paired_junction_peptides_presence.png"), width = 1200, height = 600, res = 100)
    pheatmap(
      paired_presence_numeric,
      main = "Presence of Junction Peptides in Paired 2CV-3CV Samples",
      color = c("white", "steelblue"),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      # FIX: Use "%.0f" instead of "%d" for numeric format
      number_format = "%.0f",
      annotation_col = column_ann,
      annotation_colors = ann_colors
    )
    dev.off()
    
    pdf(file.path(results_dir, "paired_junction_peptides_intensity.pdf"), width = 14, height = 8)
    pheatmap(
      log_paired_intensity,
      main = "Intensity of Junction Peptides in Paired 2CV-3CV Samples (log10)",
      color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      number_format = "%.1f",
      annotation_col = column_ann,
      annotation_colors = ann_colors
    )
    dev.off()
    
    png(file.path(results_dir, "paired_junction_peptides_intensity.png"), width = 1200, height = 600, res = 100)
    pheatmap(
      log_paired_intensity,
      main = "Intensity of Junction Peptides in Paired 2CV-3CV Samples (log10)",
      color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      number_format = "%.1f",
      annotation_col = column_ann,
      annotation_colors = ann_colors
    )
    dev.off()
  }
  
  #===============================#
  # Create Excel Output
  #===============================#
  
  # Sheet with theoretical junction-spanning peptides
  excel_theoretical <- theoretical_junction_peptides %>%
    arrange(Length, Start_Position)
  
  # Sheet with detected junction peptides details
  excel_detected <- detected_details %>%
    arrange(Length, Start_Position)
  
  # Sheet with sample-level detection of junction peptides
  excel_sample_detection <- combined_data %>%
    filter(Spans_Junction == TRUE) %>%
    select(Sample_ID, CV_Type, Peptide, Peptide_Length, Intensity, Spectral.Count, Spans_Junction) %>%
    arrange(Sample_ID, CV_Type, Peptide)
  
  # Sheet with 2CV vs 3CV detection statistics
  excel_cv_comparison <- detected_details %>%
    select(
      Peptide,
      Length,
      DNAJB1_Residues,
      PRKACA_Residues,
      Visualization,
      Total_Detected,
      Detected_2CV,
      Detected_3CV,
      Avg_Intensity_2CV,
      Avg_Intensity_3CV
    ) %>%
    mutate(
      `2CV Detection %` = Detected_2CV / sum(sample_info$cv_type == "2CV") * 100,
      `3CV Detection %` = Detected_3CV / sum(sample_info$cv_type == "3CV") * 100,
      `3CV/2CV Detection Ratio` = (Detected_3CV + 0.001) / (Detected_2CV + 0.001),
      `3CV/2CV Intensity Ratio` = (Avg_Intensity_3CV + 1) / (Avg_Intensity_2CV + 1)
    ) %>%
    arrange(Length, Peptide)
  
  # Create a list of sheets for the Excel file
  excel_sheets <- list(
    "Theoretical_Junction_Peptides" = excel_theoretical,
    "Detected_Junction_Peptides" = excel_detected,
    "Sample_Level_Detection" = excel_sample_detection,
    "2CV_vs_3CV_Comparison" = excel_cv_comparison
  )
  
  # Write Excel file with multiple sheets
  write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
} else {
  cat("No junction-spanning peptides were detected in the dataset.\n")
  
  # Create Excel with just theoretical peptides
  excel_sheets <- list(
    "Theoretical_Junction_Peptides" = theoretical_junction_peptides
  )
  
  write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
}

# Print summary information
cat("\nAnalysis complete! Results saved to:", results_dir, "\n")
cat("\nThe following files were generated:\n")
cat("1. Fusion_Junction_Peptides_Analysis.xlsx - Excel file with detailed analysis\n")

if (nrow(detected_junction_peptides) > 0) {
  cat("2. Individual peptide heatmaps for each detected junction peptide\n")
  cat("3. all_junction_peptides_presence_2CV.pdf/png - Heatmap showing presence of all junction peptides in 2CV samples\n")
  cat("4. all_junction_peptides_presence_3CV.pdf/png - Heatmap showing presence of all junction peptides in 3CV samples\n")
  cat("5. all_junction_peptides_intensity_2CV.pdf/png - Heatmap showing intensity of all junction peptides in 2CV samples\n")
  cat("6. all_junction_peptides_intensity_3CV.pdf/png - Heatmap showing intensity of all junction peptides in 3CV samples\n")
  cat("7. junction_peptide_cv_comparison.pdf/png - Plot comparing detection in 2CV vs 3CV\n")
  cat("8. junction_peptide_coverage.pdf/png - Plot showing coverage of fusion junction by detected peptides\n")
  cat("9. junction_peptide_intensity_comparison.pdf/png - Plot comparing peptide intensities between 2CV and 3CV\n")
  
  if (length(paired_samples) > 0) {
    cat("10. paired_junction_peptides_presence.pdf/png - Heatmap showing presence in paired 2CV-3CV samples\n")
    cat("11. paired_junction_peptides_intensity.pdf/png - Heatmap showing intensity in paired 2CV-3CV samples\n")
  }
}