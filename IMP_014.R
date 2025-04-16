# Script to analyze peptides derived from DNAJB1_PRKACA fusion protein
# across different samples from immunopeptidomics experiment

# Load required packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Define the path to your data files
data_path <- "/Users/dinarabadi/Desktop/2704DR"

# Define the sequences of DNAJB1 and PRKACA parts of the fusion protein
dnajb1_seq <- "GKDYYQTLGLARGASDEEIKRAYRRQALRYHPDKNKEPGAEEKFKEIAEAYDVLSDPRKREIFDRYGEE"
prkaca_seq <- "VKEFLAKAKEDFLKKWESPAQNTAHLDQFERIKTLGTGSFGRVMLVKHKETGNHYAMKILDKQKVVKLKQIEHTLNEKRILQAVNFPFLVKLEFSFKDNSNLYMVMEYVPGGEMFSHLRRIGRFSEPHARFYAAQIVLTFEYLHSLDLIYRDLKPENLLIDQQGYIQVTDFGFAKRVKGRTWTLCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFATTDWIAIYQRKVEAPFIPKFKGPGDTSNFDDYEEEEIRVSINEKCGKEFSEF"
fusion_protein <- paste0(dnajb1_seq, prkaca_seq)

# Define the synthetic peptides used for spiking
synthetic_peptides <- c(
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
  
  # Add to our list
  all_data[[sample_id]] <- data
  sample_ids <- c(sample_ids, sample_id)
}

# Get unique sample IDs (in case multiple files had the same sample ID)
sample_ids <- unique(sample_ids)
cat("Found data for", length(sample_ids), "unique sample IDs:", paste(sample_ids, collapse = ", "), "\n")

# Combine all data frames
combined_data <- bind_rows(all_data)

# Create a function to check if a peptide spans the fusion junction with more precise detection
is_fusion_junction_peptide <- function(peptide_seq) {
  # Define the end of DNAJB1 and start of PRKACA for the fusion
  # The fusion junction is between "IFDRYGEE" (DNAJB1) and "VKEFLAK" (PRKACA)
  dnajb1_end <- "IFDRYGEE"
  prkaca_start <- "VKEFLAK"
  
  # Check if the peptide spans the fusion junction
  spans_junction <- FALSE
  junction_position <- NA
  junction_details <- "Not a junction peptide"
  
  # Check if peptide is one of the synthetic peptides
  is_synthetic <- peptide_seq %in% synthetic_peptides
  
  # For the synthetic peptides, we know they span the junction
  if (is_synthetic) {
    spans_junction <- TRUE
    
    # For each synthetic peptide, identify the junction position
    # Pattern: match the end of DNAJB1 with start of PRKACA
    if (peptide_seq == "EEVKEFLAK") {
      junction_position <- 2  # EE|VKEFLAK
      junction_details <- "Junction at position 2 (EE|VKEFLAK)"
    } else if (peptide_seq == "YGEEVKEFL") {
      junction_position <- 4  # YGEE|VKEFL
      junction_details <- "Junction at position 4 (YGEE|VKEFL)"
    } else if (peptide_seq == "RYGEEVKEF") {
      junction_position <- 5  # RYGEE|VKEF
      junction_details <- "Junction at position 5 (RYGEE|VKEF)"
    } else if (peptide_seq == "RYGEEVKEFL") {
      junction_position <- 5  # RYGEE|VKEFL
      junction_details <- "Junction at position 5 (RYGEE|VKEFL)"
    } else if (peptide_seq == "EIFDRYGEEV") {
      junction_position <- 9  # EIFDRYGEE|V
      junction_details <- "Junction at position 9 (EIFDRYGEE|V)"
    } else if (peptide_seq == "IFDRYGEEV") {
      junction_position <- 8  # IFDRYGEE|V
      junction_details <- "Junction at position 8 (IFDRYGEE|V)"
    }
  } else {
    # For naturally occurring peptides, scan for junction
    if (nchar(peptide_seq) >= 6) { # Require at least 6 AA to be meaningful
      for (i in 3:(nchar(peptide_seq) - 3)) { # At least 3 AA from each side
        left_part <- substr(peptide_seq, 1, i)
        right_part <- substr(peptide_seq, i + 1, nchar(peptide_seq))
        
        # Check if left part ends with part of DNAJB1 end and right part starts with part of PRKACA start
        if (grepl(left_part, dnajb1_seq, fixed = TRUE) && 
            grepl(right_part, prkaca_seq, fixed = TRUE)) {
          
          # Additional check to ensure left part aligns with end of DNAJB1
          left_pos <- gregexpr(left_part, dnajb1_seq, fixed = TRUE)[[1]]
          if (length(left_pos) > 0 && any(left_pos + nchar(left_part) - 1 >= nchar(dnajb1_seq) - 10)) {
            
            # Additional check to ensure right part aligns with start of PRKACA
            right_pos <- gregexpr(right_part, prkaca_seq, fixed = TRUE)[[1]]
            if (length(right_pos) > 0 && any(right_pos <= 10)) {
              spans_junction <- TRUE
              junction_position <- i
              junction_details <- paste0("Junction at position ", i, " (", left_part, "|", right_part, ")")
              break
            }
          }
        }
      }
    }
  }
  
  return(list(
    spans_junction = spans_junction,
    is_synthetic = is_synthetic,
    junction_position = junction_position,
    junction_details = junction_details,
    junction_type = if(spans_junction) "Fusion Junction" else "Not Junction",
    peptide_class = if(spans_junction) paste("Fusion Junction", if(is_synthetic) "- Synthetic" else "- Natural") else "Not Junction"
  ))
}

# Extract peptides that span the fusion junction or are synthetic
junction_peptides <- combined_data %>%
  # Create a new column with peptide classification
  rowwise() %>%
  mutate(
    junction_result = list(is_fusion_junction_peptide(Peptide)),
    junction_type = junction_result$junction_type,
    spans_junction = junction_result$spans_junction,
    junction_position = junction_result$junction_position,
    junction_details = junction_result$junction_details,
    is_synthetic = junction_result$is_synthetic,
    peptide_class = junction_result$peptide_class
  ) %>%
  filter(spans_junction | is_synthetic) %>%
  select(-junction_result)

# Summarize peptides by sample
peptide_summary <- junction_peptides %>%
  group_by(Sample_ID, Peptide, junction_type, junction_details, is_synthetic, peptide_class) %>%
  summarize(
    spectral_count = sum(Spectral.Count),
    total_intensity = sum(Intensity),
    .groups = "drop"
  ) %>%
  arrange(Sample_ID, desc(total_intensity))

# Create a presence/absence matrix for heatmap visualization
peptide_matrix <- junction_peptides %>%
  group_by(Sample_ID, Peptide) %>%
  summarize(
    intensity = sum(Intensity),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Sample_ID,
    values_from = intensity,
    values_fill = 0
  ) %>%
  column_to_rownames("Peptide")

# Log transform the intensity values for better visualization (adding a small value to avoid log(0))
log_peptide_matrix <- log10(peptide_matrix + 1)

# Create a matrix showing presence (1) or absence (0)
presence_matrix <- (peptide_matrix > 0) * 1

# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("white", "steelblue", "darkblue"))(100)

# Create a summary of junction peptides with their details
junction_peptide_summary <- junction_peptides %>%
  distinct(Peptide, junction_type, junction_details, is_synthetic, peptide_class) %>%
  arrange(junction_type, Peptide)

# Create a summary of synthetic peptide detection across samples
synthetic_summary <- junction_peptides %>%
  filter(is_synthetic) %>%
  group_by(Sample_ID, Peptide) %>%
  summarize(
    detected = n() > 0,
    spectral_count = sum(Spectral.Count),
    total_intensity = sum(Intensity),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Sample_ID,
    values_from = c(detected, spectral_count, total_intensity),
    values_fill = list(detected = FALSE, spectral_count = 0, total_intensity = 0)
  )

# Write results to files
write.csv(junction_peptides, file.path(data_path, "junction_peptides_details.csv"), row.names = FALSE)
write.csv(peptide_summary, file.path(data_path, "junction_peptides_summary.csv"), row.names = FALSE)
write.csv(junction_peptide_summary, file.path(data_path, "junction_peptide_summary.csv"), row.names = FALSE)
write.csv(synthetic_summary, file.path(data_path, "synthetic_peptides_summary.csv"), row.names = FALSE)

# Create visualizations
# 1. Heatmap of peptide presence/absence across samples
pdf(file.path(data_path, "junction_peptides_heatmap.pdf"), width = 12, height = 10)
pheatmap(
  presence_matrix,
  main = "Presence of DNAJB1_PRKACA Junction Peptides Across Samples",
  color = c("white", "darkblue"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 10
)
dev.off()

# 2. Heatmap of peptide intensities across samples
pdf(file.path(data_path, "junction_peptides_intensity_heatmap.pdf"), width = 12, height = 10)
pheatmap(
  log_peptide_matrix,
  main = "Intensity of DNAJB1_PRKACA Junction Peptides Across Samples (log10)",
  color = color_palette,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 10
)
dev.off()

# 3. Bar plot of number of junction peptides per sample
junction_counts <- junction_peptides %>%
  group_by(Sample_ID, peptide_class) %>%
  summarize(count = n_distinct(Peptide), .groups = "drop")

pdf(file.path(data_path, "junction_peptides_count_by_sample.pdf"), width = 10, height = 6)
ggplot(junction_counts, aes(x = Sample_ID, y = count, fill = peptide_class)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(
    title = "Number of DNAJB1_PRKACA Junction Peptides Detected per Sample",
    x = "Sample ID",
    y = "Number of Unique Peptides",
    fill = "Peptide Class"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 4. Synthetic peptide recovery in sample 51S
synthetic_51s_samples <- grep("51S", sample_ids, value = TRUE)
if (length(synthetic_51s_samples) > 0) {
  for (sample_51s in synthetic_51s_samples) {
    synthetic_in_sample <- junction_peptides %>%
      filter(Sample_ID == sample_51s, is_synthetic) %>%
      group_by(Peptide) %>%
      summarize(
        detected = TRUE,
        spectral_count = sum(Spectral.Count),
        total_intensity = sum(Intensity),
        .groups = "drop"
      )
    
    # Add any missing synthetic peptides with FALSE for detected
    missing_peptides <- setdiff(synthetic_peptides, synthetic_in_sample$Peptide)
    if (length(missing_peptides) > 0) {
      missing_df <- data.frame(
        Peptide = missing_peptides,
        detected = FALSE,
        spectral_count = 0,
        total_intensity = 0
      )
      synthetic_in_sample <- bind_rows(synthetic_in_sample, missing_df)
    }
    
    # Create a plot for this sample
    pdf(file.path(data_path, paste0("synthetic_peptide_recovery_", sample_51s, ".pdf")), width = 10, height = 6)
    ggplot(synthetic_in_sample, aes(x = Peptide, y = total_intensity, fill = detected)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(
        title = paste("Recovery of Synthetic Peptides in Sample", sample_51s),
        x = "Peptide Sequence",
        y = "Total Intensity",
        fill = "Detected"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    dev.off()
    
    # Print detection information
    cat(paste("\nSynthetic peptide recovery in sample", sample_51s, ":\n"))
    synthetic_recovery <- synthetic_in_sample %>%
      arrange(Peptide) %>%
      select(Peptide, detected, spectral_count, total_intensity)
    print(synthetic_recovery)
  }
}

# 5. Create a detailed visualization of the fusion junction and the peptides that span it
junction_viz_data <- junction_peptide_summary %>%
  filter(junction_type == "Fusion Junction") %>%
  mutate(
    junction_pos = as.numeric(gsub("Junction at position ([0-9]+) .*", "\\1", junction_details)),
    left_part = gsub("Junction at position [0-9]+ \\((.*)\\|.*\\)", "\\1", junction_details),
    right_part = gsub("Junction at position [0-9]+ \\(.*\\|(.*)\\)", "\\1", junction_details)
  )

# Create a visualization showing where each peptide spans the junction
if (nrow(junction_viz_data) > 0) {
  pdf(file.path(data_path, "junction_peptide_map.pdf"), width = 12, height = 8)
  
  # Calculate the maximum peptide length for visualization
  max_peptide_length <- max(nchar(junction_viz_data$Peptide))
  
  # Create a plot showing the junction and how each peptide spans it
  par(mar = c(5, 10, 4, 2))
  plot(
    NULL, 
    xlim = c(1, max_peptide_length), 
    ylim = c(0, nrow(junction_viz_data) + 1), 
    xlab = "Amino Acid Position", 
    ylab = "", 
    main = "DNAJB1-PRKACA Junction Peptides",
    yaxt = "n"
  )
  
  # Add a line representing the junction
  abline(v = junction_viz_data$junction_pos[1], col = "red", lty = 2)
  text(junction_viz_data$junction_pos[1], nrow(junction_viz_data) + 0.5, "Junction", pos = 3, col = "red")
  
  # Plot each peptide
  for (i in 1:nrow(junction_viz_data)) {
    peptide <- junction_viz_data$Peptide[i]
    junction_pos <- junction_viz_data$junction_pos[i]
    
    # Get left and right parts
    left_part <- junction_viz_data$left_part[i]
    right_part <- junction_viz_data$right_part[i]
    
    # Draw the peptide
    lines(c(1, nchar(peptide)), c(i, i), lwd = 2, col = "black")
    
    # Color the parts differently
    lines(c(1, junction_pos), c(i, i), lwd = 2, col = "blue")
    lines(c(junction_pos + 1, nchar(peptide)), c(i, i), lwd = 2, col = "green")
    
    # Add peptide name
    axis(2, at = i, labels = peptide, las = 2, cex.axis = 0.7)
  }
  
  # Add a legend
  legend("topright", 
         legend = c("DNAJB1 Part", "PRKACA Part", "Junction"),
         col = c("blue", "green", "red"),
         lty = c(1, 1, 2),
         lwd = c(2, 2, 1))
  
  dev.off()
}

# 6. Create additional visualization: Compare the same peptide across different samples
for (peptide in unique(junction_peptides$Peptide)) {
  peptide_data <- junction_peptides %>%
    filter(Peptide == peptide) %>%
    group_by(Sample_ID) %>%
    summarize(
      spectral_count = sum(Spectral.Count),
      total_intensity = sum(Intensity),
      .groups = "drop"
    )
  
  # Add missing samples with zero values
  missing_samples <- setdiff(sample_ids, peptide_data$Sample_ID)
  if (length(missing_samples) > 0) {
    missing_df <- data.frame(
      Sample_ID = missing_samples,
      spectral_count = 0,
      total_intensity = 0
    )
    peptide_data <- bind_rows(peptide_data, missing_df)
  }
  
  # Sort by the order of sample_ids to maintain consistency
  peptide_data$Sample_ID <- factor(peptide_data$Sample_ID, levels = sample_ids)
  peptide_data <- peptide_data %>% arrange(Sample_ID)
  
  # Skip if all intensities are 0
  if (sum(peptide_data$total_intensity) == 0) next
  
  # Create a plot for this peptide
  pdf(file.path(data_path, paste0("peptide_comparison_", gsub("[^A-Za-z0-9]", "_", peptide), ".pdf")), width = 10, height = 6)
  p <- ggplot(peptide_data, aes(x = Sample_ID, y = total_intensity, fill = Sample_ID)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(
      title = paste("Intensity of Peptide", peptide, "Across Samples"),
      x = "Sample ID",
      y = "Total Intensity"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  dev.off()
}

# Print summary information
cat("\nSummary of DNAJB1_PRKACA junction peptides detected across samples:\n")
print(junction_counts)

cat("\nDetails of detected junction peptides:\n")
print(junction_peptide_summary %>% select(Peptide, junction_details, is_synthetic))

cat("\nAnalysis complete! Results saved to:", data_path, "\n")
cat("\nThe following files were generated:\n")
cat("1. junction_peptides_details.csv - Complete details of all junction-spanning and synthetic peptides\n")
cat("2. junction_peptides_summary.csv - Summary of junction peptides by sample\n")
cat("3. junction_peptide_summary.csv - Details about each unique junction peptide\n")
cat("4. synthetic_peptides_summary.csv - Summary of synthetic peptide detection\n")
cat("5. junction_peptides_heatmap.pdf - Heatmap showing presence/absence\n")
cat("6. junction_peptides_intensity_heatmap.pdf - Heatmap showing intensities\n")
cat("7. junction_peptides_count_by_sample.pdf - Bar plot of peptide counts\n")
cat("8. synthetic_peptide_recovery_51S.pdf - Recovery of synthetic peptides in sample 51S\n")
cat("9. junction_peptide_map.pdf - Visual map of junction-spanning peptides\n")
cat("10. Individual peptide comparison charts - One per detected peptide\n")