# DNAJB1-PRKACA Fusion Peptide Search Script
# This script searches for peptides that span the DNAJB1-PRKACA fusion junction
# Usage: Rscript search_fusion_peptides.R <tsv_file_path>

# Load required packages
library(tidyverse)
library(readr)
library(writexl)
library(ggplot2)

# Set working directory
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis")

# Define fusion protein sequence
# DNAJB1: KREIFDRYGEE | PRKACA: VKEFLAKAKED
FUSION_SEQUENCE <- "KREIFDRYGEEVKEFLAKAKED"
JUNCTION_POSITION <- 12  # Position where DNAJB1 ends and PRKACA begins (1-indexed)

# Function to create output directory structure
create_output_directory <- function(search_name) {
  output_dir <- file.path("peptide_search", search_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", output_dir, "\n")
  return(output_dir)
}

# Function to read the TSV file with proper column types
read_peptide_data <- function(tsv_file) {
  # Create a column specification to handle problematic columns
  col_spec <- cols(
    `Assigned Modifications` = col_character(),
    `Observed Modifications` = col_character(),
    Intensity = col_double(),
    .default = col_guess()
  )
  
  # Read the TSV file
  peptide_data <- read_tsv(tsv_file, col_types = col_spec, show_col_types = FALSE)
  
  # Check for parsing problems
  parsing_problems <- problems(peptide_data)
  if (nrow(parsing_problems) > 0) {
    cat("\nFound", nrow(parsing_problems), "parsing issues:\n")
    print(head(parsing_problems, 10))
    write_csv(parsing_problems, "fusion_parsing_problems.csv")
    cat("All parsing issues saved to fusion_parsing_problems.csv\n")
  }
  
  return(peptide_data)
}

# Function to generate all possible spanning peptides from fusion sequence
generate_spanning_peptides <- function(sequence, junction_pos, min_length = 8, max_length = 12) {
  spanning_peptides <- c()
  sequence_length <- nchar(sequence)
  
  # Generate all possible peptides of specified lengths
  for (length in min_length:max_length) {
    # Generate all possible start positions
    for (start_pos in 1:(sequence_length - length + 1)) {
      end_pos <- start_pos + length - 1
      
      # Check if this peptide spans the junction
      # It must include at least one AA before junction AND one AA after junction
      if (start_pos < junction_pos && end_pos >= junction_pos) {
        peptide <- substr(sequence, start_pos, end_pos)
        spanning_peptides <- c(spanning_peptides, peptide)
      }
    }
  }
  
  return(unique(spanning_peptides))
}

# Enhanced peptide summary with sample information
create_enhanced_peptide_summary <- function(filtered_data, output_dir, search_name) {
  # Get list of samples per peptide
  peptide_samples <- filtered_data %>%
    group_by(Peptide, SampleID) %>%
    summarise(
      Sample_Intensity = sum(Intensity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Peptide, desc(Sample_Intensity))
  
  # Create sample list as string
  sample_list <- peptide_samples %>%
    group_by(Peptide) %>%
    summarise(
      Samples = paste(SampleID, collapse = ", "),
      Sample_Count = n_distinct(SampleID),
      .groups = "drop"
    )
  
  # Create detailed peptide summary
  detailed_peptide_summary <- filtered_data %>%
    group_by(Peptide, `Peptide Length`) %>%
    summarise(
      Count = n(),
      Total_Intensity = sum(Intensity, na.rm = TRUE),
      Mean_Intensity = mean(Intensity, na.rm = TRUE),
      Max_Intensity = max(Intensity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(sample_list, by = "Peptide") %>%
    arrange(desc(Sample_Count), desc(Total_Intensity))
  
  # Write to Excel
  write_xlsx(detailed_peptide_summary, file.path(output_dir, paste0("peptide_summary_", search_name, ".xlsx")))
  cat("Enhanced peptide summary saved to", file.path(output_dir, paste0("peptide_summary_", search_name, ".xlsx")), "\n")
  
  return(detailed_peptide_summary)
}

# Create visualizations
create_visualizations <- function(filtered_data, peptide_summary, output_dir, search_name) {
  # 1. Peptide intensity across samples (bar chart)
  # Get all peptides (since we likely have few fusion peptides)
  all_peptides <- unique(filtered_data$Peptide)
  
  # Prepare data for visualization
  vis_data <- filtered_data %>%
    group_by(Peptide, SampleID) %>%
    summarise(Intensity = sum(Intensity, na.rm = TRUE), .groups = "drop")
  
  # Create the bar chart
  p1 <- ggplot(vis_data, aes(x = reorder(Peptide, Intensity, sum), y = Intensity, fill = SampleID)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(title = "Fusion Peptide Intensities Across Samples",
         x = "Peptide", y = "Intensity") +
    scale_y_log10()
  
  # Save the plot
  ggsave(file.path(output_dir, paste0("fusion_peptide_intensity_by_sample_", search_name, ".png")), 
         p1, width = 12, height = 7)
  cat("Created visualization: Fusion peptide intensities across samples\n")
  
  # 2. Sample comparison showing peptide counts
  sample_counts <- filtered_data %>%
    group_by(SampleID) %>%
    summarise(
      Peptide_Count = n_distinct(Peptide),
      Total_Intensity = sum(Intensity, na.rm = TRUE)
    )
  
  p2 <- ggplot(sample_counts, aes(x = reorder(SampleID, Peptide_Count), y = Peptide_Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = "Number of Fusion Peptides per Sample",
         x = "Sample ID", y = "Peptide Count")
  
  ggsave(file.path(output_dir, paste0("fusion_peptide_counts_by_sample_", search_name, ".png")), 
         p2, width = 10, height = 6)
  cat("Created visualization: Fusion peptide counts by sample\n")
  
  # 3. Intensity distribution
  p3 <- ggplot(filtered_data, aes(x = log10(Intensity))) +
    geom_histogram(bins = 20, fill = "darkgreen") +
    facet_wrap(~SampleID) +
    theme_minimal() +
    labs(title = "Fusion Peptide Intensity Distribution by Sample",
         x = "Log10(Intensity)", y = "Count")
  
  ggsave(file.path(output_dir, paste0("fusion_intensity_distribution_", search_name, ".png")), 
         p3, width = 10, height = 8)
  cat("Created visualization: Fusion peptide intensity distribution by sample\n")
  
  # 4. Peptide position in fusion sequence visualization
  if (nrow(filtered_data) > 0) {
    # Calculate position of each peptide in the fusion sequence
    peptide_positions <- filtered_data %>%
      distinct(Peptide) %>%
      mutate(
        Position_in_Fusion = map_dbl(Peptide, ~{
          pos <- str_locate(FUSION_SEQUENCE, fixed(.x))
          if (!is.na(pos[1])) pos[1] else NA
        }),
        Spans_Junction = map_lgl(Peptide, ~{
          pos <- str_locate(FUSION_SEQUENCE, fixed(.x))
          if (!is.na(pos[1])) {
            start_pos <- pos[1]
            end_pos <- pos[2]
            start_pos < JUNCTION_POSITION && end_pos >= JUNCTION_POSITION
          } else FALSE
        })
      ) %>%
      filter(!is.na(Position_in_Fusion))
    
    if (nrow(peptide_positions) > 0) {
      p4 <- ggplot(peptide_positions, aes(x = Position_in_Fusion, y = Peptide, color = Spans_Junction)) +
        geom_point(size = 3) +
        geom_vline(xintercept = JUNCTION_POSITION - 0.5, linetype = "dashed", color = "red", size = 1) +
        theme_minimal() +
        labs(title = "Fusion Peptide Positions Relative to Junction",
             x = "Position in Fusion Sequence", 
             y = "Peptide",
             color = "Spans Junction") +
        annotate("text", x = JUNCTION_POSITION - 0.5, y = Inf, label = "Junction", 
                 vjust = 2, color = "red", size = 3)
      
      ggsave(file.path(output_dir, paste0("fusion_peptide_positions_", search_name, ".png")), 
             p4, width = 12, height = max(6, nrow(peptide_positions) * 0.3))
      cat("Created visualization: Fusion peptide positions relative to junction\n")
    }
  }
}

# Main function to search for fusion peptides
search_fusion_peptides <- function(tsv_file) {
  cat("DNAJB1-PRKACA Fusion Peptide Search\n")
  cat("===================================\n")
  cat("Fusion sequence:", FUSION_SEQUENCE, "\n")
  cat("Junction position:", JUNCTION_POSITION, "\n")
  cat("Reading peptide data from:", tsv_file, "\n\n")
  
  # Read the peptide data
  peptide_data <- read_peptide_data(tsv_file)
  cat("Total entries in dataset:", nrow(peptide_data), "\n")
  
  # Generate all possible spanning peptides
  cat("Generating all possible 8-12mer spanning peptides...\n")
  spanning_peptides <- generate_spanning_peptides(FUSION_SEQUENCE, JUNCTION_POSITION)
  
  cat("Generated", length(spanning_peptides), "possible spanning peptides:\n")
  for (i in seq_along(spanning_peptides)) {
    cat(sprintf("  %2d. %s (length: %d)\n", i, spanning_peptides[i], nchar(spanning_peptides[i])))
  }
  cat("\n")
  
  # Search for these peptides in the data
  cat("Searching for fusion peptides in dataset...\n")
  found_peptides <- peptide_data %>%
    filter(Peptide %in% spanning_peptides)
  
  cat("Found", nrow(found_peptides), "entries matching fusion peptides\n")
  
  if (nrow(found_peptides) == 0) {
    cat("No fusion peptides found in the dataset.\n")
    
    # Still create output directory and save the theoretical peptides
    output_dir <- create_output_directory("DNAJB1_PRKACA_fusion")
    
    # Save theoretical peptides to file
    theoretical_peptides <- data.frame(
      Theoretical_Peptide = spanning_peptides,
      Length = nchar(spanning_peptides),
      Position_in_Fusion = sapply(spanning_peptides, function(p) {
        pos <- str_locate(FUSION_SEQUENCE, fixed(p))
        pos[1]
      }),
      Spans_Junction = TRUE
    )
    
    write_xlsx(theoretical_peptides, file.path(output_dir, "theoretical_fusion_peptides.xlsx"))
    cat("Theoretical fusion peptides saved to", file.path(output_dir, "theoretical_fusion_peptides.xlsx"), "\n")
    
    return(invisible())
  }
  
  # Create output directory
  output_dir <- create_output_directory("DNAJB1_PRKACA_fusion")
  
  # Show which specific peptides were found
  found_unique_peptides <- unique(found_peptides$Peptide)
  cat("Found the following fusion peptides:\n")
  for (i in seq_along(found_unique_peptides)) {
    peptide <- found_unique_peptides[i]
    count <- sum(found_peptides$Peptide == peptide)
    cat(sprintf("  %s (found in %d entries)\n", peptide, count))
  }
  cat("\n")
  
  # Select columns for output
  columns_to_include <- c(
    "Peptide", "Peptide Length", "Charges", "Probability", "Spectral Count", 
    "Intensity", "Protein", "Protein ID", "Entry Name", "Gene", 
    "Protein Description", "Mapped Genes", "Mapped Proteins", "SampleID",
    "detected_2cv", "detected_3cv", "detected_both", "CVType"
  )
  
  # Check which columns actually exist
  available_columns <- names(found_peptides)
  columns_to_include <- columns_to_include[columns_to_include %in% available_columns]
  
  # Add other important columns
  other_important_cols <- c("SourceFile", "Match Type", "Start", "End", "Prev AA", "Next AA")
  for (col in other_important_cols) {
    if (col %in% available_columns && !(col %in% columns_to_include)) {
      columns_to_include <- c(columns_to_include, col)
    }
  }
  
  # Create selected data
  selected_data <- found_peptides %>% select(all_of(columns_to_include))
  
  # Save main results
  write_xlsx(selected_data, file.path(output_dir, "fusion_peptide_results.xlsx"))
  cat("Fusion peptide results saved to", file.path(output_dir, "fusion_peptide_results.xlsx"), "\n")
  
  # Create intensity summary by sample
  intensity_summary <- found_peptides %>%
    group_by(SampleID) %>%
    summarise(
      Peptide_Count = n_distinct(Peptide),
      Total_Intensity = sum(Intensity, na.rm = TRUE),
      Mean_Intensity = mean(Intensity, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("Intensity summary by sample:\n")
  print(intensity_summary)
  
  write_xlsx(intensity_summary, file.path(output_dir, "fusion_intensity_summary.xlsx"))
  cat("Intensity summary saved to", file.path(output_dir, "fusion_intensity_summary.xlsx"), "\n")
  
  # Create enhanced peptide summary
  peptide_summary <- create_enhanced_peptide_summary(found_peptides, output_dir, "fusion")
  
  # Create visualizations
  create_visualizations(found_peptides, peptide_summary, output_dir, "fusion")
  
  # Save theoretical vs found comparison
  theoretical_peptides <- data.frame(
    Theoretical_Peptide = spanning_peptides,
    Length = nchar(spanning_peptides),
    Found_in_Data = spanning_peptides %in% found_unique_peptides,
    stringsAsFactors = FALSE
  )
  
  write_xlsx(theoretical_peptides, file.path(output_dir, "theoretical_vs_found_peptides.xlsx"))
  cat("Theoretical vs found peptides comparison saved to", file.path(output_dir, "theoretical_vs_found_peptides.xlsx"), "\n")
  
  # Summary statistics
  cat("\n=== SUMMARY ===\n")
  cat("Fusion sequence:", FUSION_SEQUENCE, "\n")
  cat("Total possible spanning peptides (8-12mers):", length(spanning_peptides), "\n")
  cat("Peptides found in data:", length(found_unique_peptides), "\n")
  cat("Total entries found:", nrow(found_peptides), "\n")
  cat("Samples with fusion peptides:", length(unique(found_peptides$SampleID)), "\n")
  cat("Coverage: ", round(length(found_unique_peptides)/length(spanning_peptides)*100, 1), "%\n")
}

# Function to show usage
show_usage <- function() {
  cat("DNAJB1-PRKACA Fusion Peptide Search Tool\n")
  cat("========================================\n")
  cat("This script searches for peptides that span the DNAJB1-PRKACA fusion junction.\n")
  cat("Fusion sequence: KREIFDRYGEEVKEFLAKAKED\n")
  cat("Junction position: between E and V (positions 11-12)\n\n")
  cat("Usage:\n")
  cat("  Rscript search_fusion_peptides.R <tsv_file_path>\n")
  cat("  Example: Rscript search_fusion_peptides.R unique_peptides_unmodified.tsv\n\n")
  cat("Output:\n")
  cat("  - Excel file with found fusion peptides\n")
  cat("  - Intensity summary by sample\n")
  cat("  - Enhanced peptide summary with sample information\n")
  cat("  - Theoretical vs found peptides comparison\n")
  cat("  - Visualizations (bar charts, distributions, position plots)\n")
  cat("  - All outputs saved in 'peptide_search/DNAJB1_PRKACA_fusion/' directory\n")
}

# Main execution
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    show_usage()
    return(invisible())
  }
  
  tsv_file <- args[1]
  
  if (!file.exists(tsv_file)) {
    cat("Error: TSV file not found:", tsv_file, "\n")
    return(invisible())
  }
  
  search_fusion_peptides(tsv_file)
}

# Run the main function
main()