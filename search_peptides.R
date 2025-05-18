# Enhanced script to search for specific gene or protein names in concatenated peptide data
# Usage: 
#   Rscript search_peptides.R <combined_file_path> <gene_name> [protein_name]
#   Rscript search_peptides.R --list-genes <combined_file_path>
#   Rscript search_peptides.R --search-partial <combined_file_path> <search_pattern>

# Complete Search Options: usage within Rstudio console
# Now you have four ways to search your data:
# Search by gene or protein name (exact match):
#   rsystem("Rscript search_peptides.R combined_peptides.tsv AKAP12")
#   rsystem("Rscript search_peptides.R combined_peptides.tsv NA Q9Y3F4-2")
# List all genes in your dataset:
#   rsystem("Rscript search_peptides.R --list-genes combined_peptides.tsv")
# Search across all fields (partial match):
#   rsystem("Rscript search_peptides.R --search-partial combined_peptides.tsv PHOX")
# Search for exact peptide sequences:
#   rsystem("Rscript search_peptides.R --search-peptide combined_peptides.tsv AAADFTAKV")

# Load required packages
library(tidyverse)
library(readr)
library(writexl)
library(ggplot2)

# Function to create output directory structure
create_output_directory <- function(search_name) {
  # Create main output directory
  output_dir <- file.path("peptide_search", search_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", output_dir, "\n")
  return(output_dir)
}

# Function to read the combined TSV file with proper column types
read_combined_data <- function(combined_file) {
  # First read the header to get column names
  header_row <- read_tsv(combined_file, n_max = 0, show_col_types = FALSE)
  column_names <- names(header_row)
  
  # Create a column specification to handle problematic columns
  col_spec <- cols(
    # Specify problematic columns as character
    `Assigned Modifications` = col_character(),
    `Observed Modifications` = col_character(),
    # Make sure intensity is a number
    Intensity = col_double(),
    # For all other columns, use default type detection
    .default = col_guess()
  )
  
  # Read the combined TSV file with specified column types
  combined_data <- read_tsv(combined_file, col_types = col_spec, show_col_types = FALSE)
  
  # Check for parsing problems
  parsing_problems <- problems(combined_data)
  if (nrow(parsing_problems) > 0) {
    cat("\nFound", nrow(parsing_problems), "parsing issues:\n")
    print(head(parsing_problems, 10))  # Print first 10 problems
    
    # Save the problems to a file for further investigation
    write_csv(parsing_problems, "parsing_problems.csv")
    cat("All parsing issues saved to parsing_problems.csv\n")
  }
  
  return(combined_data)
}

# Function to list all unique genes in the dataset
list_unique_genes <- function(combined_file) {
  cat("Reading combined peptide data...\n")
  combined_data <- read_combined_data(combined_file)
  
  cat("Total entries in combined data:", nrow(combined_data), "\n")
  
  # Create output directory
  output_dir <- create_output_directory("gene_list")
  
  # Extract unique gene names
  unique_genes <- unique(combined_data$Gene)
  unique_genes <- unique_genes[!is.na(unique_genes)]  # Remove NA values
  unique_genes <- sort(unique_genes)  # Sort alphabetically
  
  cat("Found", length(unique_genes), "unique gene names\n")
  
  # Save to file
  writeLines(unique_genes, file.path(output_dir, "unique_genes.txt"))
  cat("Unique gene names saved to", file.path(output_dir, "unique_genes.txt"), "\n")
  
  # Write to CSV with counts
  gene_counts <- combined_data %>%
    group_by(Gene) %>%
    summarise(
      Count = n(),
      Unique_Peptides = n_distinct(Peptide),
      .groups = "drop"
    ) %>%
    arrange(desc(Count))
  
  write_csv(gene_counts, file.path(output_dir, "gene_counts.csv"))
  cat("Gene counts saved to", file.path(output_dir, "gene_counts.csv"), "\n")
  
  # Also create a summary of protein IDs
  protein_counts <- combined_data %>%
    group_by(`Protein ID`) %>%
    summarise(
      Count = n(),
      Unique_Peptides = n_distinct(Peptide),
      Gene = first(Gene),
      .groups = "drop"
    ) %>%
    arrange(desc(Count))
  
  write_csv(protein_counts, file.path(output_dir, "protein_counts.csv"))
  cat("Protein counts saved to", file.path(output_dir, "protein_counts.csv"), "\n")
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
  # Get top peptides by total intensity
  top_peptides <- peptide_summary %>%
    arrange(desc(Total_Intensity)) %>%
    head(15) %>%
    pull(Peptide)
  
  # Prepare data for visualization
  vis_data <- filtered_data %>%
    filter(Peptide %in% top_peptides) %>%
    group_by(Peptide, SampleID) %>%
    summarise(Intensity = sum(Intensity, na.rm = TRUE), .groups = "drop")
  
  # Create the bar chart
  p1 <- ggplot(vis_data, aes(x = reorder(Peptide, Intensity, sum), y = Intensity, fill = SampleID)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Peptide Intensities Across Samples",
         x = "Peptide", y = "Intensity") +
    scale_y_log10() # Log scale for better visualization of intensities
  
  # Save the plot
  ggsave(file.path(output_dir, paste0("peptide_intensity_by_sample_", search_name, ".png")), 
         p1, width = 12, height = 7)
  cat("Created visualization: Peptide intensities across samples\n")
  
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
    labs(title = "Number of Unique Peptides per Sample",
         x = "Sample ID", y = "Peptide Count")
  
  ggsave(file.path(output_dir, paste0("peptide_counts_by_sample_", search_name, ".png")), 
         p2, width = 10, height = 6)
  cat("Created visualization: Peptide counts by sample\n")
  
  # 3. Intensity distribution
  p3 <- ggplot(filtered_data, aes(x = log10(Intensity))) +
    geom_histogram(bins = 30, fill = "darkgreen") +
    facet_wrap(~SampleID) +
    theme_minimal() +
    labs(title = "Intensity Distribution by Sample",
         x = "Log10(Intensity)", y = "Count")
  
  ggsave(file.path(output_dir, paste0("intensity_distribution_", search_name, ".png")), 
         p3, width = 10, height = 8)
  cat("Created visualization: Intensity distribution by sample\n")
  
  # 4. Heatmap of peptides across samples
  # Create a matrix of peptide intensities
  heat_data <- filtered_data %>%
    group_by(Peptide, SampleID) %>%
    summarise(Intensity = sum(Intensity, na.rm = TRUE), .groups = "drop") %>%
    mutate(log_intensity = log10(Intensity + 1)) # Add 1 to handle zeros
  
  # Only include peptides found in at least 2 samples for readability
  peptides_in_multiple_samples <- heat_data %>%
    group_by(Peptide) %>%
    summarise(sample_count = n_distinct(SampleID)) %>%
    filter(sample_count >= 2) %>%
    pull(Peptide)
  
  if(length(peptides_in_multiple_samples) > 0) {
    heat_data_filtered <- heat_data %>% 
      filter(Peptide %in% peptides_in_multiple_samples)
    
    p4 <- ggplot(heat_data_filtered, aes(x = SampleID, y = Peptide, fill = log_intensity)) +
      geom_tile() +
      scale_fill_viridis_c() +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8)) +
      labs(title = "Peptide Intensity Heatmap Across Samples",
           fill = "Log10(Intensity)")
    
    ggsave(file.path(output_dir, paste0("peptide_heatmap_", search_name, ".png")), 
           p4, width = 10, height = max(7, length(peptides_in_multiple_samples) * 0.3))
    cat("Created visualization: Peptide intensity heatmap\n")
  } else {
    cat("Not enough peptides found in multiple samples for heatmap visualization\n")
  }
}

# Function to search for partial matches across multiple fields
search_partial_gene <- function(combined_file, search_pattern) {
  cat("Searching for pattern across all relevant fields:", search_pattern, "\n")
  cat("Reading combined peptide data...\n")
  combined_data <- read_combined_data(combined_file)
  
  cat("Total entries in combined data:", nrow(combined_data), "\n")
  
  # Create output directory
  output_dir <- create_output_directory(paste0("partial_search_", search_pattern))
  
  # Search for the pattern in all relevant columns
  filtered_data <- combined_data %>%
    filter(
      grepl(search_pattern, Peptide, ignore.case = TRUE) |
        grepl(search_pattern, Gene, ignore.case = TRUE) |
        grepl(search_pattern, `Mapped Genes`, ignore.case = TRUE) |
        grepl(search_pattern, `Protein Description`, ignore.case = TRUE) |
        grepl(search_pattern, Protein, ignore.case = TRUE) |
        grepl(search_pattern, `Protein ID`, ignore.case = TRUE) |
        grepl(search_pattern, `Entry Name`, ignore.case = TRUE) |
        grepl(search_pattern, `Mapped Proteins`, ignore.case = TRUE)
    )
  
  cat("Found", nrow(filtered_data), "entries containing pattern:", search_pattern, "\n")
  
  if (nrow(filtered_data) == 0) {
    cat("No matches found for the specified pattern\n")
    return(invisible())
  }
  
  # Generate output file name
  output_file <- file.path(output_dir, paste0("partial_search_", search_pattern, ".xlsx"))
  
  # Select specific columns for the Excel output (if all columns exist)
  columns_to_include <- c(
    "Peptide", "Peptide Length", "Charges", "Probability", "Spectral Count", 
    "Intensity", "Protein", "Protein ID", "Entry Name", "Gene", 
    "Protein Description", "Mapped Genes", "Mapped Proteins", "SampleID",
    "detected_2cv", "detected_3cv", "detected_both", "CVType"
  )
  
  # Check which columns actually exist in the data
  available_columns <- names(filtered_data)
  columns_to_include <- columns_to_include[columns_to_include %in% available_columns]
  
  # Add any other important columns that may be in the data but not in our list
  other_important_cols <- c("SourceFile", "Match Type", "Start", "End", "Prev AA", "Next AA")
  for (col in other_important_cols) {
    if (col %in% available_columns && !(col %in% columns_to_include)) {
      columns_to_include <- c(columns_to_include, col)
    }
  }
  
  # Reorder columns to prioritize the most important ones
  selected_data <- filtered_data %>% select(all_of(columns_to_include))
  
  # Write to Excel file
  write_xlsx(selected_data, output_file)
  cat("Results saved to", output_file, "\n")
  
  # Show summary of matching genes and proteins
  gene_summary <- filtered_data %>%
    group_by(Gene) %>%
    summarise(
      Count = n(),
      Unique_Peptides = n_distinct(Peptide),
      .groups = "drop"
    ) %>%
    arrange(desc(Count))
  
  cat("\nMatching genes summary:\n")
  print(gene_summary)
  
  protein_summary <- filtered_data %>%
    group_by(`Protein ID`, `Entry Name`) %>%
    summarise(
      Gene = first(Gene),
      Count = n(),
      Unique_Peptides = n_distinct(Peptide),
      .groups = "drop"
    ) %>%
    arrange(desc(Count))
  
  cat("\nMatching proteins summary:\n")
  print(protein_summary)
  
  # Create enhanced peptide summary
  peptide_summary <- create_enhanced_peptide_summary(filtered_data, output_dir, paste0("partial_", search_pattern))
  
  # Create visualizations
  create_visualizations(filtered_data, peptide_summary, output_dir, paste0("partial_", search_pattern))
  
  # Create intensity summary by sample
  intensity_summary <- filtered_data %>%
    group_by(SampleID) %>%
    summarise(
      Peptide_Count = n_distinct(Peptide),
      Total_Intensity = sum(Intensity, na.rm = TRUE),
      Mean_Intensity = mean(Intensity, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_xlsx(intensity_summary, file.path(output_dir, paste0("intensity_summary_partial_", search_pattern, ".xlsx")))
  cat("Intensity summary saved to", file.path(output_dir, paste0("intensity_summary_partial_", search_pattern, ".xlsx")), "\n")
}

# Function to search for gene or protein in combined data
search_peptide_data <- function(combined_file, gene_name = NULL, protein_name = NULL) {
  cat("Searching in file:", combined_file, "\n")
  
  # Read the combined TSV file
  cat("Reading combined peptide data...\n")
  combined_data <- read_combined_data(combined_file)
  
  cat("Total entries in combined data:", nrow(combined_data), "\n")
  
  # Initialize filtered data
  filtered_data <- combined_data
  search_terms <- c()
  
  # Filter by gene name if provided
  if (!is.null(gene_name) && gene_name != "NA") {
    cat("Filtering for gene name:", gene_name, "\n")
    filtered_data <- filtered_data %>% 
      filter(Gene == gene_name)
    search_terms <- c(search_terms, gene_name)
  }
  
  # Filter by protein name if provided
  if (!is.null(protein_name) && protein_name != "NA") {
    cat("Filtering for protein name:", protein_name, "\n")
    filtered_data <- filtered_data %>% 
      filter(Protein == protein_name | `Protein ID` == protein_name | `Entry Name` == protein_name)
    search_terms <- c(search_terms, protein_name)
  }
  
  cat("Found", nrow(filtered_data), "matching entries\n")
  
  # If no matches found
  if (nrow(filtered_data) == 0) {
    cat("No matches found for the specified search criteria\n")
    return(invisible())
  }
  
  # Create a descriptive name for the search
  search_name <- paste(search_terms, collapse = "_and_")
  
  # Create output directory
  output_dir <- create_output_directory(search_name)
  
  # Generate output file name
  output_file <- file.path(output_dir, paste0("search_results_", search_name, ".xlsx"))
  
  # Select specific columns for the Excel output (if all columns exist)
  columns_to_include <- c(
    "Peptide", "Peptide Length", "Charges", "Probability", "Spectral Count", 
    "Intensity", "Protein", "Protein ID", "Entry Name", "Gene", 
    "Protein Description", "Mapped Genes", "Mapped Proteins", "SampleID",
    "detected_2cv", "detected_3cv", "detected_both", "CVType"
  )
  
  # Check which columns actually exist in the data
  available_columns <- names(filtered_data)
  columns_to_include <- columns_to_include[columns_to_include %in% available_columns]
  
  # Add any other important columns that may be in the data but not in our list
  other_important_cols <- c("SourceFile", "Match Type", "Start", "End", "Prev AA", "Next AA")
  for (col in other_important_cols) {
    if (col %in% available_columns && !(col %in% columns_to_include)) {
      columns_to_include <- c(columns_to_include, col)
    }
  }
  
  # Reorder columns to prioritize the most important ones
  selected_data <- filtered_data %>% select(all_of(columns_to_include))
  
  # Write to Excel file
  write_xlsx(selected_data, output_file)
  cat("Results saved to", output_file, "\n")
  
  # Show summary statistics
  cat("\nSummary of results:\n")
  cat("Total matches:", nrow(filtered_data), "\n")
  cat("Unique peptides:", length(unique(filtered_data$Peptide)), "\n")
  cat("Samples represented:", length(unique(filtered_data$SampleID)), "\n")
  
  # Create an intensity summary by sample
  intensity_summary <- filtered_data %>%
    group_by(SampleID) %>%
    summarise(
      Peptide_Count = n_distinct(Peptide),
      Total_Intensity = sum(Intensity, na.rm = TRUE),
      Mean_Intensity = mean(Intensity, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\nIntensity summary by sample:\n")
  print(intensity_summary)
  
  # Also save the summary to Excel
  write_xlsx(intensity_summary, file.path(output_dir, paste0("intensity_summary_", search_name, ".xlsx")))
  cat("Intensity summary saved to", file.path(output_dir, paste0("intensity_summary_", search_name, ".xlsx")), "\n")
  
  # Create enhanced peptide summary
  peptide_summary <- create_enhanced_peptide_summary(filtered_data, output_dir, search_name)
  
  # Create visualizations
  create_visualizations(filtered_data, peptide_summary, output_dir, search_name)
}

# Function to show usage information
show_usage <- function() {
  cat("Enhanced Peptide Search Tool\n")
  cat("---------------------------\n")
  cat("Usage:\n")
  cat("  1. Search by exact gene or protein name:\n")
  cat("     Rscript search_peptides.R <combined_file_path> <gene_name> [protein_name]\n")
  cat("     Example: Rscript search_peptides.R combined_peptides.tsv STRAP NA\n")
  cat("     Example: Rscript search_peptides.R combined_peptides.tsv NA Q9Y3F4\n\n")
  cat("  2. List all unique genes in the dataset:\n")
  cat("     Rscript search_peptides.R --list-genes <combined_file_path>\n")
  cat("     Example: Rscript search_peptides.R --list-genes combined_peptides.tsv\n\n")
  cat("  3. Search for partial matches across all fields:\n")
  cat("     Rscript search_peptides.R --search-partial <combined_file_path> <search_pattern>\n")
  cat("     Example: Rscript search_peptides.R --search-partial combined_peptides.tsv PHOX\n\n")
  cat("Notes:\n")
  cat("  - Use 'NA' if you don't want to specify either gene_name or protein_name\n")
  cat("  - The --search-partial option will search for the pattern across multiple fields including:\n")
  cat("    Peptide, Gene, Mapped Genes, Protein Description, Protein, Protein ID, Entry Name, and Mapped Proteins\n")
  cat("  - The --list-genes option will generate a list of all unique genes and their counts\n")
  cat("  - All outputs are now organized in a 'peptide_search' directory\n")
  cat("  - Enhanced visualization outputs include bar charts, intensity distributions, and heatmaps\n")
}

# Main function
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check if no arguments provided
  if (length(args) == 0) {
    show_usage()
    return(invisible())
  }
  
  # Handle special commands
  if (args[1] == "--list-genes") {
    if (length(args) < 2) {
      cat("Error: Please provide the combined file path\n")
      show_usage()
      return(invisible())
    }
    list_unique_genes(args[2])
    return(invisible())
  } else if (args[1] == "--search-partial") {
    if (length(args) < 3) {
      cat("Error: Please provide the combined file path and search pattern\n")
      show_usage()
      return(invisible())
    }
    search_partial_gene(args[2], args[3])
    return(invisible())
  }
  
  # Regular search mode
  combined_file <- args[1]
  gene_name <- if (length(args) >= 2) args[2] else NULL
  protein_name <- if (length(args) >= 3) args[3] else NULL
  
  # Validate input
  if ((is.null(gene_name) || gene_name == "NA") && 
      (is.null(protein_name) || protein_name == "NA")) {
    cat("Error: At least one of gene_name or protein_name must be provided\n")
    show_usage()
    return(invisible())
  }
  
  # Check if file exists
  if (!file.exists(combined_file)) {
    cat("Error: Combined file not found:", combined_file, "\n")
    return(invisible())
  }
  
  # Search for the specified gene or protein
  search_peptide_data(combined_file, gene_name, protein_name)
}

# Call the main function
main()

# # Enhanced script to search for specific gene or protein names in concatenated peptide data
# # Usage: 
# #   Rscript search_peptides.R <combined_file_path> <gene_name> [protein_name]
# #   Rscript search_peptides.R --list-genes <combined_file_path>
# #   Rscript search_peptides.R --search-partial <combined_file_path> <search_pattern>
# 
# # Complete Search Options: usage within Rstudio console
# # Now you have four ways to search your data:
# # Search by gene or protein name (exact match):
# #   rsystem("Rscript search_peptides.R combined_peptides.tsv AKAP12")
# #   rsystem("Rscript search_peptides.R combined_peptides.tsv NA Q9Y3F4-2")
# # List all genes in your dataset:
# #   rsystem("Rscript search_peptides.R --list-genes combined_peptides.tsv")
# # Search across all fields (partial match):
# #   rsystem("Rscript search_peptides.R --search-partial combined_peptides.tsv PHOX")
# # Search for exact peptide sequences:
# #   rsystem("Rscript search_peptides.R --search-peptide combined_peptides.tsv AAADFTAKV")
# 
# # Load required packages
# library(tidyverse)
# library(readr)
# library(writexl)
# 
# # Function to read the combined TSV file with proper column types
# read_combined_data <- function(combined_file) {
#   # First read the header to get column names
#   header_row <- read_tsv(combined_file, n_max = 0, show_col_types = FALSE)
#   column_names <- names(header_row)
#   
#   # Create a column specification to handle problematic columns
#   col_spec <- cols(
#     # Specify problematic columns as character
#     `Assigned Modifications` = col_character(),
#     `Observed Modifications` = col_character(),
#     # Make sure intensity is a number
#     Intensity = col_double(),
#     # For all other columns, use default type detection
#     .default = col_guess()
#   )
#   
#   # Read the combined TSV file with specified column types
#   combined_data <- read_tsv(combined_file, col_types = col_spec, show_col_types = FALSE)
#   
#   # Check for parsing problems
#   parsing_problems <- problems(combined_data)
#   if (nrow(parsing_problems) > 0) {
#     cat("\nFound", nrow(parsing_problems), "parsing issues:\n")
#     print(head(parsing_problems, 10))  # Print first 10 problems
#     
#     # Save the problems to a file for further investigation
#     write_csv(parsing_problems, "parsing_problems.csv")
#     cat("All parsing issues saved to parsing_problems.csv\n")
#   }
#   
#   return(combined_data)
# }
# 
# # Function to list all unique genes in the dataset
# list_unique_genes <- function(combined_file) {
#   cat("Reading combined peptide data...\n")
#   combined_data <- read_combined_data(combined_file)
#   
#   cat("Total entries in combined data:", nrow(combined_data), "\n")
#   
#   # Extract unique gene names
#   unique_genes <- unique(combined_data$Gene)
#   unique_genes <- unique_genes[!is.na(unique_genes)]  # Remove NA values
#   unique_genes <- sort(unique_genes)  # Sort alphabetically
#   
#   cat("Found", length(unique_genes), "unique gene names\n")
#   
#   # Save to file
#   writeLines(unique_genes, "unique_genes.txt")
#   cat("Unique gene names saved to unique_genes.txt\n")
#   
#   # Write to CSV with counts
#   gene_counts <- combined_data %>%
#     group_by(Gene) %>%
#     summarise(
#       Count = n(),
#       Unique_Peptides = n_distinct(Peptide),
#       .groups = "drop"
#     ) %>%
#     arrange(desc(Count))
#   
#   write_csv(gene_counts, "gene_counts.csv")
#   cat("Gene counts saved to gene_counts.csv\n")
#   
#   # Also create a summary of protein IDs
#   protein_counts <- combined_data %>%
#     group_by(`Protein ID`) %>%
#     summarise(
#       Count = n(),
#       Unique_Peptides = n_distinct(Peptide),
#       Gene = first(Gene),
#       .groups = "drop"
#     ) %>%
#     arrange(desc(Count))
#   
#   write_csv(protein_counts, "protein_counts.csv")
#   cat("Protein counts saved to protein_counts.csv\n")
# }
# 
# # Function to search for partial matches across multiple fields
# search_partial_gene <- function(combined_file, search_pattern) {
#   cat("Searching for pattern across all relevant fields:", search_pattern, "\n")
#   cat("Reading combined peptide data...\n")
#   combined_data <- read_combined_data(combined_file)
#   
#   cat("Total entries in combined data:", nrow(combined_data), "\n")
#   
#   # Search for the pattern in all relevant columns
#   filtered_data <- combined_data %>%
#     filter(
#       grepl(search_pattern, Peptide, ignore.case = TRUE) |
#         grepl(search_pattern, Gene, ignore.case = TRUE) |
#         grepl(search_pattern, `Mapped Genes`, ignore.case = TRUE) |
#         grepl(search_pattern, `Protein Description`, ignore.case = TRUE) |
#         grepl(search_pattern, Protein, ignore.case = TRUE) |
#         grepl(search_pattern, `Protein ID`, ignore.case = TRUE) |
#         grepl(search_pattern, `Entry Name`, ignore.case = TRUE) |
#         grepl(search_pattern, `Mapped Proteins`, ignore.case = TRUE)
#     )
#   
#   cat("Found", nrow(filtered_data), "entries containing pattern:", search_pattern, "\n")
#   
#   if (nrow(filtered_data) == 0) {
#     cat("No matches found for the specified pattern\n")
#     return(invisible())
#   }
#   
#   # Generate output file name
#   output_file <- paste0("partial_search_", search_pattern, ".xlsx")
#   
#   # Select specific columns for the Excel output (if all columns exist)
#   columns_to_include <- c(
#     "Peptide", "Peptide Length", "Charges", "Probability", "Spectral Count", 
#     "Intensity", "Protein", "Protein ID", "Entry Name", "Gene", 
#     "Protein Description", "Mapped Genes", "Mapped Proteins", "SampleID",
#     "detected_2cv", "detected_3cv", "detected_both", "CVType"
#   )
#   
#   # Check which columns actually exist in the data
#   available_columns <- names(filtered_data)
#   columns_to_include <- columns_to_include[columns_to_include %in% available_columns]
#   
#   # Add any other important columns that may be in the data but not in our list
#   other_important_cols <- c("SourceFile", "Match Type", "Start", "End", "Prev AA", "Next AA")
#   for (col in other_important_cols) {
#     if (col %in% available_columns && !(col %in% columns_to_include)) {
#       columns_to_include <- c(columns_to_include, col)
#     }
#   }
#   
#   # Reorder columns to prioritize the most important ones
#   selected_data <- filtered_data %>% select(all_of(columns_to_include))
#   
#   # Write to Excel file
#   write_xlsx(selected_data, output_file)
#   cat("Results saved to", output_file, "\n")
#   
#   # Show summary of matching genes and proteins
#   gene_summary <- filtered_data %>%
#     group_by(Gene) %>%
#     summarise(
#       Count = n(),
#       Unique_Peptides = n_distinct(Peptide),
#       .groups = "drop"
#     ) %>%
#     arrange(desc(Count))
#   
#   cat("\nMatching genes summary:\n")
#   print(gene_summary)
#   
#   protein_summary <- filtered_data %>%
#     group_by(`Protein ID`, `Entry Name`) %>%
#     summarise(
#       Gene = first(Gene),
#       Count = n(),
#       Unique_Peptides = n_distinct(Peptide),
#       .groups = "drop"
#     ) %>%
#     arrange(desc(Count))
#   
#   cat("\nMatching proteins summary:\n")
#   print(protein_summary)
#   
#   # Create intensity summary by sample
#   intensity_summary <- filtered_data %>%
#     group_by(SampleID) %>%
#     summarise(
#       Peptide_Count = n_distinct(Peptide),
#       Total_Intensity = sum(Intensity, na.rm = TRUE),
#       Mean_Intensity = mean(Intensity, na.rm = TRUE),
#       .groups = "drop"
#     )
#   
#   write_xlsx(intensity_summary, paste0("intensity_summary_partial_", search_pattern, ".xlsx"))
#   cat("Intensity summary saved to", paste0("intensity_summary_partial_", search_pattern, ".xlsx"), "\n")
# }
# 
# # Function to search for gene or protein in combined data
# search_peptide_data <- function(combined_file, gene_name = NULL, protein_name = NULL) {
#   cat("Searching in file:", combined_file, "\n")
#   
#   # Read the combined TSV file
#   cat("Reading combined peptide data...\n")
#   combined_data <- read_combined_data(combined_file)
#   
#   cat("Total entries in combined data:", nrow(combined_data), "\n")
#   
#   # Initialize filtered data
#   filtered_data <- combined_data
#   search_terms <- c()
#   
#   # Filter by gene name if provided
#   if (!is.null(gene_name) && gene_name != "NA") {
#     cat("Filtering for gene name:", gene_name, "\n")
#     filtered_data <- filtered_data %>% 
#       filter(Gene == gene_name)
#     search_terms <- c(search_terms, gene_name)
#   }
#   
#   # Filter by protein name if provided
#   if (!is.null(protein_name) && protein_name != "NA") {
#     cat("Filtering for protein name:", protein_name, "\n")
#     filtered_data <- filtered_data %>% 
#       filter(Protein == protein_name | `Protein ID` == protein_name | `Entry Name` == protein_name)
#     search_terms <- c(search_terms, protein_name)
#   }
#   
#   cat("Found", nrow(filtered_data), "matching entries\n")
#   
#   # If no matches found
#   if (nrow(filtered_data) == 0) {
#     cat("No matches found for the specified search criteria\n")
#     return(invisible())
#   }
#   
#   # Create a descriptive name for the search
#   search_name <- paste(search_terms, collapse = "_and_")
#   
#   # Generate output file name
#   output_file <- paste0("search_results_", search_name, ".xlsx")
#   
#   # Select specific columns for the Excel output (if all columns exist)
#   columns_to_include <- c(
#     "Peptide", "Peptide Length", "Charges", "Probability", "Spectral Count", 
#     "Intensity", "Protein", "Protein ID", "Entry Name", "Gene", 
#     "Protein Description", "Mapped Genes", "Mapped Proteins", "SampleID",
#     "detected_2cv", "detected_3cv", "detected_both", "CVType"
#   )
#   
#   # Check which columns actually exist in the data
#   available_columns <- names(filtered_data)
#   columns_to_include <- columns_to_include[columns_to_include %in% available_columns]
#   
#   # Add any other important columns that may be in the data but not in our list
#   other_important_cols <- c("SourceFile", "Match Type", "Start", "End", "Prev AA", "Next AA")
#   for (col in other_important_cols) {
#     if (col %in% available_columns && !(col %in% columns_to_include)) {
#       columns_to_include <- c(columns_to_include, col)
#     }
#   }
#   
#   # Reorder columns to prioritize the most important ones
#   selected_data <- filtered_data %>% select(all_of(columns_to_include))
#   
#   # Write to Excel file
#   write_xlsx(selected_data, output_file)
#   cat("Results saved to", output_file, "\n")
#   
#   # Show summary statistics
#   cat("\nSummary of results:\n")
#   cat("Total matches:", nrow(filtered_data), "\n")
#   cat("Unique peptides:", length(unique(filtered_data$Peptide)), "\n")
#   cat("Samples represented:", length(unique(filtered_data$SampleID)), "\n")
#   
#   # Create an intensity summary by sample
#   intensity_summary <- filtered_data %>%
#     group_by(SampleID) %>%
#     summarise(
#       Peptide_Count = n_distinct(Peptide),
#       Total_Intensity = sum(Intensity, na.rm = TRUE),
#       Mean_Intensity = mean(Intensity, na.rm = TRUE),
#       .groups = "drop"
#     )
#   
#   cat("\nIntensity summary by sample:\n")
#   print(intensity_summary)
#   
#   # Also save the summary to Excel
#   write_xlsx(intensity_summary, paste0("intensity_summary_", search_name, ".xlsx"))
#   cat("Intensity summary saved to", paste0("intensity_summary_", search_name, ".xlsx"), "\n")
#   
#   # Create a summary of peptide metrics
#   peptide_summary <- filtered_data %>%
#     group_by(Peptide, `Peptide Length`) %>%
#     summarise(
#       Count = n(),
#       Total_Intensity = sum(Intensity, na.rm = TRUE),
#       Mean_Intensity = mean(Intensity, na.rm = TRUE),
#       Max_Intensity = max(Intensity, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     arrange(desc(Total_Intensity))
#   
#   write_xlsx(peptide_summary, paste0("peptide_summary_", search_name, ".xlsx"))
#   cat("Peptide summary saved to", paste0("peptide_summary_", search_name, ".xlsx"), "\n")
# }
# 
# # Function to show usage information
# show_usage <- function() {
#   cat("Enhanced Peptide Search Tool\n")
#   cat("---------------------------\n")
#   cat("Usage:\n")
#   cat("  1. Search by exact gene or protein name:\n")
#   cat("     Rscript search_peptides.R <combined_file_path> <gene_name> [protein_name]\n")
#   cat("     Example: Rscript search_peptides.R combined_peptides.tsv STRAP NA\n")
#   cat("     Example: Rscript search_peptides.R combined_peptides.tsv NA Q9Y3F4\n\n")
#   cat("  2. List all unique genes in the dataset:\n")
#   cat("     Rscript search_peptides.R --list-genes <combined_file_path>\n")
#   cat("     Example: Rscript search_peptides.R --list-genes combined_peptides.tsv\n\n")
#   cat("  3. Search for partial matches across all fields:\n")
#   cat("     Rscript search_peptides.R --search-partial <combined_file_path> <search_pattern>\n")
#   cat("     Example: Rscript search_peptides.R --search-partial combined_peptides.tsv PHOX\n\n")
#   cat("Notes:\n")
#   cat("  - Use 'NA' if you don't want to specify either gene_name or protein_name\n")
#   cat("  - The --search-partial option will search for the pattern across multiple fields including:\n")
#   cat("    Peptide, Gene, Mapped Genes, Protein Description, Protein, Protein ID, Entry Name, and Mapped Proteins\n")
#   cat("  - The --list-genes option will generate a list of all unique genes and their counts\n")
# }
# 
# # Main function
# main <- function() {
#   # Parse command line arguments
#   args <- commandArgs(trailingOnly = TRUE)
#   
#   # Check if no arguments provided
#   if (length(args) == 0) {
#     show_usage()
#     return(invisible())
#   }
#   
#   # Handle special commands
#   if (args[1] == "--list-genes") {
#     if (length(args) < 2) {
#       cat("Error: Please provide the combined file path\n")
#       show_usage()
#       return(invisible())
#     }
#     list_unique_genes(args[2])
#     return(invisible())
#   } else if (args[1] == "--search-partial") {
#     if (length(args) < 3) {
#       cat("Error: Please provide the combined file path and search pattern\n")
#       show_usage()
#       return(invisible())
#     }
#     search_partial_gene(args[2], args[3])
#     return(invisible())
#   }
#   
#   # Regular search mode
#   combined_file <- args[1]
#   gene_name <- if (length(args) >= 2) args[2] else NULL
#   protein_name <- if (length(args) >= 3) args[3] else NULL
#   
#   # Validate input
#   if ((is.null(gene_name) || gene_name == "NA") && 
#       (is.null(protein_name) || protein_name == "NA")) {
#     cat("Error: At least one of gene_name or protein_name must be provided\n")
#     show_usage()
#     return(invisible())
#   }
#   
#   # Check if file exists
#   if (!file.exists(combined_file)) {
#     cat("Error: Combined file not found:", combined_file, "\n")
#     return(invisible())
#   }
#   
#   # Search for the specified gene or protein
#   search_peptide_data(combined_file, gene_name, protein_name)
# }
# 
# # Call the main function
# main()
