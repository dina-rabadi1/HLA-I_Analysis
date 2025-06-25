# ===============================
# FUSION PEPTIDE SEARCH FUNCTIONS
# ===============================
# Functions for searching DNAJB1-PRKACA fusion peptides in immunopeptidomics data
# immunopeptidomics_fusion_search_functions.R
# Author: Generated for immunopeptidomics analysis

# Define fusion protein sequence constants
FUSION_SEQUENCE <- "KREIFDRYGEEVKEFLAKAKED"
JUNCTION_POSITION <- 12  # Position where DNAJB1 ends and PRKACA begins (1-indexed)

# Function to generate all possible spanning peptides from fusion sequence
generate_spanning_peptides <- function(sequence = FUSION_SEQUENCE, junction_pos = JUNCTION_POSITION, min_length = 8, max_length = 12) {
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

# Function to prepare fusion search data based on data source
prepare_fusion_search_data <- function(data_source, combined_data = NULL, file_path = NULL) {
  
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
    
    # Create a column specification to handle problematic columns
    col_spec <- cols(
      `Assigned Modifications` = col_character(),
      `Observed Modifications` = col_character(),
      Intensity = col_double(),
      .default = col_guess()
    )
    
    # Read the TSV file
    data <- read_tsv(file_path, col_types = col_spec, show_col_types = FALSE)
    
    # Check for parsing problems
    parsing_problems <- problems(data)
    if (nrow(parsing_problems) > 0) {
      cat("Found", nrow(parsing_problems), "parsing issues (this is normal)\n")
    }
    
  } else {
    stop("ERROR: Invalid data_source. Must be 'combined_data' or 'file'")
  }
  
  # Check required columns
  required_cols <- c("Peptide", "SampleID")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("ERROR: Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  cat("Fusion search data prepared:", nrow(data), "entries from", length(unique(data$SampleID)), "samples\n")
  
  return(data)
}

# Function to perform fusion peptide search - FIXED VERSION
perform_fusion_search <- function(data, output_dir, dataset_name, timestamp, intensity_percentiles = NULL) {
  
  cat("\n=== FUSION PEPTIDE SEARCH ===\n")
  cat("Fusion sequence:", FUSION_SEQUENCE, "\n")
  cat("Junction position:", JUNCTION_POSITION, "\n")
  
  # Create analysis subdirectory
  fusion_dir <- file.path(output_dir, "fusion_search")
  plots_dir <- file.path(fusion_dir, "plots")
  tables_dir <- file.path(fusion_dir, "tables")
  
  dir.create(fusion_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Generate all possible spanning peptides
  cat("Generating all possible 8-12mer spanning peptides...\n")
  spanning_peptides <- generate_spanning_peptides()
  
  cat("Generated", length(spanning_peptides), "possible spanning peptides:\n")
  for (i in seq_along(spanning_peptides)) {
    cat(sprintf("  %2d. %s (length: %d)\n", i, spanning_peptides[i], nchar(spanning_peptides[i])))
  }
  cat("\n")
  
  # Search for these peptides in the data
  cat("Searching for fusion peptides in dataset...\n")
  found_peptides <- data %>%
    filter(Peptide %in% spanning_peptides)
  
  cat("Found", nrow(found_peptides), "entries matching fusion peptides\n")
  
  # Create theoretical peptides table
  theoretical_peptides <- data.frame(
    Theoretical_Peptide = spanning_peptides,
    Length = nchar(spanning_peptides),
    Position_in_Fusion = sapply(spanning_peptides, function(p) {
      pos <- str_locate(FUSION_SEQUENCE, fixed(p))
      pos[1]
    }),
    Spans_Junction = TRUE,
    stringsAsFactors = FALSE
  )
  
  if (nrow(found_peptides) == 0) {
    cat("No fusion peptides found in the dataset.\n")
    
    # Save theoretical peptides only
    write_xlsx(theoretical_peptides, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_01_theoretical_fusion_peptides.xlsx")))
    
    # Create summary
    summary_data <- data.frame(
      Metric = c("Fusion Sequence", "Junction Position", "Total Theoretical Peptides", "Peptides Found", "Coverage %"),
      Value = c(FUSION_SEQUENCE, JUNCTION_POSITION, length(spanning_peptides), 0, "0.0%")
    )
    
    write_csv(summary_data, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_02_fusion_search_summary.csv")))
    
    cat("✓ Theoretical fusion peptides saved (no peptides found in data)\n")
    
    return(list(
      found_peptides = data.frame(),
      theoretical_peptides = theoretical_peptides,
      summary = summary_data,
      plots_dir = plots_dir,
      tables_dir = tables_dir
    ))
  }
  
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
    "Intensity", "final_intensity", "Protein", "Protein ID", "Entry Name", "Gene", 
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
  
  # Create intensity summary by sample (use appropriate intensity column)
  intensity_col <- if ("final_intensity" %in% colnames(found_peptides)) "final_intensity" else "Intensity"
  
  intensity_summary <- found_peptides %>%
    group_by(SampleID) %>%
    summarise(
      Peptide_Count = n_distinct(Peptide),
      Total_Intensity = sum(.data[[intensity_col]], na.rm = TRUE),
      Mean_Intensity = mean(.data[[intensity_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create enhanced peptide summary
  peptide_samples <- found_peptides %>%
    group_by(Peptide, SampleID) %>%
    summarise(
      Sample_Intensity = sum(.data[[intensity_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Peptide, desc(Sample_Intensity))
  
  sample_list <- peptide_samples %>%
    group_by(Peptide) %>%
    summarise(
      Samples = paste(SampleID, collapse = ", "),
      Sample_Count = n_distinct(SampleID),
      .groups = "drop"
    )
  
  detailed_peptide_summary <- found_peptides %>%
    group_by(Peptide) %>%
    summarise(
      Count = n(),
      Length = first(`Peptide Length`),
      Total_Intensity = sum(.data[[intensity_col]], na.rm = TRUE),
      Mean_Intensity = mean(.data[[intensity_col]], na.rm = TRUE),
      Max_Intensity = max(.data[[intensity_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(sample_list, by = "Peptide") %>%
    arrange(desc(Sample_Count), desc(Total_Intensity))
  
  # Theoretical vs found comparison
  theoretical_vs_found <- theoretical_peptides %>%
    mutate(
      Found_in_Data = Theoretical_Peptide %in% found_unique_peptides,
      Entry_Count = sapply(Theoretical_Peptide, function(p) sum(found_peptides$Peptide == p))
    )
  
  # Create summary
  summary_data <- data.frame(
    Metric = c("Fusion Sequence", "Junction Position", "Total Theoretical Peptides", 
               "Peptides Found", "Total Entries Found", "Samples with Fusion Peptides", "Coverage %"),
    Value = c(FUSION_SEQUENCE, JUNCTION_POSITION, length(spanning_peptides), 
              length(found_unique_peptides), nrow(found_peptides), 
              length(unique(found_peptides$SampleID)), 
              paste0(round(length(found_unique_peptides)/length(spanning_peptides)*100, 1), "%"))
  )
  
  # Save tables
  cat("Saving fusion search tables...\n")
  write_xlsx(selected_data, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_01_fusion_peptide_results.xlsx")))
  write_xlsx(intensity_summary, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_02_fusion_intensity_summary.xlsx")))
  write_xlsx(detailed_peptide_summary, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_03_peptide_summary.xlsx")))
  write_xlsx(theoretical_vs_found, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_04_theoretical_vs_found.xlsx")))
  write_csv(summary_data, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_05_fusion_search_summary.csv")))
  
  cat("✓ Fusion search tables saved!\n")
  
  # Print summary
  cat("\n=== FUSION SEARCH SUMMARY ===\n")
  cat("Fusion sequence:", FUSION_SEQUENCE, "\n")
  cat("Total possible spanning peptides (8-12mers):", length(spanning_peptides), "\n")
  cat("Peptides found in data:", length(found_unique_peptides), "\n")
  cat("Total entries found:", nrow(found_peptides), "\n")
  cat("Samples with fusion peptides:", length(unique(found_peptides$SampleID)), "\n")
  cat("Coverage:", round(length(found_unique_peptides)/length(spanning_peptides)*100, 1), "%\n")
  
  # FIXED: Initialize results list FIRST, then add percentile analysis if needed
  results <- list(
    found_peptides = found_peptides,
    theoretical_peptides = theoretical_peptides,
    theoretical_vs_found = theoretical_vs_found,
    intensity_summary = intensity_summary,
    detailed_peptide_summary = detailed_peptide_summary,
    summary = summary_data,
    intensity_col = intensity_col,
    plots_dir = plots_dir,
    tables_dir = tables_dir
  )
  
  # Create intensity percentile analysis if percentiles are provided
  if (!is.null(intensity_percentiles)) {
    percentile_analysis <- create_fusion_intensity_percentile_analysis(
      found_peptides = found_peptides,
      intensity_percentiles = intensity_percentiles,
      intensity_col = intensity_col,
      tables_dir = tables_dir,
      timestamp = timestamp,
      dataset_name = dataset_name
    )
    
    # Add to results (now that results exists!)
    results$percentile_analysis <- percentile_analysis
  }
  
  return(results)
}

# Function to create fusion search visualizations
create_fusion_search_plots <- function(search_results, dataset_name, timestamp) {
  
  cat("Creating fusion search visualizations...\n")
  
  # Extract data from results
  found_peptides <- search_results$found_peptides
  plots_dir <- search_results$plots_dir
  intensity_col <- search_results$intensity_col
  
  if (nrow(found_peptides) == 0) {
    cat("No fusion peptides found - skipping plot creation\n")
    return(list())
  }
  
  plots <- list()
  
  # 1. Peptide intensity across samples (bar chart)
  vis_data <- found_peptides %>%
    group_by(Peptide, SampleID) %>%
    summarise(Intensity = sum(.data[[intensity_col]], na.rm = TRUE), .groups = "drop")
  
  plots$intensity_by_sample <- ggplot(vis_data, aes(x = reorder(Peptide, Intensity, sum), y = Intensity, fill = SampleID)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(title = "Fusion Peptide Intensities Across Samples",
         x = "Peptide", y = "Intensity") +
    scale_y_log10()
  
  # 2. Sample comparison showing peptide counts
  sample_counts <- found_peptides %>%
    group_by(SampleID) %>%
    summarise(
      Peptide_Count = n_distinct(Peptide),
      Total_Intensity = sum(.data[[intensity_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  plots$peptide_counts_by_sample <- ggplot(sample_counts, aes(x = reorder(SampleID, Peptide_Count), y = Peptide_Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = "Number of Fusion Peptides per Sample",
         x = "Sample ID", y = "Peptide Count")
  
  # 3. Intensity distribution
  plots$intensity_distribution <- ggplot(found_peptides, aes(x = log10(.data[[intensity_col]]))) +
    geom_histogram(bins = 20, fill = "darkgreen") +
    facet_wrap(~SampleID) +
    theme_minimal() +
    labs(title = "Fusion Peptide Intensity Distribution by Sample",
         x = paste0("Log10(", intensity_col, ")"), y = "Count")
  
  # 4. Peptide position in fusion sequence visualization
  peptide_positions <- found_peptides %>%
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
    plots$peptide_positions <- ggplot(peptide_positions, aes(x = Position_in_Fusion, y = Peptide, color = Spans_Junction)) +
      geom_point(size = 3) +
      geom_vline(xintercept = JUNCTION_POSITION - 0.5, linetype = "dashed", color = "red", size = 1) +
      theme_minimal() +
      labs(title = "Fusion Peptide Positions Relative to Junction",
           x = "Position in Fusion Sequence", 
           y = "Peptide",
           color = "Spans Junction") +
      annotate("text", x = JUNCTION_POSITION - 0.5, y = Inf, label = "Junction", 
               vjust = 2, color = "red", size = 3)
  }
  
  # Save plots
  cat("Saving fusion search plots...\n")
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_01_fusion_peptide_intensity_by_sample.png")), 
         plots$intensity_by_sample, width = 12, height = 7, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_02_fusion_peptide_counts_by_sample.png")), 
         plots$peptide_counts_by_sample, width = 10, height = 6, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_03_fusion_intensity_distribution.png")), 
         plots$intensity_distribution, width = 10, height = 8, dpi = 300)
  
  if ("peptide_positions" %in% names(plots)) {
    ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_04_fusion_peptide_positions.png")), 
           plots$peptide_positions, width = 12, height = max(6, nrow(peptide_positions) * 0.3), dpi = 300)
  }
  
  cat("✓ Fusion search plots saved!\n")
  
  return(plots)
}

# Function to create fusion peptide intensity percentile analysis
create_fusion_intensity_percentile_analysis <- function(found_peptides, intensity_percentiles, intensity_col, tables_dir, timestamp, dataset_name) {
  
  if (nrow(found_peptides) == 0) {
    cat("No fusion peptides found - skipping percentile analysis\n")
    return(NULL)
  }
  
  cat("Creating fusion peptide intensity percentile analysis...\n")
  
  # Create percentile lookup table
  percentile_values <- setNames(intensity_percentiles$Value, intensity_percentiles$Percentile)
  
  # Function to determine percentile category
  determine_percentile <- function(intensity) {
    if (is.na(intensity) || intensity == 0) return("Below 1st percentile")
    
    # Check each percentile threshold
    if (intensity >= percentile_values["99%"]) return("≥99th percentile")
    if (intensity >= percentile_values["95%"]) return("95th-99th percentile")
    if (intensity >= percentile_values["90%"]) return("90th-95th percentile")
    if (intensity >= percentile_values["75%"]) return("75th-90th percentile")
    if (intensity >= percentile_values["50%"]) return("50th-75th percentile")
    if (intensity >= percentile_values["25%"]) return("25th-50th percentile")
    if (intensity >= percentile_values["10%"]) return("10th-25th percentile")
    if (intensity >= percentile_values["5%"]) return("5th-10th percentile")
    if (intensity >= percentile_values["1%"]) return("1st-5th percentile")
    
    return("Below 1st percentile")
  }
  
  # Create detailed fusion peptide analysis with percentiles
  fusion_peptide_intensity_analysis <- found_peptides %>%
    select(Peptide, SampleID, all_of(intensity_col), 
           any_of(c("Peptide Length", "Charges", "Probability", "Spectral Count", "Protein", "Gene"))) %>%
    rename(Intensity = all_of(intensity_col)) %>%
    mutate(
      Intensity_Percentile_Category = map_chr(Intensity, determine_percentile),
      Log10_Intensity = ifelse(Intensity > 0, log10(Intensity), NA)
    ) %>%
    arrange(desc(Intensity))
  
  # Add percentile reference values for context
  percentile_reference <- data.frame(
    Percentile = names(percentile_values),
    Intensity_Threshold = percentile_values,
    stringsAsFactors = FALSE
  )
  
  # Create summary by peptide and percentile category
  peptide_percentile_summary <- fusion_peptide_intensity_analysis %>%
    group_by(Peptide, Intensity_Percentile_Category) %>%
    summarise(
      Sample_Count = n(),
      Mean_Intensity = mean(Intensity, na.rm = TRUE),
      Max_Intensity = max(Intensity, na.rm = TRUE),
      Samples = paste(SampleID, collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(Peptide, desc(Mean_Intensity))
  
  # Save the detailed analysis
  filename_detailed <- paste0(timestamp, "_", dataset_name, "_06_fusion_peptide_intensity_percentiles_detailed.csv")
  write_csv(fusion_peptide_intensity_analysis, file.path(tables_dir, filename_detailed))
  
  # Save the percentile reference
  filename_reference <- paste0(timestamp, "_", dataset_name, "_07_intensity_percentile_reference.csv")
  write_csv(percentile_reference, file.path(tables_dir, filename_reference))
  
  # Save the summary
  filename_summary <- paste0(timestamp, "_", dataset_name, "_08_fusion_peptide_percentile_summary.csv")
  write_csv(peptide_percentile_summary, file.path(tables_dir, filename_summary))
  
  cat("✓ Fusion peptide intensity percentile analysis saved!\n")
  cat("  - Detailed analysis:", filename_detailed, "\n")
  cat("  - Percentile reference:", filename_reference, "\n")
  cat("  - Summary by peptide:", filename_summary, "\n")
  
  return(list(
    detailed_analysis = fusion_peptide_intensity_analysis,
    percentile_reference = percentile_reference,
    peptide_summary = peptide_percentile_summary
  ))
}

# # ===============================
# # FUSION PEPTIDE SEARCH FUNCTIONS
# # ===============================
# # Functions for searching DNAJB1-PRKACA fusion peptides in immunopeptidomics data
# # immunopeptidomics_fusion_search_functions.R
# # Author: Generated for immunopeptidomics analysis
# 
# # Define fusion protein sequence constants
# FUSION_SEQUENCE <- "KREIFDRYGEEVKEFLAKAKED"
# JUNCTION_POSITION <- 12  # Position where DNAJB1 ends and PRKACA begins (1-indexed)
# 
# # Function to generate all possible spanning peptides from fusion sequence
# generate_spanning_peptides <- function(sequence = FUSION_SEQUENCE, junction_pos = JUNCTION_POSITION, min_length = 8, max_length = 12) {
#   spanning_peptides <- c()
#   sequence_length <- nchar(sequence)
#   
#   # Generate all possible peptides of specified lengths
#   for (length in min_length:max_length) {
#     # Generate all possible start positions
#     for (start_pos in 1:(sequence_length - length + 1)) {
#       end_pos <- start_pos + length - 1
#       
#       # Check if this peptide spans the junction
#       # It must include at least one AA before junction AND one AA after junction
#       if (start_pos < junction_pos && end_pos >= junction_pos) {
#         peptide <- substr(sequence, start_pos, end_pos)
#         spanning_peptides <- c(spanning_peptides, peptide)
#       }
#     }
#   }
#   
#   return(unique(spanning_peptides))
# }
# 
# # Function to prepare fusion search data based on data source
# prepare_fusion_search_data <- function(data_source, combined_data = NULL, file_path = NULL) {
#   
#   if (data_source == "combined_data") {
#     if (is.null(combined_data)) {
#       stop("ERROR: combined_data is NULL but data_source is set to 'combined_data'")
#     }
#     data <- combined_data
#     cat("Using combined_data from pipeline (", nrow(data), "entries)\n")
#     
#   } else if (data_source == "file") {
#     if (is.null(file_path) || !file.exists(file_path)) {
#       stop("ERROR: File path not provided or file not found: ", file_path)
#     }
#     cat("Loading data from file:", file_path, "\n")
#     
#     # Create a column specification to handle problematic columns
#     col_spec <- cols(
#       `Assigned Modifications` = col_character(),
#       `Observed Modifications` = col_character(),
#       Intensity = col_double(),
#       .default = col_guess()
#     )
#     
#     # Read the TSV file
#     data <- read_tsv(file_path, col_types = col_spec, show_col_types = FALSE)
#     
#     # Check for parsing problems
#     parsing_problems <- problems(data)
#     if (nrow(parsing_problems) > 0) {
#       cat("Found", nrow(parsing_problems), "parsing issues (this is normal)\n")
#     }
#     
#   } else {
#     stop("ERROR: Invalid data_source. Must be 'combined_data' or 'file'")
#   }
#   
#   # Check required columns
#   required_cols <- c("Peptide", "SampleID")
#   missing_cols <- setdiff(required_cols, colnames(data))
#   if (length(missing_cols) > 0) {
#     stop("ERROR: Missing required columns: ", paste(missing_cols, collapse = ", "))
#   }
#   
#   cat("Fusion search data prepared:", nrow(data), "entries from", length(unique(data$SampleID)), "samples\n")
#   
#   return(data)
# }
# 
# # Function to perform fusion peptide search
# perform_fusion_search <- function(data, output_dir, dataset_name, timestamp, intensity_percentiles = NULL) {
#   
#   cat("\n=== FUSION PEPTIDE SEARCH ===\n")
#   cat("Fusion sequence:", FUSION_SEQUENCE, "\n")
#   cat("Junction position:", JUNCTION_POSITION, "\n")
#   
#   # Create analysis subdirectory
#   fusion_dir <- file.path(output_dir, "fusion_search")
#   plots_dir <- file.path(fusion_dir, "plots")
#   tables_dir <- file.path(fusion_dir, "tables")
#   
#   dir.create(fusion_dir, recursive = TRUE, showWarnings = FALSE)
#   dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
#   dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
#   
#   # Generate all possible spanning peptides
#   cat("Generating all possible 8-12mer spanning peptides...\n")
#   spanning_peptides <- generate_spanning_peptides()
#   
#   cat("Generated", length(spanning_peptides), "possible spanning peptides:\n")
#   for (i in seq_along(spanning_peptides)) {
#     cat(sprintf("  %2d. %s (length: %d)\n", i, spanning_peptides[i], nchar(spanning_peptides[i])))
#   }
#   cat("\n")
#   
#   # Search for these peptides in the data
#   cat("Searching for fusion peptides in dataset...\n")
#   found_peptides <- data %>%
#     filter(Peptide %in% spanning_peptides)
#   
#   cat("Found", nrow(found_peptides), "entries matching fusion peptides\n")
#   
#   # Create theoretical peptides table
#   theoretical_peptides <- data.frame(
#     Theoretical_Peptide = spanning_peptides,
#     Length = nchar(spanning_peptides),
#     Position_in_Fusion = sapply(spanning_peptides, function(p) {
#       pos <- str_locate(FUSION_SEQUENCE, fixed(p))
#       pos[1]
#     }),
#     Spans_Junction = TRUE,
#     stringsAsFactors = FALSE
#   )
#   
#   if (nrow(found_peptides) == 0) {
#     cat("No fusion peptides found in the dataset.\n")
#     
#     # Save theoretical peptides only
#     write_xlsx(theoretical_peptides, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_01_theoretical_fusion_peptides.xlsx")))
#     
#     # Create summary
#     summary_data <- data.frame(
#       Metric = c("Fusion Sequence", "Junction Position", "Total Theoretical Peptides", "Peptides Found", "Coverage %"),
#       Value = c(FUSION_SEQUENCE, JUNCTION_POSITION, length(spanning_peptides), 0, "0.0%")
#     )
#     
#     write_csv(summary_data, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_02_fusion_search_summary.csv")))
#     
#     cat("✓ Theoretical fusion peptides saved (no peptides found in data)\n")
#     
#     return(list(
#       found_peptides = data.frame(),
#       theoretical_peptides = theoretical_peptides,
#       summary = summary_data,
#       plots_dir = plots_dir,
#       tables_dir = tables_dir
#     ))
#   }
#   
#   # Show which specific peptides were found
#   found_unique_peptides <- unique(found_peptides$Peptide)
#   cat("Found the following fusion peptides:\n")
#   for (i in seq_along(found_unique_peptides)) {
#     peptide <- found_unique_peptides[i]
#     count <- sum(found_peptides$Peptide == peptide)
#     cat(sprintf("  %s (found in %d entries)\n", peptide, count))
#   }
#   cat("\n")
#   
#   # Select columns for output
#   columns_to_include <- c(
#     "Peptide", "Peptide Length", "Charges", "Probability", "Spectral Count", 
#     "Intensity", "final_intensity", "Protein", "Protein ID", "Entry Name", "Gene", 
#     "Protein Description", "Mapped Genes", "Mapped Proteins", "SampleID",
#     "detected_2cv", "detected_3cv", "detected_both", "CVType"
#   )
#   
#   # Check which columns actually exist
#   available_columns <- names(found_peptides)
#   columns_to_include <- columns_to_include[columns_to_include %in% available_columns]
#   
#   # Add other important columns
#   other_important_cols <- c("SourceFile", "Match Type", "Start", "End", "Prev AA", "Next AA")
#   for (col in other_important_cols) {
#     if (col %in% available_columns && !(col %in% columns_to_include)) {
#       columns_to_include <- c(columns_to_include, col)
#     }
#   }
#   
#   # Create selected data
#   selected_data <- found_peptides %>% select(all_of(columns_to_include))
#   
#   # Create intensity summary by sample (use appropriate intensity column)
#   intensity_col <- if ("final_intensity" %in% colnames(found_peptides)) "final_intensity" else "Intensity"
#   
#   intensity_summary <- found_peptides %>%
#     group_by(SampleID) %>%
#     summarise(
#       Peptide_Count = n_distinct(Peptide),
#       Total_Intensity = sum(.data[[intensity_col]], na.rm = TRUE),
#       Mean_Intensity = mean(.data[[intensity_col]], na.rm = TRUE),
#       .groups = "drop"
#     )
#   
#   # Create enhanced peptide summary
#   peptide_samples <- found_peptides %>%
#     group_by(Peptide, SampleID) %>%
#     summarise(
#       Sample_Intensity = sum(.data[[intensity_col]], na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     arrange(Peptide, desc(Sample_Intensity))
#   
#   sample_list <- peptide_samples %>%
#     group_by(Peptide) %>%
#     summarise(
#       Samples = paste(SampleID, collapse = ", "),
#       Sample_Count = n_distinct(SampleID),
#       .groups = "drop"
#     )
#   
#   detailed_peptide_summary <- found_peptides %>%
#     group_by(Peptide) %>%
#     summarise(
#       Count = n(),
#       Length = first(`Peptide Length`),
#       Total_Intensity = sum(.data[[intensity_col]], na.rm = TRUE),
#       Mean_Intensity = mean(.data[[intensity_col]], na.rm = TRUE),
#       Max_Intensity = max(.data[[intensity_col]], na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     left_join(sample_list, by = "Peptide") %>%
#     arrange(desc(Sample_Count), desc(Total_Intensity))
#   
#   # Theoretical vs found comparison
#   theoretical_vs_found <- theoretical_peptides %>%
#     mutate(
#       Found_in_Data = Theoretical_Peptide %in% found_unique_peptides,
#       Entry_Count = sapply(Theoretical_Peptide, function(p) sum(found_peptides$Peptide == p))
#     )
#   
#   # Create summary
#   summary_data <- data.frame(
#     Metric = c("Fusion Sequence", "Junction Position", "Total Theoretical Peptides", 
#                "Peptides Found", "Total Entries Found", "Samples with Fusion Peptides", "Coverage %"),
#     Value = c(FUSION_SEQUENCE, JUNCTION_POSITION, length(spanning_peptides), 
#               length(found_unique_peptides), nrow(found_peptides), 
#               length(unique(found_peptides$SampleID)), 
#               paste0(round(length(found_unique_peptides)/length(spanning_peptides)*100, 1), "%"))
#   )
#   
#   #new
#   ## Function to create fusion peptide intensity percentile analysis
#   create_fusion_intensity_percentile_analysis <- function(found_peptides, intensity_percentiles, intensity_col, tables_dir, timestamp, dataset_name) {
#     
#     if (nrow(found_peptides) == 0) {
#       cat("No fusion peptides found - skipping percentile analysis\n")
#       return(NULL)
#     }
#     
#     cat("Creating fusion peptide intensity percentile analysis...\n")
#     
#     # Create percentile lookup table
#     percentile_values <- setNames(intensity_percentiles$Value, intensity_percentiles$Percentile)
#     
#     # Function to determine percentile category
#     determine_percentile <- function(intensity) {
#       if (is.na(intensity) || intensity == 0) return("Below 1st percentile")
#       
#       # Check each percentile threshold
#       if (intensity >= percentile_values["99%"]) return("≥99th percentile")
#       if (intensity >= percentile_values["95%"]) return("95th-99th percentile")
#       if (intensity >= percentile_values["90%"]) return("90th-95th percentile")
#       if (intensity >= percentile_values["75%"]) return("75th-90th percentile")
#       if (intensity >= percentile_values["50%"]) return("50th-75th percentile")
#       if (intensity >= percentile_values["25%"]) return("25th-50th percentile")
#       if (intensity >= percentile_values["10%"]) return("10th-25th percentile")
#       if (intensity >= percentile_values["5%"]) return("5th-10th percentile")
#       if (intensity >= percentile_values["1%"]) return("1st-5th percentile")
#       
#       return("Below 1st percentile")
#     }
#     
#     # Create detailed fusion peptide analysis with percentiles
#     fusion_peptide_intensity_analysis <- found_peptides %>%
#       select(Peptide, SampleID, all_of(intensity_col), 
#              any_of(c("Peptide Length", "Charges", "Probability", "Spectral Count", "Protein", "Gene"))) %>%
#       rename(Intensity = all_of(intensity_col)) %>%
#       mutate(
#         Intensity_Percentile_Category = map_chr(Intensity, determine_percentile),
#         Log10_Intensity = ifelse(Intensity > 0, log10(Intensity), NA)
#       ) %>%
#       arrange(desc(Intensity))
#     
#     # Add percentile reference values for context
#     percentile_reference <- data.frame(
#       Percentile = names(percentile_values),
#       Intensity_Threshold = percentile_values,
#       stringsAsFactors = FALSE
#     )
#     
#     # Create summary by peptide and percentile category
#     peptide_percentile_summary <- fusion_peptide_intensity_analysis %>%
#       group_by(Peptide, Intensity_Percentile_Category) %>%
#       summarise(
#         Sample_Count = n(),
#         Mean_Intensity = mean(Intensity, na.rm = TRUE),
#         Max_Intensity = max(Intensity, na.rm = TRUE),
#         Samples = paste(SampleID, collapse = ", "),
#         .groups = "drop"
#       ) %>%
#       arrange(Peptide, desc(Mean_Intensity))
#     
#     # Save the detailed analysis
#     filename_detailed <- paste0(timestamp, "_", dataset_name, "_06_fusion_peptide_intensity_percentiles_detailed.csv")
#     write_csv(fusion_peptide_intensity_analysis, file.path(tables_dir, filename_detailed))
#     
#     # Save the percentile reference
#     filename_reference <- paste0(timestamp, "_", dataset_name, "_07_intensity_percentile_reference.csv")
#     write_csv(percentile_reference, file.path(tables_dir, filename_reference))
#     
#     # Save the summary
#     filename_summary <- paste0(timestamp, "_", dataset_name, "_08_fusion_peptide_percentile_summary.csv")
#     write_csv(peptide_percentile_summary, file.path(tables_dir, filename_summary))
#     
#     cat("✓ Fusion peptide intensity percentile analysis saved!\n")
#     cat("  - Detailed analysis:", filename_detailed, "\n")
#     cat("  - Percentile reference:", filename_reference, "\n")
#     cat("  - Summary by peptide:", filename_summary, "\n")
#     
#     
#     #new
#     ## Add this to your perform_fusion_search function before the final return(results):
#     # Create intensity percentile analysis (requires intensity_percentiles from main analysis)
#     if (exists("intensity_percentiles")) {
#       percentile_analysis <- create_fusion_intensity_percentile_analysis(
#         found_peptides = found_peptides,
#         intensity_percentiles = intensity_percentiles,
#         intensity_col = intensity_col,
#         tables_dir = tables_dir,
#         timestamp = timestamp,
#         dataset_name = dataset_name
#       )
#       
#       # Add to results
#       results$percentile_analysis <- percentile_analysis
#     }
#     
#     # Function to create fusion peptide intensity percentile analysis
#     create_fusion_intensity_percentile_analysis <- function(found_peptides, intensity_percentiles, intensity_col, tables_dir, timestamp, dataset_name) {
#       
#       if (nrow(found_peptides) == 0) {
#         cat("No fusion peptides found - skipping percentile analysis\n")
#         return(NULL)
#       }
#       
#       cat("Creating fusion peptide intensity percentile analysis...\n")
#       
#       # Create percentile lookup table
#       percentile_values <- setNames(intensity_percentiles$Value, intensity_percentiles$Percentile)
#       
#       # Function to determine percentile category
#       determine_percentile <- function(intensity) {
#         if (is.na(intensity) || intensity == 0) return("Below 1st percentile")
#         
#         # Check each percentile threshold
#         if (intensity >= percentile_values["99%"]) return("≥99th percentile")
#         if (intensity >= percentile_values["95%"]) return("95th-99th percentile")
#         if (intensity >= percentile_values["90%"]) return("90th-95th percentile")
#         if (intensity >= percentile_values["75%"]) return("75th-90th percentile")
#         if (intensity >= percentile_values["50%"]) return("50th-75th percentile")
#         if (intensity >= percentile_values["25%"]) return("25th-50th percentile")
#         if (intensity >= percentile_values["10%"]) return("10th-25th percentile")
#         if (intensity >= percentile_values["5%"]) return("5th-10th percentile")
#         if (intensity >= percentile_values["1%"]) return("1st-5th percentile")
#         
#         return("Below 1st percentile")
#       }
#       
#       # Create detailed fusion peptide analysis with percentiles
#       fusion_peptide_intensity_analysis <- found_peptides %>%
#         select(Peptide, SampleID, all_of(intensity_col), 
#                any_of(c("Peptide Length", "Charges", "Probability", "Spectral Count", "Protein", "Gene"))) %>%
#         rename(Intensity = all_of(intensity_col)) %>%
#         mutate(
#           Intensity_Percentile_Category = map_chr(Intensity, determine_percentile),
#           Log10_Intensity = ifelse(Intensity > 0, log10(Intensity), NA)
#         ) %>%
#         arrange(desc(Intensity))
#       
#       # Add percentile reference values for context
#       percentile_reference <- data.frame(
#         Percentile = names(percentile_values),
#         Intensity_Threshold = percentile_values,
#         stringsAsFactors = FALSE
#       )
#       
#       # Create summary by peptide and percentile category
#       peptide_percentile_summary <- fusion_peptide_intensity_analysis %>%
#         group_by(Peptide, Intensity_Percentile_Category) %>%
#         summarise(
#           Sample_Count = n(),
#           Mean_Intensity = mean(Intensity, na.rm = TRUE),
#           Max_Intensity = max(Intensity, na.rm = TRUE),
#           Samples = paste(SampleID, collapse = ", "),
#           .groups = "drop"
#         ) %>%
#         arrange(Peptide, desc(Mean_Intensity))
#       
#       # Save the detailed analysis
#       filename_detailed <- paste0(timestamp, "_", dataset_name, "_06_fusion_peptide_intensity_percentiles_detailed.csv")
#       write_csv(fusion_peptide_intensity_analysis, file.path(tables_dir, filename_detailed))
#       
#       # Save the percentile reference
#       filename_reference <- paste0(timestamp, "_", dataset_name, "_07_intensity_percentile_reference.csv")
#       write_csv(percentile_reference, file.path(tables_dir, filename_reference))
#       
#       # Save the summary
#       filename_summary <- paste0(timestamp, "_", dataset_name, "_08_fusion_peptide_percentile_summary.csv")
#       write_csv(peptide_percentile_summary, file.path(tables_dir, filename_summary))
#       
#       cat("✓ Fusion peptide intensity percentile analysis saved!\n")
#       cat("  - Detailed analysis:", filename_detailed, "\n")
#       cat("  - Percentile reference:", filename_reference, "\n")
#       cat("  - Summary by peptide:", filename_summary, "\n")
#       
#       return(list(
#         detailed_analysis = fusion_peptide_intensity_analysis,
#         percentile_reference = percentile_reference,
#         peptide_summary = peptide_percentile_summary
#       ))
#     }
#     
#     return(list(
#       detailed_analysis = fusion_peptide_intensity_analysis,
#       percentile_reference = percentile_reference,
#       peptide_summary = peptide_percentile_summary
#     ))
#   }
#   
#   # Save tables
#   cat("Saving fusion search tables...\n")
#   write_xlsx(selected_data, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_01_fusion_peptide_results.xlsx")))
#   write_xlsx(intensity_summary, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_02_fusion_intensity_summary.xlsx")))
#   write_xlsx(detailed_peptide_summary, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_03_peptide_summary.xlsx")))
#   write_xlsx(theoretical_vs_found, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_04_theoretical_vs_found.xlsx")))
#   write_csv(summary_data, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_05_fusion_search_summary.csv")))
#   
#   cat("✓ Fusion search tables saved!\n")
#   
#   # Print summary
#   cat("\n=== FUSION SEARCH SUMMARY ===\n")
#   cat("Fusion sequence:", FUSION_SEQUENCE, "\n")
#   cat("Total possible spanning peptides (8-12mers):", length(spanning_peptides), "\n")
#   cat("Peptides found in data:", length(found_unique_peptides), "\n")
#   cat("Total entries found:", nrow(found_peptides), "\n")
#   cat("Samples with fusion peptides:", length(unique(found_peptides$SampleID)), "\n")
#   cat("Coverage:", round(length(found_unique_peptides)/length(spanning_peptides)*100, 1), "%\n")
#   
#   # Return results for plotting
#   results <- list(
#     found_peptides = found_peptides,
#     theoretical_peptides = theoretical_peptides,
#     theoretical_vs_found = theoretical_vs_found,
#     intensity_summary = intensity_summary,
#     detailed_peptide_summary = detailed_peptide_summary,
#     summary = summary_data,
#     intensity_col = intensity_col,
#     plots_dir = plots_dir,
#     tables_dir = tables_dir
#   )
#   
#   return(results)
# }
# 
# # Function to create fusion search visualizations
# create_fusion_search_plots <- function(search_results, dataset_name, timestamp) {
#   
#   cat("Creating fusion search visualizations...\n")
#   
#   # Extract data from results
#   found_peptides <- search_results$found_peptides
#   plots_dir <- search_results$plots_dir
#   intensity_col <- search_results$intensity_col
#   
#   if (nrow(found_peptides) == 0) {
#     cat("No fusion peptides found - skipping plot creation\n")
#     return(list())
#   }
#   
#   plots <- list()
#   
#   # 1. Peptide intensity across samples (bar chart)
#   vis_data <- found_peptides %>%
#     group_by(Peptide, SampleID) %>%
#     summarise(Intensity = sum(.data[[intensity_col]], na.rm = TRUE), .groups = "drop")
#   
#   plots$intensity_by_sample <- ggplot(vis_data, aes(x = reorder(Peptide, Intensity, sum), y = Intensity, fill = SampleID)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
#     labs(title = "Fusion Peptide Intensities Across Samples",
#          x = "Peptide", y = "Intensity") +
#     scale_y_log10()
#   
#   # 2. Sample comparison showing peptide counts
#   sample_counts <- found_peptides %>%
#     group_by(SampleID) %>%
#     summarise(
#       Peptide_Count = n_distinct(Peptide),
#       Total_Intensity = sum(.data[[intensity_col]], na.rm = TRUE),
#       .groups = "drop"
#     )
#   
#   plots$peptide_counts_by_sample <- ggplot(sample_counts, aes(x = reorder(SampleID, Peptide_Count), y = Peptide_Count)) +
#     geom_bar(stat = "identity", fill = "steelblue") +
#     theme_minimal() +
#     labs(title = "Number of Fusion Peptides per Sample",
#          x = "Sample ID", y = "Peptide Count")
#   
#   # 3. Intensity distribution
#   plots$intensity_distribution <- ggplot(found_peptides, aes(x = log10(.data[[intensity_col]]))) +
#     geom_histogram(bins = 20, fill = "darkgreen") +
#     facet_wrap(~SampleID) +
#     theme_minimal() +
#     labs(title = "Fusion Peptide Intensity Distribution by Sample",
#          x = paste0("Log10(", intensity_col, ")"), y = "Count")
#   
#   # 4. Peptide position in fusion sequence visualization
#   peptide_positions <- found_peptides %>%
#     distinct(Peptide) %>%
#     mutate(
#       Position_in_Fusion = map_dbl(Peptide, ~{
#         pos <- str_locate(FUSION_SEQUENCE, fixed(.x))
#         if (!is.na(pos[1])) pos[1] else NA
#       }),
#       Spans_Junction = map_lgl(Peptide, ~{
#         pos <- str_locate(FUSION_SEQUENCE, fixed(.x))
#         if (!is.na(pos[1])) {
#           start_pos <- pos[1]
#           end_pos <- pos[2]
#           start_pos < JUNCTION_POSITION && end_pos >= JUNCTION_POSITION
#         } else FALSE
#       })
#     ) %>%
#     filter(!is.na(Position_in_Fusion))
#   
#   if (nrow(peptide_positions) > 0) {
#     plots$peptide_positions <- ggplot(peptide_positions, aes(x = Position_in_Fusion, y = Peptide, color = Spans_Junction)) +
#       geom_point(size = 3) +
#       geom_vline(xintercept = JUNCTION_POSITION - 0.5, linetype = "dashed", color = "red", size = 1) +
#       theme_minimal() +
#       labs(title = "Fusion Peptide Positions Relative to Junction",
#            x = "Position in Fusion Sequence", 
#            y = "Peptide",
#            color = "Spans Junction") +
#       annotate("text", x = JUNCTION_POSITION - 0.5, y = Inf, label = "Junction", 
#                vjust = 2, color = "red", size = 3)
#   }
#   
#   # Save plots
#   cat("Saving fusion search plots...\n")
#   
#   ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_01_fusion_peptide_intensity_by_sample.png")), 
#          plots$intensity_by_sample, width = 12, height = 7, dpi = 300)
#   
#   ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_02_fusion_peptide_counts_by_sample.png")), 
#          plots$peptide_counts_by_sample, width = 10, height = 6, dpi = 300)
#   
#   ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_03_fusion_intensity_distribution.png")), 
#          plots$intensity_distribution, width = 10, height = 8, dpi = 300)
#   
#   if ("peptide_positions" %in% names(plots)) {
#     ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_04_fusion_peptide_positions.png")), 
#            plots$peptide_positions, width = 12, height = max(6, nrow(peptide_positions) * 0.3), dpi = 300)
#   }
#   
#   cat("✓ Fusion search plots saved!\n")
#   
#   return(plots)
# }