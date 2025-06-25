# Multi-Omics Neoantigen Analysis Function - FIXED VERSION
# Created to be sourced and called from RMD files
# Identifies potential targetable neoantigens from immunopeptidome data

# Load required libraries
if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if (!require(readxl)) install.packages("readxl"); library(readxl)
if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if (!require(pheatmap)) install.packages("pheatmap"); library(pheatmap)
if (!require(VennDiagram)) install.packages("VennDiagram"); library(VennDiagram)
if (!require(RColorBrewer)) install.packages("RColorBrewer"); library(RColorBrewer)
if (!require(gridExtra)) install.packages("gridExtra"); library(gridExtra)
if (!require(rlang)) install.packages("rlang"); library(rlang)

#================================================
# HELPER FUNCTIONS
#================================================

# Function to process gene information for multi-gene peptides
process_gene_info <- function(gene_string) {
  if (is.na(gene_string) || gene_string == "" || gene_string == "NA") {
    return(list(
      primary_gene = NA,
      all_genes = NA,
      is_multi_gene = FALSE,
      gene_count = 0
    ))
  }
  
  genes <- trimws(unlist(strsplit(as.character(gene_string), ";")))
  genes <- genes[genes != "" & genes != "NA"]
  
  if (length(genes) == 0) {
    return(list(
      primary_gene = NA,
      all_genes = NA,
      is_multi_gene = FALSE,
      gene_count = 0
    ))
  }
  
  return(list(
    primary_gene = genes[1],
    all_genes = paste(genes, collapse = "; "),
    is_multi_gene = length(genes) > 1,
    gene_count = length(genes)
  ))
}

# Function to determine expression status
get_expression_status <- function(log2fc, up_thresh = 1, down_thresh = -1) {
  if (is.na(log2fc)) return("MISSING")
  if (log2fc > up_thresh) return("UP")
  if (log2fc < down_thresh) return("DOWN")
  return("NO_CHANGE")
}

# Function to combine LFQ and TMT proteome data
combine_proteome_data <- function(lfq_fc, tmt_fc) {
  lfq_status <- get_expression_status(lfq_fc)
  tmt_status <- get_expression_status(tmt_fc)
  
  # If both missing
  if (lfq_status == "MISSING" && tmt_status == "MISSING") {
    return(list(combined_fc = NA, combined_status = "MISSING", agreement = "BOTH_MISSING"))
  }
  
  # If one missing, use the other
  if (lfq_status == "MISSING") {
    return(list(combined_fc = tmt_fc, combined_status = tmt_status, agreement = "TMT_ONLY"))
  }
  if (tmt_status == "MISSING") {
    return(list(combined_fc = lfq_fc, combined_status = lfq_status, agreement = "LFQ_ONLY"))
  }
  
  # If both present, check agreement
  if (lfq_status == tmt_status) {
    # Agreement - take average
    avg_fc <- mean(c(lfq_fc, tmt_fc), na.rm = TRUE)
    return(list(combined_fc = avg_fc, combined_status = lfq_status, agreement = "AGREE"))
  } else {
    # Disagreement - flag as conflicting, use the one with larger magnitude
    if (abs(lfq_fc) > abs(tmt_fc)) {
      return(list(combined_fc = lfq_fc, combined_status = paste0("CONFLICT_", lfq_status), agreement = "CONFLICT_LFQ"))
    } else {
      return(list(combined_fc = tmt_fc, combined_status = paste0("CONFLICT_", tmt_status), agreement = "CONFLICT_TMT"))
    }
  }
}

# Function to assign tiers
assign_tier <- function(transcriptome_status, proteome_status) {
  tier_map <- list(
    "UP_UP" = "Tier_1_Trans_UP_Prot_UP",
    "UP_DOWN" = "Tier_2_Trans_UP_Prot_DOWN", 
    "DOWN_UP" = "Tier_3_Trans_DOWN_Prot_UP",
    "UP_NO_CHANGE" = "Tier_4_Trans_UP_Prot_NOCHANGE",
    "NO_CHANGE_UP" = "Tier_5_Trans_NOCHANGE_Prot_UP",
    "DOWN_NO_CHANGE" = "Tier_6_Trans_DOWN_Prot_NOCHANGE",
    "NO_CHANGE_DOWN" = "Tier_7_Trans_NOCHANGE_Prot_DOWN",
    "NO_CHANGE_NO_CHANGE" = "Tier_8_Trans_NOCHANGE_Prot_NOCHANGE",
    "DOWN_DOWN" = "Tier_9_Trans_DOWN_Prot_DOWN",
    "UP_MISSING" = "Tier_10_Trans_UP_Prot_MISSING",
    "DOWN_MISSING" = "Tier_11_Trans_DOWN_Prot_MISSING",
    "NO_CHANGE_MISSING" = "Tier_12_Trans_NOCHANGE_Prot_MISSING",
    "MISSING_UP" = "Tier_13_Trans_MISSING_Prot_UP",
    "MISSING_DOWN" = "Tier_13_Trans_MISSING_Prot_DOWN",
    "MISSING_NO_CHANGE" = "Tier_13_Trans_MISSING_Prot_NOCHANGE",
    "MISSING_MISSING" = "Tier_14_Trans_MISSING_Prot_MISSING"
  )
  
  # Handle conflict cases
  if (grepl("CONFLICT", proteome_status)) {
    base_prot_status <- gsub("CONFLICT_", "", proteome_status)
    key <- paste(transcriptome_status, base_prot_status, sep = "_")
    tier <- tier_map[[key]]
    if (!is.null(tier)) {
      return(paste0(tier, "_CONFLICT"))
    }
  }
  
  key <- paste(transcriptome_status, proteome_status, sep = "_")
  tier <- tier_map[[key]]
  
  if (is.null(tier)) {
    return("Tier_Unknown")
  }
  
  return(tier)
}

#' Extract RU148 peptide ratios from analysis results
#' 
#' Reads RU148 analysis CSV files and extracts log2_fold_change_immuno values
#' for use in neoantigen filtering
#' 
#' @param ru148_output_dir Path to RU148 analysis output directory
#' @param verbose Print detailed information
#' 
#' @return Data frame with Peptide and log2_fold_change_immuno columns
extract_ru148_peptide_ratios <- function(ru148_output_dir, verbose = TRUE) {
  
  if (verbose) {
    cat("=== EXTRACTING RU148 PEPTIDE RATIOS ===\n")
    cat("RU148 output directory:", ru148_output_dir, "\n")
  }
  
  # Check if directory exists
  if (!dir.exists(ru148_output_dir)) {
    if (verbose) cat("❌ RU148 output directory not found\n")
    return(data.frame(Peptide = character(0), log2_fold_change_immuno = numeric(0)))
  }
  
  # Look for RU148 CSV files in tables subdirectory
  tables_dir <- file.path(ru148_output_dir, "tables")
  
  if (!dir.exists(tables_dir)) {
    if (verbose) cat("❌ RU148 tables directory not found\n")
    return(data.frame(Peptide = character(0), log2_fold_change_immuno = numeric(0)))
  }
  
  # Find tumor_upregulated and tumor_exclusive CSV files
  csv_files <- list.files(tables_dir, pattern = "tumor_(upregulated|exclusive).*combined\\.csv$", 
                          full.names = TRUE, recursive = TRUE)
  
  if (length(csv_files) == 0) {
    if (verbose) cat("❌ No RU148 tumor CSV files found\n")
    return(data.frame(Peptide = character(0), log2_fold_change_immuno = numeric(0)))
  }
  
  if (verbose) {
    cat("Found RU148 CSV files:\n")
    for (file in csv_files) {
      cat("  -", basename(file), "\n")
    }
  }
  
  # Read and combine all RU148 files
  all_ru148_data <- data.frame()
  
  for (file in csv_files) {
    if (verbose) cat("Reading:", basename(file), "\n")
    
    tryCatch({
      # Read the CSV file
      file_data <- read_csv(file, show_col_types = FALSE)
      
      # Check for required columns
      if (!"Peptide" %in% colnames(file_data)) {
        if (verbose) cat("  ⚠ No 'Peptide' column found, skipping\n")
        next
      }
      
      if (!"log2_fold_change_immuno" %in% colnames(file_data)) {
        if (verbose) cat("  ⚠ No 'log2_fold_change_immuno' column found, skipping\n")
        next
      }
      
      # Extract relevant columns
      peptide_ratios <- file_data %>%
        select(Peptide, log2_fold_change_immuno) %>%
        filter(!is.na(log2_fold_change_immuno))
      
      if (verbose) cat("  ✓ Extracted", nrow(peptide_ratios), "peptide ratios\n")
      
      # Add to combined data
      all_ru148_data <- bind_rows(all_ru148_data, peptide_ratios)
      
    }, error = function(e) {
      if (verbose) cat("  ❌ Error reading file:", e$message, "\n")
    })
  }
  
  # Remove duplicates (keep the first occurrence)
  if (nrow(all_ru148_data) > 0) {
    all_ru148_data <- all_ru148_data %>%
      distinct(Peptide, .keep_all = TRUE)
    
    if (verbose) {
      cat("\n✓ Combined RU148 data:\n")
      cat("  - Total unique peptides:", nrow(all_ru148_data), "\n")
      cat("  - log2FC range:", round(range(all_ru148_data$log2_fold_change_immuno, na.rm = TRUE), 2), "\n")
      cat("  - Mean log2FC:", round(mean(all_ru148_data$log2_fold_change_immuno, na.rm = TRUE), 2), "\n")
    }
  } else {
    if (verbose) cat("❌ No valid RU148 data extracted\n")
  }
  
  return(all_ru148_data)
}

#' Apply Tumor vs Normal Immunopeptidome Filter - FIXED VERSION
#' 
#' Filters peptides based on tumor/normal tissue expression ratios to reduce
#' potential off-target effects on healthy tissue. KEEPS peptides not found
#' in the specific tumor/normal pair for evaluation in other samples.
#' 
#' @param data Combined immunopeptidome data
#' @param normal_sample_id Sample ID for normal tissue (e.g., "148N")
#' @param tumor_sample_id Sample ID for tumor tissue (e.g., "148T")
#' @param log2fc_threshold Minimum log2 fold change (tumor/normal) to keep peptide
#' @param keep_normal_exclusive Keep peptides found only in normal tissue
#' @param keep_tumor_exclusive Keep peptides found only in tumor tissue
#' @param verbose Print detailed filtering information
#' 
#' @return List containing filtered data and filtering statistics
apply_tumor_normal_immunopeptidome_filter <- function(
    data,
    normal_sample_id,
    tumor_sample_id, 
    log2fc_threshold = 1.0,
    keep_normal_exclusive = FALSE,
    keep_tumor_exclusive = TRUE,
    verbose = TRUE
) {
  
  if (verbose) {
    cat("\n=== IMMUNOPEPTIDOME NORMAL TISSUE FILTERING ===\n")
    cat("Normal sample ID:", normal_sample_id, "\n")
    cat("Tumor sample ID:", tumor_sample_id, "\n")
    cat("Log2FC threshold:", log2fc_threshold, "(", round(2^log2fc_threshold, 2), "-fold)\n")
    cat("Keep normal-exclusive peptides:", keep_normal_exclusive, "\n")
    cat("Keep tumor-exclusive peptides:", keep_tumor_exclusive, "\n")
    cat("Strategy: Keep peptides not found in T/N pair for evaluation in other samples\n\n")
  }
  
  # Check if required samples exist in data
  available_samples <- unique(data$SampleID)
  normal_exists <- normal_sample_id %in% available_samples
  tumor_exists <- tumor_sample_id %in% available_samples
  
  if (!normal_exists && !tumor_exists) {
    if (verbose) {
      cat("WARNING: Neither normal (", normal_sample_id, ") nor tumor (", tumor_sample_id, ") samples found in data.\n")
      cat("Available samples:", paste(available_samples, collapse = ", "), "\n")
      cat("Skipping immunopeptidome normal tissue filtering.\n\n")
    }
    return(list(
      filtered_data = data,
      filtering_stats = data.frame(
        category = "No filtering applied",
        count = nrow(data),
        percentage = 100
      ),
      peptides_removed = data.frame(),
      filter_applied = FALSE
    ))
  }
  
  if (!normal_exists) {
    if (verbose) {
      cat("WARNING: Normal sample (", normal_sample_id, ") not found in data.\n")
      cat("Cannot perform tumor/normal filtering. Keeping all peptides.\n\n")
    }
    return(list(
      filtered_data = data,
      filtering_stats = data.frame(
        category = "Normal sample missing",
        count = nrow(data),
        percentage = 100
      ),
      peptides_removed = data.frame(),
      filter_applied = FALSE
    ))
  }
  
  if (!tumor_exists) {
    if (verbose) {
      cat("WARNING: Tumor sample (", tumor_sample_id, ") not found in data.\n")
      cat("Cannot perform tumor/normal filtering. Keeping all peptides.\n\n")
    }
    return(list(
      filtered_data = data,
      filtering_stats = data.frame(
        category = "Tumor sample missing",
        count = nrow(data),
        percentage = 100
      ),
      peptides_removed = data.frame(),
      filter_applied = FALSE
    ))
  }
  
  # Get peptides from normal and tumor samples
  normal_peptides <- data %>%
    filter(SampleID == normal_sample_id) %>%
    select(Peptide, final_intensity) %>%
    rename(normal_intensity = final_intensity)
  
  tumor_peptides <- data %>%
    filter(SampleID == tumor_sample_id) %>%
    select(Peptide, final_intensity) %>%
    rename(tumor_intensity = final_intensity)
  
  # Create comprehensive peptide comparison
  peptide_comparison <- data %>%
    select(Peptide) %>%
    distinct() %>%
    left_join(normal_peptides, by = "Peptide") %>%
    left_join(tumor_peptides, by = "Peptide") %>%
    mutate(
      # Classify peptides by presence
      in_normal = !is.na(normal_intensity) & normal_intensity > 0,
      in_tumor = !is.na(tumor_intensity) & tumor_intensity > 0,
      
      # Handle zero/NA intensities for log2FC calculation
      normal_intensity_adj = ifelse(is.na(normal_intensity) | normal_intensity <= 0, 0.1, normal_intensity),
      tumor_intensity_adj = ifelse(is.na(tumor_intensity) | tumor_intensity <= 0, 0.1, tumor_intensity),
      
      # Calculate log2 fold change (tumor/normal)
      log2fc_tumor_normal = log2(tumor_intensity_adj / normal_intensity_adj),
      
      # Classify peptides
      peptide_category = case_when(
        in_tumor & in_normal ~ "Both_tumor_and_normal",
        in_tumor & !in_normal ~ "Tumor_exclusive",
        !in_tumor & in_normal ~ "Normal_exclusive",
        TRUE ~ "Neither_detected"
      ),
      
      # FIXED: Determine if peptide should be kept
      keep_peptide = case_when(
        peptide_category == "Tumor_exclusive" ~ keep_tumor_exclusive,
        peptide_category == "Normal_exclusive" ~ keep_normal_exclusive,
        peptide_category == "Both_tumor_and_normal" ~ log2fc_tumor_normal >= log2fc_threshold,
        peptide_category == "Neither_detected" ~ TRUE,  # KEEP for other samples!
        TRUE ~ FALSE
      )
    )
  
  # Filter the original data
  peptides_to_keep <- peptide_comparison %>%
    filter(keep_peptide) %>%
    pull(Peptide)
  
  filtered_data <- data %>%
    filter(Peptide %in% peptides_to_keep)
  
  # Create filtering statistics with updated language
  filtering_stats <- peptide_comparison %>%
    group_by(peptide_category) %>%
    summarise(
      total_count = n(),
      kept_count = sum(keep_peptide),
      removed_count = sum(!keep_peptide),
      .groups = "drop"
    ) %>%
    mutate(
      kept_percentage = round(100 * kept_count / total_count, 1),
      removed_percentage = round(100 * removed_count / total_count, 1)
    )
  
  # Peptides that were removed for review
  peptides_removed <- peptide_comparison %>%
    filter(!keep_peptide) %>%
    arrange(peptide_category, desc(abs(log2fc_tumor_normal)))
  
  # Detailed reporting with updated language
  if (verbose) {
    cat("FILTERING RESULTS:\n")
    cat("Total unique peptides analyzed:", nrow(peptide_comparison), "\n")
    cat("Peptides kept:", length(peptides_to_keep), "\n")
    cat("Peptides removed:", nrow(peptide_comparison) - length(peptides_to_keep), "\n\n")
    
    cat("BREAKDOWN BY CATEGORY:\n")
    for (i in 1:nrow(filtering_stats)) {
      category <- filtering_stats$peptide_category[i]
      action_description <- case_when(
        category == "Neither_detected" ~ "(kept for other PDX samples)",
        category == "Tumor_exclusive" ~ "(tumor-only is good)",
        category == "Normal_exclusive" ~ "(normal-only removed)",
        category == "Both_tumor_and_normal" ~ paste0("(requires ≥", log2fc_threshold, " log2FC)")
      )
      
      cat(sprintf("%-25s: %4d total, %4d kept (%5.1f%%), %4d removed (%5.1f%%) %s\n",
                  filtering_stats$peptide_category[i],
                  filtering_stats$total_count[i],
                  filtering_stats$kept_count[i],
                  filtering_stats$kept_percentage[i],
                  filtering_stats$removed_count[i],
                  filtering_stats$removed_percentage[i],
                  action_description))
    }
    
    # Show some examples of removed peptides (only those actually removed)
    actually_removed <- peptides_removed %>%
      filter(peptide_category != "Neither_detected")  # Don't show "Neither" as removed
    
    if (nrow(actually_removed) > 0) {
      cat("\nEXAMPLES OF REMOVED PEPTIDES:\n")
      examples <- actually_removed %>%
        group_by(peptide_category) %>%
        slice_head(n = 3) %>%
        ungroup()
      
      for (i in 1:min(10, nrow(examples))) {
        cat(sprintf("  %s (%s): Tumor=%.1f, Normal=%.1f, Log2FC=%.2f\n",
                    examples$Peptide[i],
                    examples$peptide_category[i],
                    examples$tumor_intensity_adj[i],
                    examples$normal_intensity_adj[i],
                    examples$log2fc_tumor_normal[i]))
      }
    }
    cat("\n")
  }
  
  return(list(
    filtered_data = filtered_data,
    filtering_stats = filtering_stats,
    peptide_comparison = peptide_comparison,
    peptides_removed = peptides_removed,
    filter_applied = TRUE,
    filter_parameters = list(
      normal_sample_id = normal_sample_id,
      tumor_sample_id = tumor_sample_id,
      log2fc_threshold = log2fc_threshold,
      keep_normal_exclusive = keep_normal_exclusive,
      keep_tumor_exclusive = keep_tumor_exclusive
    )
  ))
}

#' Load and process transcriptome data with flexible path handling
load_transcriptome_data <- function(data_directory, verbose = TRUE) {
  
  # Try multiple possible paths for transcriptome data
  # Handle both absolute and relative paths
  if (grepl("^~/", data_directory)) {
    # If it starts with ~/, expand it
    data_directory <- path.expand(data_directory)
  }
  
  possible_paths <- c(
    # Direct paths from data_directory
    file.path(data_directory, "transcriptome/fix_transcriptome/20250605_160122_Processed_Normalized_Gene_counts_FLCdb_Panel_1.csv"),
    file.path(data_directory, "transcriptome/Processed_Normalized_Gene_counts_FLCdb_Panel_1.csv"),
    file.path(data_directory, "transcriptome/Normalized_Gene_counts_FLCdb_Panel_1.xlsx"),
    
    # Try parent directory
    file.path(dirname(data_directory), "transcriptome/Normalized_Gene_counts_FLCdb_Panel_1.xlsx"),
    file.path(dirname(data_directory), "transcriptome/fix_transcriptome/20250605_160122_Processed_Normalized_Gene_counts_FLCdb_Panel_1.csv"),
    
    # Try relative to current working directory
    file.path("data", "transcriptome/fix_transcriptome/20250605_160122_Processed_Normalized_Gene_counts_FLCdb_Panel_1.csv"),
    file.path("data", "transcriptome/Normalized_Gene_counts_FLCdb_Panel_1.xlsx"),
    
    # Try without data prefix
    "transcriptome/fix_transcriptome/20250605_160122_Processed_Normalized_Gene_counts_FLCdb_Panel_1.csv",
    "transcriptome/Normalized_Gene_counts_FLCdb_Panel_1.xlsx"
  )
  
  transcriptome_data <- NULL
  
  for (path in possible_paths) {
    if (file.exists(path)) {
      if (verbose) cat("Found transcriptome file:", path, "\n")
      
      if (grepl("\\.csv$", path)) {
        transcriptome_data <- read_csv(path, show_col_types = FALSE)
      } else if (grepl("\\.xlsx$", path)) {
        transcriptome_data <- read_excel(path)
      }
      
      if (verbose) cat("Loaded transcriptome data with", nrow(transcriptome_data), "genes\n")
      break
    }
  }
  
  if (is.null(transcriptome_data)) {
    if (verbose) {
      cat("WARNING: No transcriptome file found. Tried:\n")
      for (path in possible_paths) {
        cat("  -", path, "\n")
      }
    }
    return(data.frame(
      symbol = character(0),
      log2FC = numeric(0)
    ))
  }
  
  return(transcriptome_data)
}

#' Load and process proteome data with flexible path handling
load_proteome_data <- function(data_directory, verbose = TRUE) {
  
  # Handle both absolute and relative paths
  if (grepl("^~/", data_directory)) {
    data_directory <- path.expand(data_directory)
  }
  
  # Try multiple possible paths for LFQ data
  lfq_possible_paths <- c(
    # Direct paths from data_directory
    file.path(data_directory, "proteome/Levin2023/adg7038_Table_S2_LFQ.xlsx"),
    
    # Try parent directory
    file.path(dirname(data_directory), "proteome/Levin2023/adg7038_Table_S2_LFQ.xlsx"),
    
    # Try relative to current working directory
    file.path("data", "proteome/Levin2023/adg7038_Table_S2_LFQ.xlsx"),
    
    # Try without data prefix
    "proteome/Levin2023/adg7038_Table_S2_LFQ.xlsx"
  )
  
  # Try multiple possible paths for TMT data
  tmt_possible_paths <- c(
    # Direct paths from data_directory
    file.path(data_directory, "proteome/Levin2023/adg7038_Table_S1_TMT.xlsx"),
    
    # Try parent directory
    file.path(dirname(data_directory), "proteome/Levin2023/adg7038_Table_S1_TMT.xlsx"),
    
    # Try relative to current working directory
    file.path("data", "proteome/Levin2023/adg7038_Table_S1_TMT.xlsx"),
    
    # Try without data prefix
    "proteome/Levin2023/adg7038_Table_S1_TMT.xlsx"
  )
  
  lfq_data <- NULL
  tmt_data <- NULL
  
  # Load LFQ data
  for (path in lfq_possible_paths) {
    if (file.exists(path)) {
      if (verbose) cat("Found LFQ file:", path, "\n")
      tryCatch({
        # Read from the specific sheet with the processed data
        lfq_data_raw <- read_excel(path, sheet = "Significant and 1.5x_2")
        
        # Process the data to standardize column names and extract relevant info
        lfq_data <- lfq_data_raw %>%
          rename_with(~ gsub(" ", "_", .), everything()) %>%
          mutate(
            Gene_Name = if("Gene_Name" %in% colnames(.)) trimws(Gene_Name) else NA_character_,
            Log2_Difference = if("Log2_Difference" %in% colnames(.)) Log2_Difference else NA_real_,
            P.value = if("P.value" %in% colnames(.)) P.value else 
              if("P_value" %in% colnames(.)) P_value else NA_real_
          )
        
        if (verbose) cat("Loaded LFQ data with", nrow(lfq_data), "entries\n")
        break
      }, error = function(e) {
        if (verbose) cat("Error reading LFQ file:", e$message, "\n")
        # Try reading without specifying sheet
        tryCatch({
          lfq_data <- read_excel(path)
          if (verbose) cat("Loaded LFQ data (default sheet) with", nrow(lfq_data), "entries\n")
          break
        }, error = function(e2) {
          if (verbose) cat("Error reading LFQ file (default sheet):", e2$message, "\n")
        })
      })
    }
  }
  
  # Load TMT data
  for (path in tmt_possible_paths) {
    if (file.exists(path)) {
      if (verbose) cat("Found TMT file:", path, "\n")
      tryCatch({
        # Read from the specific sheet with the processed data
        tmt_data_raw <- read_excel(path, sheet = "Significant and 1.5x_2")
        
        # Process the data to standardize column names and extract relevant info
        tmt_data <- tmt_data_raw %>%
          rename_with(~ gsub(" ", "_", .), everything()) %>%
          mutate(
            Gene_Name = if("Gene_Name" %in% colnames(.)) trimws(Gene_Name) else NA_character_,
            Log2_Difference = if("Log2_Difference" %in% colnames(.)) Log2_Difference else NA_real_,
            P.value = if("P.value" %in% colnames(.)) P.value else 
              if("P_value" %in% colnames(.)) P_value else NA_real_
          )
        
        if (verbose) cat("Loaded TMT data with", nrow(tmt_data), "entries\n")
        break
      }, error = function(e) {
        if (verbose) cat("Error reading TMT file:", e$message, "\n")
        # Try reading without specifying sheet
        tryCatch({
          tmt_data <- read_excel(path)
          if (verbose) cat("Loaded TMT data (default sheet) with", nrow(tmt_data), "entries\n")
          break
        }, error = function(e2) {
          if (verbose) cat("Error reading TMT file (default sheet):", e2$message, "\n")
        })
      })
    }
  }
  
  # Create empty data frames if files not found
  if (is.null(lfq_data)) {
    if (verbose) cat("WARNING: No LFQ proteome file found\n")
    lfq_data <- data.frame(
      Gene_Name = character(0),
      Log2_Difference = numeric(0)
    )
  }
  
  if (is.null(tmt_data)) {
    if (verbose) cat("WARNING: No TMT proteome file found\n")
    tmt_data <- data.frame(
      Gene_Name = character(0),
      Log2_Difference = numeric(0)
    )
  }
  
  return(list(lfq = lfq_data, tmt = tmt_data))
}

#================================================
# MAIN ANALYSIS FUNCTION - FIXED VERSION
#================================================

analyze_neoantigens <- function(
    combined_data,
    min_samples_threshold = 5,
    samples_to_exclude = c(),
    negative_control_samples = c(),
    peptide_length_min = 8,
    peptide_length_max = 12,
    log2fc_up_threshold = 1,
    log2fc_down_threshold = -1,
    log2fc_immuno_threshold = 2,  # Existing parameter
    
    # NEW: Immunopeptidome normal tissue filtering parameters
    apply_immuno_normal_filter = FALSE,
    immuno_normal_log2fc_threshold = 1.0,
    immuno_normal_sample_id = "148N",
    immuno_tumor_sample_id = "148T",
    keep_normal_exclusive_peptides = FALSE,
    keep_tumor_exclusive_peptides = TRUE,
    
    data_directory = "directory/data",
    output_dir = "neoantigen_analysis_results"
) {
  
  cat("Starting multi-omics neoantigen analysis...\n")
  
  # Create descriptive output directory name based on parameters
  excl_count <- length(samples_to_exclude)
  negctrl_count <- length(negative_control_samples)
  
  # Create parameter-based folder name (UPDATED with immunofilter info)
  param_folder <- paste0(
    "neoantigen_analysis_sample", min_samples_threshold,
    "_logfc", log2fc_up_threshold,
    "_excl", excl_count,
    "_negctrl", negctrl_count,
    if (apply_immuno_normal_filter) paste0("_immunofilter", immuno_normal_log2fc_threshold) else "_noImmunofilter"
  )
  
  # Update output directory to include parameter folder
  if (basename(output_dir) == "neoantigen_analysis_results") {
    # If using default name, replace with parameter-based name
    output_dir <- file.path(dirname(output_dir), param_folder)
  } else {
    # If custom path provided, append parameter folder
    output_dir <- file.path(output_dir, param_folder)
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  #================================================
  # STEP 1: Process immunopeptidome data
  #================================================
  
  cat("Processing immunopeptidome data...\n")
  
  # Filter out excluded samples
  if (length(samples_to_exclude) > 0) {
    combined_data <- combined_data %>%
      filter(!SampleID %in% samples_to_exclude)
    cat("Excluded", length(samples_to_exclude), "samples\n")
  }
  
  # Filter out peptides found in negative control samples
  if (length(negative_control_samples) > 0) {
    peptides_in_controls <- combined_data %>%
      filter(SampleID %in% negative_control_samples) %>%
      pull(Peptide) %>%
      unique()
    
    combined_data <- combined_data %>%
      filter(!Peptide %in% peptides_in_controls)
    cat("Filtered out", length(peptides_in_controls), "peptides found in negative controls\n")
  }
  
  # NEW: Apply immunopeptidome normal tissue filtering (EARLY FILTERING)
  immunofilter_results <- NULL
  if (apply_immuno_normal_filter) {
    cat("\nApplying immunopeptidome normal tissue filtering...\n")
    
    immunofilter_results <- apply_tumor_normal_immunopeptidome_filter(
      data = combined_data,
      normal_sample_id = immuno_normal_sample_id,
      tumor_sample_id = immuno_tumor_sample_id,
      log2fc_threshold = immuno_normal_log2fc_threshold,
      keep_normal_exclusive = keep_normal_exclusive_peptides,
      keep_tumor_exclusive = keep_tumor_exclusive_peptides,
      verbose = TRUE
    )
    
    # Update combined_data with filtered results
    if (immunofilter_results$filter_applied) {
      combined_data <- immunofilter_results$filtered_data
      cat("✓ Immunopeptidome normal tissue filtering applied successfully\n")
    } else {
      cat("⚠ Immunopeptidome normal tissue filtering was not applied\n")
    }
  } else {
    cat("Immunopeptidome normal tissue filtering: DISABLED\n")
  }
  
  # Save filtered peptides
  write_csv(combined_data, file.path(output_dir, "filtered_peptides_log2fcImmuno.csv"))
  
  # Summary output
  cat("Peptides after all immunopeptidome filtering:", nrow(combined_data), "\n")
  if (exists("filtered_out") && !is.null(filtered_out)) {
    cat("Peptides filtered out by log2fc_immuno threshold:", nrow(filtered_out), "\n")
  }
  
  # Save immunofilter results if applied
  if (!is.null(immunofilter_results) && immunofilter_results$filter_applied) {
    # Save filtering statistics
    write_csv(immunofilter_results$filtering_stats, 
              file.path(output_dir, "immunofilter_statistics.csv"))
    
    # Save detailed peptide comparison
    write_csv(immunofilter_results$peptide_comparison, 
              file.path(output_dir, "immunofilter_peptide_comparison.csv"))
    
    # Save removed peptides for review
    if (nrow(immunofilter_results$peptides_removed) > 0) {
      write_csv(immunofilter_results$peptides_removed, 
                file.path(output_dir, "immunofilter_removed_peptides.csv"))
    }
    
    # Save filter parameters
    filter_params_df <- data.frame(
      Parameter = c("normal_sample_id", "tumor_sample_id", "log2fc_threshold", 
                    "keep_normal_exclusive", "keep_tumor_exclusive"),
      Value = c(immunofilter_results$filter_parameters$normal_sample_id,
                immunofilter_results$filter_parameters$tumor_sample_id,
                immunofilter_results$filter_parameters$log2fc_threshold,
                immunofilter_results$filter_parameters$keep_normal_exclusive,
                immunofilter_results$filter_parameters$keep_tumor_exclusive)
    )
    write_csv(filter_params_df, file.path(output_dir, "immunofilter_parameters.csv"))
    
    cat("✓ Immunofilter results saved:\n")
    cat("  - immunofilter_statistics.csv\n")
    cat("  - immunofilter_peptide_comparison.csv\n")
    cat("  - immunofilter_removed_peptides.csv\n")
    cat("  - immunofilter_parameters.csv\n")
  }
  
  # Filter by peptide length
  combined_data <- combined_data %>%
    mutate(peptide_length = nchar(Peptide)) %>%
    filter(peptide_length >= peptide_length_min & peptide_length <= peptide_length_max)
  
  # Check available columns and adapt
  available_cols <- colnames(combined_data)
  cat("Available columns:", paste(available_cols, collapse = ", "), "\n")
  
  # Find protein column (could be Protein.ID, ProteinID, Protein, etc.)
  protein_col <- NULL
  protein_patterns <- c("Protein ID", "Protein.ID", "ProteinID", "Protein_ID", "Protein", "protein")
  for (pattern in protein_patterns) {
    if (pattern %in% available_cols) {
      protein_col <- pattern
      break
    }
  }
  
  # Find gene column
  gene_col <- NULL
  gene_patterns <- c("Gene", "gene", "Gene.Name", "GeneName", "Gene_Name")
  for (pattern in gene_patterns) {
    if (pattern %in% available_cols) {
      gene_col <- pattern
      break
    }
  }
  
  # Find intensity column
  intensity_col <- NULL
  intensity_patterns <- c("final_intensity", "Intensity", "intensity", "Intensity.Value", "IntensityValue")
  for (pattern in intensity_patterns) {
    if (pattern %in% available_cols) {
      intensity_col <- pattern
      break
    }
  }
  
  cat("Using columns - Gene:", gene_col, "Protein:", protein_col, "Intensity:", intensity_col, "\n")
  
  # Create peptide summary with sample presence
  peptide_summary <- combined_data %>%
    group_by(Peptide) %>%
    summarise(
      samples_present = n_distinct(SampleID),
      sample_list = paste(sort(unique(SampleID)), collapse = ", "),
      total_intensity = if (!is.null(intensity_col)) sum(.data[[intensity_col]], na.rm = TRUE) else 0,
      mean_intensity = if (!is.null(intensity_col)) mean(.data[[intensity_col]], na.rm = TRUE) else 0,
      genes = if (!is.null(gene_col)) paste(unique(.data[[gene_col]][!is.na(.data[[gene_col]])]), collapse = "; ") else "Unknown",
      proteins = if (!is.null(protein_col)) paste(unique(.data[[protein_col]][!is.na(.data[[protein_col]])]), collapse = "; ") else "Unknown",
      peptide_length = first(peptide_length),
      .groups = "drop"
    )
  
  # Filter for majority presence
  majority_peptides <- peptide_summary %>%
    filter(samples_present >= min_samples_threshold)
  
  cat("Found", nrow(majority_peptides), "peptides present in ≥", min_samples_threshold, "samples\n")
  
  # Process gene information
  gene_info <- map_dfr(majority_peptides$genes, process_gene_info)
  majority_peptides <- bind_cols(majority_peptides, gene_info)
  
  #================================================
  # STEP 2: Load and process omics data
  #================================================
  
  cat("Loading multi-omics data...\n")
  
  # Load transcriptome data
  transcriptome_data <- load_transcriptome_data(data_directory, verbose = TRUE)
  
  # Load proteome data
  proteome_data <- load_proteome_data(data_directory, verbose = TRUE)
  lfq_data <- proteome_data$lfq
  tmt_data <- proteome_data$tmt
  
  #================================================
  # STEP 3: Integrate omics data
  #================================================
  
  cat("Integrating multi-omics data...\n")
  
  # Start with majority peptides
  integrated_data <- majority_peptides
  
  # Add transcriptome data if available
  if (nrow(transcriptome_data) > 0) {
    cat("Integrating transcriptome data...\n")
    
    # Check column names and adapt
    trans_cols <- colnames(transcriptome_data)
    cat("Transcriptome columns:", paste(trans_cols, collapse = ", "), "\n")
    
    # Find gene symbol and log2FC columns
    symbol_col <- NULL
    fc_col <- NULL
    
    symbol_patterns <- c("symbol", "Symbol", "gene_name", "Gene_Name", "gene", "Gene")
    for (pattern in symbol_patterns) {
      if (pattern %in% trans_cols) {
        symbol_col <- pattern
        break
      }
    }
    
    fc_patterns <- c("log2FC", "Log2FC", "log2_fold_change", "logFC", "LogFC")
    for (pattern in fc_patterns) {
      if (pattern %in% trans_cols) {
        fc_col <- pattern
        break
      }
    }
    
    if (!is.null(symbol_col) && !is.null(fc_col)) {
      trans_summary <- transcriptome_data %>%
        select(all_of(c(symbol_col, fc_col))) %>%
        rename(gene_symbol = !!symbol_col, transcriptome_log2fc = !!fc_col)
      
      integrated_data <- integrated_data %>%
        left_join(trans_summary, by = c("primary_gene" = "gene_symbol"), relationship = "many-to-many")
      
      cat("Successfully integrated transcriptome data\n")
    } else {
      cat("WARNING: Could not find required columns in transcriptome data\n")
      integrated_data$transcriptome_log2fc <- NA_real_
    }
  } else {
    cat("No transcriptome data available\n")
    integrated_data$transcriptome_log2fc <- NA_real_
  }
  
  # Add LFQ data if available
  if (nrow(lfq_data) > 0) {
    cat("Integrating LFQ proteome data...\n")
    
    # Check column names and adapt
    lfq_cols <- colnames(lfq_data)
    cat("LFQ columns:", paste(lfq_cols, collapse = ", "), "\n")
    
    # Find gene name and log2FC columns with better pattern matching
    gene_col <- NULL
    fc_col <- NULL
    
    gene_patterns <- c("Gene_Name", "gene_name", "Gene", "gene", "symbol", "Symbol")
    for (pattern in gene_patterns) {
      if (pattern %in% lfq_cols) {
        gene_col <- pattern
        break
      }
    }
    
    fc_patterns <- c("Log2_Difference", "log2_difference", "log2FC", "Log2FC", "logFC", "Log2_fold_change")
    for (pattern in fc_patterns) {
      if (pattern %in% lfq_cols) {
        fc_col <- pattern
        break
      }
    }
    
    if (!is.null(gene_col) && !is.null(fc_col)) {
      lfq_summary <- lfq_data %>%
        select(all_of(c(gene_col, fc_col))) %>%
        rename(gene_name = !!gene_col, lfq_log2fc = !!fc_col) %>%
        filter(!is.na(gene_name) & !is.na(lfq_log2fc))
      
      integrated_data <- integrated_data %>%
        left_join(lfq_summary, by = c("primary_gene" = "gene_name"), relationship = "many-to-many")
      
      cat("Successfully integrated LFQ data with", nrow(lfq_summary), "valid entries\n")
    } else {
      cat("WARNING: Could not find required columns in LFQ data\n")
      cat("  Available columns:", paste(lfq_cols, collapse = ", "), "\n")
      cat("  Looking for gene column (tried:", paste(gene_patterns, collapse = ", "), ")\n")
      cat("  Looking for FC column (tried:", paste(fc_patterns, collapse = ", "), ")\n")
      integrated_data$lfq_log2fc <- NA_real_
    }
  } else {
    integrated_data$lfq_log2fc <- NA_real_
  }
  
  # Add TMT data if available
  if (nrow(tmt_data) > 0) {
    cat("Integrating TMT proteome data...\n")
    
    # Check column names and adapt
    tmt_cols <- colnames(tmt_data)
    cat("TMT columns:", paste(tmt_cols, collapse = ", "), "\n")
    
    # Find gene name and log2FC columns with better pattern matching
    gene_col <- NULL
    fc_col <- NULL
    
    gene_patterns <- c("Gene_Name", "gene_name", "Gene", "gene", "symbol", "Symbol")
    for (pattern in gene_patterns) {
      if (pattern %in% tmt_cols) {
        gene_col <- pattern
        break
      }
    }
    
    fc_patterns <- c("Log2_Difference", "log2_difference", "log2FC", "Log2FC", "logFC", "Log2_fold_change")
    for (pattern in fc_patterns) {
      if (pattern %in% tmt_cols) {
        fc_col <- pattern
        break
      }
    }
    
    if (!is.null(gene_col) && !is.null(fc_col)) {
      tmt_summary <- tmt_data %>%
        select(all_of(c(gene_col, fc_col))) %>%
        rename(gene_name = !!gene_col, tmt_log2fc = !!fc_col) %>%
        filter(!is.na(gene_name) & !is.na(tmt_log2fc))
      
      integrated_data <- integrated_data %>%
        left_join(tmt_summary, by = c("primary_gene" = "gene_name"), relationship = "many-to-many")
      
      cat("Successfully integrated TMT data with", nrow(tmt_summary), "valid entries\n")
    } else {
      cat("WARNING: Could not find required columns in TMT data\n")
      cat("  Available columns:", paste(tmt_cols, collapse = ", "), "\n")
      cat("  Looking for gene column (tried:", paste(gene_patterns, collapse = ", "), ")\n")
      cat("  Looking for FC column (tried:", paste(fc_patterns, collapse = ", "), ")\n")
      integrated_data$tmt_log2fc <- NA_real_
    }
  } else {
    integrated_data$tmt_log2fc <- NA_real_
  }
  
  #================================================
  # STEP 4: Calculate expression statuses and tiers
  #================================================
  
  cat("Calculating expression statuses and tier assignments...\n")
  
  # Ensure all required columns exist
  if (!"transcriptome_log2fc" %in% colnames(integrated_data)) {
    integrated_data$transcriptome_log2fc <- NA_real_
  }
  if (!"lfq_log2fc" %in% colnames(integrated_data)) {
    integrated_data$lfq_log2fc <- NA_real_
  }
  if (!"tmt_log2fc" %in% colnames(integrated_data)) {
    integrated_data$tmt_log2fc <- NA_real_
  }
  
  # Calculate expression statuses
  integrated_data <- integrated_data %>%
    mutate(
      transcriptome_status = map_chr(transcriptome_log2fc, get_expression_status),
      lfq_status = map_chr(lfq_log2fc, get_expression_status),
      tmt_status = map_chr(tmt_log2fc, get_expression_status)
    )
  
  # Combine proteome data
  proteome_combined <- map2(integrated_data$lfq_log2fc, integrated_data$tmt_log2fc, combine_proteome_data)
  
  integrated_data <- integrated_data %>%
    mutate(
      combined_proteome_log2fc = map_dbl(proteome_combined, ~ .x$combined_fc),
      combined_proteome_status = map_chr(proteome_combined, ~ .x$combined_status),
      proteome_agreement = map_chr(proteome_combined, ~ .x$agreement)
    )
  
  # Assign tiers
  integrated_data <- integrated_data %>%
    mutate(
      tier = map2_chr(transcriptome_status, combined_proteome_status, assign_tier)
    )
  
  #================================================
  # STEP 5: Generate outputs
  #================================================
  
  cat("Generating outputs...\n")
  
  # Save main results
  write_csv(integrated_data, file.path(output_dir, "all_majority_peptides_integrated.csv"))
  
  # Save tier-specific files
  tiers <- unique(integrated_data$tier)
  for (tier in tiers) {
    tier_data <- integrated_data %>% filter(tier == !!tier)
    if (nrow(tier_data) > 0) {
      filename <- paste0(gsub("[^A-Za-z0-9_]", "_", tier), ".csv")
      write_csv(tier_data, file.path(output_dir, filename))
    }
  }
  
  # Generate summary statistics
  summary_stats <- integrated_data %>%
    group_by(tier) %>%
    summarise(
      peptide_count = n(),
      multi_gene_count = sum(is_multi_gene),
      single_gene_count = sum(!is_multi_gene),
      multi_gene_percentage = round(100 * sum(is_multi_gene) / n(), 1),
      single_gene_percentage = round(100 * sum(!is_multi_gene) / n(), 1),
      mean_samples_present = mean(samples_present),
      mean_intensity = mean(mean_intensity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(peptide_count))
  
  # Add overall dataset statistics
  overall_stats <- data.frame(
    tier = "OVERALL_DATASET",
    peptide_count = nrow(integrated_data),
    multi_gene_count = sum(integrated_data$is_multi_gene),
    single_gene_count = sum(!integrated_data$is_multi_gene),
    multi_gene_percentage = round(100 * sum(integrated_data$is_multi_gene) / nrow(integrated_data), 1),
    single_gene_percentage = round(100 * sum(!integrated_data$is_multi_gene) / nrow(integrated_data), 1),
    mean_samples_present = mean(integrated_data$samples_present),
    mean_intensity = mean(integrated_data$mean_intensity, na.rm = TRUE)
  )
  
  # Combine tier stats with overall stats
  summary_stats_enhanced <- rbind(summary_stats, overall_stats)
  
  write_csv(summary_stats_enhanced, file.path(output_dir, "tier_summary_statistics.csv"))
  
  # Extract and save single gene peptides
  single_gene_peptides <- integrated_data %>%
    filter(!is_multi_gene) %>%
    arrange(tier, desc(mean_intensity))
  
  write_csv(single_gene_peptides, file.path(output_dir, "single_gene_peptides_only.csv"))
  
  cat("✓ Saved single gene peptides (", nrow(single_gene_peptides), " peptides) to single_gene_peptides_only.csv\n")
  
  #================================================
  # STEP 6: Generate visualizations (with error handling)
  #================================================
  
  cat("Generating visualizations...\n")
  
  tryCatch({
    # 1. Tier distribution bar plot
    if (nrow(summary_stats) > 0) {
      p1 <- ggplot(summary_stats, aes(x = reorder(tier, -peptide_count), y = peptide_count, fill = tier)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = peptide_count), vjust = -0.5, size = 3) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
        labs(title = "Distribution of Peptides by Tier", 
             x = "Tier", y = "Number of Peptides")
      
      ggsave(file.path(output_dir, "tier_distribution.png"), p1, width = 12, height = 8, dpi = 300)
      cat("Created tier distribution plot\n")
    }
  }, error = function(e) {
    cat("Error creating tier distribution plot:", e$message, "\n")
  })
  
  tryCatch({
    # 2. Sample presence histogram
    p2 <- ggplot(integrated_data, aes(x = samples_present, fill = is_multi_gene)) +
      geom_histogram(binwidth = 1, alpha = 0.7, position = "stack") +
      theme_minimal() +
      labs(title = "Distribution of Peptides by Sample Presence", 
           x = "Number of Samples Present", y = "Number of Peptides",
           fill = "Multi-gene Peptide")
    
    ggsave(file.path(output_dir, "sample_presence_distribution.png"), p2, width = 10, height = 6, dpi = 300)
    cat("Created sample presence distribution plot\n")
  }, error = function(e) {
    cat("Error creating sample presence plot:", e$message, "\n")
  })
  
  tryCatch({
    # 3. Multi-omics correlation scatter plot (if data available)
    if (sum(!is.na(integrated_data$transcriptome_log2fc)) > 0 && 
        sum(!is.na(integrated_data$combined_proteome_log2fc)) > 0) {
      
      plot_data <- integrated_data %>%
        filter(!is.na(transcriptome_log2fc) & !is.na(combined_proteome_log2fc))
      
      if (nrow(plot_data) > 0) {
        p3 <- ggplot(plot_data, aes(x = transcriptome_log2fc, y = combined_proteome_log2fc, 
                                    color = tier, shape = is_multi_gene)) +
          geom_point(alpha = 0.7, size = 2) +
          geom_hline(yintercept = c(log2fc_up_threshold, log2fc_down_threshold), 
                     linetype = "dashed", alpha = 0.5) +
          geom_vline(xintercept = c(log2fc_up_threshold, log2fc_down_threshold), 
                     linetype = "dashed", alpha = 0.5) +
          theme_minimal() +
          labs(title = "Transcriptome vs Proteome Expression",
               x = "Transcriptome log2FC", y = "Combined Proteome log2FC",
               color = "Tier", shape = "Multi-gene")
        
        ggsave(file.path(output_dir, "transcriptome_vs_proteome_scatter.png"), p3, width = 12, height = 8, dpi = 300)
        cat("Created transcriptome vs proteome scatter plot\n")
      }
    } else {
      cat("Insufficient data for transcriptome vs proteome scatter plot\n")
    }
  }, error = function(e) {
    cat("Error creating transcriptome vs proteome scatter plot:", e$message, "\n")
  })
  
  tryCatch({
    # 4. Heatmap for top peptides (with better validation)
    top_peptides <- integrated_data %>%
      arrange(desc(samples_present), desc(mean_intensity)) %>%
      head(min(50, nrow(integrated_data)))
    
    if (nrow(top_peptides) > 5 && nrow(combined_data) > 0) {
      # Create presence matrix
      all_samples <- unique(combined_data$SampleID)
      all_samples <- all_samples[!is.na(all_samples)]  # Remove any NA sample IDs
      
      if (length(all_samples) > 0) {
        presence_matrix <- matrix(0, nrow = nrow(top_peptides), ncol = length(all_samples))
        rownames(presence_matrix) <- make.unique(as.character(top_peptides$Peptide))
        colnames(presence_matrix) <- make.unique(as.character(all_samples))
        
        for (i in 1:nrow(top_peptides)) {
          peptide <- top_peptides$Peptide[i]
          samples_with_peptide <- combined_data %>%
            filter(Peptide == peptide) %>%
            pull(SampleID) %>%
            unique()
          samples_with_peptide <- samples_with_peptide[!is.na(samples_with_peptide)]
          
          # Only set to 1 for samples that exist in our matrix
          valid_samples <- intersect(samples_with_peptide, all_samples)
          if (length(valid_samples) > 0) {
            presence_matrix[i, valid_samples] <- 1
          }
        }
        
        # Only create heatmap if we have actual data and valid matrix
        if (sum(presence_matrix) > 0 && nrow(presence_matrix) > 1 && ncol(presence_matrix) > 1) {
          # Create annotation
          row_annotation <- data.frame(
            Tier = top_peptides$tier,
            Multi_Gene = top_peptides$is_multi_gene,
            row.names = rownames(presence_matrix)
          )
          
          png(file.path(output_dir, "peptide_presence_heatmap.png"), width = 1200, height = 800, res = 150)
          pheatmap(presence_matrix,
                   main = "Peptide Presence Across Samples (Top 50)",
                   color = c("white", "darkblue"),
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   annotation_row = row_annotation,
                   fontsize_row = 6,
                   fontsize_col = 8,
                   show_rownames = TRUE,
                   show_colnames = TRUE)
          dev.off()
          cat("Created peptide presence heatmap\n")
        } else {
          cat("No peptide presence data available for heatmap\n")
        }
      } else {
        cat("No valid sample IDs found for heatmap\n")
      }
    } else {
      cat("Insufficient data for peptide presence heatmap\n")
    }
  }, error = function(e) {
    cat("Error creating peptide presence heatmap:", e$message, "\n")
  })
  
  tryCatch({
    # 5. Venn diagram of omics overlaps (if multiple omics available)
    overlap_lists <- list()
    
    # Add immunopeptidome genes
    immuno_genes <- integrated_data$primary_gene[!is.na(integrated_data$primary_gene)]
    if (length(immuno_genes) > 0) {
      overlap_lists[["Immunopeptidome"]] <- immuno_genes
    }
    
    # Add transcriptome genes if available
    if (nrow(transcriptome_data) > 0) {
      trans_gene_col <- NULL
      symbol_patterns <- c("symbol", "Symbol", "gene_name", "Gene_Name", "gene", "Gene")
      for (pattern in symbol_patterns) {
        if (pattern %in% colnames(transcriptome_data)) {
          trans_gene_col <- pattern
          break
        }
      }
      if (!is.null(trans_gene_col)) {
        trans_genes <- transcriptome_data[[trans_gene_col]][!is.na(transcriptome_data[[trans_gene_col]])]
        if (length(trans_genes) > 0) {
          overlap_lists[["Transcriptome"]] <- trans_genes
        }
      }
    }
    
    # Add LFQ genes if available
    if (nrow(lfq_data) > 0) {
      lfq_gene_col <- NULL
      gene_patterns <- c("Gene_Name", "gene_name", "Gene", "gene", "symbol", "Symbol")
      for (pattern in gene_patterns) {
        if (pattern %in% colnames(lfq_data)) {
          lfq_gene_col <- pattern
          break
        }
      }
      if (!is.null(lfq_gene_col)) {
        lfq_genes <- lfq_data[[lfq_gene_col]][!is.na(lfq_data[[lfq_gene_col]])]
        if (length(lfq_genes) > 0) {
          overlap_lists[["LFQ_Proteome"]] <- lfq_genes
        }
      }
    }
    
    # Add TMT genes if available
    if (nrow(tmt_data) > 0) {
      tmt_gene_col <- NULL
      gene_patterns <- c("Gene_Name", "gene_name", "Gene", "gene", "symbol", "Symbol")
      for (pattern in gene_patterns) {
        if (pattern %in% colnames(tmt_data)) {
          tmt_gene_col <- pattern
          break
        }
      }
      if (!is.null(tmt_gene_col)) {
        tmt_genes <- tmt_data[[tmt_gene_col]][!is.na(tmt_data[[tmt_gene_col]])]
        if (length(tmt_genes) > 0) {
          overlap_lists[["TMT_Proteome"]] <- tmt_genes
        }
      }
    }
    
    # Create Venn diagram if we have at least 2 datasets
    if (length(overlap_lists) >= 2) {
      png(file.path(output_dir, "omics_overlap_venn.png"), width = 800, height = 800, res = 150)
      
      if (length(overlap_lists) == 2) {
        draw.pairwise.venn(
          area1 = length(overlap_lists[[1]]),
          area2 = length(overlap_lists[[2]]),
          cross.area = length(intersect(overlap_lists[[1]], overlap_lists[[2]])),
          category = names(overlap_lists),
          fill = c("lightblue", "lightcoral"),
          alpha = 0.5,
          lty = "blank"
        )
      } else if (length(overlap_lists) == 3) {
        draw.triple.venn(
          area1 = length(overlap_lists[[1]]),
          area2 = length(overlap_lists[[2]]),
          area3 = length(overlap_lists[[3]]),
          n12 = length(intersect(overlap_lists[[1]], overlap_lists[[2]])),
          n23 = length(intersect(overlap_lists[[2]], overlap_lists[[3]])),
          n13 = length(intersect(overlap_lists[[1]], overlap_lists[[3]])),
          n123 = length(Reduce(intersect, overlap_lists)),
          category = names(overlap_lists),
          fill = c("lightblue", "lightcoral", "lightgreen"),
          alpha = 0.5,
          lty = "blank"
        )
      }
      
      # After the draw.pairwise.venn() or draw.triple.venn() calls, add:
      grid::grid.newpage()  # Clear any previous plots
      
      dev.off()
      cat("Created omics overlap Venn diagram\n")
    } else {
      cat("Insufficient datasets for Venn diagram (need at least 2)\n")
    }
  }, error = function(e) {
    cat("Error creating Venn diagram:", e$message, "\n")
  })
  
  #================================================
  # STEP 7: Generate comprehensive analysis summary
  #================================================
  
  cat("Generating comprehensive analysis summary...\n")
  
  # Create analysis summary
  analysis_summary <- create_neoantigen_analysis_summary(
    integrated_data = integrated_data,
    summary_stats = summary_stats,
    analysis_params = list(
      min_samples_threshold = min_samples_threshold,
      samples_to_exclude = samples_to_exclude,
      negative_control_samples = negative_control_samples,
      peptide_length_min = peptide_length_min,
      peptide_length_max = peptide_length_max,
      log2fc_up_threshold = log2fc_up_threshold,
      log2fc_down_threshold = log2fc_down_threshold,
      log2fc_immuno_threshold = log2fc_immuno_threshold,
      apply_immuno_normal_filter = apply_immuno_normal_filter,
      immuno_normal_log2fc_threshold = immuno_normal_log2fc_threshold,
      immuno_normal_sample_id = immuno_normal_sample_id,
      immuno_tumor_sample_id = immuno_tumor_sample_id,
      keep_normal_exclusive_peptides = keep_normal_exclusive_peptides,
      keep_tumor_exclusive_peptides = keep_tumor_exclusive_peptides,
      data_directory = data_directory
    ),
    data_availability = list(
      transcriptome = nrow(transcriptome_data) > 0,
      lfq = nrow(lfq_data) > 0,
      tmt = nrow(tmt_data) > 0
    ),
    omics_support_counts = list(
      transcriptome = sum(!is.na(integrated_data$transcriptome_log2fc)),
      lfq = sum(!is.na(integrated_data$lfq_log2fc)),
      tmt = sum(!is.na(integrated_data$tmt_log2fc)),
      total_peptides = nrow(integrated_data)
    ),
    immunofilter_results = immunofilter_results,
    output_dir = output_dir
  )
  
  #================================================
  # STEP 8: Return results
  #================================================
  
  cat("Analysis complete!\n")
  cat("Results saved to:", output_dir, "\n")
  cat("Total peptides analyzed:", nrow(integrated_data), "\n")
  cat("Tier distribution:\n")
  print(summary_stats)
  
  # Provide data availability summary
  cat("\nData availability summary:\n")
  cat("- Transcriptome data:", ifelse(nrow(transcriptome_data) > 0, "Available", "Not available"), "\n")
  cat("- LFQ proteome data:", ifelse(nrow(lfq_data) > 0, "Available", "Not available"), "\n")
  cat("- TMT proteome data:", ifelse(nrow(tmt_data) > 0, "Available", "Not available"), "\n")
  
  # Count peptides with omics support
  peptides_with_transcriptome <- sum(!is.na(integrated_data$transcriptome_log2fc))
  peptides_with_lfq <- sum(!is.na(integrated_data$lfq_log2fc))
  peptides_with_tmt <- sum(!is.na(integrated_data$tmt_log2fc))
  
  cat("Peptides with omics support:\n")
  cat("- With transcriptome data:", peptides_with_transcriptome, "/", nrow(integrated_data), "\n")
  cat("- With LFQ proteome data:", peptides_with_lfq, "/", nrow(integrated_data), "\n")
  cat("- With TMT proteome data:", peptides_with_tmt, "/", nrow(integrated_data), "\n")
  
  # Immunofilter summary
  if (!is.null(immunofilter_results) && immunofilter_results$filter_applied) {
    cat("\nImmunopeptidome normal tissue filtering summary:\n")
    cat("- Filter applied: YES\n")
    cat("- Normal sample:", immuno_normal_sample_id, "\n")
    cat("- Tumor sample:", immuno_tumor_sample_id, "\n")
    cat("- Log2FC threshold:", immuno_normal_log2fc_threshold, "(", round(2^immuno_normal_log2fc_threshold, 1), "-fold)\n")
    if (nrow(immunofilter_results$filtering_stats) > 0) {
      total_before <- sum(immunofilter_results$filtering_stats$total_count)
      total_after <- sum(immunofilter_results$filtering_stats$kept_count)
      cat("- Peptides before filtering:", total_before, "\n")
      cat("- Peptides after filtering:", total_after, "\n")
      cat("- Peptides removed:", total_before - total_after, "\n")
    }
  } else {
    cat("\nImmunopeptidome normal tissue filtering: DISABLED\n")
  }
  
  return(list(
    integrated_data = integrated_data,
    summary_stats = summary_stats,
    output_dir = output_dir,
    immunofilter_results = immunofilter_results,  # NEW: Include filter results
    data_availability = list(
      transcriptome = nrow(transcriptome_data) > 0,
      lfq = nrow(lfq_data) > 0,
      tmt = nrow(tmt_data) > 0
    ),
    omics_support_counts = list(
      transcriptome = peptides_with_transcriptome,
      lfq = peptides_with_lfq,
      tmt = peptides_with_tmt,
      total_peptides = nrow(integrated_data)
    ),
    # NEW: Add immunofilter summary
    immunofilter_summary = if (!is.null(immunofilter_results) && immunofilter_results$filter_applied) {
      list(
        filter_applied = TRUE,
        normal_sample = immuno_normal_sample_id,
        tumor_sample = immuno_tumor_sample_id,
        log2fc_threshold = immuno_normal_log2fc_threshold,
        total_peptides_before = sum(immunofilter_results$filtering_stats$total_count),
        total_peptides_after = sum(immunofilter_results$filtering_stats$kept_count),
        peptides_removed = sum(immunofilter_results$filtering_stats$removed_count),
        filtering_stats = immunofilter_results$filtering_stats
      )
    } else {
      list(filter_applied = FALSE)
    }
  ))
}

#' Create Comprehensive Neoantigen Analysis Summary - ENHANCED VERSION
#' 
#' Generates detailed text and CSV summaries of the analysis including immunofilter breakdown
#' 
create_neoantigen_analysis_summary <- function(
    integrated_data, 
    summary_stats, 
    analysis_params, 
    data_availability, 
    omics_support_counts,
    immunofilter_results = NULL,  # NEW: Include immunofilter results
    output_dir
) {
  
  # Create analysis timestamp
  analysis_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  #================================================
  # 1. CREATE TEXT SUMMARY
  #================================================
  
  text_summary <- paste0(
    "=================================================================================\n",
    "                    MULTI-OMICS NEOANTIGEN ANALYSIS SUMMARY\n",
    "=================================================================================\n\n",
    
    "Analysis Date: ", analysis_timestamp, "\n",
    "Output Directory: ", basename(output_dir), "\n\n",
    
    "=== ANALYSIS PARAMETERS ===\n",
    "Minimum samples threshold: ", analysis_params$min_samples_threshold, " samples\n",
    "Peptide length range: ", analysis_params$peptide_length_min, "-", analysis_params$peptide_length_max, " amino acids\n",
    "Log2FC thresholds: UP > ", analysis_params$log2fc_up_threshold, ", DOWN < ", analysis_params$log2fc_down_threshold, "\n",
    "Excluded samples: ", ifelse(length(analysis_params$samples_to_exclude) > 0, 
                                 paste(analysis_params$samples_to_exclude, collapse = ", "), "None"), "\n",
    "Negative control samples: ", ifelse(length(analysis_params$negative_control_samples) > 0, 
                                         paste(analysis_params$negative_control_samples, collapse = ", "), "None"), "\n\n"
  )
  
  # NEW: Add detailed immunofilter section
  if (!is.null(analysis_params$apply_immuno_normal_filter) && analysis_params$apply_immuno_normal_filter && 
      !is.null(immunofilter_results) && immunofilter_results$filter_applied) {
    
    text_summary <- paste0(text_summary,
                           "=== IMMUNOPEPTIDOME NORMAL TISSUE FILTERING (DETAILED) ===\n",
                           "Strategy: Keep peptides not found in T/N pair for evaluation in other PDX samples\n",
                           "Normal sample ID: ", analysis_params$immuno_normal_sample_id, "\n",
                           "Tumor sample ID: ", analysis_params$immuno_tumor_sample_id, "\n",
                           "Log2FC threshold: ", analysis_params$immuno_normal_log2fc_threshold, " (", 
                           round(2^analysis_params$immuno_normal_log2fc_threshold, 1), "-fold)\n",
                           "Keep normal-exclusive peptides: ", analysis_params$keep_normal_exclusive_peptides, "\n",
                           "Keep tumor-exclusive peptides: ", analysis_params$keep_tumor_exclusive_peptides, "\n\n",
                           
                           "FILTERING RESULTS:\n",
                           "Total unique peptides analyzed: ", sum(immunofilter_results$filtering_stats$total_count), "\n",
                           "Peptides kept: ", sum(immunofilter_results$filtering_stats$kept_count), "\n",
                           "Peptides removed: ", sum(immunofilter_results$filtering_stats$removed_count), "\n",
                           "Retention rate: ", round(100 * sum(immunofilter_results$filtering_stats$kept_count) / 
                                                       sum(immunofilter_results$filtering_stats$total_count), 1), "%\n\n",
                           
                           "BREAKDOWN BY CATEGORY:\n"
    )
    
    # Add detailed breakdown for each category
    for (i in 1:nrow(immunofilter_results$filtering_stats)) {
      category <- immunofilter_results$filtering_stats$peptide_category[i]
      total <- immunofilter_results$filtering_stats$total_count[i]
      kept <- immunofilter_results$filtering_stats$kept_count[i]
      removed <- immunofilter_results$filtering_stats$removed_count[i]
      kept_pct <- immunofilter_results$filtering_stats$kept_percentage[i]
      removed_pct <- immunofilter_results$filtering_stats$removed_percentage[i]
      
      # Add biological interpretation
      interpretation <- case_when(
        category == "Both_tumor_and_normal" ~ paste0("(requires ≥", analysis_params$immuno_normal_log2fc_threshold, " log2FC)"),
        category == "Neither_detected" ~ "(kept for other PDX samples)",
        category == "Normal_exclusive" ~ "(normal-only removed - reduces off-target risk)",
        category == "Tumor_exclusive" ~ "(tumor-only is ideal - kept all)",
        TRUE ~ ""
      )
      
      text_summary <- paste0(text_summary,
                             sprintf("%-25s: %5d total, %5d kept (%5.1f%%), %5d removed (%5.1f%%) %s\n",
                                     category, total, kept, kept_pct, removed, removed_pct, interpretation)
      )
    }
    text_summary <- paste0(text_summary, "\n")
    
  } else if (!is.null(analysis_params$apply_immuno_normal_filter) && analysis_params$apply_immuno_normal_filter) {
    text_summary <- paste0(text_summary,
                           "=== IMMUNOPEPTIDOME NORMAL TISSUE FILTERING ===\n",
                           "Applied: YES (but filtering was not successful - see analysis log)\n",
                           "Normal sample: ", analysis_params$immuno_normal_sample_id, "\n",
                           "Tumor sample: ", analysis_params$immuno_tumor_sample_id, "\n\n"
    )
  } else {
    text_summary <- paste0(text_summary,
                           "=== IMMUNOPEPTIDOME NORMAL TISSUE FILTERING ===\n",
                           "Applied: NO (disabled in parameters)\n\n"
    )
  }
  
  text_summary <- paste0(text_summary,
                         "=== DATA AVAILABILITY ===\n",
                         "Transcriptome data: ", ifelse(data_availability$transcriptome, "✓ Available", "✗ Not available"), "\n",
                         "LFQ proteome data: ", ifelse(data_availability$lfq, "✓ Available", "✗ Not available"), "\n",
                         "TMT proteome data: ", ifelse(data_availability$tmt, "✓ Available", "✗ Not available"), "\n\n",
                         
                         "=== FINAL ANALYSIS RESULTS ===\n",
                         "Total peptides analyzed: ", omics_support_counts$total_peptides, "\n",
                         "Peptides with transcriptome support: ", omics_support_counts$transcriptome, " (", 
                         round(100 * omics_support_counts$transcriptome / omics_support_counts$total_peptides, 1), "%)\n",
                         "Peptides with LFQ proteome support: ", omics_support_counts$lfq, " (", 
                         round(100 * omics_support_counts$lfq / omics_support_counts$total_peptides, 1), "%)\n",
                         "Peptides with TMT proteome support: ", omics_support_counts$tmt, " (", 
                         round(100 * omics_support_counts$tmt / omics_support_counts$total_peptides, 1), "%)\n\n"
  )
  
  # Add tier-specific results
  if (nrow(summary_stats) > 0) {
    text_summary <- paste0(text_summary,
                           "=== TIER DISTRIBUTION IN YOUR DATA ===\n"
    )
    
    for (i in 1:nrow(summary_stats)) {
      tier_name <- summary_stats$tier[i]
      count <- summary_stats$peptide_count[i]
      percentage <- round(100 * count / sum(summary_stats$peptide_count), 1)
      
      text_summary <- paste0(text_summary,
                             sprintf("%-40s: %4d peptides (%5.1f%%)\n", tier_name, count, percentage)
      )
    }
    text_summary <- paste0(text_summary, "\n")
  }
  
  # Add filtering efficiency summary
  if (!is.null(immunofilter_results) && immunofilter_results$filter_applied) {
    total_before <- sum(immunofilter_results$filtering_stats$total_count)
    total_after <- sum(immunofilter_results$filtering_stats$kept_count)
    retention_rate <- round(100 * total_after / total_before, 1)
    
    text_summary <- paste0(text_summary,
                           "=== FILTERING EFFICIENCY SUMMARY ===\n",
                           "This analysis used a conservative filtering approach to minimize off-target risks:\n",
                           "• Started with ", format(total_before, big.mark = ","), " unique peptides\n",
                           "• Retained ", format(total_after, big.mark = ","), " peptides (", retention_rate, "%) for analysis\n",
                           "• Removed only normal-exclusive peptides and tumor-depleted peptides\n",
                           "• Preserved peptides not found in 148T/148N for evaluation in other PDX models\n",
                           "• Final analysis set: ", format(omics_support_counts$total_peptides, big.mark = ","), " peptides meeting presence criteria\n\n"
    )
  }
  
  # Save text summary
  writeLines(text_summary, file.path(output_dir, "ANALYSIS_SUMMARY.txt"))
  
  #================================================
  # 2. CREATE ENHANCED CSV SUMMARY
  #================================================
  
  # Base metrics
  base_metrics <- c(
    "Analysis_Date",
    "Output_Directory",
    "Min_Samples_Threshold",
    "Peptide_Length_Range",
    "Log2FC_Up_Threshold",
    "Log2FC_Down_Threshold",
    "Excluded_Samples_Count",
    "Negative_Control_Samples_Count",
    "Transcriptome_Data_Available",
    "LFQ_Proteome_Data_Available", 
    "TMT_Proteome_Data_Available",
    "Total_Peptides_Analyzed",
    "Peptides_With_Transcriptome_Support",
    "Peptides_With_LFQ_Support",
    "Peptides_With_TMT_Support",
    "Transcriptome_Support_Percentage",
    "LFQ_Support_Percentage",
    "TMT_Support_Percentage",
    "Total_Tiers_Present",
    "Most_Common_Tier",
    "Tier_1_Peptides_Count",
    "High_Priority_Tiers_Count"
  )
  
  base_values <- c(
    analysis_timestamp,
    basename(output_dir),
    analysis_params$min_samples_threshold,
    paste0(analysis_params$peptide_length_min, "-", analysis_params$peptide_length_max, "aa"),
    analysis_params$log2fc_up_threshold,
    analysis_params$log2fc_down_threshold,
    length(analysis_params$samples_to_exclude),
    length(analysis_params$negative_control_samples),
    data_availability$transcriptome,
    data_availability$lfq,
    data_availability$tmt,
    omics_support_counts$total_peptides,
    omics_support_counts$transcriptome,
    omics_support_counts$lfq,
    omics_support_counts$tmt,
    round(100 * omics_support_counts$transcriptome / omics_support_counts$total_peptides, 1),
    round(100 * omics_support_counts$lfq / omics_support_counts$total_peptides, 1),
    round(100 * omics_support_counts$tmt / omics_support_counts$total_peptides, 1),
    nrow(summary_stats),
    ifelse(nrow(summary_stats) > 0, summary_stats$tier[1], "None"),
    ifelse(nrow(summary_stats) > 0, 
           sum(summary_stats$peptide_count[grepl("Tier_1_", summary_stats$tier)]), 0),
    ifelse(nrow(summary_stats) > 0,
           sum(summary_stats$peptide_count[grepl("Tier_[1-5]_", summary_stats$tier)]), 0)
  )
  
  # Add immunofilter metrics if available
  if (!is.null(immunofilter_results) && immunofilter_results$filter_applied) {
    immunofilter_metrics <- c(
      "Immunofilter_Applied",
      "Immunofilter_Log2FC_Threshold",
      "Immunofilter_Normal_Sample",
      "Immunofilter_Tumor_Sample",
      "Immunofilter_Strategy",
      "Immunofilter_Total_Peptides_Analyzed",
      "Immunofilter_Peptides_Kept",
      "Immunofilter_Peptides_Removed",
      "Immunofilter_Retention_Percentage",
      "Immunofilter_Both_TumorNormal_Total",
      "Immunofilter_Both_TumorNormal_Kept",
      "Immunofilter_Both_TumorNormal_Removed",
      "Immunofilter_Neither_Detected_Total",
      "Immunofilter_Neither_Detected_Kept",
      "Immunofilter_Normal_Exclusive_Total",
      "Immunofilter_Normal_Exclusive_Removed",
      "Immunofilter_Tumor_Exclusive_Total",
      "Immunofilter_Tumor_Exclusive_Kept"
    )
    
    total_analyzed <- sum(immunofilter_results$filtering_stats$total_count)
    total_kept <- sum(immunofilter_results$filtering_stats$kept_count)
    total_removed <- sum(immunofilter_results$filtering_stats$removed_count)
    
    # Extract category-specific stats
    both_stats <- immunofilter_results$filtering_stats[immunofilter_results$filtering_stats$peptide_category == "Both_tumor_and_normal", ]
    neither_stats <- immunofilter_results$filtering_stats[immunofilter_results$filtering_stats$peptide_category == "Neither_detected", ]
    normal_stats <- immunofilter_results$filtering_stats[immunofilter_results$filtering_stats$peptide_category == "Normal_exclusive", ]
    tumor_stats <- immunofilter_results$filtering_stats[immunofilter_results$filtering_stats$peptide_category == "Tumor_exclusive", ]
    
    immunofilter_values <- c(
      TRUE,
      analysis_params$immuno_normal_log2fc_threshold,
      analysis_params$immuno_normal_sample_id,
      analysis_params$immuno_tumor_sample_id,
      "Keep_peptides_not_in_TN_pair_for_other_PDX",
      total_analyzed,
      total_kept,
      total_removed,
      round(100 * total_kept / total_analyzed, 1),
      ifelse(nrow(both_stats) > 0, both_stats$total_count, 0),
      ifelse(nrow(both_stats) > 0, both_stats$kept_count, 0),
      ifelse(nrow(both_stats) > 0, both_stats$removed_count, 0),
      ifelse(nrow(neither_stats) > 0, neither_stats$total_count, 0),
      ifelse(nrow(neither_stats) > 0, neither_stats$kept_count, 0),
      ifelse(nrow(normal_stats) > 0, normal_stats$total_count, 0),
      ifelse(nrow(normal_stats) > 0, normal_stats$removed_count, 0),
      ifelse(nrow(tumor_stats) > 0, tumor_stats$total_count, 0),
      ifelse(nrow(tumor_stats) > 0, tumor_stats$kept_count, 0)
    )
    
    # Combine metrics
    all_metrics <- c(base_metrics, immunofilter_metrics)
    all_values <- c(base_values, immunofilter_values)
    
  } else {
    # Add basic immunofilter status
    immunofilter_metrics <- c(
      "Immunofilter_Applied",
      "Immunofilter_Status"
    )
    
    immunofilter_values <- c(
      ifelse(!is.null(analysis_params$apply_immuno_normal_filter), analysis_params$apply_immuno_normal_filter, FALSE),
      ifelse(!is.null(analysis_params$apply_immuno_normal_filter) && analysis_params$apply_immuno_normal_filter, "Enabled_but_failed", "Disabled")
    )
    
    # Combine metrics
    all_metrics <- c(base_metrics, immunofilter_metrics)
    all_values <- c(base_values, immunofilter_values)
  }
  
  csv_summary <- data.frame(
    Metric = all_metrics,
    Value = all_values
  )
  
  # Save CSV summary
  write_csv(csv_summary, file.path(output_dir, "ANALYSIS_SUMMARY.csv"))
  
  cat("✓ Enhanced analysis summary files created:\n")
  cat("  - ANALYSIS_SUMMARY.txt (comprehensive guide with immunofilter details)\n")
  cat("  - ANALYSIS_SUMMARY.csv (key metrics including immunofilter breakdown)\n")
  
  return(list(
    text_summary = text_summary,
    csv_summary = csv_summary
  ))
}