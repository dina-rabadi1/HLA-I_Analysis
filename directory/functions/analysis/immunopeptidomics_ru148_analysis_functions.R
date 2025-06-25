#' RU148 Multi-Omics Analysis Functions
#' 
#' Functions for performing comprehensive multi-omics analysis of RU148 tumor vs normal
#' focusing on tumor-exclusive and tumor-upregulated peptides
#' Author: Generated for Dina Rabadi's HLA-I Analysis Pipeline
#' immunopeptidomics_ru148_analysis.R

#' RU148 Multi-Omics Analysis Functions
#' 
#' Functions for performing comprehensive multi-omics analysis of RU148 tumor vs normal
#' focusing on tumor-exclusive and tumor-upregulated peptides
#' Author: Generated for Dina Rabadi's HLA-I Analysis Pipeline
#' 

# Load required libraries
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(readxl)) install.packages("readxl")
if (!require(writexl)) install.packages("writexl")
if (!require(openxlsx)) install.packages("openxlsx")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(pheatmap)) install.packages("pheatmap")
if (!require(VennDiagram)) install.packages("VennDiagram")
if (!require(RColorBrewer)) install.packages("RColorBrewer")
if (!require(gridExtra)) install.packages("gridExtra")

library(tidyverse)
library(readxl)
library(writexl)
library(openxlsx)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(gridExtra)

#' Process RU148 Transcriptome Data
#' 
#' Reads original transcriptome file and creates focused RU148 dataset
#' 
#' @param transcriptome_path Path to original transcriptome Excel file
#' @param verbose Print progress messages
#' 
#' @return Processed RU148 transcriptome data
#' 
process_ru148_transcriptome <- function(
    transcriptome_path = "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/directory/data/transcriptome/Normalized_Gene_counts_FLCdb_Panel_1.xlsx",
    verbose = TRUE
) {
  
  if (verbose) cat("Processing RU148-specific transcriptome data...\n")
  
  # Read the original transcriptome file
  if (!file.exists(transcriptome_path)) {
    warning("Transcriptome file not found: ", transcriptome_path)
    return(data.frame(
      symbol = character(0),
      RU148_N = numeric(0),
      RU148_T8 = numeric(0),
      RU148_T11 = numeric(0),
      RU148_T_avg = numeric(0),
      log2_fold_change_transcriptome_ru148 = numeric(0),
      geneID = character(0)
    ))
  }
  
  transcriptome_raw <- read_excel(transcriptome_path)
  
  # Check if RU148 columns exist
  ru148_cols <- c("RU148_N", "RU148_T8", "RU148_T11")
  missing_cols <- ru148_cols[!ru148_cols %in% colnames(transcriptome_raw)]
  
  if (length(missing_cols) > 0) {
    warning("Missing RU148 columns in transcriptome: ", paste(missing_cols, collapse = ", "))
    # Create empty columns for missing ones
    for (col in missing_cols) {
      transcriptome_raw[[col]] <- NA_real_
    }
  }
  
  # Process RU148 data
  ru148_transcriptome <- transcriptome_raw %>%
    select(all_of(c("symbol", "geneID", ru148_cols))) %>%
    rowwise() %>%
    mutate(
      # Average the tumor samples
      RU148_T_avg = mean(c(RU148_T8, RU148_T11), na.rm = TRUE),
      # Calculate log2 fold change
      log2_fold_change_transcriptome_ru148 = log2(
        ifelse(RU148_T_avg == 0 | is.na(RU148_T_avg), 0.1, RU148_T_avg) / 
          ifelse(RU148_N == 0 | is.na(RU148_N), 0.1, RU148_N)
      )
    ) %>%
    ungroup() %>%
    # Remove rows with all NA values
    filter(!is.na(symbol))
  
  if (verbose) {
    cat("  Processed", nrow(ru148_transcriptome), "genes for RU148\n")
    cat("  RU148_N range:", round(min(ru148_transcriptome$RU148_N, na.rm = TRUE), 1), 
        "to", round(max(ru148_transcriptome$RU148_N, na.rm = TRUE), 1), "\n")
    cat("  RU148_T_avg range:", round(min(ru148_transcriptome$RU148_T_avg, na.rm = TRUE), 1), 
        "to", round(max(ru148_transcriptome$RU148_T_avg, na.rm = TRUE), 1), "\n")
  }
  
  return(ru148_transcriptome)
}

#' Prepare RU148 Multi-Omics Data
#' 
#' Loads and processes all omics datasets for RU148 analysis
#' 
#' @param combined_data Combined immunopeptidome data from main pipeline
#' @param lfq_path Path to LFQ proteome data
#' @param tmt_path Path to TMT proteome data
#' @param transcriptome_path Path to original transcriptome Excel file
#' @param verbose Print progress messages
#' 
#' @return List containing all processed omics datasets
#' 
prepare_ru148_omics_data <- function(
    combined_data,
    lfq_path = "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/directory/data/proteome/Levin2023/adg7038_Table_S2_LFQ.xlsx",
    tmt_path = "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/directory/data/proteome/Levin2023/adg7038_Table_S1_TMT.xlsx",
    transcriptome_path = "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/directory/data/transcriptome/Normalized_Gene_counts_FLCdb_Panel_1.xlsx",
    verbose = TRUE
) {
  
  if (verbose) cat("=== PREPARING RU148 MULTI-OMICS DATA ===\n")
  
  # 1. Process immunopeptidome data for RU148 from combined_data
  if (verbose) cat("Processing RU148 immunopeptidome data from combined dataset...\n")
  
  # Filter combined_data for RU148 samples
  ru148_immuno <- combined_data %>%
    filter(grepl("148", SampleID)) %>%
    mutate(
      # Create sample type indicator
      sample_type = ifelse(grepl("148T", SourceFile) | grepl("148T", SampleID), "148T", "148N"),
      peptide_length = nchar(Peptide)
    ) %>%
    # Filter for peptide length 8-12
    filter(peptide_length >= 8 & peptide_length <= 12) %>%
    # Aggregate by peptide and sample type
    group_by(Peptide, sample_type) %>%
    summarize(
      peptide_length = first(peptide_length),
      spectral_count = sum(`Spectral Count`, na.rm = TRUE),
      total_intensity = sum(Intensity, na.rm = TRUE),
      protein_ids = paste(unique(`Protein ID`), collapse = "; "),
      genes = paste(unique(Gene), collapse = "; "),
      .groups = "drop"
    ) %>%
    # Pivot to wide format
    pivot_wider(
      names_from = sample_type,
      values_from = c(spectral_count, total_intensity, protein_ids, genes),
      values_fill = list(spectral_count = 0, total_intensity = 0, protein_ids = NA, genes = NA)
    ) %>%
    # Calculate fold changes and detection status
    mutate(
      # Calculate adjusted intensities for log calculation (handle zeros)
      total_intensity_148N_adj = ifelse(total_intensity_148N == 0, 0.1, total_intensity_148N),
      total_intensity_148T_adj = ifelse(total_intensity_148T == 0, 0.1, total_intensity_148T),
      
      # Calculate log2 fold change
      log2_fold_change_immuno = log2(total_intensity_148T_adj / total_intensity_148N_adj),
      
      # Determine detection status
      detection_status = case_when(
        total_intensity_148T > 0 & total_intensity_148N == 0 ~ "Tumor-exclusive",
        total_intensity_148N > 0 & total_intensity_148T == 0 ~ "Normal-exclusive",
        total_intensity_148T > 0 & total_intensity_148N > 0 ~ "Detected in both",
        TRUE ~ "Not detected"
      ),
      
      # Extract primary gene
      genes_combined = coalesce(genes_148T, genes_148N),
      proteins_combined = coalesce(protein_ids_148T, protein_ids_148N),
      primary_gene = sapply(strsplit(genes_combined, ";\\s*"), function(x) trimws(x[1]))
    )
  
  if (verbose) {
    cat("  Processed", nrow(ru148_immuno), "unique RU148 peptides\n")
    
    # Add detection status summary
    detection_summary <- ru148_immuno %>% count(detection_status, sort = TRUE)
    cat("  Detection status breakdown:\n")
    for (i in 1:nrow(detection_summary)) {
      cat("    ", detection_summary$detection_status[i], ":", detection_summary$n[i], "\n")
    }
  }
  
  # 2. Process RU148 transcriptome data
  if (verbose) cat("Processing RU148 transcriptome data...\n")
  ru148_transcriptome <- process_ru148_transcriptome(transcriptome_path, verbose = verbose)
  
  # 3. Process LFQ proteome data
  if (verbose) cat("Processing LFQ proteome data...\n")
  
  if (file.exists(lfq_path)) {
    lfq_data <- read_excel(lfq_path, sheet = "Significant and 1.5x_2")
    lfq_processed <- lfq_data %>%
      rename_with(~ gsub(" ", "_", .), everything()) %>%
      mutate(
        Gene_Name = trimws(Gene_Name),
        log2_fold_change_lfq = Log2_Difference,
        p_value_lfq = P.value
      ) %>%
      select(Gene_Name, log2_fold_change_lfq, p_value_lfq, Protein_Name)
    
    if (verbose) cat("  Processed", nrow(lfq_processed), "entries from LFQ proteome\n")
  } else {
    warning("LFQ file not found: ", lfq_path)
    lfq_processed <- data.frame(
      Gene_Name = character(0),
      log2_fold_change_lfq = numeric(0),
      p_value_lfq = numeric(0),
      Protein_Name = character(0)
    )
  }
  
  # 4. Process TMT proteome data
  if (verbose) cat("Processing TMT proteome data...\n")
  
  if (file.exists(tmt_path)) {
    tmt_data <- read_excel(tmt_path, sheet = "Significant and 1.5x_2")
    tmt_processed <- tmt_data %>%
      rename_with(~ gsub(" ", "_", .), everything()) %>%
      mutate(
        Gene_Name = trimws(Gene_Name),
        log2_fold_change_tmt = Log2_Difference,
        p_value_tmt = P.value
      ) %>%
      select(Gene_Name, log2_fold_change_tmt, p_value_tmt, Protein_Name)
    
    if (verbose) cat("  Processed", nrow(tmt_processed), "entries from TMT proteome\n")
  } else {
    warning("TMT file not found: ", tmt_path)
    tmt_processed <- data.frame(
      Gene_Name = character(0),
      log2_fold_change_tmt = numeric(0),
      p_value_tmt = numeric(0),
      Protein_Name = character(0)
    )
  }
  
  if (verbose) cat("✓ RU148 multi-omics data preparation completed\n\n")
  
  return(list(
    immunopeptidome = ru148_immuno,
    transcriptome = ru148_transcriptome,
    lfq_proteome = lfq_processed,
    tmt_proteome = tmt_processed
  ))
}

#' Perform RU148 Tumor-Exclusive Analysis
#' 
#' Analyzes peptides found exclusively in tumor samples
#' 
#' @param omics_data List of omics datasets from prepare_ru148_omics_data
#' @param verbose Print progress messages
#' 
#' @return Processed tumor-exclusive analysis results
#' 
perform_ru148_tumor_exclusive_analysis <- function(
    omics_data,
    verbose = TRUE
) {
  
  if (verbose) cat("=== PERFORMING TUMOR-EXCLUSIVE ANALYSIS ===\n")
  
  # Extract tumor-exclusive peptides
  tumor_exclusive <- omics_data$immunopeptidome %>%
    filter(detection_status == "Tumor-exclusive") %>%
    arrange(desc(total_intensity_148T))
  
  if (verbose) cat("Found", nrow(tumor_exclusive), "tumor-exclusive peptides\n")
  
  if (nrow(tumor_exclusive) == 0) {
    warning("No tumor-exclusive peptides found")
    return(list(
      peptides = tumor_exclusive,
      combined_analysis = data.frame(),
      regulation_categories = data.frame()
    ))
  }
  
  # Join with transcriptome data
  combined_analysis <- tumor_exclusive %>%
    left_join(
      omics_data$transcriptome %>% select(symbol, log2_fold_change_transcriptome_ru148, geneID),
      by = c("primary_gene" = "symbol"),
      relationship = "many-to-many"
    ) %>%
    # Join with LFQ data
    left_join(
      omics_data$lfq_proteome,
      by = c("primary_gene" = "Gene_Name"),
      relationship = "many-to-many"
    ) %>%
    # Join with TMT data
    left_join(
      omics_data$tmt_proteome,
      by = c("primary_gene" = "Gene_Name"),
      relationship = "many-to-many"
    ) %>%
    # Create regulation categories
    mutate(
      # Transcriptome regulation
      transcriptome_regulation = case_when(
        is.na(log2_fold_change_transcriptome_ru148) ~ "No data",
        log2_fold_change_transcriptome_ru148 > 1 ~ "Upregulated",
        log2_fold_change_transcriptome_ru148 < -1 ~ "Downregulated",
        TRUE ~ "Unchanged"
      ),
      
      # LFQ regulation
      lfq_regulation = case_when(
        is.na(log2_fold_change_lfq) ~ "No data",
        log2_fold_change_lfq > 1 ~ "Upregulated",
        log2_fold_change_lfq < -1 ~ "Downregulated",
        TRUE ~ "Unchanged"
      ),
      
      # TMT regulation
      tmt_regulation = case_when(
        is.na(log2_fold_change_tmt) ~ "No data",
        log2_fold_change_tmt > 1 ~ "Upregulated",
        log2_fold_change_tmt < -1 ~ "Downregulated",
        TRUE ~ "Unchanged"
      ),
      
      # Overall regulation pattern
      regulation_pattern = paste(
        "T:", transcriptome_regulation,
        "| L:", lfq_regulation,
        "| M:", tmt_regulation
      ),
      
      # Simplified category
      regulation_category = case_when(
        transcriptome_regulation == "Upregulated" & 
          lfq_regulation == "Upregulated" & 
          tmt_regulation == "Upregulated" ~ "Consistent upregulation",
        
        transcriptome_regulation == "Upregulated" & 
          (lfq_regulation == "Upregulated" | tmt_regulation == "Upregulated") ~ "Mostly upregulated",
        
        transcriptome_regulation == "Upregulated" ~ "Transcriptome-specific upregulation",
        lfq_regulation == "Upregulated" ~ "LFQ-specific upregulation",
        tmt_regulation == "Upregulated" ~ "TMT-specific upregulation",
        
        (transcriptome_regulation == "Downregulated" | 
           lfq_regulation == "Downregulated" | 
           tmt_regulation == "Downregulated") ~ "Discordant regulation",
        
        TRUE ~ "Peptide-specific (no omics support)"
      )
    )
  
  # Create summary of regulation categories
  regulation_summary <- combined_analysis %>%
    group_by(regulation_category) %>%
    summarize(
      count = n(),
      avg_intensity = mean(total_intensity_148T, na.rm = TRUE),
      genes = paste(unique(primary_gene[!is.na(primary_gene)]), collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(desc(count))
  
  if (verbose) {
    cat("✓ Tumor-exclusive analysis completed\n")
    cat("Regulation categories:\n")
    print(regulation_summary)
  }
  
  return(list(
    peptides = tumor_exclusive,
    combined_analysis = combined_analysis,
    regulation_categories = regulation_summary
  ))
}

#' Perform RU148 Tumor-Upregulated Analysis
#' 
#' Analyzes peptides that are shared but highly upregulated in tumor (>4-fold)
#' 
#' @param omics_data List of omics datasets from prepare_ru148_omics_data
#' @param fold_change_threshold log2 fold change threshold (default: 2, meaning 4-fold)
#' @param verbose Print progress messages
#' 
#' @return Processed tumor-upregulated analysis results
#' 
perform_ru148_tumor_upregulated_analysis <- function(
    omics_data,
    fold_change_threshold = 2,
    verbose = TRUE
) {
  
  if (verbose) cat("=== PERFORMING TUMOR-UPREGULATED ANALYSIS ===\n")
  if (verbose) cat("Using fold change threshold: log2FC >", fold_change_threshold, "(", 2^fold_change_threshold, "-fold)\n")
  
  # Extract tumor-upregulated peptides (shared but highly upregulated)
  tumor_upregulated <- omics_data$immunopeptidome %>%
    filter(
      detection_status == "Detected in both",
      log2_fold_change_immuno > fold_change_threshold
    ) %>%
    arrange(desc(log2_fold_change_immuno))
  
  if (verbose) cat("Found", nrow(tumor_upregulated), "tumor-upregulated peptides\n")
  
  if (nrow(tumor_upregulated) == 0) {
    warning("No tumor-upregulated peptides found")
    return(list(
      peptides = tumor_upregulated,
      combined_analysis = data.frame(),
      regulation_categories = data.frame()
    ))
  }
  
  # Join with transcriptome data
  combined_analysis <- tumor_upregulated %>%
    left_join(
      omics_data$transcriptome %>% select(symbol, log2_fold_change_transcriptome_ru148, geneID),
      by = c("primary_gene" = "symbol"),
      relationship = "many-to-many"
    ) %>%
    # Join with LFQ data
    left_join(
      omics_data$lfq_proteome,
      by = c("primary_gene" = "Gene_Name"),
      relationship = "many-to-many"
    ) %>%
    # Join with TMT data
    left_join(
      omics_data$tmt_proteome,
      by = c("primary_gene" = "Gene_Name"),
      relationship = "many-to-many"
    ) %>%
    # Create regulation categories (same logic as tumor-exclusive)
    mutate(
      # Transcriptome regulation
      transcriptome_regulation = case_when(
        is.na(log2_fold_change_transcriptome_ru148) ~ "No data",
        log2_fold_change_transcriptome_ru148 > 1 ~ "Upregulated",
        log2_fold_change_transcriptome_ru148 < -1 ~ "Downregulated",
        TRUE ~ "Unchanged"
      ),
      
      # LFQ regulation
      lfq_regulation = case_when(
        is.na(log2_fold_change_lfq) ~ "No data",
        log2_fold_change_lfq > 1 ~ "Upregulated",
        log2_fold_change_lfq < -1 ~ "Downregulated",
        TRUE ~ "Unchanged"
      ),
      
      # TMT regulation
      tmt_regulation = case_when(
        is.na(log2_fold_change_tmt) ~ "No data",
        log2_fold_change_tmt > 1 ~ "Upregulated",
        log2_fold_change_tmt < -1 ~ "Downregulated",
        TRUE ~ "Unchanged"
      ),
      
      # Overall regulation pattern
      regulation_pattern = paste(
        "T:", transcriptome_regulation,
        "| L:", lfq_regulation,
        "| M:", tmt_regulation
      ),
      
      # Simplified category
      regulation_category = case_when(
        transcriptome_regulation == "Upregulated" & 
          lfq_regulation == "Upregulated" & 
          tmt_regulation == "Upregulated" ~ "Consistent upregulation",
        
        transcriptome_regulation == "Upregulated" & 
          (lfq_regulation == "Upregulated" | tmt_regulation == "Upregulated") ~ "Mostly upregulated",
        
        transcriptome_regulation == "Upregulated" ~ "Transcriptome-specific upregulation",
        lfq_regulation == "Upregulated" ~ "LFQ-specific upregulation",
        tmt_regulation == "Upregulated" ~ "TMT-specific upregulation",
        
        (transcriptome_regulation == "Downregulated" | 
           lfq_regulation == "Downregulated" | 
           tmt_regulation == "Downregulated") ~ "Discordant regulation",
        
        TRUE ~ "Peptide-specific (no omics support)"
      )
    )
  
  # Create summary of regulation categories
  regulation_summary <- combined_analysis %>%
    group_by(regulation_category) %>%
    summarize(
      count = n(),
      avg_log2fc_immuno = mean(log2_fold_change_immuno, na.rm = TRUE),
      avg_intensity_tumor = mean(total_intensity_148T, na.rm = TRUE),
      genes = paste(unique(primary_gene[!is.na(primary_gene)]), collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(desc(count))
  
  if (verbose) {
    cat("✓ Tumor-upregulated analysis completed\n")
    cat("Regulation categories:\n")
    print(regulation_summary)
  }
  
  return(list(
    peptides = tumor_upregulated,
    combined_analysis = combined_analysis,
    regulation_categories = regulation_summary
  ))
}

#' Create Enhanced RU148 Analysis Visualizations
#' 
#' Creates comprehensive visualizations including functions from original script
#' 
#' @param tumor_exclusive_results Results from perform_ru148_tumor_exclusive_analysis
#' @param tumor_upregulated_results Results from perform_ru148_tumor_upregulated_analysis
#' @param output_dir Output directory for plots
#' @param dataset_name Dataset name for file naming
#' @param timestamp Timestamp for file naming
#' @param verbose Print progress messages
#' 
#' @return List of created plots
#' 
create_ru148_analysis_plots <- function(
    tumor_exclusive_results,
    tumor_upregulated_results,
    output_dir,
    dataset_name,
    timestamp,
    verbose = TRUE
) {
  
  if (verbose) cat("=== CREATING RU148 ANALYSIS VISUALIZATIONS ===\n")
  
  # Create plot directories
  plot_dir_exclusive <- file.path(output_dir, "plots", "tumor_exclusive")
  plot_dir_upregulated <- file.path(output_dir, "plots", "tumor_upregulated")
  
  if (!dir.exists(plot_dir_exclusive)) dir.create(plot_dir_exclusive, recursive = TRUE)
  if (!dir.exists(plot_dir_upregulated)) dir.create(plot_dir_upregulated, recursive = TRUE)
  
  created_plots <- list()
  
  # 1. Regulation category bar plots
  if (nrow(tumor_exclusive_results$regulation_categories) > 0) {
    if (verbose) cat("Creating tumor-exclusive regulation category plot...\n")
    
    p1 <- ggplot(tumor_exclusive_results$regulation_categories, 
                 aes(x = reorder(regulation_category, count), y = count, fill = regulation_category)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = count), hjust = -0.1) +
      coord_flip() +
      theme_minimal() +
      scale_fill_brewer(palette = "Set3") +
      labs(
        title = "RU148 Tumor-Exclusive Peptides: Regulation Categories",
        subtitle = paste("Total peptides:", sum(tumor_exclusive_results$regulation_categories$count)),
        x = "Regulation Category",
        y = "Number of Peptides",
        fill = "Category"
      ) +
      theme(legend.position = "none")
    
    plot_file <- file.path(plot_dir_exclusive, paste0(timestamp, "_", dataset_name, "_01_regulation_categories.png"))
    ggsave(plot_file, p1, width = 12, height = 8, dpi = 300)
    created_plots[["tumor_exclusive_regulation"]] <- plot_file
  }
  
  if (nrow(tumor_upregulated_results$regulation_categories) > 0) {
    if (verbose) cat("Creating tumor-upregulated regulation category plot...\n")
    
    p2 <- ggplot(tumor_upregulated_results$regulation_categories, 
                 aes(x = reorder(regulation_category, count), y = count, fill = regulation_category)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = count), hjust = -0.1) +
      coord_flip() +
      theme_minimal() +
      scale_fill_brewer(palette = "Set3") +
      labs(
        title = "RU148 Tumor-Upregulated Peptides: Regulation Categories",
        subtitle = paste("Total peptides:", sum(tumor_upregulated_results$regulation_categories$count)),
        x = "Regulation Category",
        y = "Number of Peptides",
        fill = "Category"
      ) +
      theme(legend.position = "none")
    
    plot_file <- file.path(plot_dir_upregulated, paste0(timestamp, "_", dataset_name, "_01_regulation_categories.png"))
    ggsave(plot_file, p2, width = 12, height = 8, dpi = 300)
    created_plots[["tumor_upregulated_regulation"]] <- plot_file
  }
  
  # 2. Volcano plots
  if (nrow(tumor_upregulated_results$combined_analysis) > 0) {
    if (verbose) cat("Creating volcano plot...\n")
    
    volcano_data <- tumor_upregulated_results$combined_analysis %>%
      filter(!is.na(log2_fold_change_immuno))
    
    if (nrow(volcano_data) > 0) {
      p_volcano <- ggplot(volcano_data, aes(x = log2_fold_change_immuno, y = -log10(0.05), 
                                            color = regulation_category)) +
        geom_point(alpha = 0.7, size = 2) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
        theme_minimal() +
        scale_color_brewer(palette = "Set2") +
        labs(
          title = "RU148 Tumor-Exclusive Peptides: Transcriptome vs Immunopeptidome",
          x = "Log2 Fold Change Transcriptome (Tumor/Normal)",
          y = "Log10(Peptide Intensity in Tumor + 1)",
          color = "Regulation Category"
        ) +
        theme(legend.position = "bottom")
      
      plot_file <- file.path(plot_dir_exclusive, paste0(timestamp, "_", dataset_name, "_02_transcriptome_scatter.png"))
      ggsave(plot_file, p_volcano, width = 12, height = 8, dpi = 300)
      created_plots[["tumor_exclusive_transcriptome"]] <- plot_file
    }
  }
  
  if (nrow(tumor_upregulated_results$combined_analysis) > 0) {
    if (verbose) cat("Creating tumor-upregulated multi-omics scatter plots...\n")
    
    plot_data <- tumor_upregulated_results$combined_analysis %>%
      filter(!is.na(log2_fold_change_transcriptome_ru148))
    
    if (nrow(plot_data) > 0) {
      p4 <- ggplot(plot_data, aes(x = log2_fold_change_transcriptome_ru148, y = log2_fold_change_immuno, 
                                  color = regulation_category)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_hline(yintercept = 2, linetype = "dashed", alpha = 0.5, color = "red") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
        theme_minimal() +
        scale_color_brewer(palette = "Set2") +
        labs(
          title = "RU148 Tumor-Upregulated Peptides: Transcriptome vs Immunopeptidome",
          x = "Log2 Fold Change Transcriptome (Tumor/Normal)",
          y = "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
          color = "Regulation Category"
        ) +
        theme(legend.position = "bottom")
      
      plot_file <- file.path(plot_dir_upregulated, paste0(timestamp, "_", dataset_name, "_02_transcriptome_scatter.png"))
      ggsave(plot_file, p4, width = 12, height = 8, dpi = 300)
      created_plots[["tumor_upregulated_transcriptome"]] <- plot_file
    }
  }
  
  # 4. Heatmaps for multi-omics comparison
  # Tumor-exclusive heatmap (transcriptome/LFQ/TMT only)
  if (nrow(tumor_exclusive_results$combined_analysis) > 0 && "primary_gene" %in% colnames(tumor_exclusive_results$combined_analysis)) {
    if (verbose) cat("Creating tumor-exclusive heatmap...\n")
    
    heatmap_data <- tumor_exclusive_results$combined_analysis %>%
      filter(!is.na(primary_gene)) %>%
      select(primary_gene, Peptide, log2_fold_change_transcriptome_ru148, 
             log2_fold_change_lfq, log2_fold_change_tmt, regulation_category) %>%
      distinct()
    
    if (nrow(heatmap_data) > 5 && any(!is.na(heatmap_data$log2_fold_change_transcriptome_ru148))) {
      # Create matrix (3 columns only: transcriptome, LFQ, TMT)
      heatmap_matrix <- heatmap_data %>%
        select(log2_fold_change_transcriptome_ru148, log2_fold_change_lfq, log2_fold_change_tmt) %>%
        as.matrix()
      
      # Replace NA with 0 for visualization
      heatmap_matrix[is.na(heatmap_matrix)] <- 0
      
      # Create unique row names
      row_labels <- paste0(heatmap_data$primary_gene, " (", substr(heatmap_data$Peptide, 1, 8), "...)")
      rownames(heatmap_matrix) <- make.unique(row_labels)
      
      # Create annotation
      row_annotation <- data.frame(
        Regulation = heatmap_data$regulation_category,
        row.names = rownames(heatmap_matrix)
      )
      
      # Create heatmap
      plot_file <- file.path(plot_dir_exclusive, paste0(timestamp, "_", dataset_name, "_03_multiomics_heatmap.png"))
      
      png(plot_file, width = 1000, height = max(800, nrow(heatmap_matrix)*30), res = 150)
      pheatmap(
        heatmap_matrix,
        main = "RU148 Tumor-Exclusive Peptides: Multi-Omics Fold Changes",
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-3, 3, length.out = 101),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        annotation_row = row_annotation,
        display_numbers = TRUE,
        number_format = "%.1f",
        fontsize_row = 8,
        fontsize_col = 10,
        labels_col = c("Transcriptome", "LFQ Proteome", "TMT Proteome")
      )
      dev.off()
      
      created_plots[["tumor_exclusive_heatmap"]] <- plot_file
    }
  }
  
  # Tumor-upregulated heatmap (immunopeptidome + transcriptome/LFQ/TMT)
  if (nrow(tumor_upregulated_results$combined_analysis) > 0 && "primary_gene" %in% colnames(tumor_upregulated_results$combined_analysis)) {
    if (verbose) cat("Creating tumor-upregulated heatmap...\n")
    
    heatmap_data <- tumor_upregulated_results$combined_analysis %>%
      filter(!is.na(primary_gene)) %>%
      select(primary_gene, Peptide, log2_fold_change_immuno, log2_fold_change_transcriptome_ru148, 
             log2_fold_change_lfq, log2_fold_change_tmt, regulation_category) %>%
      distinct()
    
    if (nrow(heatmap_data) > 5) {
      # Create matrix
      heatmap_matrix <- heatmap_data %>%
        select(log2_fold_change_immuno, log2_fold_change_transcriptome_ru148, 
               log2_fold_change_lfq, log2_fold_change_tmt) %>%
        as.matrix()
      
      # Replace NA with 0 for visualization
      heatmap_matrix[is.na(heatmap_matrix)] <- 0
      
      # Create unique row names
      row_labels <- paste0(heatmap_data$primary_gene, " (", substr(heatmap_data$Peptide, 1, 8), "...)")
      rownames(heatmap_matrix) <- make.unique(row_labels)
      
      # Create annotation
      row_annotation <- data.frame(
        Regulation = heatmap_data$regulation_category,
        row.names = rownames(heatmap_matrix)
      )
      
      # Create heatmap
      plot_file <- file.path(plot_dir_upregulated, paste0(timestamp, "_", dataset_name, "_03_multiomics_heatmap.png"))
      
      png(plot_file, width = 1200, height = max(800, nrow(heatmap_matrix)*30), res = 150)
      pheatmap(
        heatmap_matrix,
        main = "RU148 Tumor-Upregulated Peptides: Multi-Omics Fold Changes",
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-3, 3, length.out = 101),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        annotation_row = row_annotation,
        display_numbers = TRUE,
        number_format = "%.1f",
        fontsize_row = 8,
        fontsize_col = 10,
        labels_col = c("Immunopeptidome", "Transcriptome", "LFQ Proteome", "TMT Proteome")
      )
      dev.off()
      
      created_plots[["tumor_upregulated_heatmap"]] <- plot_file
    }
  }
  
  # 5. Venn diagram of regulation overlap
  if (verbose) cat("Creating regulation overlap visualization...\n")
  
  # Check if we have data and required columns before proceeding
  exclusive_has_data <- nrow(tumor_exclusive_results$combined_analysis) > 0 && 
    "primary_gene" %in% colnames(tumor_exclusive_results$combined_analysis)
  upregulated_has_data <- nrow(tumor_upregulated_results$combined_analysis) > 0 && 
    "primary_gene" %in% colnames(tumor_upregulated_results$combined_analysis)
  
  # Combine data for overlap analysis
  all_genes_exclusive <- data.frame()
  all_genes_upregulated <- data.frame()
  
  if (exclusive_has_data) {
    all_genes_exclusive <- tumor_exclusive_results$combined_analysis %>%
      filter(!is.na(primary_gene)) %>%
      select(primary_gene, transcriptome_regulation, lfq_regulation, tmt_regulation) %>%
      distinct()
  }
  
  if (upregulated_has_data) {
    all_genes_upregulated <- tumor_upregulated_results$combined_analysis %>%
      filter(!is.na(primary_gene)) %>%
      select(primary_gene, transcriptome_regulation, lfq_regulation, tmt_regulation) %>%
      distinct()
  }
  
  # Create overlap summary for both analyses
  for (analysis_type in c("exclusive", "upregulated")) {
    data_to_use <- if (analysis_type == "exclusive") all_genes_exclusive else all_genes_upregulated
    plot_dir_to_use <- if (analysis_type == "exclusive") plot_dir_exclusive else plot_dir_upregulated
    
    if (nrow(data_to_use) > 0) {
      # Count genes upregulated in each dataset
      genes_trans_up <- data_to_use %>% filter(transcriptome_regulation == "Upregulated") %>% pull(primary_gene)
      genes_lfq_up <- data_to_use %>% filter(lfq_regulation == "Upregulated") %>% pull(primary_gene)
      genes_tmt_up <- data_to_use %>% filter(tmt_regulation == "Upregulated") %>% pull(primary_gene)
      
      if (length(genes_trans_up) > 0 || length(genes_lfq_up) > 0 || length(genes_tmt_up) > 0) {
        # Create Venn diagram
        if (length(genes_trans_up) > 0 && length(genes_lfq_up) > 0 && length(genes_tmt_up) > 0) {
          plot_file <- file.path(plot_dir_to_use, paste0(timestamp, "_", dataset_name, "_06_regulation_venn.png"))
          
          venn_colors <- brewer.pal(3, "Set1")
          
          png(plot_file, width = 800, height = 800, res = 120)
          draw.triple.venn(
            area1 = length(genes_trans_up),
            area2 = length(genes_lfq_up),
            area3 = length(genes_tmt_up),
            n12 = length(intersect(genes_trans_up, genes_lfq_up)),
            n23 = length(intersect(genes_lfq_up, genes_tmt_up)),
            n13 = length(intersect(genes_trans_up, genes_tmt_up)),
            n123 = length(intersect(intersect(genes_trans_up, genes_lfq_up), genes_tmt_up)),
            category = c("Transcriptome", "LFQ Proteome", "TMT Proteome"),
            fill = venn_colors,
            alpha = 0.5,
            lty = "blank",
            cex = 1.5,
            cat.cex = 1.2,
            cat.col = venn_colors
          )
          dev.off()
          
          created_plots[[paste0("tumor_", analysis_type, "_venn")]] <- plot_file
        }
        
        # Create overlap summary bar plot
        overlap_summary <- data.frame(
          Category = c(
            "Transcriptome only",
            "LFQ only", 
            "TMT only",
            "Trans + LFQ",
            "Trans + TMT",
            "LFQ + TMT",
            "All three"
          ),
          Count = c(
            length(setdiff(genes_trans_up, union(genes_lfq_up, genes_tmt_up))),
            length(setdiff(genes_lfq_up, union(genes_trans_up, genes_tmt_up))),
            length(setdiff(genes_tmt_up, union(genes_trans_up, genes_lfq_up))),
            length(setdiff(intersect(genes_trans_up, genes_lfq_up), genes_tmt_up)),
            length(setdiff(intersect(genes_trans_up, genes_tmt_up), genes_lfq_up)),
            length(setdiff(intersect(genes_lfq_up, genes_tmt_up), genes_trans_up)),
            length(intersect(intersect(genes_trans_up, genes_lfq_up), genes_tmt_up))
          )
        ) %>%
          filter(Count > 0) %>%
          arrange(desc(Count))
        
        if (nrow(overlap_summary) > 0) {
          p_overlap <- ggplot(overlap_summary, aes(x = reorder(Category, Count), y = Count, fill = Category)) +
            geom_bar(stat = "identity") +
            geom_text(aes(label = Count), hjust = -0.1) +
            coord_flip() +
            theme_minimal() +
            scale_fill_brewer(palette = "Set2") +
            labs(
              title = paste("RU148 Tumor-", str_to_title(analysis_type), "Genes: Upregulation Overlap"),
              x = "Overlap Category",
              y = "Number of Genes",
              fill = "Category"
            ) +
            theme(legend.position = "none")
          
          plot_file <- file.path(plot_dir_to_use, paste0(timestamp, "_", dataset_name, "_04_regulation_overlap.png"))
          ggsave(plot_file, p_overlap, width = 10, height = 6, dpi = 300)
          created_plots[[paste0("tumor_", analysis_type, "_overlap")]] <- plot_file
        }
      }
    }
  }
  
  # 6. Additional scatter plots (LFQ vs TMT, etc.)
  if (nrow(tumor_upregulated_results$combined_analysis) > 0 && "log2_fold_change_lfq" %in% colnames(tumor_upregulated_results$combined_analysis)) {
    if (verbose) cat("Creating additional scatter plots...\n")
    
    # LFQ vs TMT scatter
    plot_data_lfq_tmt <- tumor_upregulated_results$combined_analysis %>%
      filter(!is.na(log2_fold_change_lfq) & !is.na(log2_fold_change_tmt))
    
    if (nrow(plot_data_lfq_tmt) > 0) {
      p_lfq_tmt <- ggplot(plot_data_lfq_tmt, aes(x = log2_fold_change_lfq, y = log2_fold_change_tmt, 
                                                 color = regulation_category)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        geom_abline(intercept = 0, slope = 1, linetype = "dotted", alpha = 0.5) +
        theme_minimal() +
        scale_color_brewer(palette = "Set2") +
        labs(
          title = "RU148 Tumor-Upregulated Peptides: LFQ vs TMT Proteome",
          x = "Log2 Fold Change LFQ Proteome (Tumor/Normal)",
          y = "Log2 Fold Change TMT Proteome (Tumor/Normal)",
          color = "Regulation Category"
        ) +
        theme(legend.position = "bottom")
      
      plot_file <- file.path(plot_dir_upregulated, paste0(timestamp, "_", dataset_name, "_07_lfq_vs_tmt_scatter.png"))
      ggsave(plot_file, p_lfq_tmt, width = 10, height = 8, dpi = 300)
      created_plots[["lfq_vs_tmt_scatter"]] <- plot_file
    }
  }
  
  if (verbose) cat("✓ RU148 analysis visualizations completed\n")
  cat("Created", length(created_plots), "visualization files\n")
  
  return(created_plots)
}

#' Save RU148 Analysis Results
#' 
#' Saves all RU148 analysis results to Excel and CSV files
#' 
#' @param tumor_exclusive_results Results from tumor-exclusive analysis
#' @param tumor_upregulated_results Results from tumor-upregulated analysis
#' @param output_dir Output directory
#' @param dataset_name Dataset name for file naming
#' @param timestamp Timestamp for file naming
#' @param verbose Print progress messages
#' 
#' @return List of saved files
#' 
save_ru148_analysis_results <- function(
    tumor_exclusive_results,
    tumor_upregulated_results,
    output_dir,
    dataset_name,
    timestamp,
    verbose = TRUE
) {
  
  if (verbose) cat("=== SAVING RU148 ANALYSIS RESULTS ===\n")
  
  # Create table directories
  table_dir_exclusive <- file.path(output_dir, "tables", "tumor_exclusive")
  table_dir_upregulated <- file.path(output_dir, "tables", "tumor_upregulated")
  
  if (!dir.exists(table_dir_exclusive)) dir.create(table_dir_exclusive, recursive = TRUE)
  if (!dir.exists(table_dir_upregulated)) dir.create(table_dir_upregulated, recursive = TRUE)
  
  saved_files <- list()
  
  # 1. Save tumor-exclusive results as CSV files
  if (verbose) cat("Saving tumor-exclusive results...\n")
  
  # Combined analysis CSV
  csv_file_exclusive_combined <- file.path(table_dir_exclusive, 
                                           paste0(timestamp, "_", dataset_name, "_tumor_exclusive_combined.csv"))
  write_csv(tumor_exclusive_results$combined_analysis, csv_file_exclusive_combined)
  saved_files[["tumor_exclusive_combined_csv"]] <- csv_file_exclusive_combined
  
  # Regulation categories CSV
  csv_file_exclusive_categories <- file.path(table_dir_exclusive, 
                                             paste0(timestamp, "_", dataset_name, "_tumor_exclusive_regulation_categories.csv"))
  write_csv(tumor_exclusive_results$regulation_categories, csv_file_exclusive_categories)
  saved_files[["tumor_exclusive_categories_csv"]] <- csv_file_exclusive_categories
  
  # Peptides only CSV
  csv_file_exclusive_peptides <- file.path(table_dir_exclusive, 
                                           paste0(timestamp, "_", dataset_name, "_tumor_exclusive_peptides.csv"))
  write_csv(tumor_exclusive_results$peptides, csv_file_exclusive_peptides)
  saved_files[["tumor_exclusive_peptides_csv"]] <- csv_file_exclusive_peptides
  
  # Save regulation-specific files
  if (nrow(tumor_exclusive_results$combined_analysis) > 0 && "regulation_category" %in% colnames(tumor_exclusive_results$combined_analysis)) {
    for (category in unique(tumor_exclusive_results$combined_analysis$regulation_category)) {
      if (!is.na(category)) {
        sheet_data <- tumor_exclusive_results$combined_analysis %>%
          filter(regulation_category == category)
        
        if (nrow(sheet_data) > 0) {
          safe_name <- gsub("[^A-Za-z0-9_]", "_", category)
          csv_file_category <- file.path(table_dir_exclusive, 
                                         paste0(timestamp, "_", dataset_name, "_tumor_exclusive_", safe_name, ".csv"))
          write_csv(sheet_data, csv_file_category)
          saved_files[[paste0("tumor_exclusive_", safe_name, "_csv")]] <- csv_file_category
        }
      }
    }
  }
  
  # 2. Save tumor-upregulated results as CSV files
  if (verbose) cat("Saving tumor-upregulated results...\n")
  
  # Combined analysis CSV
  csv_file_upregulated_combined <- file.path(table_dir_upregulated, 
                                             paste0(timestamp, "_", dataset_name, "_tumor_upregulated_combined.csv"))
  write_csv(tumor_upregulated_results$combined_analysis, csv_file_upregulated_combined)
  saved_files[["tumor_upregulated_combined_csv"]] <- csv_file_upregulated_combined
  
  # Regulation categories CSV
  csv_file_upregulated_categories <- file.path(table_dir_upregulated, 
                                               paste0(timestamp, "_", dataset_name, "_tumor_upregulated_regulation_categories.csv"))
  write_csv(tumor_upregulated_results$regulation_categories, csv_file_upregulated_categories)
  saved_files[["tumor_upregulated_categories_csv"]] <- csv_file_upregulated_categories
  
  # Peptides only CSV
  csv_file_upregulated_peptides <- file.path(table_dir_upregulated, 
                                             paste0(timestamp, "_", dataset_name, "_tumor_upregulated_peptides.csv"))
  write_csv(tumor_upregulated_results$peptides, csv_file_upregulated_peptides)
  saved_files[["tumor_upregulated_peptides_csv"]] <- csv_file_upregulated_peptides
  
  # Save regulation-specific files
  if (nrow(tumor_upregulated_results$combined_analysis) > 0 && "regulation_category" %in% colnames(tumor_upregulated_results$combined_analysis)) {
    for (category in unique(tumor_upregulated_results$combined_analysis$regulation_category)) {
      if (!is.na(category)) {
        sheet_data <- tumor_upregulated_results$combined_analysis %>%
          filter(regulation_category == category)
        
        if (nrow(sheet_data) > 0) {
          safe_name <- gsub("[^A-Za-z0-9_]", "_", category)
          csv_file_category <- file.path(table_dir_upregulated, 
                                         paste0(timestamp, "_", dataset_name, "_tumor_upregulated_", safe_name, ".csv"))
          write_csv(sheet_data, csv_file_category)
          saved_files[[paste0("tumor_upregulated_", safe_name, "_csv")]] <- csv_file_category
        }
      }
    }
  }
  
  # 3. Create summary file
  if (verbose) cat("Creating analysis summary...\n")
  
  summary_data <- data.frame(
    Metric = c(
      "Analysis Type",
      "Dataset Name",
      "Timestamp",
      "Tumor-Exclusive Peptides",
      "Tumor-Upregulated Peptides (>4-fold)",
      "Exclusive: Consistent Upregulation",
      "Exclusive: Mostly Upregulated", 
      "Exclusive: Transcriptome-Specific",
      "Exclusive: Peptide-Specific",
      "Upregulated: Consistent Upregulation",
      "Upregulated: Mostly Upregulated",
      "Upregulated: Transcriptome-Specific", 
      "Upregulated: Peptide-Specific"
    ),
    Value = c(
      "RU148 Multi-Omics Analysis",
      dataset_name,
      timestamp,
      nrow(tumor_exclusive_results$peptides),
      nrow(tumor_upregulated_results$peptides),
      
      # Exclusive categories
      ifelse(nrow(tumor_exclusive_results$regulation_categories) > 0,
             tumor_exclusive_results$regulation_categories %>% 
               filter(regulation_category == "Consistent upregulation") %>% 
               pull(count) %>% ifelse(length(.) > 0, ., 0), 0),
      
      ifelse(nrow(tumor_exclusive_results$regulation_categories) > 0,
             tumor_exclusive_results$regulation_categories %>% 
               filter(regulation_category == "Mostly upregulated") %>% 
               pull(count) %>% ifelse(length(.) > 0, ., 0), 0),
      
      ifelse(nrow(tumor_exclusive_results$regulation_categories) > 0,
             tumor_exclusive_results$regulation_categories %>% 
               filter(regulation_category == "Transcriptome-specific upregulation") %>% 
               pull(count) %>% ifelse(length(.) > 0, ., 0), 0),
      
      ifelse(nrow(tumor_exclusive_results$regulation_categories) > 0,
             tumor_exclusive_results$regulation_categories %>% 
               filter(regulation_category == "Peptide-specific (no omics support)") %>% 
               pull(count) %>% ifelse(length(.) > 0, ., 0), 0),
      
      # Upregulated categories
      ifelse(nrow(tumor_upregulated_results$regulation_categories) > 0,
             tumor_upregulated_results$regulation_categories %>% 
               filter(regulation_category == "Consistent upregulation") %>% 
               pull(count) %>% ifelse(length(.) > 0, ., 0), 0),
      
      ifelse(nrow(tumor_upregulated_results$regulation_categories) > 0,
             tumor_upregulated_results$regulation_categories %>% 
               filter(regulation_category == "Mostly upregulated") %>% 
               pull(count) %>% ifelse(length(.) > 0, ., 0), 0),
      
      ifelse(nrow(tumor_upregulated_results$regulation_categories) > 0,
             tumor_upregulated_results$regulation_categories %>% 
               filter(regulation_category == "Transcriptome-specific upregulation") %>% 
               pull(count) %>% ifelse(length(.) > 0, ., 0), 0),
      
      ifelse(nrow(tumor_upregulated_results$regulation_categories) > 0,
             tumor_upregulated_results$regulation_categories %>% 
               filter(regulation_category == "Peptide-specific (no omics support)") %>% 
               pull(count) %>% ifelse(length(.) > 0, ., 0), 0)
    )
  )
  
  summary_file <- file.path(output_dir, "tables", 
                            paste0(timestamp, "_", dataset_name, "_ru148_analysis_summary.csv"))
  write_csv(summary_data, summary_file)
  saved_files[["summary"]] <- summary_file
  
  if (verbose) {
    cat("✓ RU148 analysis results saved\n")
    cat("Files created:\n")
    for (name in names(saved_files)) {
      cat("  ", name, ":", basename(saved_files[[name]]), "\n")
    }
  }
  
  return(saved_files)
}

#' Main RU148 Analysis Function
#' 
#' Performs complete RU148 multi-omics analysis
#' 
#' @param combined_data Combined immunopeptidome data from main pipeline
#' @param data_directory Directory containing input data
#' @param output_dir Output directory for results
#' @param dataset_name Dataset name
#' @param timestamp Timestamp for file naming
#' @param verbose Print progress messages
#' 
#' @return List containing all analysis results
#' 
perform_ru148_analysis <- function(
    combined_data,
    data_directory,
    output_dir,
    dataset_name,
    timestamp,
    verbose = TRUE
) {
  
  if (verbose) cat("\n=== STARTING RU148 MULTI-OMICS ANALYSIS ===\n")
  
  # Create output directory structure
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Set up file paths relative to data directory
  transcriptome_path <- file.path(dirname(data_directory), "transcriptome", "Normalized_Gene_counts_FLCdb_Panel_1.xlsx")
  lfq_path <- file.path(dirname(data_directory), "proteome", "Levin2023", "adg7038_Table_S2_LFQ.xlsx")
  tmt_path <- file.path(dirname(data_directory), "proteome", "Levin2023", "adg7038_Table_S1_TMT.xlsx")
  
  # 1. Prepare omics data
  omics_data <- prepare_ru148_omics_data(
    combined_data = combined_data,
    lfq_path = lfq_path,
    tmt_path = tmt_path,
    transcriptome_path = transcriptome_path,
    verbose = verbose
  )
  
  # 2. Perform tumor-exclusive analysis
  tumor_exclusive_results <- perform_ru148_tumor_exclusive_analysis(
    omics_data = omics_data,
    verbose = verbose
  )
  
  # 3. Perform tumor-upregulated analysis
  tumor_upregulated_results <- perform_ru148_tumor_upregulated_analysis(
    omics_data = omics_data,
    fold_change_threshold = 1,  # 2-fold change
    verbose = verbose
  )
  
  # 4. Create visualizations
  created_plots <- create_ru148_analysis_plots(
    tumor_exclusive_results = tumor_exclusive_results,
    tumor_upregulated_results = tumor_upregulated_results,
    output_dir = output_dir,
    dataset_name = dataset_name,
    timestamp = timestamp,
    verbose = verbose
  )
  
  # 5. Save results
  saved_files <- save_ru148_analysis_results(
    tumor_exclusive_results = tumor_exclusive_results,
    tumor_upregulated_results = tumor_upregulated_results,
    output_dir = output_dir,
    dataset_name = dataset_name,
    timestamp = timestamp,
    verbose = verbose
  )
  
  if (verbose) {
    cat("\n=== RU148 ANALYSIS COMPLETED ===\n")
    cat("Results saved to:", output_dir, "\n")
    cat("Total plots created:", length(created_plots), "\n")
    cat("Total files saved:", length(saved_files), "\n")
  }
  
  return(list(
    omics_data = omics_data,
    tumor_exclusive_results = tumor_exclusive_results,
    tumor_upregulated_results = tumor_upregulated_results,
    created_plots = created_plots,
    saved_files = saved_files
  ))
}