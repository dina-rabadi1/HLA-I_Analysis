#' HLA-I Peptide Analysis Integration and Report Generation
#' peptide_integration.R
#' Functions for multi-omics data integration and report generation
#' @author Your Name
#' @version 1.0

source("peptide_core_utils.R")
source("peptide_data_processing.R")

# Combine immunopeptidome data with transcriptome data
integrate_immuno_transcriptome <- function(immuno_data, transcriptome_data,
                                           gene_col_immuno = "primary_gene",
                                           gene_col_transcriptome = "symbol") {
  
  # Ensure correct column names
  if(!gene_col_immuno %in% colnames(immuno_data)) {
    stop("Gene column '", gene_col_immuno, "' not found in immunopeptidome data")
  }
  
  if(!gene_col_transcriptome %in% colnames(transcriptome_data)) {
    stop("Gene column '", gene_col_transcriptome, "' not found in transcriptome data")
  }
  
  # Join the immunopeptidome and transcriptome data based on gene symbol
  combined_data <- immuno_data %>%
    left_join(
      transcriptome_data,
      by = setNames(gene_col_transcriptome, gene_col_immuno),
      relationship = "many-to-many"  # Explicitly acknowledge many-to-many relationship
    )
  
  # Add correlation statistics if fold change columns exist
  if("log2_fold_change" %in% colnames(immuno_data) && 
     "log2_fold_change_transcriptome" %in% colnames(combined_data)) {
    
    combined_data <- combined_data %>%
      mutate(
        immuno_trans_correlation = log2_fold_change * log2_fold_change_transcriptome,
        expression_category = case_when(
          is.na(log2_fold_change_transcriptome) ~ "No transcriptome data",
          log2_fold_change_transcriptome > 1 & log2_fold_change > 1 ~ "Up in both",
          log2_fold_change_transcriptome < -1 & log2_fold_change < -1 ~ "Down in both",
          log2_fold_change_transcriptome > 1 & log2_fold_change < -1 ~ "Up in transcriptome, down in immunopeptidome",
          log2_fold_change_transcriptome < -1 & log2_fold_change > 1 ~ "Down in transcriptome, up in immunopeptidome",
          TRUE ~ "No significant change"
        )
      )
  }
  
  return(combined_data)
}

# Combine immunopeptidome data with proteomics data (LFQ or TMT)
integrate_immuno_proteomics <- function(immuno_data, proteomics_data, proteomics_type = "LFQ",
                                        gene_col_immuno = "primary_gene",
                                        gene_col_proteomics = "Gene_Name") {
  
  # Set suffix for column naming
  suffix <- tolower(proteomics_type)
  
  # Ensure correct column names
  if(!gene_col_immuno %in% colnames(immuno_data)) {
    stop("Gene column '", gene_col_immuno, "' not found in immunopeptidome data")
  }
  
  if(!gene_col_proteomics %in% colnames(proteomics_data)) {
    stop("Gene column '", gene_col_proteomics, "' not found in proteomics data")
  }
  
  # Join the immunopeptidome and proteomics data based on gene symbol
  combined_data <- immuno_data %>%
    left_join(
      proteomics_data,
      by = setNames(gene_col_proteomics, gene_col_immuno),
      relationship = "many-to-many"  # Explicitly acknowledge many-to-many relationship
    )
  
  # Add correlation statistics if fold change columns exist
  fc_col_proteomics <- paste0("log2_fold_change_", suffix)
  
  if("log2_fold_change" %in% colnames(immuno_data) && 
     fc_col_proteomics %in% colnames(combined_data)) {
    
    correlation_col <- paste0("immuno_", suffix, "_correlation")
    category_col <- paste0("expression_category_", suffix)
    
    combined_data <- combined_data %>%
      mutate(
        !!correlation_col := case_when(
          !is.na(!!sym(fc_col_proteomics)) ~ log2_fold_change * !!sym(fc_col_proteomics),
          TRUE ~ NA_real_
        ),
        !!category_col := case_when(
          is.na(!!sym(fc_col_proteomics)) ~ paste0("No ", proteomics_type, " proteome data"),
          !!sym(fc_col_proteomics) > 1 & log2_fold_change > 1 ~ "Up in both",
          !!sym(fc_col_proteomics) < -1 & log2_fold_change < -1 ~ "Down in both",
          !!sym(fc_col_proteomics) > 1 & log2_fold_change < -1 ~ paste0("Up in ", proteomics_type, ", down in immunopeptidome"),
          !!sym(fc_col_proteomics) < -1 & log2_fold_change > 1 ~ paste0("Down in ", proteomics_type, ", up in immunopeptidome"),
          TRUE ~ "No significant change"
        )
      )
  }
  
  return(combined_data)
}

# Perform multi-omics integration (4-way)
integrate_multi_omics <- function(immuno_data, transcriptome_data = NULL, 
                                  lfq_data = NULL, tmt_data = NULL) {
  
  # Start with immunopeptidome data
  combined_data <- immuno_data
  
  # Add transcriptome data if provided
  if(!is.null(transcriptome_data)) {
    combined_data <- integrate_immuno_transcriptome(combined_data, transcriptome_data)
  }
  
  # Add LFQ proteomics data if provided
  if(!is.null(lfq_data)) {
    combined_data <- integrate_immuno_proteomics(combined_data, lfq_data, "LFQ")
  }
  
  # Add TMT proteomics data if provided
  if(!is.null(tmt_data)) {
    combined_data <- integrate_immuno_proteomics(combined_data, tmt_data, "TMT")
  }
  
  return(combined_data)
}

# Define public neoantigen criteria for multi-omics data
identify_public_neoantigens <- function(multi_omics_data, fc_threshold = 1) {
  
  # Check if we have the necessary data
  has_transcriptome <- "log2_fold_change_transcriptome" %in% colnames(multi_omics_data)
  has_lfq <- "log2_fold_change_lfq" %in% colnames(multi_omics_data)
  has_tmt <- "log2_fold_change_tmt" %in% colnames(multi_omics_data)
  
  # Define criteria based on available data
  multi_omics_data <- multi_omics_data %>%
    mutate(
      # Criteria 1: Upregulated in all available datasets (log2FC > 1)
      upregulated_in_all = case_when(
        # If we have all four datasets
        has_transcriptome & has_lfq & has_tmt ~ 
          log2_fold_change > fc_threshold & 
          !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold &
          !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > fc_threshold &
          !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > fc_threshold,
        
        # If we have three datasets (immunopeptidome, transcriptome, and one proteomics)
        has_transcriptome & has_lfq & !has_tmt ~
          log2_fold_change > fc_threshold & 
          !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold &
          !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > fc_threshold,
        
        has_transcriptome & !has_lfq & has_tmt ~
          log2_fold_change > fc_threshold & 
          !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold &
          !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > fc_threshold,
        
        # If we have two datasets (immunopeptidome and one other)
        has_transcriptome & !has_lfq & !has_tmt ~
          log2_fold_change > fc_threshold & 
          !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold,
        
        !has_transcriptome & has_lfq & !has_tmt ~
          log2_fold_change > fc_threshold & 
          !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > fc_threshold,
        
        !has_transcriptome & !has_lfq & has_tmt ~
          log2_fold_change > fc_threshold & 
          !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > fc_threshold,
        
        # Default
        TRUE ~ FALSE
      ),
      
      # Criteria 2: Upregulated in immunopeptidome and at least X other datasets
      # where X depends on how many datasets we have
      upregulated_in_multiple = case_when(
        log2_fold_change <= fc_threshold ~ FALSE,  # Must be upregulated in immunopeptidome
        
        # If we have all four datasets, require at least 2 others
        has_transcriptome & has_lfq & has_tmt ~ 
          sum(
            !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold,
            !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > fc_threshold,
            !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > fc_threshold
          ) >= 2,
        
        # If we have three datasets, require at least 1 other
        has_transcriptome & has_lfq & !has_tmt ~ 
          sum(
            !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold,
            !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > fc_threshold
          ) >= 1,
        
        has_transcriptome & !has_lfq & has_tmt ~ 
          sum(
            !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold,
            !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > fc_threshold
          ) >= 1,
        
        !has_transcriptome & has_lfq & has_tmt ~ 
          sum(
            !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > fc_threshold,
            !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > fc_threshold
          ) >= 1,
        
        # If we have two datasets total, require the other one to be upregulated
        has_transcriptome & !has_lfq & !has_tmt ~ 
          !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold,
        
        !has_transcriptome & has_lfq & !has_tmt ~ 
          !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > fc_threshold,
        
        !has_transcriptome & !has_lfq & has_tmt ~ 
          !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > fc_threshold,
        
        # Default
        TRUE ~ FALSE
      ),
      
      # Criteria 3: Tumor-specific in immunopeptidome and upregulated in other datasets
      tumor_specific_upregulated = case_when(
        !("detection_status" %in% colnames(multi_omics_data)) ~ FALSE,
        detection_status != "Tumor-specific" ~ FALSE,
        
        # If we have all four datasets, require at least 2 others
        has_transcriptome & has_lfq & has_tmt ~ 
          sum(
            !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold,
            !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > fc_threshold,
            !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > fc_threshold
          ) >= 2,
        
        # If we have three datasets, require at least 1 other
        has_transcriptome & has_lfq & !has_tmt ~ 
          sum(
            !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold,
            !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > fc_threshold
          ) >= 1,
        
        has_transcriptome & !has_lfq & has_tmt ~ 
          sum(
            !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold,
            !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > fc_threshold
          ) >= 1,
        
        !has_transcriptome & has_lfq & has_tmt ~ 
          sum(
            !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > fc_threshold,
            !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > fc_threshold
          ) >= 1,
        
        # If we have two datasets total, require the other one to be upregulated
        has_transcriptome & !has_lfq & !has_tmt ~ 
          !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > fc_threshold,
        
        !has_transcriptome & has_lfq & !has_tmt ~ 
          !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > fc_threshold,
        
        !has_transcriptome & !has_lfq & has_tmt ~ 
          !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > fc_threshold,
        
        # Default
        TRUE ~ FALSE
      )
    )
  
  # Calculate combined score and overall classification
  multi_omics_data <- multi_omics_data %>%
    mutate(
      # Combined public neoantigen score (sum of criteria)
      public_neoantigen_score = as.integer(upregulated_in_all) * 3 + 
        as.integer(upregulated_in_multiple) * 2 + 
        as.integer(tumor_specific_upregulated) * 1,
      
      # Classify as potential public neoantigen if any criteria are met
      potential_public_neoantigen = public_neoantigen_score > 0,
      
      # Classification label
      public_neoantigen_classification = case_when(
        upregulated_in_all ~ paste0("Tier 1 (Upregulated in all ", 
                                    sum(1, has_transcriptome, has_lfq, has_tmt), " datasets)"),
        upregulated_in_multiple ~ paste0("Tier 2 (Upregulated in immunopeptidome + others)"),
        tumor_specific_upregulated ~ "Tier 3 (Tumor-specific + upregulated in others)",
        TRUE ~ "Not a potential public neoantigen"
      )
    )
  
  # Extract only public neoantigens
  public_neoantigens <- multi_omics_data %>%
    filter(potential_public_neoantigen) %>%
    arrange(desc(public_neoantigen_score), desc(log2_fold_change))
  
  # Return both datasets
  return(list(
    all_data = multi_omics_data,
    public_neoantigens = public_neoantigens
  ))
}

# Generate comprehensive Excel report
generate_excel_report <- function(data_list, output_file, header_style = NULL) {
  # Create a default header style if not provided
  if (is.null(header_style)) {
    header_style <- createStyle(textDecoration = "bold", fgFill = "#D9D9D9")
  }
  
  # Create workbook
  wb <- createWorkbook()
  
  # Add each data sheet with formatting
  for (name in names(data_list)) {
    # Skip empty datasets
    if (is.null(data_list[[name]]) || nrow(data_list[[name]]) == 0) {
      next
    }
    
    # Sanitize sheet name (max 31 chars, no special chars)
    sheet_name <- substr(gsub("[^A-Za-z0-9_]", "_", name), 1, 31)
    
    # Add worksheet
    addWorksheet(wb, sheet_name)
    
    # Write data
    writeData(wb, sheet_name, data_list[[name]], headerStyle = header_style)
    
    # Freeze header row
    freezePane(wb, sheet_name, firstRow = TRUE)
    
    # Auto-adjust column widths
    setColWidths(wb, sheet_name, cols = 1:ncol(data_list[[name]]), widths = "auto")
    
    # Add conditional formatting for important columns
    # Log2 fold change columns - highlight >1 (red) and <-1 (blue)
    fc_cols <- grep("fold_change|FC", colnames(data_list[[name]]), value = TRUE)
    for (col in fc_cols) {
      col_idx <- which(colnames(data_list[[name]]) == col)
      conditionalFormatting(wb, sheet_name, cols = col_idx, 
                            rows = 2:(nrow(data_list[[name]]) + 1),
                            rule = ">1", 
                            style = createStyle(bgFill = "#FFCCCC"))
      conditionalFormatting(wb, sheet_name, cols = col_idx, 
                            rows = 2:(nrow(data_list[[name]]) + 1),
                            rule = "<-1", 
                            style = createStyle(bgFill = "#CCCCFF"))
    }
    
    # P-value columns - highlight <0.05 (green)
    pval_cols <- grep("p_value|p_val|p\\.value|pval", colnames(data_list[[name]]), value = TRUE, ignore.case = TRUE)
    for (col in pval_cols) {
      col_idx <- which(colnames(data_list[[name]]) == col)
      conditionalFormatting(wb, sheet_name, cols = col_idx, 
                            rows = 2:(nrow(data_list[[name]]) + 1),
                            rule = "<0.05", 
                            style = createStyle(bgFill = "#E2EFDA"))
    }
    
    # Public neoantigen columns - highlight TRUE (yellow)
    neoant_cols <- grep("public_neoantigen|neoantigen", colnames(data_list[[name]]), value = TRUE, ignore.case = TRUE)
    for (col in neoant_cols) {
      if (!is.numeric(data_list[[name]][[col]])) {
        col_idx <- which(colnames(data_list[[name]]) == col)
        conditionalFormatting(wb, sheet_name, cols = col_idx, 
                              rows = 2:(nrow(data_list[[name]]) + 1),
                              rule = "==TRUE", 
                              style = createStyle(bgFill = "#FFEB9C"))
      }
    }
  }
  
  # Save workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat("Excel report saved to:", output_file, "\n")
}

# Function to prepare Excel data for different analysis types
prepare_excel_data <- function(data, analysis_type = "tumor_normal") {
  excel_sheets <- list()
  
  if (analysis_type == "tumor_normal") {
    # Basic sheets for tumor-normal analysis
    excel_sheets[["Combined_Analysis"]] <- data
    
    # Extract only immunopeptidome data if it exists
    if ("log2_fold_change" %in% colnames(data)) {
      immuno_cols <- c("Peptide", "peptide_length", "primary_gene", "genes_combined", "proteins_combined",
                       grep("spectral_count|intensity|detection_status|fold_change$|expression_category$", 
                            colnames(data), value = TRUE))
      
      excel_sheets[["Immunopeptidome_Only"]] <- data %>%
        select(all_of(intersect(immuno_cols, colnames(data))))
    }
    
    # Extract only transcriptome data if it exists
    if ("log2_fold_change_transcriptome" %in% colnames(data)) {
      trans_cols <- c("primary_gene", 
                      grep("transcriptome", colnames(data), value = TRUE))
      
      excel_sheets[["Transcriptome_Only"]] <- data %>%
        select(all_of(intersect(trans_cols, colnames(data)))) %>%
        distinct()
    }
    
    # Extract only LFQ proteome data if it exists
    if ("log2_fold_change_lfq" %in% colnames(data)) {
      lfq_cols <- c("primary_gene", 
                    grep("lfq", colnames(data), value = TRUE, ignore.case = TRUE))
      
      excel_sheets[["LFQ_Proteome_Only"]] <- data %>%
        select(all_of(intersect(lfq_cols, colnames(data)))) %>%
        distinct()
    }
    
    # Extract only TMT proteome data if it exists
    if ("log2_fold_change_tmt" %in% colnames(data)) {
      tmt_cols <- c("primary_gene", 
                    grep("tmt", colnames(data), value = TRUE, ignore.case = TRUE))
      
      excel_sheets[["TMT_Proteome_Only"]] <- data %>%
        select(all_of(intersect(tmt_cols, colnames(data)))) %>%
        distinct()
    }
    
    # Add pairwise comparison sheets
    # Immunopeptidome vs Transcriptome
    if (all(c("log2_fold_change", "log2_fold_change_transcriptome") %in% colnames(data))) {
      excel_sheets[["Immuno_vs_Transcriptome"]] <- data %>%
        filter(!is.na(log2_fold_change_transcriptome)) %>%
        select(
          Peptide, primary_gene, peptide_length,
          log2_fold_change, log2_fold_change_transcriptome,
          any_of(c("expression_category", "expression_category_transcriptome", "detection_status"))
        ) %>%
        arrange(desc(log2_fold_change))
    }
    
    # Immunopeptidome vs LFQ
    if (all(c("log2_fold_change", "log2_fold_change_lfq") %in% colnames(data))) {
      excel_sheets[["Immuno_vs_LFQ"]] <- data %>%
        filter(!is.na(log2_fold_change_lfq)) %>%
        select(
          Peptide, primary_gene, peptide_length,
          log2_fold_change, log2_fold_change_lfq,
          any_of(c("p_value_lfq", "expression_category_lfq", "detection_status"))
        ) %>%
        arrange(desc(log2_fold_change))
    }
    
    # Immunopeptidome vs TMT
    if (all(c("log2_fold_change", "log2_fold_change_tmt") %in% colnames(data))) {
      excel_sheets[["Immuno_vs_TMT"]] <- data %>%
        filter(!is.na(log2_fold_change_tmt)) %>%
        select(
          Peptide, primary_gene, peptide_length,
          log2_fold_change, log2_fold_change_tmt,
          any_of(c("p_value_tmt", "expression_category_tmt", "detection_status"))
        ) %>%
        arrange(desc(log2_fold_change))
    }
    
    # Add status-based sheets
    if ("detection_status" %in% colnames(data)) {
      # Tumor-specific peptides
      excel_sheets[["Tumor_Specific_Peptides"]] <- data %>%
        filter(detection_status == "Tumor-specific")
      
      # Normal-specific peptides
      excel_sheets[["Normal_Specific_Peptides"]] <- data %>%
        filter(detection_status == "Normal-specific")
    }
    
    # Add fusion protein sheets if available
    if ("from_fusion" %in% colnames(data)) {
      excel_sheets[["Fusion_Peptides"]] <- data %>%
        filter(from_fusion == TRUE)
      
      # Junction-spanning fusion peptides
      if ("spans_junction" %in% colnames(data)) {
        excel_sheets[["Junction_Spanning_Peptides"]] <- data %>%
          filter(spans_junction == TRUE)
      }
    }
    
    # Add public neoantigen sheets if available
    if ("potential_public_neoantigen" %in% colnames(data)) {
      excel_sheets[["Public_Neoantigens"]] <- data %>%
        filter(potential_public_neoantigen == TRUE) %>%
        arrange(desc(public_neoantigen_score))
      
      # Add tier-specific sheets
      if ("public_neoantigen_classification" %in% colnames(data)) {
        tier_values <- unique(data$public_neoantigen_classification)
        
        for (tier in tier_values) {
          if (!grepl("Not a potential", tier)) {
            sheet_name <- paste0("Tier_", substr(tier, 5, 5), "_Neoantigens")
            excel_sheets[[sheet_name]] <- data %>%
              filter(public_neoantigen_classification == tier)
          }
        }
      }
    }
  } else if (analysis_type == "multi_sample") {
    # Basic sheets for multi-sample analysis
    excel_sheets[["All_Peptides"]] <- data
    
    # Private vs shared peptides
    if ("is_private" %in% colnames(data)) {
      excel_sheets[["Private_Peptides"]] <- data %>%
        filter(is_private == TRUE)
      
      excel_sheets[["Shared_Peptides"]] <- data %>%
        filter(is_shared == TRUE)
    }
    
    # Different sharing categories
    if ("sharing_category" %in% colnames(data)) {
      sharing_categories <- unique(data$sharing_category)
      
      for (category in sharing_categories) {
        if (category != "Not detected") {
          sheet_name <- gsub(" ", "_", category)
          sheet_name <- gsub("[()]", "", sheet_name)
          
          excel_sheets[[sheet_name]] <- data %>%
            filter(sharing_category == category)
        }
      }
    }
    
    # Add fusion protein sheets if available
    if ("from_fusion" %in% colnames(data)) {
      excel_sheets[["Fusion_Peptides"]] <- data %>%
        filter(from_fusion == TRUE)
      
      # Junction-spanning fusion peptides
      if ("spans_junction" %in% colnames(data)) {
        excel_sheets[["Junction_Spanning_Peptides"]] <- data %>%
          filter(spans_junction == TRUE)
      }
    }
  }
  
  return(excel_sheets)
}