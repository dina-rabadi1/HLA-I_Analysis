#' HLA-I Analysis Data Processing Functions
#' peptide_data_processing.R
#' Functions for processing peptide data in different experimental contexts
#' @author Your Name
#' @version 1.0

source("peptide_core_utils.R")

# Process immunopeptidome data with flexible sample types
process_immunopeptidome_data <- function(data_files, 
                                         experiment_type = "tumor_normal",
                                         sample_pattern = NULL,
                                         filter_peptide_length = c(8, 12)) {
  
  # Read peptide files
  if (experiment_type == "tumor_normal") {
    # For tumor/normal pairs, use specific pattern extraction
    if (is.null(sample_pattern)) {
      # Default pattern for tumor/normal experiment (e.g., "148T" or "148N")
      sample_pattern <- ".*_([0-9]+[TN])_.*"
    }
  } else if (experiment_type == "multi_sample") {
    # For multiple samples, use more flexible extraction
    if (is.null(sample_pattern)) {
      # Default pattern for multi-sample (e.g., extract sample ID with number)
      sample_pattern <- ".*_([0-9]+[A-Z]?)_.*"
    }
  }
  
  # Read the peptide files
  peptide_data <- read_peptide_files(data_files, sample_id_pattern = sample_pattern)
  
  # Filter by peptide length if specified
  if (!is.null(filter_peptide_length) && length(filter_peptide_length) == 2) {
    peptide_data <- peptide_data %>%
      mutate(Peptide_Length = nchar(Peptide)) %>%
      filter(Peptide_Length >= filter_peptide_length[1] & 
               Peptide_Length <= filter_peptide_length[2])
    
    cat("Filtered peptides to length", filter_peptide_length[1], "to", 
        filter_peptide_length[2], "amino acids\n")
  }
  
  return(peptide_data)
}

# Create a peptide-sample matrix for any experiment type
create_peptide_sample_matrix <- function(peptide_data, value_col = "Intensity", 
                                         id_col = "Sample_ID", peptide_col = "Peptide") {
  
  # Create a summary of peptide detection for each sample
  peptide_summary <- peptide_data %>%
    group_by(!!sym(id_col), !!sym(peptide_col)) %>%
    summarize(
      peptide_length = first(nchar(!!sym(peptide_col))),
      spectral_count = sum(Spectral.Count, na.rm = TRUE),
      total_intensity = sum(!!sym(value_col), na.rm = TRUE),
      protein_ids = paste(unique(na.omit(Protein.ID)), collapse = "; "),
      genes = paste(unique(na.omit(Gene)), collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(!!sym(peptide_col), !!sym(id_col))
  
  # Create a wide format table with all samples
  peptide_matrix <- peptide_summary %>%
    select(!!sym(id_col), !!sym(peptide_col), peptide_length, spectral_count, total_intensity, protein_ids, genes) %>%
    pivot_wider(
      names_from = !!sym(id_col),
      values_from = c(spectral_count, total_intensity, protein_ids, genes),
      values_fill = list(spectral_count = 0, total_intensity = 0, protein_ids = NA, genes = NA)
    )
  
  return(list(
    summary = peptide_summary,
    matrix = peptide_matrix
  ))
}

# Process tumor-normal comparison data
process_tumor_normal_data <- function(peptide_matrix, tumor_id, normal_id, 
                                      intensity_prefix = "total_intensity_") {
  
  # Extract column names for intensity values
  tumor_col <- paste0(intensity_prefix, tumor_id)
  normal_col <- paste0(intensity_prefix, normal_id)
  
  # Check if columns exist
  if (!(tumor_col %in% colnames(peptide_matrix) && normal_col %in% colnames(peptide_matrix))) {
    stop("Tumor or normal columns not found in data matrix. Check sample IDs and prefix.")
  }
  
  # Calculate fold changes and create differential expression analysis
  diff_expr <- peptide_matrix %>%
    mutate(
      # Replace zero with small value to prevent division by zero or Inf
      normal_intensity_adj = ifelse(!!sym(normal_col) == 0, 0.1, !!sym(normal_col)),
      tumor_intensity_adj = ifelse(!!sym(tumor_col) == 0, 0.1, !!sym(tumor_col)),
      
      # Calculate fold changes (log2)
      log2_fold_change = log2(tumor_intensity_adj / normal_intensity_adj),
      
      # Determine if peptide is specific to tumor or normal
      detection_status = case_when(
        !!sym(tumor_col) > 0 & !!sym(normal_col) == 0 ~ "Tumor-specific",
        !!sym(normal_col) > 0 & !!sym(tumor_col) == 0 ~ "Normal-specific",
        !!sym(tumor_col) > 0 & !!sym(normal_col) > 0 ~ "Detected in both",
        TRUE ~ "Not detected"
      ),
      
      # Simplified category for plotting
      expression_category = case_when(
        log2_fold_change > 1 ~ "Up in Tumor (FC > 2)",
        log2_fold_change < -1 ~ "Down in Tumor (FC < 0.5)",
        TRUE ~ "Similar (-1 < log2FC < 1)"
      )
    ) %>%
    # Clean up protein and gene info
    mutate(
      genes_combined = coalesce(!!sym(paste0("genes_", tumor_id)), 
                                !!sym(paste0("genes_", normal_id))),
      proteins_combined = coalesce(!!sym(paste0("protein_ids_", tumor_id)), 
                                   !!sym(paste0("protein_ids_", normal_id)))
    ) %>%
    # Extract primary gene for later comparison with transcriptome
    mutate(
      primary_gene = sapply(strsplit(genes_combined, ";\\s*"), function(x) trimws(x[1]))
    ) %>%
    # Sort by fold change for easier viewing
    arrange(desc(log2_fold_change))
  
  return(diff_expr)
}

# Process multi-sample analysis
process_multi_sample_data <- function(peptide_matrix, samples = NULL, 
                                      intensity_prefix = "total_intensity_") {
  
  # If no samples provided, extract all samples
  if (is.null(samples)) {
    # Extract all sample IDs from column names
    intensity_cols <- grep(paste0("^", intensity_prefix), colnames(peptide_matrix), value = TRUE)
    samples <- gsub(intensity_prefix, "", intensity_cols)
  }
  
  # Verify all samples exist in the dataset
  for (sample in samples) {
    if (!(paste0(intensity_prefix, sample) %in% colnames(peptide_matrix))) {
      stop("Sample '", sample, "' not found in data matrix. Check sample IDs and prefix.")
    }
  }
  
  # Calculate shared vs. private peptides across samples
  multi_sample_analysis <- peptide_matrix %>%
    # Add info about which samples each peptide is found in
    rowwise() %>%
    mutate(
      detected_samples = sum(sapply(samples, function(s) {
        col_name <- paste0(intensity_prefix, s)
        ifelse(!!sym(col_name) > 0, 1, 0)
      })),
      sample_list = paste(sort(samples[sapply(samples, function(s) {
        col_name <- paste0(intensity_prefix, s)
        !!sym(col_name) > 0
      })]), collapse = ", "),
      is_private = detected_samples == 1,
      is_shared = detected_samples > 1,
      sharing_category = case_when(
        detected_samples == 0 ~ "Not detected",
        detected_samples == 1 ~ "Private (1 sample)",
        detected_samples <= 3 ~ paste0("Shared (", detected_samples, " samples)"),
        detected_samples <= length(samples)/2 ~ "Moderately shared",
        detected_samples <= length(samples)-1 ~ "Highly shared",
        detected_samples == length(samples) ~ "Universal",
        TRUE ~ "Unknown"
      )
    )
  
  return(multi_sample_analysis)
}

# Process transcriptome data for tumor-normal comparison
process_transcriptome_data <- function(transcriptome_file, tumor_col = NULL, normal_col = NULL,
                                       gene_col = "symbol") {
  
  # Read transcriptome data
  transcriptome_data <- read_excel(transcriptome_file)
  
  # Print column names to help with debugging
  cat("Transcriptome data columns:", paste(colnames(transcriptome_data), collapse=", "), "\n")
  
  # Process data
  transcriptome_processed <- transcriptome_data
  
  # Handle gene column name variations
  if (!(gene_col %in% colnames(transcriptome_processed))) {
    # Try alternative names
    if ("gene_symbol" %in% colnames(transcriptome_processed)) {
      transcriptome_processed <- transcriptome_processed %>%
        rename(!!gene_col := gene_symbol)
    } else if ("Symbol" %in% colnames(transcriptome_processed)) {
      transcriptome_processed <- transcriptome_processed %>%
        rename(!!gene_col := Symbol)
    } else {
      # Create an empty symbol column if needed
      transcriptome_processed[[gene_col]] <- NA_character_
      cat("Warning: No suitable gene symbol column found\n")
    }
  }
  
  # Auto-detect tumor and normal columns if not provided
  if (is.null(tumor_col) && is.null(normal_col)) {
    # Look for common patterns in column names
    tumor_cols <- grep("_T[0-9]*$|T_Average|Tumor$", colnames(transcriptome_processed), value = TRUE)
    normal_cols <- grep("_N$|Normal$", colnames(transcriptome_processed), value = TRUE)
    
    if (length(tumor_cols) > 0 && length(normal_cols) > 0) {
      # If multiple matches, use the first
      if (length(tumor_cols) > 1) {
        # Prefer average if available
        avg_col <- grep("Average", tumor_cols, value = TRUE)
        if (length(avg_col) > 0) {
          tumor_col <- avg_col[1]
        } else {
          tumor_col <- tumor_cols[1]
        }
      } else {
        tumor_col <- tumor_cols[1]
      }
      
      normal_col <- normal_cols[1]
      
      cat("Auto-detected tumor column:", tumor_col, "\n")
      cat("Auto-detected normal column:", normal_col, "\n")
    } else {
      stop("Could not auto-detect tumor and normal columns. Please specify them.")
    }
  }
  
  # Calculate log2 fold change
  if (tumor_col %in% colnames(transcriptome_processed) && 
      normal_col %in% colnames(transcriptome_processed)) {
    
    transcriptome_processed <- transcriptome_processed %>%
      mutate(
        log2_fold_change_transcriptome = log2(
          ifelse(!!sym(tumor_col) == 0, 0.1, !!sym(tumor_col)) / 
            ifelse(!!sym(normal_col) == 0, 0.1, !!sym(normal_col))
        )
      )
  } else {
    transcriptome_processed$log2_fold_change_transcriptome <- NA_real_
    cat("Warning: Unable to calculate transcriptome log2 fold change\n")
  }
  
  return(transcriptome_processed)
}

# Process proteomics data
process_proteomics_data <- function(proteomics_file, data_type = "LFQ", sheet_name = NULL,
                                    p_value_cutoff = 0.05, gene_col = "Gene_Name",
                                    fold_change_col = "Log2_Difference", p_value_col = "P.value") {
  
  # Set suffix for output columns
  suffix <- tolower(data_type)
  
  # Read the proteomics data
  if (!is.null(sheet_name)) {
    proteomics_data <- read_excel(proteomics_file, sheet = sheet_name)
  } else {
    proteomics_data <- read_excel(proteomics_file)
  }
  
  # Print column names to help with debugging
  cat(data_type, "proteomics data columns:", 
      paste(colnames(proteomics_data), collapse=", "), "\n")
  
  # Process data
  proteomics_processed <- proteomics_data %>%
    # Ensure consistent column names
    rename_with(~ gsub(" ", "_", .), everything())
  
  # Adjust column names if they don't match expected names
  if (!(gene_col %in% colnames(proteomics_processed))) {
    # Try alternative names
    gene_alternatives <- c("Gene", "gene", "Gene_name", "gene_name", "Symbol", "symbol")
    for (alt in gene_alternatives) {
      if (alt %in% colnames(proteomics_processed)) {
        proteomics_processed <- proteomics_processed %>%
          rename(!!gene_col := !!alt)
        break
      }
    }
  }
  
  if (!(fold_change_col %in% colnames(proteomics_processed))) {
    # Try alternative names
    fc_alternatives <- c("log2FC", "Log2FC", "log2_FC", "Log2_FC", "Fold_Change", "log2FoldChange")
    for (alt in fc_alternatives) {
      if (alt %in% colnames(proteomics_processed)) {
        proteomics_processed <- proteomics_processed %>%
          rename(!!fold_change_col := !!alt)
        break
      }
    }
  }
  
  if (!(p_value_col %in% colnames(proteomics_processed))) {
    # Try alternative names
    p_alternatives <- c("pvalue", "p_value", "PValue", "P_Value", "Pval", "padj", "adj_P_Val")
    for (alt in p_alternatives) {
      if (alt %in% colnames(proteomics_processed)) {
        proteomics_processed <- proteomics_processed %>%
          rename(!!p_value_col := !!alt)
        break
      }
    }
  }
  
  # Filter and standardize
  if (all(c(gene_col, fold_change_col, p_value_col) %in% colnames(proteomics_processed))) {
    proteomics_processed <- proteomics_processed %>%
      # Filter for significant entries
      filter(!!sym(p_value_col) <= p_value_cutoff) %>%
      # Clean up and standardize
      mutate(
        !!gene_col := trimws(!!sym(gene_col)),
        !!paste0("log2_fold_change_", suffix) := !!sym(fold_change_col),
        !!paste0("p_value_", suffix) := !!sym(p_value_col),
        !!paste0("protein_category_", suffix) := case_when(
          !!sym(fold_change_col) > 1 ~ "Up in Tumor (FC > 2)",
          !!sym(fold_change_col) < -1 ~ "Down in Tumor (FC < 0.5)",
          TRUE ~ "Similar (-1 < log2FC < 1)"
        )
      )
    
    # Select relevant columns
    result_cols <- c(
      gene_col,
      paste0("log2_fold_change_", suffix),
      paste0("p_value_", suffix),
      paste0("protein_category_", suffix)
    )
    
    # Add protein name column if available
    protein_name_cols <- c("Protein_Name", "Protein_name", "protein_name", "Protein")
    for (col in protein_name_cols) {
      if (col %in% colnames(proteomics_processed)) {
        result_cols <- c(result_cols, col)
        break
      }
    }
    
    proteomics_processed <- proteomics_processed %>%
      select(all_of(result_cols))
    
    cat("Processed", nrow(proteomics_processed), "significant entries from", 
        data_type, "proteome data\n")
  } else {
    stop("Required columns not found in proteomics data. Please check your file.")
  }
  
  return(proteomics_processed)
}

# Identify fusion peptides
identify_fusion_peptides <- function(peptide_data, fusion_sequence, 
                                     fusion_parts = NULL, junction_position = NULL) {
  
  # If fusion parts and junction position not provided, ask for them
  if (is.null(fusion_parts) || is.null(junction_position)) {
    stop("Please provide fusion_parts (list with 'part1' and 'part2') and junction_position")
  }
  
  # Extract fusion parts
  dnajb1_part <- fusion_parts$part1
  prkaca_part <- fusion_parts$part2
  
  # Define function to check if a peptide spans the fusion junction
  is_fusion_junction_peptide <- function(peptide_seq) {
    # Check if peptide is entirely within the fusion protein
    from_fusion <- grepl(peptide_seq, fusion_sequence, fixed = TRUE)
    
    # Check if the peptide spans the fusion junction
    spans_junction <- FALSE
    
    if (from_fusion && nchar(peptide_seq) >= 6) { # Require at least 6 AA to be meaningful
      for (i in 3:(nchar(peptide_seq) - 3)) { # At least 3 AA from each side
        left_part <- substr(peptide_seq, 1, i)
        right_part <- substr(peptide_seq, i + 1, nchar(peptide_seq))
        
        # Check if left part is in part1 and right part in part2
        if (grepl(left_part, dnajb1_part, fixed = TRUE) && 
            grepl(right_part, prkaca_part, fixed = TRUE)) {
          
          # Additional check to ensure left part aligns with end of part1
          left_pos <- gregexpr(left_part, dnajb1_part, fixed = TRUE)[[1]]
          if (length(left_pos) > 0 && any(left_pos + nchar(left_part) - 1 >= nchar(dnajb1_part) - 5)) {
            
            # Additional check to ensure right part aligns with start of part2
            right_pos <- gregexpr(right_part, prkaca_part, fixed = TRUE)[[1]]
            if (length(right_pos) > 0 && any(right_pos <= 5)) {
              spans_junction <- TRUE
              break
            }
          }
        }
      }
    }
    
    # Check if peptide is from either part of the fusion protein
    from_part1 <- grepl(peptide_seq, dnajb1_part, fixed = TRUE)
    from_part2 <- grepl(peptide_seq, prkaca_part, fixed = TRUE)
    
    return(list(
      spans_junction = spans_junction,
      from_part1 = from_part1,
      from_part2 = from_part2,
      from_fusion = from_fusion | spans_junction | from_part1 | from_part2
    ))
  }
  
  # Apply to all peptides
  peptide_data_with_fusion <- peptide_data %>%
    rowwise() %>%
    mutate(
      fusion_info = list(is_fusion_junction_peptide(Peptide)),
      from_fusion = fusion_info$from_fusion,
      spans_junction = fusion_info$spans_junction,
      from_part1 = fusion_info$from_part1,
      from_part2 = fusion_info$from_part2,
      fusion_peptide_type = case_when(
        spans_junction ~ "Junction-spanning",
        from_part1 ~ "Part1",
        from_part2 ~ "Part2",
        TRUE ~ "Not from fusion"
      )
    ) %>%
    select(-fusion_info)
  
  # Extract only fusion peptides
  fusion_peptides <- peptide_data_with_fusion %>%
    filter(from_fusion) %>%
    arrange(desc(spans_junction))
  
  return(list(
    all_data_with_fusion = peptide_data_with_fusion,
    fusion_peptides = fusion_peptides
  ))
}