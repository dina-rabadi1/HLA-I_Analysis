# Multi-omics integration module for shared peptide analysis
# src/modules/multi_omics_integration.R

#' Run comprehensive multi-omics integration analysis for shared peptides
#' 
#' @param peptide_data Data frame containing peptide data
#' @param transcriptome_data Processed transcriptome data frame (optional)
#' @param lfq_data Processed LFQ proteome data frame (optional)
#' @param tmt_data Processed TMT proteome data frame (optional)
#' @param config Configuration list with settings
#' @return List containing integrated analysis results
run_multi_omics_integration <- function(peptide_data, transcriptome_data = NULL, 
                                        lfq_data = NULL, tmt_data = NULL, config) {
  cat("\n--- Running Multi-Omics Integration for Shared Peptides ---\n")
  
  # Extract relevant config parameters
  sample_ids <- config$samples$include
  exclude_samples <- config$samples$exclude_from_shared
  
  # Remove excluded samples
  sample_ids <- setdiff(sample_ids, exclude_samples)
  
  # Filter data to only selected samples
  filtered_data <- peptide_data %>%
    dplyr::filter(SampleID %in% sample_ids)
  
  if (nrow(filtered_data) == 0) {
    stop("No peptide data found for the selected samples")
  }
  
  # Create peptide sharing report
  peptide_report <- filtered_data %>%
    dplyr::select(Peptide, SampleID, `Peptide Length`, Gene, Protein) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Peptide, `Peptide Length`) %>%
    dplyr::summarize(
      sample_list = paste(sort(SampleID), collapse = ", "),
      sample_count = dplyr::n_distinct(SampleID),
      genes = paste(unique(na.omit(Gene)), collapse = "; "),
      proteins = paste(unique(na.omit(Protein)), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(sample_count), Peptide)
  
  # Add primary gene for each peptide (first gene in list)
  peptide_report <- peptide_report %>%
    dplyr::mutate(
      primary_gene = sapply(strsplit(genes, ";"), function(x) trimws(x[1]))
    )
  
  # Define thresholds for shared peptides
  if (!is.null(config$analysis$transcriptome$sharing_thresholds)) {
    sharing_thresholds <- config$analysis$transcriptome$sharing_thresholds
  } else {
    # Default thresholds based on number of samples
    n_samples <- length(sample_ids)
    sharing_thresholds <- c(
      2,  # Minimum sharing
      ceiling(n_samples * 0.5),  # Shared in half of samples
      ceiling(n_samples * 0.75),  # Shared in 75% of samples
      n_samples  # Shared in all samples
    )
    sharing_thresholds <- unique(sharing_thresholds)  # Remove duplicates
  }
  
  cat("Analyzing peptides at sharing thresholds:", paste(sharing_thresholds, collapse = ", "), "\n")
  
  # Process shared peptides at each threshold
  threshold_results <- list()
  
  for (threshold in sharing_thresholds) {
    cat("Processing peptides shared in >=", threshold, "samples...\n")
    
    # Get peptides shared at or above this threshold
    shared_peptides <- peptide_report %>%
      dplyr::filter(sample_count >= threshold)
    
    if (nrow(shared_peptides) == 0) {
      cat("No peptides shared in >=", threshold, "samples\n")
      next
    }
    
    cat("Found", nrow(shared_peptides), "peptides shared in >=", threshold, "samples\n")
    
    # Initialize integrated data with shared peptides
    integrated_data <- shared_peptides %>%
      dplyr::select(Peptide, `Peptide Length`, sample_count, primary_gene, genes)
    
    # Add immunopeptidome detection status
    integrated_data <- integrated_data %>%
      dplyr::mutate(
        detection_status = "Detected"  # All are detected since we're looking at shared peptides
      )
    
    # Integrate with transcriptome data if available
    if (!is.null(transcriptome_data) && nrow(transcriptome_data) > 0) {
      integrated_data <- integrate_transcriptome_data(integrated_data, transcriptome_data, config)
    }
    
    # Integrate with LFQ proteome data if available
    if (!is.null(lfq_data) && nrow(lfq_data) > 0) {
      integrated_data <- integrate_lfq_data(integrated_data, lfq_data, config)
    }
    
    # Integrate with TMT proteome data if available
    if (!is.null(tmt_data) && nrow(tmt_data) > 0) {
      integrated_data <- integrate_tmt_data(integrated_data, tmt_data, config)
    }
    
    # Add public neoantigen classification
    integrated_data <- add_neoantigen_classification(integrated_data)
    
    # Store results for this threshold
    threshold_results[[paste0("shared_", threshold)]] <- integrated_data
    
    # Print summary 
    cat("Integrated shared peptides (threshold >=", threshold, "samples) with:\n")
    if ("transcriptome_status" %in% colnames(integrated_data)) {
      cat("  - Transcriptome:", sum(!is.na(integrated_data$transcriptome_status)), "matches\n")
    }
    if ("lfq_status" %in% colnames(integrated_data)) {
      cat("  - LFQ Proteome:", sum(!is.na(integrated_data$lfq_status)), "matches\n")
    }
    if ("tmt_status" %in% colnames(integrated_data)) {
      cat("  - TMT Proteome:", sum(!is.na(integrated_data$tmt_status)), "matches\n")
    }
    
    # Print neoantigen stats if available
    if ("public_neoantigen_classification" %in% colnames(integrated_data)) {
      neoantigens <- integrated_data %>%
        dplyr::filter(public_neoantigen_classification != "Not a public neoantigen")
      
      if (nrow(neoantigens) > 0) {
        cat("  - Identified", nrow(neoantigens), "potential public neoantigens\n")
      }
    }
  }
  
  # Return combined results
  return(list(
    peptide_report = peptide_report,
    threshold_results = threshold_results,
    sharing_thresholds = sharing_thresholds
  ))
}

#' Integrate transcriptome data with shared peptides
#' 
#' @param integrated_data Data frame with shared peptide information
#' @param transcriptome_data Transcriptome data frame
#' @param config Configuration settings
#' @return Integrated data frame with transcriptome information
integrate_transcriptome_data <- function(integrated_data, transcriptome_data, config) {
  cat("Integrating transcriptome data...\n")
  
  # Try to find tumor-normal fold changes
  # Check if we already have a log2_fold_change column in transcriptome data
  if ("log2_fold_change" %in% colnames(transcriptome_data)) {
    # Join with transcriptome data to get fold changes
    result <- integrated_data %>%
      dplyr::left_join(
        transcriptome_data %>%
          dplyr::select(gene_symbol, log2_fold_change),
        by = c("primary_gene" = "gene_symbol"),
        relationship = "many-to-many"
      ) %>%
      # Rename to be more specific
      dplyr::rename(log2_fold_change_transcriptome = log2_fold_change)
  } else {
    # Try to calculate fold changes from tumor/normal samples
    # First, identify tumor and normal samples
    sample_ids <- config$samples$include
    tumor_samples <- grep("T$|[0-9]$", sample_ids, value = TRUE)
    normal_samples <- grep("N$", sample_ids, value = TRUE)
    
    if (length(tumor_samples) > 0 && length(normal_samples) > 0) {
      # Find expression columns for tumor and normal samples
      tumor_expr_cols <- NULL
      normal_expr_cols <- NULL
      
      # Try different patterns to find matching columns
      for (tumor_id in tumor_samples) {
        # Try direct ID match first
        tumor_pattern <- paste0("RU", tumor_id)
        matched_cols <- grep(tumor_pattern, colnames(transcriptome_data), value = TRUE)
        
        # Filter to exclude normal samples (_N suffix)
        matched_cols <- matched_cols[!grepl("_N", matched_cols)]
        
        if (length(matched_cols) > 0) {
          tumor_expr_cols <- c(tumor_expr_cols, matched_cols)
        }
      }
      
      for (normal_id in normal_samples) {
        # Try direct ID match first
        normal_pattern <- paste0("RU", normal_id)
        matched_cols <- grep(normal_pattern, colnames(transcriptome_data), value = TRUE)
        
        if (length(matched_cols) > 0) {
          normal_expr_cols <- c(normal_expr_cols, matched_cols)
        }
      }
      
      # If we found matching columns, calculate fold changes
      if (length(tumor_expr_cols) > 0 && length(normal_expr_cols) > 0) {
        cat("Calculating tumor vs normal fold changes using", 
            length(tumor_expr_cols), "tumor columns and", 
            length(normal_expr_cols), "normal columns\n")
        
        # Calculate mean expression for tumor and normal samples
        transcriptome_data$tumor_avg <- rowMeans(
          transcriptome_data[, tumor_expr_cols, drop = FALSE], 
          na.rm = TRUE
        )
        
        transcriptome_data$normal_avg <- rowMeans(
          transcriptome_data[, normal_expr_cols, drop = FALSE], 
          na.rm = TRUE
        )
        
        # Calculate log2 fold change with small offset to avoid division by zero
        transcriptome_data <- transcriptome_data %>%
          dplyr::mutate(
            tumor_avg_adj = ifelse(tumor_avg == 0, 0.1, tumor_avg),
            normal_avg_adj = ifelse(normal_avg == 0, 0.1, normal_avg),
            log2_fold_change_transcriptome = log2(tumor_avg_adj / normal_avg_adj)
          )
        
        # Join with integrated data
        result <- integrated_data %>%
          dplyr::left_join(
            transcriptome_data %>%
              dplyr::select(gene_symbol, log2_fold_change_transcriptome),
            by = c("primary_gene" = "gene_symbol")
          )
      } else {
        # If we couldn't find matching columns, just use the gene symbols for joining
        # to get expression levels without fold changes
        result <- integrated_data %>%
          dplyr::left_join(
            transcriptome_data %>%
              dplyr::select(gene_symbol, dplyr::starts_with("RU")),
            by = c("primary_gene" = "gene_symbol"),
            relationship = "many-to-many"
          )
        
        # Calculate average expression across all samples
        expr_cols <- grep("^RU", colnames(result), value = TRUE)
        if (length(expr_cols) > 0) {
          result$avg_expression <- rowMeans(
            result[, expr_cols, drop = FALSE], 
            na.rm = TRUE
          )
        }
      }
    } else {
      # If we don't have tumor/normal samples, just join to get expression values
      result <- integrated_data %>%
        dplyr::left_join(
          transcriptome_data %>%
            dplyr::select(gene_symbol, dplyr::starts_with("RU")),
          by = c("primary_gene" = "gene_symbol"),
          relationship = "many-to-many"
        )
      
      # Calculate average expression across all samples
      expr_cols <- grep("^RU", colnames(result), value = TRUE)
      if (length(expr_cols) > 0) {
        result$avg_expression <- rowMeans(
          result[, expr_cols, drop = FALSE], 
          na.rm = TRUE
        )
      }
    }
  }
  
  # Add transcriptome status if fold change is available
  if ("log2_fold_change_transcriptome" %in% colnames(result)) {
    result <- result %>%
      dplyr::mutate(
        transcriptome_status = dplyr::case_when(
          is.na(log2_fold_change_transcriptome) ~ NA_character_,
          log2_fold_change_transcriptome > 1 ~ "Up in transcriptome",
          log2_fold_change_transcriptome < -1 ~ "Down in transcriptome",
          TRUE ~ "Unchanged in transcriptome"
        )
      )
  }
  
  return(result)
}

#' Integrate LFQ proteomics data with shared peptides
#' 
#' @param integrated_data Data frame with shared peptide information
#' @param lfq_data LFQ proteomics data frame
#' @param config Configuration settings
#' @return Integrated data frame with LFQ proteome information
integrate_lfq_data <- function(integrated_data, lfq_data, config) {
  cat("Integrating LFQ proteome data...\n")
  
  # Try to find log2 fold change column
  fc_col <- NULL
  for (col in c("log2_fold_change_lfq", "Log2_Difference", "log2fc", "log2FC", "Log2FC", "Log2_Fold_Change")) {
    if (col %in% colnames(lfq_data)) {
      fc_col <- col
      break
    }
  }
  
  if (!is.null(fc_col)) {
    # Ensure consistent column naming
    if (fc_col != "log2_fold_change_lfq") {
      lfq_data <- lfq_data %>%
        dplyr::rename(log2_fold_change_lfq = !!fc_col)
    }
    
    # Join with integrated data
    result <- integrated_data %>%
      dplyr::left_join(
        lfq_data %>%
          dplyr::select(gene_symbol, log2_fold_change_lfq),
        by = c("primary_gene" = "gene_symbol")
      )
    
    # Add LFQ status based on fold change
    result <- result %>%
      dplyr::mutate(
        lfq_status = dplyr::case_when(
          is.na(log2_fold_change_lfq) ~ NA_character_,
          log2_fold_change_lfq > 1 ~ "Up in LFQ proteome",
          log2_fold_change_lfq < -1 ~ "Down in LFQ proteome",
          TRUE ~ "Unchanged in LFQ proteome"
        )
      )
  } else {
    # If no fold change column is available, just join with gene information
    result <- integrated_data %>%
      dplyr::left_join(
        lfq_data %>%
          dplyr::select(gene_symbol),
        by = c("primary_gene" = "gene_symbol"),
        relationship = "many-to-many"
      )
    
    # Add a flag to indicate presence in LFQ data
    result <- result %>%
      dplyr::mutate(
        in_lfq = !is.na(gene_symbol)
      )
    
    # Add NA for lfq_status since we don't have fold change info
    result$lfq_status <- NA_character_
  }
  
  return(result)
}

#' Integrate TMT proteomics data with shared peptides
#' 
#' @param integrated_data Data frame with shared peptide information
#' @param tmt_data TMT proteomics data frame
#' @param config Configuration settings
#' @return Integrated data frame with TMT proteome information
integrate_tmt_data <- function(integrated_data, tmt_data, config) {
  cat("Integrating TMT proteome data...\n")
  
  # Try to find log2 fold change column
  fc_col <- NULL
  for (col in c("log2_fold_change_tmt", "Log2_Difference", "log2fc", "log2FC", "Log2FC", "Log2_Fold_Change")) {
    if (col %in% colnames(tmt_data)) {
      fc_col <- col
      break
    }
  }
  
  if (!is.null(fc_col)) {
    # Ensure consistent column naming
    if (fc_col != "log2_fold_change_tmt") {
      tmt_data <- tmt_data %>%
        dplyr::rename(log2_fold_change_tmt = !!fc_col)
    }
    
    # Join with integrated data
    result <- integrated_data %>%
      dplyr::left_join(
        tmt_data %>%
          dplyr::select(gene_symbol, log2_fold_change_tmt),
        by = c("primary_gene" = "gene_symbol")
      )
    
    # Add TMT status based on fold change
    result <- result %>%
      dplyr::mutate(
        tmt_status = dplyr::case_when(
          is.na(log2_fold_change_tmt) ~ NA_character_,
          log2_fold_change_tmt > 1 ~ "Up in TMT proteome",
          log2_fold_change_tmt < -1 ~ "Down in TMT proteome",
          TRUE ~ "Unchanged in TMT proteome"
        )
      )
  } else {
    # If no fold change column is available, just join with gene information
    result <- integrated_data %>%
      dplyr::left_join(
        tmt_data %>%
          dplyr::select(gene_symbol),
        by = c("primary_gene" = "gene_symbol"),
        relationship = "many-to-many"
      )
    
    # Add a flag to indicate presence in TMT data
    result <- result %>%
      dplyr::mutate(
        in_tmt = !is.na(gene_symbol)
      )
    
    # Add NA for tmt_status since we don't have fold change info
    result$tmt_status <- NA_character_
  }
  
  return(result)
}

#' Add public neoantigen classification based on multi-omics data
#' 
#' @param integrated_data Data frame with integrated multi-omics data
#' @return Data frame with added neoantigen classification
add_neoantigen_classification <- function(integrated_data) {
  # Check which omics data sources we have
  has_transcriptome <- "transcriptome_status" %in% colnames(integrated_data)
  has_lfq <- "lfq_status" %in% colnames(integrated_data)
  has_tmt <- "tmt_status" %in% colnames(integrated_data)
  
  # Only proceed if we have at least one omics source besides immunopeptidome
  if (!has_transcriptome && !has_lfq && !has_tmt) {
    return(integrated_data)
  }
  
  # Count available data sources for each peptide
  integrated_data <- integrated_data %>%
    dplyr::mutate(
      data_sources_count = 1 +  # Immunopeptidome is always there
        (has_transcriptome & !is.na(transcriptome_status)) +
        (has_lfq & !is.na(lfq_status)) +
        (has_tmt & !is.na(tmt_status))
    )
  
  # Categorize based on upregulation status in available datasets
  integrated_data <- integrated_data %>%
    dplyr::mutate(
      # Check if upregulated in transcriptome
      up_in_transcriptome = has_transcriptome & 
        !is.na(transcriptome_status) & 
        transcriptome_status == "Up in transcriptome",
      
      # Check if upregulated in LFQ
      up_in_lfq = has_lfq & 
        !is.na(lfq_status) & 
        lfq_status == "Up in LFQ proteome",
      
      # Check if upregulated in TMT
      up_in_tmt = has_tmt & 
        !is.na(tmt_status) & 
        tmt_status == "Up in TMT proteome",
      
      # Count upregulated sources
      upregulated_sources = sum(up_in_transcriptome, up_in_lfq, up_in_tmt, na.rm = TRUE),
      
      # Define neoantigen tiers
      public_neoantigen_classification = dplyr::case_when(
        # Tier 1: Upregulated in all available omics datasets
        upregulated_sources == data_sources_count - 1 ~ "Tier 1 (Up in all datasets)",
        
        # Tier 2: Upregulated in majority of available datasets
        upregulated_sources >= (data_sources_count - 1) / 2 ~ "Tier 2 (Up in majority of datasets)",
        
        # Tier 3: Upregulated in at least one dataset
        upregulated_sources > 0 ~ "Tier 3 (Up in at least one dataset)",
        
        # Not a neoantigen
        TRUE ~ "Not a public neoantigen"
      ),
      
      # Define a numeric score for sorting
      neoantigen_score = dplyr::case_when(
        public_neoantigen_classification == "Tier 1 (Up in all datasets)" ~ 3,
        public_neoantigen_classification == "Tier 2 (Up in majority of datasets)" ~ 2,
        public_neoantigen_classification == "Tier 3 (Up in at least one dataset)" ~ 1,
        TRUE ~ 0
      )
    )
  
  return(integrated_data)
}

#' Generate comprehensive multi-omics integration report
#' 
#' @param multi_omics_results Results from run_multi_omics_integration
#' @param output_dir Directory to save report files
#' @param config Configuration settings
#' @return List of file paths for generated reports
generate_multi_omics_report <- function(multi_omics_results, output_dir, config) {
  cat("\nGenerating multi-omics integration report...\n")
  
  # Make sure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize report files list
  report_files <- list()
  
  # Generate reports for each threshold level
  for (threshold_name in names(multi_omics_results$threshold_results)) {
    threshold_data <- multi_omics_results$threshold_results[[threshold_name]]
    threshold_value <- as.numeric(gsub("shared_", "", threshold_name))
    
    # Export integrated data
    integrated_file <- file.path(output_dir, paste0("integrated_", threshold_name, ".csv"))
    write.csv(threshold_data, integrated_file, row.names = FALSE)
    report_files[[paste0("integrated_", threshold_name)]] <- integrated_file
    
    # If we have neoantigen classifications, export those separately
    if ("public_neoantigen_classification" %in% colnames(threshold_data)) {
      neoantigen_data <- threshold_data %>%
        dplyr::filter(public_neoantigen_classification != "Not a public neoantigen") %>%
        dplyr::arrange(dplyr::desc(neoantigen_score), dplyr::desc(sample_count))
      
      if (nrow(neoantigen_data) > 0) {
        neoantigen_file <- file.path(output_dir, paste0("neoantigens_", threshold_name, ".csv"))
        write.csv(neoantigen_data, neoantigen_file, row.names = FALSE)
        report_files[[paste0("neoantigens_", threshold_name)]] <- neoantigen_file
      }
    }
  }
  
  return(report_files)
}

# ADD this function at the bottom of the file

#' Generate visualizations for multi-omics integration results
#' 
#' @param multi_omics_results Results from run_multi_omics_integration
#' @param viz_dir Directory to save visualizations
#' @param config Configuration settings
#' @return List of paths to generated visualizations
generate_multi_omics_visualizations <- function(multi_omics_results, viz_dir, config) {
  cat("Generating multi-omics visualizations...\n")
  
  # Make sure output directory exists
  if (!dir.exists(viz_dir)) {
    dir.create(viz_dir, recursive = TRUE)
  }
  
  # Initialize results
  viz_paths <- list()
  
  # Generate visualizations for each threshold
  if (!is.null(multi_omics_results$threshold_results)) {
    for (threshold_name in names(multi_omics_results$threshold_results)) {
      threshold_data <- multi_omics_results$threshold_results[[threshold_name]]
      threshold_value <- as.numeric(gsub("shared_", "", threshold_name))
      
      # Create barplot of detected categories
      barplot_file <- file.path(viz_dir, paste0(threshold_name, "_categories.png"))
      
      tryCatch({
        # Count by neoantigen category
        if ("public_neoantigen_classification" %in% colnames(threshold_data)) {
          cat_data <- threshold_data %>%
            dplyr::count(public_neoantigen_classification) %>%
            dplyr::mutate(pct = n / sum(n) * 100)
          
          # Save barplot
          png(barplot_file, width = 800, height = 600)
          barplot(cat_data$n, 
                  names.arg = cat_data$public_neoantigen_classification,
                  col = c("red", "orange", "green", "gray")[1:nrow(cat_data)],
                  main = paste0("Peptides shared in ≥", threshold_value, " samples"),
                  xlab = "Category", ylab = "Count")
          dev.off()
          
          viz_paths[[paste0(threshold_name, "_barplot")]] <- barplot_file
        }
      }, error = function(e) {
        cat("Error creating barplot for threshold", threshold_name, ":", e$message, "\n")
      })
    }
  }
  
  return(viz_paths)
}