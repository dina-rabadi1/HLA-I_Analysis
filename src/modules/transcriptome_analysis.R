# Transcriptome analysis module for peptide pipeline
# src/modules/transcriptome_analysis.R

#' Run transcriptome analysis and integration with peptide data
#' 
#' @param peptide_data Data frame containing peptide data
#' @param config Configuration list with settings
#' @return List containing transcriptome analysis results
run_transcriptome_analysis <- function(peptide_data, config) {
  # Validate configuration settings
  if (!config$analysis$do_transcriptome) {
    cat("Transcriptome analysis disabled in config\n")
    return(NULL)
  }
  
  if (!file.exists(config$input$transcriptome_file)) {
    stop("Transcriptome file not found: ", config$input$transcriptome_file)
  }
  
  cat("\nLoading transcriptome data from:", config$input$transcriptome_file, "\n")
  
  # Try to load the transcriptome data
  tryCatch({
    transcriptome_data <- readxl::read_excel(config$input$transcriptome_file)
    cat("Successfully loaded transcriptome data with", nrow(transcriptome_data), "genes and", 
        ncol(transcriptome_data), "columns\n")
  }, error = function(e) {
    stop("Error loading transcriptome data: ", e$message)
  })
  
  # Print column names to help with debugging
  cat("Transcriptome data columns:", paste(head(colnames(transcriptome_data), 10), collapse = ", "), 
      ifelse(length(colnames(transcriptome_data)) > 10, "...", ""), "\n")

  
    
  # Extract relevant columns for analysis
  # Need to identify gene symbol column
  symbol_col <- NULL
  for (possible_name in c("symbol", "gene_symbol", "Symbol", "SYMBOL", "gene")) {
    if (possible_name %in% colnames(transcriptome_data)) {
      symbol_col <- possible_name
      break
    }
  }
  
  if (is.null(symbol_col)) {
    stop("Could not identify gene symbol column in transcriptome data")
  }
  
  cat("Using", symbol_col, "as gene symbol column\n")
  
  # Prepare a cleaned version of transcriptome data
  transcriptome_processed <- transcriptome_data %>%
    dplyr::rename(gene_symbol = !!symbol_col) %>%
    dplyr::filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    dplyr::mutate(gene_symbol = as.character(gene_symbol))  # Ensure it's character type
  
  # 1. Paired tumor-normal analysis
  paired_results <- NULL
  if (!is.null(config$analysis$tumor_normal) && config$analysis$tumor_normal$do_tumor_normal) {
    paired_results <- run_paired_tumor_normal_analysis(transcriptome_processed, config)
  }
  
  # 2. Sample-specific and shared peptide analysis
  sample_results <- run_sample_transcriptome_analysis(peptide_data, transcriptome_processed, config)
  
  # 3. Combined multi-omics analysis (if enabled)
  omics_results <- NULL
  if (!is.null(config$analysis$do_multi_omics) && config$analysis$do_multi_omics) {
    omics_results <- run_multi_omics_analysis(peptide_data, transcriptome_processed, config)
  }
  
  # Combine all results
  results <- list(
    transcriptome_data = transcriptome_processed,
    paired_analysis = paired_results,
    sample_analysis = sample_results,
    multi_omics = omics_results
  )
  
  return(results)
}

#' Run paired tumor-normal transcriptome analysis
#' 
#' @param transcriptome_data Processed transcriptome data frame
#' @param config Configuration list with settings
#' @return List containing paired analysis results
run_paired_tumor_normal_analysis <- function(transcriptome_data, config) {
  cat("\n--- Running paired tumor-normal transcriptome analysis ---\n")
  
  # Get pairs from config
  pairs <- config$analysis$tumor_normal$pairs
  
  if (length(pairs) == 0) {
    cat("No tumor-normal pairs defined in config\n")
    return(NULL)
  }
  
  # Initialize results list for each pair
  pair_results <- list()
  
  # Process each pair
  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]
    tumor_id <- pair$tumor
    normal_id <- pair$normal
    
    cat("Analyzing tumor-normal pair:", tumor_id, "vs", normal_id, "\n")
    
    # Identify relevant columns in transcriptome data
    tumor_cols <- grep(paste0(gsub("[TN]$", "", tumor_id), ".*[_]?T"), 
                       colnames(transcriptome_data), value = TRUE)
    normal_cols <- grep(paste0(gsub("[TN]$", "", normal_id), ".*[_]?N"), 
                        colnames(transcriptome_data), value = TRUE)
    
    # Skip if columns not found
    if (length(tumor_cols) == 0) {
      cat("Warning: No tumor columns found for", tumor_id, "\n")
      next
    }
    if (length(normal_cols) == 0) {
      cat("Warning: No normal columns found for", normal_id, "\n")
      next
    }
    
    cat("Found", length(tumor_cols), "tumor columns:", paste(tumor_cols, collapse = ", "), "\n")
    cat("Found", length(normal_cols), "normal columns:", paste(normal_cols, collapse = ", "), "\n")
    
    # Calculate average if multiple columns per condition
    pair_results[[paste0(tumor_id, "_vs_", normal_id)]] <- calculate_differential_expression(
      transcriptome_data, 
      tumor_cols, 
      normal_cols,
      tumor_id, 
      normal_id
    )
  }
  
  return(pair_results)
}

#' Calculate differential expression between tumor and normal samples
#'
#' @param transcriptome_data Transcriptome data frame
#' @param tumor_cols Tumor column names
#' @param normal_cols Normal column names
#' @param tumor_id Tumor sample ID
#' @param normal_id Normal sample ID
#' @return Data frame with differential expression analysis
calculate_differential_expression <- function(transcriptome_data, tumor_cols, normal_cols, 
                                              tumor_id, normal_id) {
  # Calculate average expression for tumor and normal
  transcriptome_data$tumor_avg <- rowMeans(
    transcriptome_data[, tumor_cols, drop = FALSE], 
    na.rm = TRUE
  )
  
  transcriptome_data$normal_avg <- rowMeans(
    transcriptome_data[, normal_cols, drop = FALSE], 
    na.rm = TRUE
  )
  
  # Calculate log2 fold change with small offset to avoid division by zero
  diff_expr <- transcriptome_data %>%
    dplyr::mutate(
      tumor_avg_adj = ifelse(tumor_avg == 0, 0.1, tumor_avg),
      normal_avg_adj = ifelse(normal_avg == 0, 0.1, normal_avg),
      log2_fold_change = log2(tumor_avg_adj / normal_avg_adj),
      expression_status = dplyr::case_when(
        log2_fold_change > 1 ~ "Up in Tumor (FC > 2)",
        log2_fold_change < -1 ~ "Down in Tumor (FC < 0.5)",
        TRUE ~ "Similar (-1 < log2FC < 1)"
      )
    ) %>%
    dplyr::select(gene_symbol, tumor_avg, normal_avg, log2_fold_change, expression_status)
  
  cat("Calculated differential expression:", nrow(diff_expr), "genes\n")
  cat("Upregulated in tumor:", sum(diff_expr$expression_status == "Up in Tumor (FC > 2)"), "\n")
  cat("Downregulated in tumor:", sum(diff_expr$expression_status == "Down in Tumor (FC < 0.5)"), "\n")
  cat("Similar expression:", sum(diff_expr$expression_status == "Similar (-1 < log2FC < 1)"), "\n")
  
  return(diff_expr)
}

#' Run sample-specific transcriptome integration analysis
#' 
#' @param peptide_data Peptide data frame
#' @param transcriptome_data Processed transcriptome data frame
#' @param config Configuration list with settings
#' @return List containing sample-specific analysis results
run_sample_transcriptome_analysis <- function(peptide_data, transcriptome_data, config) {
  cat("\n--- Running sample-specific transcriptome integration ---\n")
  
  # Get included samples
  sample_ids <- config$samples$include
  
  # Initialize results list
  sample_results <- list()
  
  # Process each sample separately
  for (sample_id in sample_ids) {
    cat("Analyzing sample:", sample_id, "\n")
    
    # Filter peptide data to current sample
    sample_peptides <- peptide_data %>%
      dplyr::filter(SampleID == sample_id)
    
    if (nrow(sample_peptides) == 0) {
      cat("No peptides found for sample", sample_id, "\n")
      next
    }
    
    # Extract genes from peptides
    peptide_genes <- sample_peptides %>%
      dplyr::filter(!is.na(Gene) & Gene != "") %>%
      dplyr::select(Peptide, Gene) %>%
      dplyr::distinct()
  
    # Find matching transcriptome columns - simpler approach
    sample_cols <- grep(paste0("RU", sample_id), colnames(transcriptome_data), value = TRUE)
    
    # Filter out normal columns if we're looking for non-normal samples
    if (!grepl("N$", sample_id)) {
      sample_cols <- sample_cols[!grepl("_N[0-9]*$", sample_cols)]
    }
    
    if (length(sample_cols) == 0) {
      cat("Warning: No transcriptome columns found for sample", sample_id, "\n")
      next
    }
          
    # # Find matching transcriptome columns using more flexible pattern matching
    # # First, try the original pattern
    # sample_base_id <- gsub("[TN]$", "", sample_id)
    # sample_suffix <- substr(sample_id, nchar(sample_id), nchar(sample_id))
    # 
    # # Find columns with sample ID prefix, but exclude those with _N suffix (normal samples)
    # sample_cols <- grep(paste0("RU", sample_id, "[^N]*$"), colnames(transcriptome_data), value = TRUE)
    # 
    # # If no direct match, try a broader search for patterns like RU51_M, etc.
    # if (length(sample_cols) == 0) {
    #   broader_pattern <- paste0("RU", sample_id, "_[PM]")
    #   sample_cols <- grep(broader_pattern, colnames(transcriptome_data), value = TRUE)
    # }
    # 
    # # If still no match, try even broader
    # if (length(sample_cols) == 0) {
    #   even_broader_pattern <- paste0("RU", sample_id)
    #   candidate_cols <- grep(even_broader_pattern, colnames(transcriptome_data), value = TRUE)
    #   # Exclude columns ending with _N (normal samples)
    #   sample_cols <- candidate_cols[!grepl("_N[0-9]*$", candidate_cols)]
    # }
    # 
    # if (length(sample_cols) == 0) {
    #   cat("Warning: No transcriptome columns found for sample", sample_id, "\n")
    #   next
    # }
    
    cat("Found", length(sample_cols), "transcriptome columns:", paste(sample_cols, collapse = ", "), "\n")
    
    # Calculate average expression for sample
    transcriptome_data$sample_expr <- rowMeans(
      transcriptome_data[, sample_cols, drop = FALSE], 
      na.rm = TRUE
    )
    
    # Match peptide genes with transcriptome
    matched_data <- peptide_genes %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        gene_list = strsplit(Gene, ";")[[1]][1]  # Take first gene if multiple
      ) %>%
      dplyr::inner_join(
        transcriptome_data %>%
          dplyr::select(gene_symbol, sample_expr),
        by = c("gene_list" = "gene_symbol")
      ) %>%
      dplyr::ungroup()
    
    # Categorize expression levels
    matched_data <- matched_data %>%
      dplyr::mutate(
        expression_level = dplyr::case_when(
          sample_expr > quantile(transcriptome_data$sample_expr, 0.75, na.rm = TRUE) ~ "High (top 25%)",
          sample_expr > quantile(transcriptome_data$sample_expr, 0.5, na.rm = TRUE) ~ "Medium-high",
          sample_expr > quantile(transcriptome_data$sample_expr, 0.25, na.rm = TRUE) ~ "Medium-low",
          TRUE ~ "Low (bottom 25%)"
        )
      )

    # Store results with proper structure
    sample_results[[sample_id]] <- list(
      matched_data = matched_data,
      match_rate = nrow(matched_data) / nrow(peptide_genes),
      summary = list(
        total_peptides = nrow(sample_peptides),  # Add this field
        has_transcriptome = nrow(matched_data),  # Add this field
        expression_levels = matched_data %>%
          dplyr::count(expression_level) %>%
          dplyr::mutate(percentage = n / sum(n) * 100)
      )
    )
    
    cat("Matched", nrow(matched_data), "out of", nrow(peptide_genes), 
        sprintf("peptide genes (%.1f%%)\n", nrow(matched_data) / nrow(peptide_genes) * 100))
  }
  
  # Perform shared peptide analysis if enabled
  shared_results <- NULL
  if (!is.null(config$analysis$do_shared_peptide) && config$analysis$do_shared_peptide) {
    shared_results <- analyze_shared_peptides_transcriptome(peptide_data, transcriptome_data, config)
  }
  
  return(list(
    sample_specific = sample_results,
    shared_analysis = shared_results
  ))
}

#' Analyze shared peptides with transcriptome data
#' 
#' @param peptide_data Peptide data frame
#' @param transcriptome_data Processed transcriptome data frame
#' @param config Configuration list with settings
#' @return List containing shared peptide analysis results
analyze_shared_peptides_transcriptome <- function(peptide_data, transcriptome_data, config) {
  cat("\n--- Analyzing shared peptides with transcriptome data ---\n")
  
  # Extract relevant config parameters
  sample_ids <- config$samples$include
  exclude_samples <- config$samples$exclude_from_shared
  
  # Remove excluded samples
  sample_ids <- setdiff(sample_ids, exclude_samples)
  
  # First filter data to only selected samples
  filtered_data <- filter_by_samples(peptide_data, sample_ids)
  
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
  
  # Define sharing thresholds from config if available, otherwise use defaults
  sharing_thresholds <- config$analysis$transcriptome$sharing_thresholds
  if (is.null(sharing_thresholds)) {
    # Default thresholds based on number of samples
    n_samples <- length(sample_ids)
    
    # Define different thresholds based on number of samples
    if (n_samples <= 3) {
      sharing_thresholds <- c(2, n_samples)  # Just use min and max for few samples
    } else if (n_samples <= 6) {
      sharing_thresholds <- c(2, n_samples %/% 2, n_samples)  # 2, half, and all
    } else {
      # Create a more granular set of thresholds for larger sample sets
      sharing_thresholds <- c(2, 3, n_samples %/% 3, n_samples %/% 2, 
                              ceiling(n_samples * 0.75), n_samples)
      sharing_thresholds <- unique(sharing_thresholds)  # Remove duplicates
      sharing_thresholds <- sharing_thresholds[sharing_thresholds <= n_samples]  # Ensure all are valid
    }
  }
  
  sharing_thresholds <- sort(unique(sharing_thresholds))  # Sort and deduplicate
  
  cat("Analyzing peptides at sharing thresholds:", paste(sharing_thresholds, collapse = ", "), "\n")
  
  # Get the primary gene for each peptide (first gene in list)
  peptide_report$primary_gene <- sapply(strsplit(peptide_report$genes, ";"), 
                                        function(x) trimws(x[1]))
  
  # Create a list to store results for each threshold
  threshold_results <- list()
  
  # Process each threshold
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
    
    threshold_matched <- shared_peptides %>%
      dplyr::filter(!is.na(primary_gene) & primary_gene != "") %>%
      dplyr::inner_join(
        transcriptome_data %>%
          dplyr::select(gene_symbol, dplyr::starts_with("RU")),
        by = c("primary_gene" = "gene_symbol")
      )
    
    # Calculate average expression for each gene across all samples
    expr_cols <- grep("^RU", colnames(threshold_matched), value = TRUE)
    threshold_matched$avg_expression <- rowMeans(
      threshold_matched[, expr_cols, drop = FALSE], 
      na.rm = TRUE
    )
    
    # Categorize expression levels
    threshold_matched <- threshold_matched %>%
      dplyr::mutate(
        expression_level = dplyr::case_when(
          avg_expression > quantile(transcriptome_data$sample_expr, 0.75, na.rm = TRUE) ~ "High (top 25%)",
          avg_expression > quantile(transcriptome_data$sample_expr, 0.5, na.rm = TRUE) ~ "Medium-high",
          avg_expression > quantile(transcriptome_data$sample_expr, 0.25, na.rm = TRUE) ~ "Medium-low",
          TRUE ~ "Low (bottom 25%)"
        )
      )
    
    # Store results for this threshold
    threshold_results[[paste0("shared_", threshold)]] <- list(
      shared_peptides = shared_peptides,
      matched_data = threshold_matched,
      match_rate = nrow(threshold_matched) / nrow(shared_peptides),
      summary = list(
        total_peptides = nrow(shared_peptides),
        has_transcriptome = nrow(threshold_matched),
        expression_levels = threshold_matched %>%
          dplyr::count(expression_level) %>%
          dplyr::mutate(percentage = n / sum(n) * 100)
      )
    )
    
    cat("Matched", nrow(threshold_matched), "out of", nrow(shared_peptides), 
        sprintf("peptide genes (%.1f%%)\n", nrow(threshold_matched) / nrow(shared_peptides) * 100))
  }
  
  # Return combined results
  return(list(
    peptide_report = peptide_report,
    threshold_results = threshold_results,
    sharing_thresholds = sharing_thresholds
  ))
}

#' Run multi-omics integration analysis
#' 
#' @param peptide_data Peptide data frame
#' @param transcriptome_data Processed transcriptome data frame
#' @param config Configuration list with settings
#' @return List containing multi-omics analysis results
run_multi_omics_analysis <- function(peptide_data, transcriptome_data, config) {
  cat("\n--- Running multi-omics integration analysis ---\n")
  
  # Only proceed if additional omics data is available
  if (!file.exists(config$input$lfq_file) && !file.exists(config$input$tmt_file)) {
    cat("No additional omics data files available for multi-omics analysis\n")
    return(NULL)
  }
  
  # Initialize results list
  omics_results <- list(
    data_sources = c("Immunopeptidome", "Transcriptome")
  )
  
  # Add LFQ proteomics if available
  lfq_processed <- NULL
  if (file.exists(config$input$lfq_file)) {
    cat("Loading LFQ proteomics data from:", config$input$lfq_file, "\n")
    
    tryCatch({
      lfq_data <- readxl::read_excel(config$input$lfq_file, sheet = "Significant and 1.5x_2")
      
      # Process LFQ data
      lfq_processed <- process_lfq_data(lfq_data)
      
      omics_results$lfq_data <- lfq_processed
      omics_results$data_sources <- c(omics_results$data_sources, "LFQ Proteome")
      
      cat("Processed", nrow(lfq_processed), "entries from LFQ proteome data\n")
    }, error = function(e) {
      cat("Error processing LFQ data:", e$message, "\n")
    })
  }
  
  # Add TMT proteomics if available
  tmt_processed <- NULL
  if (file.exists(config$input$tmt_file)) {
    cat("Loading TMT proteomics data from:", config$input$tmt_file, "\n")
    
    tryCatch({
      tmt_data <- readxl::read_excel(config$input$tmt_file, sheet = "Significant and 1.5x_2")
      
      # Process TMT data
      tmt_processed <- process_tmt_data(tmt_data)
      
      omics_results$tmt_data <- tmt_processed
      omics_results$data_sources <- c(omics_results$data_sources, "TMT Proteome")
      
      cat("Processed", nrow(tmt_processed), "entries from TMT proteome data\n")
    }, error = function(e) {
      cat("Error processing TMT data:", e$message, "\n")
    })
  }
  
  # Create a comprehensive integration across all available data sources
  # Only proceed if we have at least 2 data sources (immunopeptidome + at least one other)
  if (length(omics_results$data_sources) > 1) {
    cat("Performing multi-omics integration with:", paste(omics_results$data_sources, collapse = ", "), "\n")
    
    # Run appropriate integration based on available data
    tryCatch({
      omics_results$integrated_analysis <- integrate_omics_data(
        peptide_data, 
        transcriptome_data, 
        lfq_processed,
        tmt_processed,
        config
      )
    }, error = function(e) {
      cat("Error in integrate_omics_data:", e$message, "\n")
    })
  }
  
  return(omics_results)
}

#' Process LFQ proteomics data
#' 
#' @param lfq_data Raw LFQ data frame
#' @return Processed LFQ data frame
process_lfq_data <- function(lfq_data) {
  # Clean up column names
  lfq_processed <- lfq_data %>%
    dplyr::rename_with(~ gsub(" ", "_", .), dplyr::everything())
  
  # Try to find key columns
  gene_col <- NULL
  for (col in c("Gene_Name", "gene", "symbol", "gene_symbol", "Gene", "Symbol")) {
    if (col %in% colnames(lfq_processed)) {
      gene_col <- col
      break
    }
  }
  
  if (is.null(gene_col)) {
    stop("Cannot find gene name column in LFQ data")
  }
  
  # Process and standardize
  lfq_processed <- lfq_processed %>%
    dplyr::filter(!is.na(!!dplyr::sym(gene_col))) %>%
    dplyr::rename(gene_symbol = !!gene_col)
  
  # Try to identify log2FC column
  fc_col <- NULL
  for (col in c("Log2_Difference", "log2fc", "log2FC", "Log2FC", "Log2_Fold_Change")) {
    if (col %in% colnames(lfq_processed)) {
      fc_col <- col
      break
    }
  }
  
  if (!is.null(fc_col)) {
    lfq_processed <- lfq_processed %>%
      dplyr::rename(log2_fold_change_lfq = !!fc_col)
  } else {
    # If no FC column found, log a warning
    cat("Warning: No log2 fold change column found in LFQ data\n")
  }
  
  # Add protein category if possible
  if ("log2_fold_change_lfq" %in% colnames(lfq_processed)) {
    lfq_processed <- lfq_processed %>%
      dplyr::mutate(
        protein_category_lfq = dplyr::case_when(
          log2_fold_change_lfq > 1 ~ "Up in Tumor (FC > 2)",
          log2_fold_change_lfq < -1 ~ "Down in Tumor (FC < 0.5)",
          TRUE ~ "Similar (-1 < log2FC < 1)"
        )
      )
  }
  
  return(lfq_processed)
}

#' Process TMT proteomics data
#' 
#' @param tmt_data Raw TMT data frame
#' @return Processed TMT data frame
process_tmt_data <- function(tmt_data) {
  # Clean up column names
  tmt_processed <- tmt_data %>%
    dplyr::rename_with(~ gsub(" ", "_", .), dplyr::everything())
  
  # Try to find key columns
  gene_col <- NULL
  for (col in c("Gene_Name", "gene", "symbol", "gene_symbol", "Gene", "Symbol")) {
    if (col %in% colnames(tmt_processed)) {
      gene_col <- col
      break
    }
  }
  
  if (is.null(gene_col)) {
    stop("Cannot find gene name column in TMT data")
  }
  
  # Process and standardize
  tmt_processed <- tmt_processed %>%
    dplyr::filter(!is.na(!!dplyr::sym(gene_col))) %>%
    dplyr::rename(gene_symbol = !!gene_col)
  
  # Try to identify log2FC column
  fc_col <- NULL
  for (col in c("Log2_Difference", "log2fc", "log2FC", "Log2FC", "Log2_Fold_Change")) {
    if (col %in% colnames(tmt_processed)) {
      fc_col <- col
      break
    }
  }
  
  if (!is.null(fc_col)) {
    tmt_processed <- tmt_processed %>%
      dplyr::rename(log2_fold_change_tmt = !!fc_col)
  } else {
    # If no FC column found, log a warning
    cat("Warning: No log2 fold change column found in TMT data\n")
  }
  
  # Add protein category if possible
  if ("log2_fold_change_tmt" %in% colnames(tmt_processed)) {
    tmt_processed <- tmt_processed %>%
      dplyr::mutate(
        protein_category_tmt = dplyr::case_when(
          log2_fold_change_tmt > 1 ~ "Up in Tumor (FC > 2)",
          log2_fold_change_tmt < -1 ~ "Down in Tumor (FC < 0.5)",
          TRUE ~ "Similar (-1 < log2FC < 1)"
        )
      )
  }
  
  return(tmt_processed)
}

#' Integrate multiple omics data sources
#' 
#' @param peptide_data Peptide data frame
#' @param transcriptome_data Transcriptome data frame
#' @param lfq_data LFQ proteomics data frame (or NULL if not available)
#' @param tmt_data TMT proteomics data frame (or NULL if not available)
#' @param config Configuration list with settings
#' @return Integrated multi-omics analysis results
integrate_omics_data <- function(peptide_data, transcriptome_data, lfq_data, tmt_data, config) {
  cat("Integrating multi-omics data sources...\n")
  
  # Get sample IDs from config
  sample_ids <- config$samples$include
  
  # Process all samples together first
  # Create a peptide-gene mapping
  peptide_genes <- peptide_data %>%
    dplyr::filter(!is.na(Gene) & Gene != "") %>%
    dplyr::select(Peptide, SampleID, Gene) %>%
    dplyr::distinct() %>%
    # Extract primary gene (first gene in list)
    dplyr::rowwise() %>%
    dplyr::mutate(
      primary_gene = strsplit(Gene, ";")[[1]][1]
    ) %>%
    dplyr::ungroup()
  
  # First, precalculate sample-specific expression values
  # Create a new column for each sample's expression
  for (sample_id in sample_ids) {
    sample_base_id <- gsub("[TN]$", "", sample_id)
    sample_suffix <- substr(sample_id, nchar(sample_id), nchar(sample_id))
    
    # Find columns matching this sample
    sample_cols <- grep(paste0(sample_base_id, ".*[_]?", sample_suffix), 
                        colnames(transcriptome_data), value = TRUE)
    
    if (length(sample_cols) > 0) {
      cat("Calculating expression for", sample_id, "using", length(sample_cols), "columns:", 
          paste(sample_cols, collapse=", "), "\n")
      
      # Calculate the mean expression for this sample
      expr_col_name <- paste0("expr_", sample_id)
      transcriptome_data[[expr_col_name]] <- rowMeans(
        transcriptome_data[, sample_cols, drop = FALSE], 
        na.rm = TRUE
      )
    } else {
      cat("Warning: No columns found for sample", sample_id, "\n")
    }
  }
  
  # Calculate the overall average expression across all samples
  expr_cols <- grep("^expr_", colnames(transcriptome_data), value = TRUE)
  if (length(expr_cols) > 0) {
    transcriptome_data$sample_expr <- rowMeans(
      transcriptome_data[, expr_cols, drop = FALSE], 
      na.rm = TRUE
    )
  } else {
    # Fallback in case no sample-specific expressions were calculated
    cat("Warning: No sample-specific expression columns were created, using a generic average\n")
    
    # Try to use all RU columns as a fallback
    ru_cols <- grep("^RU", colnames(transcriptome_data), value = TRUE)
    if (length(ru_cols) > 0) {
      transcriptome_data$sample_expr <- rowMeans(
        transcriptome_data[, ru_cols, drop = FALSE], 
        na.rm = TRUE
      )
    } else {
      # Last resort - create a dummy column to avoid errors
      cat("Warning: No RU columns found either, creating a dummy expression column\n")
      transcriptome_data$sample_expr <- 1
    }
  }
  
  # Now perform the join with the precalculated expression values
  cat("Joining peptide genes with transcriptome data...\n")
  integrated_data <- peptide_genes %>%
    dplyr::inner_join(
      transcriptome_data %>%
        dplyr::select(gene_symbol, sample_expr, dplyr::starts_with("expr_")),
      by = c("primary_gene" = "gene_symbol"),
      relationship = "many-to-many"  # Explicitly set relationship
    )
  
  cat("Successfully joined", nrow(integrated_data), "rows\n")
  
  # Add transcriptome fold change if available
  if ("log2_fold_change" %in% colnames(transcriptome_data)) {
    integrated_data <- integrated_data %>%
      dplyr::left_join(
        transcriptome_data %>%
          dplyr::select(gene_symbol, log2_fold_change),
        by = c("primary_gene" = "gene_symbol"),
        relationship = "many-to-many"
      ) %>%
      dplyr::rename(log2_fold_change_transcriptome = log2_fold_change)
  } else if (length(grep("_T", sample_ids, value = TRUE)) > 0 && 
             length(grep("_N", sample_ids, value = TRUE)) > 0) {
    # If we have tumor and normal samples but no explicit fold change,
    # calculate it on the fly
    tumor_samples <- grep("_T", sample_ids, value = TRUE)
    normal_samples <- grep("_N", sample_ids, value = TRUE)
    
    cat("Calculating tumor/normal fold changes using:", 
        paste(tumor_samples, collapse=", "), "vs", 
        paste(normal_samples, collapse=", "), "\n")
    
    # Get expression columns for tumor and normal
    tumor_expr_cols <- paste0("expr_", tumor_samples)
    normal_expr_cols <- paste0("expr_", normal_samples)
    
    # Calculate mean expression for tumor and normal
    tumor_expr_cols_present <- tumor_expr_cols[tumor_expr_cols %in% colnames(integrated_data)]
    normal_expr_cols_present <- normal_expr_cols[normal_expr_cols %in% colnames(integrated_data)]
    
    if (length(tumor_expr_cols_present) > 0 && length(normal_expr_cols_present) > 0) {
      integrated_data <- integrated_data %>%
        dplyr::mutate(
          tumor_avg = rowMeans(dplyr::across(dplyr::all_of(tumor_expr_cols_present)), na.rm = TRUE),
          normal_avg = rowMeans(dplyr::across(dplyr::all_of(normal_expr_cols_present)), na.rm = TRUE),
          # Add small value to avoid division by zero
          tumor_avg_adj = ifelse(tumor_avg == 0, 0.1, tumor_avg),
          normal_avg_adj = ifelse(normal_avg == 0, 0.1, normal_avg),
          log2_fold_change_transcriptome = log2(tumor_avg_adj / normal_avg_adj)
        )
    }
  }
  
  #' Process LFQ proteomics data
  #' 
  #' @param lfq_data Raw LFQ data frame
  #' @return Processed LFQ data frame
  process_lfq_data <- function(lfq_data) {
    # Clean up column names
    lfq_processed <- lfq_data %>%
      dplyr::rename_with(~ gsub(" ", "_", .), dplyr::everything())
    
    # Try to find key columns
    gene_col <- NULL
    for (col in c("Gene_Name", "gene", "symbol", "gene_symbol", "Gene", "Symbol")) {
      if (col %in% colnames(lfq_processed)) {
        gene_col <- col
        break
      }
    }
    
    if (is.null(gene_col)) {
      stop("Cannot find gene name column in LFQ data")
    }
    
    # Process and standardize
    lfq_processed <- lfq_processed %>%
      dplyr::filter(!is.na(!!dplyr::sym(gene_col))) %>%
      dplyr::rename(gene_symbol = !!gene_col)
    
    # Specifically look for the Log2_Difference column
    if ("Log2_Difference" %in% colnames(lfq_processed)) {
      cat("Found Log2_Difference column in LFQ data\n")
      lfq_processed <- lfq_processed %>%
        dplyr::rename(log2_fold_change_lfq = Log2_Difference) %>%
        # Add protein category
        dplyr::mutate(
          protein_category_lfq = dplyr::case_when(
            log2_fold_change_lfq > 1 ~ "Up in Tumor (FC > 2)",
            log2_fold_change_lfq < -1 ~ "Down in Tumor (FC < 0.5)",
            TRUE ~ "Similar (-1 < log2FC < 1)"
          )
        )
    } else {
      # If Log2_Difference not found, try other possibilities
      fc_col <- NULL
      for (col in c("log2fc", "log2FC", "Log2FC", "Log2_Fold_Change")) {
        if (col %in% colnames(lfq_processed)) {
          fc_col <- col
          break
        }
      }
      
      if (!is.null(fc_col)) {
        cat("Found alternative fold change column:", fc_col, "in LFQ data\n")
        lfq_processed <- lfq_processed %>%
          dplyr::rename(log2_fold_change_lfq = !!fc_col) %>%
          # Add protein category
          dplyr::mutate(
            protein_category_lfq = dplyr::case_when(
              log2_fold_change_lfq > 1 ~ "Up in Tumor (FC > 2)",
              log2_fold_change_lfq < -1 ~ "Down in Tumor (FC < 0.5)",
              TRUE ~ "Similar (-1 < log2FC < 1)"
            )
          )
      } else {
        # If no FC column found, log a warning
        cat("Warning: No log2 fold change column found in LFQ data\n")
        
        # Try to find columns for tumor and normal to calculate fold change
        tumor_cols <- grep("tumor|T$|_T_", colnames(lfq_processed), ignore.case = TRUE, value = TRUE)
        normal_cols <- grep("normal|N$|_N_", colnames(lfq_processed), ignore.case = TRUE, value = TRUE)
        
        if (length(tumor_cols) > 0 && length(normal_cols) > 0) {
          cat("Attempting to calculate fold change from", 
              paste(tumor_cols, collapse=", "), "vs", 
              paste(normal_cols, collapse=", "), "\n")
          
          tryCatch({
            # Convert to numeric if needed
            for (col in c(tumor_cols, normal_cols)) {
              if (!is.numeric(lfq_processed[[col]])) {
                lfq_processed[[col]] <- as.numeric(lfq_processed[[col]])
              }
            }
            
            # Calculate fold change
            lfq_processed <- lfq_processed %>%
              dplyr::mutate(
                tumor_avg = rowMeans(dplyr::across(dplyr::all_of(tumor_cols)), na.rm = TRUE),
                normal_avg = rowMeans(dplyr::across(dplyr::all_of(normal_cols)), na.rm = TRUE),
                # Add small value to avoid division by zero
                tumor_avg_adj = ifelse(tumor_avg == 0, 0.1, tumor_avg),
                normal_avg_adj = ifelse(normal_avg == 0, 0.1, normal_avg),
                log2_fold_change_lfq = log2(tumor_avg_adj / normal_avg_adj),
                protein_category_lfq = dplyr::case_when(
                  log2_fold_change_lfq > 1 ~ "Up in Tumor (FC > 2)",
                  log2_fold_change_lfq < -1 ~ "Down in Tumor (FC < 0.5)",
                  TRUE ~ "Similar (-1 < log2FC < 1)"
                )
              )
            
            cat("Successfully calculated LFQ fold change from expression values\n")
          }, error = function(e) {
            cat("Warning: Could not calculate fold change from columns:", e$message, "\n")
          })
        }
      }
    }
    
    return(lfq_processed)
  }
  
  #' Process TMT proteomics data
  #' 
  #' @param tmt_data Raw TMT data frame
  #' @return Processed TMT data frame
  process_tmt_data <- function(tmt_data) {
    # Clean up column names
    tmt_processed <- tmt_data %>%
      dplyr::rename_with(~ gsub(" ", "_", .), dplyr::everything())
    
    # Try to find key columns
    gene_col <- NULL
    for (col in c("Gene_Name", "gene", "symbol", "gene_symbol", "Gene", "Symbol")) {
      if (col %in% colnames(tmt_processed)) {
        gene_col <- col
        break
      }
    }
    
    if (is.null(gene_col)) {
      stop("Cannot find gene name column in TMT data")
    }
    
    # Process and standardize
    tmt_processed <- tmt_processed %>%
      dplyr::filter(!is.na(!!dplyr::sym(gene_col))) %>%
      dplyr::rename(gene_symbol = !!gene_col)
    
    # Specifically look for the Log2_Difference column
    if ("Log2_Difference" %in% colnames(tmt_processed)) {
      cat("Found Log2_Difference column in TMT data\n")
      tmt_processed <- tmt_processed %>%
        dplyr::rename(log2_fold_change_tmt = Log2_Difference) %>%
        # Add protein category
        dplyr::mutate(
          protein_category_tmt = dplyr::case_when(
            log2_fold_change_tmt > 1 ~ "Up in Tumor (FC > 2)",
            log2_fold_change_tmt < -1 ~ "Down in Tumor (FC < 0.5)",
            TRUE ~ "Similar (-1 < log2FC < 1)"
          )
        )
    } else {
      # If Log2_Difference not found, try other possibilities
      fc_col <- NULL
      for (col in c("log2fc", "log2FC", "Log2FC", "Log2_Fold_Change")) {
        if (col %in% colnames(tmt_processed)) {
          fc_col <- col
          break
        }
      }
      
      if (!is.null(fc_col)) {
        cat("Found alternative fold change column:", fc_col, "in TMT data\n")
        tmt_processed <- tmt_processed %>%
          dplyr::rename(log2_fold_change_tmt = !!fc_col) %>%
          # Add protein category
          dplyr::mutate(
            protein_category_tmt = dplyr::case_when(
              log2_fold_change_tmt > 1 ~ "Up in Tumor (FC > 2)",
              log2_fold_change_tmt < -1 ~ "Down in Tumor (FC < 0.5)",
              TRUE ~ "Similar (-1 < log2FC < 1)"
            )
          )
      } else {
        # If no FC column found, log a warning
        cat("Warning: No log2 fold change column found in TMT data\n")
        
        # Try to find columns for tumor and normal to calculate fold change
        tumor_cols <- grep("tumor|T$|_T_", colnames(tmt_processed), ignore.case = TRUE, value = TRUE)
        normal_cols <- grep("normal|N$|_N_", colnames(tmt_processed), ignore.case = TRUE, value = TRUE)
        
        if (length(tumor_cols) > 0 && length(normal_cols) > 0) {
          cat("Attempting to calculate fold change from", 
              paste(tumor_cols, collapse=", "), "vs", 
              paste(normal_cols, collapse=", "), "\n")
          
          tryCatch({
            # Convert to numeric if needed
            for (col in c(tumor_cols, normal_cols)) {
              if (!is.numeric(tmt_processed[[col]])) {
                tmt_processed[[col]] <- as.numeric(tmt_processed[[col]])
              }
            }
            
            # Calculate fold change
            tmt_processed <- tmt_processed %>%
              dplyr::mutate(
                tumor_avg = rowMeans(dplyr::across(dplyr::all_of(tumor_cols)), na.rm = TRUE),
                normal_avg = rowMeans(dplyr::across(dplyr::all_of(normal_cols)), na.rm = TRUE),
                # Add small value to avoid division by zero
                tumor_avg_adj = ifelse(tumor_avg == 0, 0.1, tumor_avg),
                normal_avg_adj = ifelse(normal_avg == 0, 0.1, normal_avg),
                log2_fold_change_tmt = log2(tumor_avg_adj / normal_avg_adj),
                protein_category_tmt = dplyr::case_when(
                  log2_fold_change_tmt > 1 ~ "Up in Tumor (FC > 2)",
                  log2_fold_change_tmt < -1 ~ "Down in Tumor (FC < 0.5)",
                  TRUE ~ "Similar (-1 < log2FC < 1)"
                )
              )
            
            cat("Successfully calculated TMT fold change from expression values\n")
          }, error = function(e) {
            cat("Warning: Could not calculate fold change from columns:", e$message, "\n")
          })
        }
      }
    }
    
    return(tmt_processed)
  }
  
  # Calculate correlations between data sources
  # We'll calculate correlations for peptides that have data in multiple sources
  if ("log2_fold_change_transcriptome" %in% colnames(integrated_data)) {
    if ("log2_fold_change_lfq" %in% colnames(integrated_data)) {
      integrated_data <- integrated_data %>%
        dplyr::mutate(
          immuno_trans_correlation = dplyr::case_when(
            !is.na(log2_fold_change_transcriptome) ~ log2_fold_change_transcriptome, 
            TRUE ~ NA_real_
          )
        )
    }
    
    if ("log2_fold_change_tmt" %in% colnames(integrated_data)) {
      integrated_data <- integrated_data %>%
        dplyr::mutate(
          immuno_tmt_correlation = dplyr::case_when(
            !is.na(log2_fold_change_tmt) ~ log2_fold_change_tmt,
            TRUE ~ NA_real_
          )
        )
    }
  }
  
  # Identify potential public neoantigens (similar to 148_4way script)
  # This is only applicable if we have fold changes from multiple data sources
  has_fc_transcriptome <- "log2_fold_change_transcriptome" %in% colnames(integrated_data)
  has_fc_lfq <- "log2_fold_change_lfq" %in% colnames(integrated_data)
  has_fc_tmt <- "log2_fold_change_tmt" %in% colnames(integrated_data)
  
  # Only attempt public neoantigen identification if we have fold changes from multiple sources
  if (sum(has_fc_transcriptome, has_fc_lfq, has_fc_tmt) >= 2) {
    integrated_data <- integrated_data %>%
      dplyr::mutate(
        # Upregulated in all available datasets
        upregulated_in_all = dplyr::case_when(
          # Different conditions based on available data
          has_fc_transcriptome & has_fc_lfq & has_fc_tmt ~ 
            (!is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > 1 &
               !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > 1 &
               !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > 1),
          has_fc_transcriptome & has_fc_lfq & !has_fc_tmt ~
            (!is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > 1 &
               !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > 1),
          has_fc_transcriptome & !has_fc_lfq & has_fc_tmt ~
            (!is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > 1 &
               !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > 1),
          !has_fc_transcriptome & has_fc_lfq & has_fc_tmt ~
            (!is.na(log2_fold_change_lfq) & log2_fold_change_lfq > 1 &
               !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > 1),
          TRUE ~ FALSE
        ),
        # Upregulated in at least 2 datasets
        upregulated_in_2_datasets = dplyr::case_when(
          # Calculate sum of TRUE values for available datasets
          TRUE ~ as.logical(
            sum(
              if(has_fc_transcriptome) !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > 1 else FALSE,
              if(has_fc_lfq) !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > 1 else FALSE,
              if(has_fc_tmt) !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > 1 else FALSE,
              na.rm = TRUE
            ) >= 2
          )
        ),
        # Combined score
        public_neoantigen_score = as.integer(upregulated_in_all) * 3 + 
          as.integer(upregulated_in_2_datasets) * 2,
        
        # Flag as potential public neoantigen if any criteria are met
        potential_public_neoantigen = public_neoantigen_score > 0,
        
        # Classification label
        public_neoantigen_classification = dplyr::case_when(
          upregulated_in_all ~ "Tier 1 (Upregulated in all datasets)",
          upregulated_in_2_datasets ~ "Tier 2 (Upregulated in at least 2 datasets)",
          TRUE ~ "Not a potential public neoantigen"
        )
      )
  }
  
  # Summarize results
  cat("Integrated data summary:\n")
  cat("Total peptide-gene pairs analyzed:", nrow(integrated_data), "\n")
  
  if ("potential_public_neoantigen" %in% colnames(integrated_data)) {
    potential_neoantigens <- integrated_data %>%
      dplyr::filter(potential_public_neoantigen)
    
    cat("Potential public neoantigens identified:", nrow(potential_neoantigens), "\n")
    
    if (nrow(potential_neoantigens) > 0) {
      # Summarize by tier
      tier_summary <- potential_neoantigens %>%
        dplyr::count(public_neoantigen_classification) %>%
        dplyr::arrange(desc(n))
      
      print(tier_summary)
    }
  }
  
  # Create sample-specific integrated results
  sample_integrated <- list()
  
  for (sample_id in sample_ids) {
    # Filter to current sample
    sample_data <- integrated_data %>%
      dplyr::filter(SampleID == sample_id)
    
    if (nrow(sample_data) == 0) {
      cat("No integrated data for sample", sample_id, "\n")
      next
    }
    
    # Create sample results
    sample_integrated[[sample_id]] <- list(
      integrated_data = sample_data,
      summary = list(
        total_peptides = nrow(sample_data),
        has_transcriptome = sum(!is.na(sample_data$sample_expr)),
        has_lfq = if("log2_fold_change_lfq" %in% colnames(sample_data)) 
          sum(!is.na(sample_data$log2_fold_change_lfq)) else 0,
        has_tmt = if("log2_fold_change_tmt" %in% colnames(sample_data))
          sum(!is.na(sample_data$log2_fold_change_tmt)) else 0
      )
    )
    
    # Add neoantigen counts if available
    if ("potential_public_neoantigen" %in% colnames(sample_data)) {
      sample_integrated[[sample_id]]$summary$potential_neoantigens <- 
        sum(sample_data$potential_public_neoantigen, na.rm = TRUE)
    }
    
    cat("Sample", sample_id, ":", nrow(sample_data), "peptide-gene pairs\n")
  }
  
  # Return results
  return(list(
    combined_data = integrated_data,
    sample_specific = sample_integrated,
    potential_neoantigens = if(exists("potential_neoantigens")) potential_neoantigens else NULL
  ))
}