# GTEx Analysis Functions
# File: functions/analysis/immunopeptidomics_gtex_analysis_functions.R

#' Perform comprehensive GTEx normal tissue expression analysis
#'
#' @param neoantigen_output_dir Directory containing neoantigen analysis results
#' @param gtex_output_dir Output directory for GTEx analysis results
#' @param dataset_name Name of the dataset
#' @param timestamp Analysis timestamp
#' @param gtex_data_path Path to GTEx median TPM file
#' @param gtex_tmp_thresholds Vector of TPM thresholds to test
#' @param exclude_immune_privileged Boolean to exclude brain, eye, testis
#' @param exclude_high_turnover Boolean to exclude skin, GI tissues
#' @param include_all_tissues Boolean to override tissue filtering
#' @param verbose Boolean for detailed progress reporting
#'
#' @return List containing analysis results and summary statistics
perform_gtex_analysis <- function(
    neoantigen_output_dir,
    gtex_output_dir,
    dataset_name,
    timestamp,
    gtex_data_path,
    gtex_tpm_thresholds = c(0.5, 1, 2, 5),
    exclude_immune_privileged = TRUE,
    exclude_high_turnover = TRUE,
    include_all_tissues = FALSE,
    verbose = TRUE
) {
  
  if (verbose) cat("Starting GTEx normal tissue expression analysis...\n")
  
  # Load required libraries
  library(data.table)
  library(dplyr)
  library(biomaRt)
  library(ggplot2)
  
  # Create output directories
  dir.create(gtex_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(gtex_output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(gtex_output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  
  if (verbose) cat("✓ GTEx analysis directories created\n")
  
  #-------------------------------
  # Step 1: Load GTEx median TPM data
  #-------------------------------
  if (verbose) cat("Loading GTEx data...\n")
  
  # Check if GTEx file exists
  if (!file.exists(gtex_data_path)) {
    stop("GTEx file not found at: ", gtex_data_path)
  }
  
  # Skip the first two lines (header info)
  gtex <- fread(gtex_data_path, skip = 2)
  
  # Strip Ensembl version numbers
  gtex$ensembl_id <- sub("\\..*", "", gtex$Name)
  
  if (verbose) cat("✓ GTEx data loaded:", nrow(gtex), "genes across", ncol(gtex)-3, "tissues\n")
  
  #-------------------------------
  # Step 2: Map Ensembl → Gene Symbol
  #-------------------------------
  if (verbose) cat("Mapping Ensembl IDs to gene symbols...\n")
  
  # Connect to Ensembl
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Get gene mappings
  gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = gtex$ensembl_id,
                    mart = mart)
  
  # Merge gene symbols into GTEx data
  gtex_annotated <- gtex %>%
    mutate(ensembl_id = sub("\\..*", "", Name)) %>%
    left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id")) %>%
    filter(hgnc_symbol != "" & !is.na(hgnc_symbol))  # Remove entries without gene symbols
  
  if (verbose) cat("✓ Gene mapping completed:", nrow(gtex_annotated), "genes with symbols\n")
  
  #-------------------------------
  # Step 3: Auto-detect input files
  #-------------------------------
  input_files <- detect_neoantigen_files(neoantigen_output_dir, verbose)
  
  if (length(input_files) == 0) {
    stop("No neoantigen analysis input files found in: ", neoantigen_output_dir)
  }
  
  #-------------------------------
  # Step 4: Process each input file with each threshold
  #-------------------------------
  all_results <- list()
  
  for (file_name in names(input_files)) {
    if (verbose) cat("\nProcessing", file_name, "...\n")
    
    # Load peptide data
    file_path <- input_files[[file_name]]
    peptides <- fread(file_path)
    
    if (verbose) cat("  ✓ Loaded", nrow(peptides), "peptides\n")
    
    # Check if primary_gene column exists
    if (!"primary_gene" %in% names(peptides)) {
      if (verbose) cat("  WARNING: No 'primary_gene' column found. Skipping this file.\n")
      next
    }
    
    # Process each TPM threshold
    for (threshold in gtex_tpm_thresholds) {
      if (verbose) cat("  Processing TPM threshold:", threshold, "\n")
      
      # Identify genes that are expressed in at least one tissue
      gtex_expressed_genes <- gtex_annotated %>%
        rowwise() %>%
        mutate(max_tpm = max(c_across(where(is.numeric)), na.rm = TRUE)) %>%
        filter(max_tpm > threshold) %>%
        pull(hgnc_symbol) %>%
        unique()
      
      # Flag peptides from genes expressed in healthy tissue
      peptides_processed <- peptides %>%
        mutate(
          expressed_in_gtex = primary_gene %in% gtex_expressed_genes,
          tpm_threshold = threshold,
          input_file = file_name
        )
      
      # Calculate summary statistics
      total_peptides <- nrow(peptides)
      expressed_peptides <- sum(peptides_processed$expressed_in_gtex)
      not_expressed_peptides <- sum(!peptides_processed$expressed_in_gtex)
      
      # Store results
      result_key <- paste(file_name, "tpm", threshold, sep = "_")
      all_results[[result_key]] <- list(
        peptides = peptides_processed,
        expressed_genes = gtex_expressed_genes,
        summary = list(
          total_peptides = total_peptides,
          peptides_expressed = expressed_peptides,
          peptides_not_expressed = not_expressed_peptides,
          percent_not_expressed = round(not_expressed_peptides / total_peptides * 100, 1),
          threshold = threshold,
          file_name = file_name
        )
      )
      
      if (verbose) {
        cat("    ✓ TPM", threshold, ":", not_expressed_peptides, "/", total_peptides, 
            "peptides NOT expressed in GTEx (",
            round(not_expressed_peptides / total_peptides * 100, 1), "%)\n")
      }
    }
  }
  
  #-------------------------------
  # Step 5: Save results
  #-------------------------------
  if (verbose) cat("\nSaving GTEx analysis results...\n")
  
  save_gtex_results_simple(all_results, gtex_output_dir, dataset_name, timestamp)
  
  #-------------------------------
  # Step 6: Create visualizations
  #-------------------------------
  if (verbose) cat("Creating GTEx analysis visualizations...\n")
  
  create_gtex_plots_simple(all_results, gtex_output_dir, dataset_name, timestamp)
  
  #-------------------------------
  # Step 7: Create summary report
  #-------------------------------
  summary_stats <- create_gtex_summary_simple(all_results, gtex_output_dir, dataset_name, timestamp)
  
  if (verbose) {
    cat("✓ GTEx analysis completed successfully!\n")
    cat("  - Analyzed", length(input_files), "input files\n")
    cat("  - Tested", length(gtex_tpm_thresholds), "TPM thresholds\n")
    cat("  - Generated", length(all_results), "result sets\n")
    cat("  - Results saved in:", basename(gtex_output_dir), "\n")
    
    print_gtex_key_findings_simple(summary_stats, input_files)
  }
  
  return(list(
    results = all_results,
    summary = summary_stats,
    input_files_processed = input_files
  ))
}

#' Auto-detect neoantigen analysis input files
detect_neoantigen_files <- function(neoantigen_output_dir, verbose) {
  
  if (verbose) cat("Auto-detecting neoantigen analysis input files...\n")
  
  input_files <- list()
  
  if (!dir.exists(neoantigen_output_dir)) {
    if (verbose) cat("Directory does not exist:", neoantigen_output_dir, "\n")
    return(input_files)
  }
  
  # Check if the provided directory directly contains the target files
  all_majority_file <- list.files(neoantigen_output_dir, pattern = "all_majority_peptides_integrated\\.csv$", full.names = TRUE, recursive = TRUE)
  single_gene_file <- list.files(neoantigen_output_dir, pattern = "single_gene_peptides_only\\.csv$", full.names = TRUE, recursive = TRUE)
  
  if (length(all_majority_file) > 0) {
    input_files[["all_majority"]] <- all_majority_file[1]
  }
  if (length(single_gene_file) > 0) {
    input_files[["single_gene"]] <- single_gene_file[1]
  }
  
  if (verbose) {
    if (length(input_files) > 0) {
      cat("✓ Found", length(input_files), "input files for GTEx analysis:\n")
      for (name in names(input_files)) {
        cat("  -", name, ":", input_files[[name]], "\n")
      }
    } else {
      cat("No target files found. Looking for:\n")
      cat("  - all_majority_peptides_integrated.csv\n")
      cat("  - single_gene_peptides_only.csv\n")
    }
  }
  
  return(input_files)
}

#' Save GTEx analysis results to files
save_gtex_results_simple <- function(all_results, gtex_output_dir, dataset_name, timestamp) {
  
  for (result_key in names(all_results)) {
    result <- all_results[[result_key]]
    
    # Create output file names
    output_name <- paste0(timestamp, "_", dataset_name, "_gtex_", result_key)
    
    # Save peptides with GTEx flags
    fwrite(result$peptides, 
           file.path(gtex_output_dir, "tables", paste0(output_name, "_peptides_with_gtex_flags.csv")))
    
    # Save filtered peptides (not expressed in GTEx)
    peptides_not_expressed <- result$peptides %>% filter(!expressed_in_gtex)
    
    fwrite(peptides_not_expressed,
           file.path(gtex_output_dir, "tables", paste0(output_name, "_peptides_not_in_gtex.csv")))
  }
}

#' Create simple GTEx visualizations
create_gtex_plots_simple <- function(all_results, gtex_output_dir, dataset_name, timestamp) {
  
  # Debug: Check the structure of the first result
  cat("Debugging plot data structure...\n")
  first_result <- all_results[[1]]
  cat("First result summary structure:\n")
  print(str(first_result$summary))
  
  # Compile threshold summary data
  summary_data <- data.frame()
  for (result_key in names(all_results)) {
    result <- all_results[[result_key]]
    cat("Processing result:", result_key, "\n")
    cat("Summary fields:", names(result$summary), "\n")
    
    new_row <- data.frame(
      File = result$summary$file_name,
      TPM_Threshold = result$summary$threshold,
      Total_Peptides = result$summary$total_peptides,
      Good_Candidates = result$summary$peptides_not_expressed_in_normal,
      Percent_Good_Candidates = result$summary$percent_good_candidates
    )
    cat("New row dimensions:", dim(new_row), "\n")
    summary_data <- rbind(summary_data, new_row)
  }
  
  cat("Final summary_data dimensions:", dim(summary_data), "\n")
  print(head(summary_data))
  
  # Create threshold comparison plot
  p1 <- ggplot(summary_data, aes(x = factor(TPM_Threshold), y = Good_Candidates, fill = File)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_text(aes(label = Good_Candidates), 
              position = position_dodge(width = 0.9), 
              vjust = -0.3, size = 3) +
    labs(title = "Good Neoantigen Candidates (NOT Expressed in Normal Tissue)",
         subtitle = "Higher TPM threshold = More candidates (riskier)",
         x = "TPM Threshold", 
         y = "Number of Candidate Peptides", 
         fill = "Input File") +
    theme_minimal() +
    scale_fill_brewer(type = "qual", palette = "Set2")
  
  ggsave(file.path(gtex_output_dir, "plots", paste0(timestamp, "_", dataset_name, "_gtex_threshold_comparison.png")),
         p1, width = 12, height = 6, dpi = 300)
  
  # Create percentage plot
  p2 <- ggplot(summary_data, aes(x = factor(TPM_Threshold), y = Percent_Good_Candidates, fill = File)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_text(aes(label = paste0(Percent_Good_Candidates, "%")), 
              position = position_dodge(width = 0.9), 
              vjust = -0.3, size = 3) +
    labs(title = "Percentage of Good Neoantigen Candidates",
         subtitle = "Higher TPM threshold = Higher percentage (riskier)",
         x = "TPM Threshold", 
         y = "Percentage of Candidate Peptides", 
         fill = "Input File") +
    theme_minimal() +
    scale_fill_brewer(type = "qual", palette = "Set2")
  
  ggsave(file.path(gtex_output_dir, "plots", paste0(timestamp, "_", dataset_name, "_gtex_percentage_comparison.png")),
         p2, width = 12, height = 6, dpi = 300)
}

#' Create summary report for GTEx analysis
create_gtex_summary_simple <- function(all_results, gtex_output_dir, dataset_name, timestamp) {
  
  # Compile summary statistics across all analyses
  summary_table <- data.frame()
  for (result_key in names(all_results)) {
    result <- all_results[[result_key]]
    summary_table <- rbind(summary_table, data.frame(
      Analysis = result_key,
      Input_File = result$summary$file_name,
      TPM_Threshold = result$summary$threshold,
      Total_Peptides = result$summary$total_peptides,
      Peptides_Expressed_in_Normal = result$summary$peptides_expressed_in_normal,
      Peptides_NOT_Expressed_in_Normal = result$summary$peptides_not_expressed_in_normal,
      Percent_Good_Candidates = result$summary$percent_good_candidates
    ))
  }
  
  # Save summary
  fwrite(summary_table, 
         file.path(gtex_output_dir, "tables", paste0(timestamp, "_", dataset_name, 
                                                     "_gtex_analysis_summary.csv")))
  
  return(summary_table)
}

#' Print key findings from GTEx analysis
print_gtex_key_findings_simple <- function(summary_stats, input_files) {
  
  cat("\n=== KEY FINDINGS ===\n")
  for (file_name in names(input_files)) {
    tpm1_rows <- summary_stats[summary_stats$Input_File == file_name & summary_stats$TPM_Threshold == 1, ]
    
    if (nrow(tmp1_rows) > 0) {
      cat("📋", file_name, "(TPM ≥ 1):\n")
      cat("  - Total peptides:", tmp1_rows$Total_Peptides, "\n")
      cat("  - Not expressed in GTEx:", tpm1_rows$Peptides_Not_Expressed,
          "(", tmp1_rows$Percent_Not_Expressed, "%)\n")
    }
  }
}

#' # GTEx Analysis Functions
#' # File: functions/analysis/immunopeptidomics_gtex_analysis_functions.R
#' 
#' #' Perform comprehensive GTEx normal tissue expression analysis
#' #'
#' #' @param neoantigen_output_dir Directory containing neoantigen analysis results
#' #' @param gtex_output_dir Output directory for GTEx analysis results
#' #' @param dataset_name Name of the dataset
#' #' @param timestamp Analysis timestamp
#' #' @param gtex_data_path Path to GTEx median TPM file
#' #' @param gtex_tpm_thresholds Vector of TPM thresholds to test
#' #' @param exclude_immune_privileged Boolean to exclude brain, eye, testis
#' #' @param exclude_high_turnover Boolean to exclude skin, GI tissues
#' #' @param include_all_tissues Boolean to override tissue filtering
#' #' @param verbose Boolean for detailed progress reporting
#' #'
#' #' @return List containing analysis results and summary statistics
#' perform_gtex_analysis <- function(
#'     neoantigen_output_dir,
#'     gtex_output_dir,
#'     dataset_name,
#'     timestamp,
#'     gtex_data_path,
#'     gtex_tpm_thresholds = c(0.5, 1, 2, 5),
#'     exclude_immune_privileged = TRUE,
#'     exclude_high_turnover = TRUE,
#'     include_all_tissues = FALSE,
#'     verbose = TRUE
#' ) {
#'   
#'   if (verbose) cat("Starting GTEx normal tissue expression analysis...\n")
#'   
#'   # Create output directories
#'   dir.create(gtex_output_dir, recursive = TRUE, showWarnings = FALSE)
#'   dir.create(file.path(gtex_output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
#'   dir.create(file.path(gtex_output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
#'   
#'   if (verbose) cat("✓ GTEx analysis directories created\n")
#'   
#'   #-------------------------------
#'   # Step 1: Load GTEx median TPM data
#'   #-------------------------------
#'   if (verbose) cat("Loading GTEx data...\n")
#'   
#'   # Check if GTEx file exists
#'   if (!file.exists(gtex_data_path)) {
#'     stop("GTEx file not found at: ", gtex_data_path)
#'   }
#'   
#'   # Load required libraries for GTEx analysis
#'   if (!require(data.table, quietly = TRUE)) install.packages("data.table")
#'   if (!require(biomaRt, quietly = TRUE)) {
#'     if (!require(BiocManager, quietly = TRUE)) install.packages("BiocManager")
#'     BiocManager::install("biomaRt")
#'   }
#'   library(data.table)
#'   library(biomaRt)
#'   
#'   # Skip the first two lines (header info)
#'   gtex <- fread(gtex_data_path, skip = 2)
#'   
#'   # Strip Ensembl version numbers
#'   gtex$ensembl_id <- sub("\\..*", "", gtex$Name)
#'   
#'   if (verbose) cat("✓ GTEx data loaded:", nrow(gtex), "genes across", ncol(gtex)-3, "tissues\n")
#'   
#'   #-------------------------------
#'   # Step 2: Map Ensembl → Gene Symbol
#'   #-------------------------------
#'   if (verbose) cat("Mapping Ensembl IDs to gene symbols...\n")
#'   
#'   # Connect to Ensembl
#'   mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#'   
#'   # Get gene mappings
#'   gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
#'                     filters = "ensembl_gene_id",
#'                     values = gtex$ensembl_id,
#'                     mart = mart)
#'   
#'   # Merge gene symbols into GTEx data
#'   gtex_annotated <- gtex %>%
#'     dplyr::mutate(ensembl_id = sub("\\..*", "", Name)) %>%
#'     dplyr::left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id")) %>%
#'     dplyr::filter(hgnc_symbol != "" & !is.na(hgnc_symbol))  # Remove entries without gene symbols
#'   
#'   if (verbose) cat("✓ Gene mapping completed:", nrow(gtex_annotated), "genes with symbols\n")
#'   
#'   #-------------------------------
#'   # Step 3: Define tissue filtering
#'   #-------------------------------
#'   tissue_info <- setup_tissue_filtering(gtex_annotated, exclude_immune_privileged, 
#'                                         exclude_high_turnover, include_all_tissues, verbose)
#'   
#'   #-------------------------------
#'   # Step 4: Auto-detect input files and perform analysis
#'   #-------------------------------
#'   input_files <- detect_neoantigen_files(neoantigen_output_dir, verbose)
#'   
#'   if (length(input_files) == 0) {
#'     stop("No neoantigen analysis input files found in: ", neoantigen_output_dir)
#'   }
#'   
#'   # Initialize results storage
#'   all_gtex_results <- list()
#'   
#'   #-------------------------------
#'   # Step 5: Process each input file
#'   #-------------------------------
#'   for (file_name in names(input_files)) {
#'     if (verbose) cat("\nProcessing", file_name, "...\n")
#'     
#'     file_results <- process_gtex_file(
#'       file_path = input_files[[file_name]],
#'       file_name = file_name,
#'       gtex_annotated = gtex_annotated,
#'       tissue_info = tissue_info,
#'       gtex_tpm_thresholds = gtex_tpm_thresholds,
#'       verbose = verbose
#'     )
#'     
#'     # Merge results
#'     all_gtex_results <- c(all_gtex_results, file_results)
#'   }
#'   
#'   #-------------------------------
#'   # Step 6: Save results and create visualizations
#'   #-------------------------------
#'   if (verbose) cat("\nSaving GTEx analysis results...\n")
#'   
#'   save_gtex_results(all_gtex_results, gtex_output_dir, dataset_name, timestamp)
#'   
#'   #-------------------------------
#'   # Step 7: Create visualizations
#'   #-------------------------------
#'   if (verbose) cat("Creating GTEx analysis visualizations...\n")
#'   
#'   plot_results <- create_gtex_visualizations(
#'     all_gtex_results, input_files, gtex_output_dir, dataset_name, timestamp, verbose
#'   )
#'   
#'   #-------------------------------
#'   # Step 8: Create summary report
#'   #-------------------------------
#'   summary_stats <- create_gtex_summary(all_gtex_results, gtex_output_dir, dataset_name, timestamp)
#'   
#'   if (verbose) {
#'     cat("✓ GTEx analysis completed successfully!\n")
#'     cat("  - Analyzed", length(input_files), "input files\n")
#'     cat("  - Tested", length(gtex_tmp_thresholds), "TPM thresholds\n")
#'     cat("  - Generated visualizations and data tables\n")
#'     cat("  - Results saved in:", basename(gtex_output_dir), "\n")
#'     
#'     print_gtex_key_findings(summary_stats, input_files)
#'   }
#'   
#'   return(list(
#'     results = all_gtex_results,
#'     summary = summary_stats,
#'     plots_created = plot_results,
#'     input_files_processed = input_files
#'   ))
#' }
#' 
#' #' Setup tissue filtering configuration
#' setup_tissue_filtering <- function(gtex_annotated, exclude_immune_privileged, 
#'                                    exclude_high_turnover, include_all_tissues, verbose) {
#'   
#'   if (verbose) cat("Setting up tissue filtering...\n")
#'   
#'   # Get all tissue column names (exclude Name, Description, ensembl_id, hgnc_symbol)
#'   tissue_columns <- setdiff(names(gtex_annotated), c("Name", "Description", "ensembl_id", "hgnc_symbol"))
#'   
#'   # Define tissue exclusion lists based on GTEx naming conventions
#'   immune_privileged_tissues <- c(
#'     "Brain - Amygdala", "Brain - Anterior cingulate cortex (BA24)", "Brain - Caudate (basal ganglia)",
#'     "Brain - Cerebellar Hemisphere", "Brain - Cerebellum", "Brain - Cortex", "Brain - Frontal Cortex (BA9)",
#'     "Brain - Hippocampus", "Brain - Hypothalamus", "Brain - Nucleus accumbens (basal ganglia)",
#'     "Brain - Putamen (basal ganglia)", "Brain - Spinal cord (cervical c-1)", "Brain - Substantia nigra",
#'     "Testis", "Eye"
#'   )
#'   
#'   high_turnover_tissues <- c(
#'     "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)",
#'     "Colon - Sigmoid", "Colon - Transverse", "Small Intestine - Terminal Ileum"
#'   )
#'   
#'   # Filter tissue columns based on availability
#'   available_immune_privileged <- intersect(immune_privileged_tissues, tissue_columns)
#'   available_high_turnover <- intersect(high_turnover_tissues, tissue_columns)
#'   
#'   # Determine which tissues to exclude
#'   tissues_to_exclude <- c()
#'   if (exclude_immune_privileged && !include_all_tissues) {
#'     tissues_to_exclude <- c(tissues_to_exclude, available_immune_privileged)
#'   }
#'   if (exclude_high_turnover && !include_all_tissues) {
#'     tissues_to_exclude <- c(tissues_to_exclude, available_high_turnover)
#'   }
#'   
#'   # Create filtered tissue lists
#'   all_analysis_tissues <- tissue_columns
#'   filtered_analysis_tissues <- setdiff(tissue_columns, tissues_to_exclude)
#'   
#'   if (verbose) {
#'     cat("✓ Tissue filtering configured:\n")
#'     cat("  - Total available tissues:", length(all_analysis_tissues), "\n")
#'     cat("  - Tissues after filtering:", length(filtered_analysis_tissues), "\n")
#'     cat("  - Excluded immune-privileged:", length(available_immune_privileged), "\n")
#'     cat("  - Excluded high-turnover:", length(available_high_turnover), "\n")
#'   }
#'   
#'   return(list(
#'     all_tissues = all_analysis_tissues,
#'     filtered_tissues = filtered_analysis_tissues,
#'     excluded_immune_privileged = available_immune_privileged,
#'     excluded_high_turnover = available_high_turnover
#'   ))
#' }
#' 
#' #' Auto-detect neoantigen analysis input files
#' detect_neoantigen_files <- function(neoantigen_output_dir, verbose) {
#'   
#'   if (verbose) cat("Auto-detecting neoantigen analysis input files...\n")
#'   
#'   input_files <- list()
#'   
#'   if (!dir.exists(neoantigen_output_dir)) {
#'     if (verbose) cat("Directory does not exist:", neoantigen_output_dir, "\n")
#'     return(input_files)
#'   }
#'   
#'   # Check if the provided directory directly contains the target files
#'   all_majority_file <- list.files(neoantigen_output_dir, pattern = "all_majority_peptides_integrated\\.csv$", full.names = TRUE)
#'   single_gene_file <- list.files(neoantigen_output_dir, pattern = "single_gene_peptides_only\\.csv$", full.names = TRUE)
#'   
#'   if (length(all_majority_file) > 0 || length(single_gene_file) > 0) {
#'     # Files found directly in the provided directory
#'     if (length(all_majority_file) > 0) {
#'       input_files[[paste0("all_majority_", basename(neoantigen_output_dir))]] <- all_majority_file[1]
#'     }
#'     if (length(single_gene_file) > 0) {
#'       input_files[[paste0("single_gene_", basename(neoantigen_output_dir))]] <- single_gene_file[1]
#'     }
#'   } else {
#'     # Look in subdirectories (original behavior)
#'     analysis_subdirs <- list.dirs(neoantigen_output_dir, recursive = FALSE)
#'     
#'     for (subdir in analysis_subdirs) {
#'       # Look for the target files
#'       all_majority_file <- list.files(subdir, pattern = "all_majority_peptides_integrated\\.csv$", full.names = TRUE)
#'       single_gene_file <- list.files(subdir, pattern = "single_gene_peptides_only\\.csv$", full.names = TRUE)
#'       
#'       if (length(all_majority_file) > 0) {
#'         input_files[[paste0("all_majority_", basename(subdir))]] <- all_majority_file[1]
#'       }
#'       if (length(single_gene_file) > 0) {
#'         input_files[[paste0("single_gene_", basename(subdir))]] <- single_gene_file[1]
#'       }
#'     }
#'   }
#'   
#'   if (verbose) {
#'     if (length(input_files) > 0) {
#'       cat("✓ Found", length(input_files), "input files for GTEx analysis:\n")
#'       for (name in names(input_files)) {
#'         cat("  -", name, ":", input_files[[name]], "\n")
#'       }
#'     } else {
#'       cat("No target files found. Looking for:\n")
#'       cat("  - all_majority_peptides_integrated.csv\n")
#'       cat("  - single_gene_peptides_only.csv\n")
#'       cat("In directory:", neoantigen_output_dir, "\n")
#'       cat("Available files:\n")
#'       available_files <- list.files(neoantigen_output_dir)
#'       for (file in head(available_files, 10)) {  # Show first 10 files
#'         cat("  -", file, "\n")
#'       }
#'       if (length(available_files) > 10) {
#'         cat("  ... and", length(available_files) - 10, "more files\n")
#'       }
#'     }
#'   }
#'   
#'   return(input_files)
#' }
#' 
#' #' Process a single file for GTEx analysis
#' process_gtex_file <- function(file_path, file_name, gtex_annotated, tissue_info, 
#'                               gtex_tpm_thresholds, verbose) {
#'   
#'   file_results <- list()
#'   
#'   if (!file.exists(file_path)) {
#'     if (verbose) cat("  ERROR: File not found:", file_path, "\n")
#'     return(file_results)
#'   }
#'   
#'   # Load peptide data
#'   peptides <- fread(file_path)
#'   if (verbose) cat("  ✓ Loaded", nrow(peptides), "peptides\n")
#'   
#'   # Check if primary_gene column exists
#'   if (!"primary_gene" %in% names(peptides)) {
#'     if (verbose) cat("  WARNING: No 'primary_gene' column found. Skipping this file.\n")
#'     return(file_results)
#'   }
#'   
#'   # Process for each TPM threshold
#'   for (threshold in gtex_tpm_thresholds) {
#'     if (verbose) cat("  Processing TPM threshold:", threshold, "\n")
#'     
#'     # Get expressed genes for both tissue sets
#'     expressed_genes <- get_expressed_genes(gtex_annotated, tissue_info, threshold)
#'     
#'     # Flag peptides based on GTEx expression
#'     peptides_processed <- peptides %>%
#'       dplyr::mutate(
#'         expressed_in_gtex_all = primary_gene %in% expressed_genes$all_tissues,
#'         expressed_in_gtex_filtered = primary_gene %in% expressed_genes$filtered_tissues,
#'         tpm_threshold = threshold,
#'         input_file = file_name
#'       )
#'     
#'     # Create detailed GTEx expression data for each gene
#'     gene_expression_summary <- create_gene_expression_summary(
#'       gtex_annotated, peptides, tissue_info, threshold
#'     )
#'     
#'     # Calculate summary statistics
#'     summary_stats <- calculate_gtex_summary_stats(peptides_processed, expressed_genes, peptides)
#'     
#'     # Store results
#'     result_key <- paste(file_name, "tpm", threshold, sep = "_")
#'     file_results[[result_key]] <- list(
#'       peptides = peptides_processed,
#'       gene_expression = gene_expression_summary,
#'       expressed_genes_all = expressed_genes$all_tissues,
#'       expressed_genes_filtered = expressed_genes$filtered_tissues,
#'       summary = summary_stats
#'     )
#'     
#'     if (verbose) {
#'       cat("    ✓ TPM", threshold, "- All tissues: ", 
#'           summary_stats$peptides_not_expressed_all, "/", summary_stats$total_peptides, "peptides not expressed\n")
#'       cat("    ✓ TPM", threshold, "- Filtered tissues: ", 
#'           summary_stats$peptides_not_expressed_filtered, "/", summary_stats$total_peptides, "peptides not expressed\n")
#'     }
#'   }
#'   
#'   return(file_results)
#' }
#' 
#' #' Get expressed genes for given threshold and tissue sets
#' get_expressed_genes <- function(gtex_annotated, tissue_info, threshold) {
#'   
#'   # Analysis with all tissues
#'   expressed_genes_all <- gtex_annotated %>%
#'     dplyr::rowwise() %>%
#'     dplyr::mutate(max_tpm_all = max(dplyr::c_across(dplyr::all_of(tissue_info$all_tissues)), na.rm = TRUE)) %>%
#'     dplyr::filter(max_tpm_all > threshold) %>%
#'     dplyr::pull(hgnc_symbol) %>%
#'     unique()
#'   
#'   # Analysis with filtered tissues
#'   expressed_genes_filtered <- gtex_annotated %>%
#'     dplyr::rowwise() %>%
#'     dplyr::mutate(max_tpm_filtered = max(dplyr::c_across(dplyr::all_of(tissue_info$filtered_tissues)), na.rm = TRUE)) %>%
#'     dplyr::filter(max_tpm_filtered > threshold) %>%
#'     dplyr::pull(hgnc_symbol) %>%
#'     unique()
#'   
#'   return(list(
#'     all_tissues = expressed_genes_all,
#'     filtered_tissues = expressed_genes_filtered
#'   ))
#' }
#' 
#' #' Create gene expression summary statistics
#' create_gene_expression_summary <- function(gtex_annotated, peptides, tissue_info, threshold) {
#'   
#'   gtex_annotated %>%
#'     dplyr::filter(hgnc_symbol %in% peptides$primary_gene) %>%
#'     dplyr::rowwise() %>%
#'     dplyr::mutate(
#'       mean_tmp_all = mean(dplyr::c_across(dplyr::all_of(tissue_info$all_tissues)), na.rm = TRUE),
#'       max_tpm_all = max(dplyr::c_across(dplyr::all_of(tissue_info$all_tissues)), na.rm = TRUE),
#'       median_tpm_all = median(dplyr::c_across(dplyr::all_of(tissue_info$all_tissues)), na.rm = TRUE),
#'       tissues_above_threshold_all = sum(dplyr::c_across(dplyr::all_of(tissue_info$all_tissues)) > threshold, na.rm = TRUE),
#'       mean_tpm_filtered = mean(dplyr::c_across(dplyr::all_of(tissue_info$filtered_tissues)), na.rm = TRUE),
#'       max_tpm_filtered = max(dplyr::c_across(dplyr::all_of(tissue_info$filtered_tissues)), na.rm = TRUE),
#'       median_tpm_filtered = median(dplyr::c_across(dplyr::all_of(tissue_info$filtered_tissues)), na.rm = TRUE),
#'       tissues_above_threshold_filtered = sum(dplyr::c_across(dplyr::all_of(tissue_info$filtered_tissues)) > threshold, na.rm = TRUE),
#'       tpm_threshold = threshold
#'     ) %>%
#'     dplyr::select(hgnc_symbol, mean_tpm_all, max_tpm_all, median_tpm_all, tissues_above_threshold_all,
#'                   mean_tpm_filtered, max_tpm_filtered, median_tpm_filtered, tissues_above_threshold_filtered, tpm_threshold)
#' }
#' 
#' #' Calculate summary statistics for GTEx analysis
#' calculate_gtex_summary_stats <- function(peptides_processed, expressed_genes, original_peptides) {
#'   
#'   list(
#'     total_peptides = nrow(original_peptides),
#'     peptides_expressed_all = sum(peptides_processed$expressed_in_gtex_all),
#'     peptides_not_expressed_all = sum(!peptides_processed$expressed_in_gtex_all),
#'     peptides_expressed_filtered = sum(peptides_processed$expressed_in_gtex_filtered),
#'     peptides_not_expressed_filtered = sum(!peptides_processed$expressed_in_gtex_filtered),
#'     unique_genes_total = length(unique(original_peptides$primary_gene)),
#'     unique_genes_expressed_all = length(intersect(unique(original_peptides$primary_gene), expressed_genes$all_tissues)),
#'     unique_genes_expressed_filtered = length(intersect(unique(original_peptides$primary_gene), expressed_genes$filtered_tissues))
#'   )
#' }
#' 
#' #' Save GTEx analysis results to files
#' save_gtex_results <- function(all_gtex_results, gtex_output_dir, dataset_name, timestamp) {
#'   
#'   for (result_key in names(all_gtex_results)) {
#'     result <- all_gtex_results[[result_key]]
#'     
#'     # Create file-specific output directory
#'     file_output_dir <- file.path(gtex_output_dir, "tables", result_key)
#'     dir.create(file_output_dir, recursive = TRUE, showWarnings = FALSE)
#'     
#'     # Create output file prefix
#'     output_name <- paste0(timestamp, "_", dataset_name, "_gtex_", result_key)
#'     
#'     # Save peptides with GTEx flags
#'     fwrite(result$peptides, 
#'            file.path(file_output_dir, paste0(output_name, "_peptides_with_gtex_flags.csv")))
#'     
#'     # Save gene expression summary
#'     fwrite(result$gene_expression,
#'            file.path(file_output_dir, paste0(output_name, "_gene_expression_summary.csv")))
#'     
#'     # Save filtered peptides (not expressed in GTEx)
#'     peptides_not_expressed_all <- result$peptides %>% dplyr::filter(!expressed_in_gtex_all)
#'     peptides_not_expressed_filtered <- result$peptides %>% dplyr::filter(!expressed_in_gtex_filtered)
#'     
#'     fwrite(peptides_not_expressed_all,
#'            file.path(file_output_dir, paste0(output_name, "_peptides_not_in_gtex_all_tissues.csv")))
#'     
#'     fwrite(peptides_not_expressed_filtered,
#'            file.path(file_output_dir, paste0(output_name, "_peptides_not_in_gtex_filtered_tissues.csv")))
#'     
#'     # Save tissue-by-tissue breakdown
#'     tissue_breakdown <- get_tissue_breakdown(result, all_gtex_results)
#'     if (!is.null(tissue_breakdown)) {
#'       fwrite(tissue_breakdown,
#'              file.path(file_output_dir, paste0(output_name, "_raw_gtex_expression_by_tissue.csv")))
#'     }
#'   }
#' }
#' 
#' #' Get tissue-by-tissue expression breakdown
#' get_tissue_breakdown <- function(result, all_gtex_results) {
#'   # This would need access to gtex_annotated - consider passing it as parameter
#'   # or restructuring to include this data in the result object
#'   return(NULL)  # Placeholder - implement based on data structure needs
#' }
#' 
#' #' Create GTEx analysis visualizations
#' create_gtex_visualizations <- function(all_gtex_results, input_files, gtex_output_dir, 
#'                                        dataset_name, timestamp, verbose) {
#'   
#'   plot_results <- list()
#'   
#'   # Create plots for each input file
#'   for (file_name in names(input_files)) {
#'     if (verbose) cat("  Creating plots for", file_name, "...\n")
#'     
#'     file_plots <- create_file_specific_plots(
#'       all_gtex_results, file_name, gtex_output_dir, dataset_name, timestamp
#'     )
#'     
#'     plot_results[[file_name]] <- file_plots
#'   }
#'   
#'   return(plot_results)
#' }
#' 
#' #' Create plots for a specific input file
#' create_file_specific_plots <- function(all_gtex_results, file_name, gtex_output_dir, 
#'                                        dataset_name, timestamp) {
#'   
#'   # Filter results for this file
#'   file_results <- all_gtex_results[grepl(paste0("^", file_name), names(all_gtex_results))]
#'   
#'   if (length(file_results) == 0) return(list())
#'   
#'   plot_output_dir <- file.path(gtex_output_dir, "plots", file_name)
#'   dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
#'   
#'   plots_created <- list()
#'   
#'   # 1. Bar chart: peptides filtered at different TPM thresholds
#'   p1 <- create_threshold_comparison_plot(file_results, file_name)
#'   plot_file_1 <- file.path(plot_output_dir, paste0(timestamp, "_", dataset_name, 
#'                                                    "_gtex_", file_name, "_01_threshold_comparison.png"))
#'   ggsave(plot_file_1, p1, width = 10, height = 6, dpi = 300)
#'   plots_created$threshold_comparison <- plot_file_1
#'   
#'   # 2. Before/after filtering: peptide count by gene (using TPM 1 as example)
#'   tpm1_result <- file_results[[paste(file_name, "tpm", "1", sep = "_")]]
#'   if (!is.null(tpm1_result)) {
#'     p2 <- create_gene_count_plot(tpm1_result, file_name)
#'     plot_file_2 <- file.path(plot_output_dir, paste0(timestamp, "_", dataset_name,
#'                                                      "_gtex_", file_name, "_02_gene_peptide_counts.png"))
#'     ggsave(plot_file_2, p2, width = 12, height = 8, dpi = 300)
#'     plots_created$gene_counts <- plot_file_2
#'     
#'     # 3. Pie chart: expressed vs. not expressed in GTEx
#'     p3 <- create_expression_pie_chart(tpm1_result, file_name)
#'     plot_file_3 <- file.path(plot_output_dir, paste0(timestamp, "_", dataset_name,
#'                                                      "_gtex_", file_name, "_03_expression_pie_chart.png"))
#'     ggsave(plot_file_3, p3, width = 8, height = 8, dpi = 300)
#'     plots_created$expression_pie <- plot_file_3
#'     
#'     # 4. Box plots: intensity comparison (if intensity data available)
#'     if ("total_intensity_148T" %in% names(tmp1_result$peptides)) {
#'       p4 <- create_intensity_comparison_plot(tmp1_result, file_name)
#'       if (!is.null(p4)) {
#'         plot_file_4 <- file.path(plot_output_dir, paste0(timestamp, "_", dataset_name,
#'                                                          "_gtex_", file_name, "_04_intensity_comparison.png"))
#'         ggsave(plot_file_4, p4, width = 10, height = 6, dpi = 300)
#'         plots_created$intensity_comparison <- plot_file_4
#'       }
#'     }
#'     
#'     # 5. Tissue expression heatmap (requires gtex_annotated - implement separately)
#'     # This would need additional data structure modifications
#'   }
#'   
#'   return(plots_created)
#' }
#' 
#' # Individual plot creation functions would go here...
#' # (create_threshold_comparison_plot, create_gene_count_plot, etc.)
#' # These would contain the specific ggplot code from the original script
#' 
#' #' Create summary report for GTEx analysis
#' create_gtex_summary <- function(all_gtex_results, gtex_output_dir, dataset_name, timestamp) {
#'   
#'   # Compile summary statistics across all analyses
#'   summary_table <- data.frame()
#'   for (result_key in names(all_gtex_results)) {
#'     result <- all_gtex_results[[result_key]]
#'     summary_table <- rbind(summary_table, data.frame(
#'       Analysis = result_key,
#'       Input_File = result$peptides$input_file[1],
#'       TPM_Threshold = result$peptides$tpm_threshold[1],
#'       Total_Peptides = result$summary$total_peptides,
#'       Peptides_Not_Expressed_All_Tissues = result$summary$peptides_not_expressed_all,
#'       Peptides_Not_Expressed_Filtered_Tissues = result$summary$peptides_not_expressed_filtered,
#'       Percent_Not_Expressed_All = round(result$summary$peptides_not_expressed_all / result$summary$total_peptides * 100, 1),
#'       Percent_Not_Expressed_Filtered = round(result$summary$peptides_not_expressed_filtered / result$summary$total_peptides * 100, 1),
#'       Unique_Genes_Total = result$summary$unique_genes_total,
#'       Unique_Genes_Expressed_All = result$summary$unique_genes_expressed_all,
#'       Unique_Genes_Expressed_Filtered = result$summary$unique_genes_expressed_filtered
#'     ))
#'   }
#'   
#'   # Save summary
#'   fwrite(summary_table, 
#'          file.path(gtex_output_dir, "tables", paste0(timestamp, "_", dataset_name, 
#'                                                      "_gtex_analysis_summary.csv")))
#'   
#'   return(summary_table)
#' }
#' 
#' #' Print key findings from GTEx analysis
#' print_gtex_key_findings <- function(summary_stats, input_files) {
#'   
#'   cat("\n=== KEY FINDINGS ===\n")
#'   for (file_name in names(input_files)) {
#'     tpm1_key <- paste(file_name, "tpm", "1", sep = "_")
#'     tmp1_row <- summary_stats[summary_stats$Analysis == tmp1_key, ]
#'     
#'     if (nrow(tpm1_row) > 0) {
#'       cat("📋", gsub("_", " ", file_name), "(TPM ≥ 1):\n")
#'       cat("  - Total peptides:", tmp1_row$Total_Peptides, "\n")
#'       cat("  - Not expressed (all tissues):", tpm1_row$Peptides_Not_Expressed_All_Tissues,
#'           "(", tpm1_row$Percent_Not_Expressed_All, "%)\n")
#'       cat("  - Not expressed (filtered):", tpm1_row$Peptides_Not_Expressed_Filtered_Tissues,
#'           "(", tpm1_row$Percent_Not_Expressed_Filtered, "%)\n")
#'     }
#'   }
#' }
#' 
#' # GTEx Analysis Plotting Functions
#' # Add these functions to: functions/analysis/immunopeptidomics_gtex_analysis_functions.R
#' 
#' #' Create threshold comparison bar plot
#' #' @param file_results List of results for a specific input file
#' #' @param file_name Name of the input file
#' #' @return ggplot object
#' create_threshold_comparison_plot <- function(file_results, file_name) {
#'   
#'   # Compile threshold summary data
#'   threshold_summary <- data.frame()
#'   for (result_key in names(file_results)) {
#'     result <- file_results[[result_key]]
#'     threshold <- result$peptides$tmp_threshold[1]
#'     threshold_summary <- rbind(threshold_summary, data.frame(
#'       TPM_Threshold = threshold,
#'       All_Tissues_Not_Expressed = result$summary$peptides_not_expressed_all,
#'       Filtered_Tissues_Not_Expressed = result$summary$peptides_not_expressed_filtered,
#'       Total_Peptides = result$summary$total_peptides
#'     ))
#'   }
#'   
#'   # Reshape for plotting
#'   threshold_long <- threshold_summary %>%
#'     pivot_longer(cols = c(All_Tissues_Not_Expressed, Filtered_Tissues_Not_Expressed),
#'                  names_to = "Analysis_Type", values_to = "Peptides_Not_Expressed") %>%
#'     mutate(Analysis_Type = gsub("_", " ", Analysis_Type))
#'   
#'   # Create the plot
#'   p <- ggplot(threshold_long, aes(x = factor(TPM_Threshold), y = Peptides_Not_Expressed, fill = Analysis_Type)) +
#'     geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
#'     geom_text(aes(label = Peptides_Not_Expressed), 
#'               position = position_dodge(width = 0.9), 
#'               vjust = -0.3, size = 3) +
#'     labs(title = paste("Peptides Not Expressed in GTEx by TPM Threshold"),
#'          subtitle = paste("Analysis:", gsub("_", " ", file_name)),
#'          x = "TPM Threshold", 
#'          y = "Number of Peptides", 
#'          fill = "Tissue Set") +
#'     theme_minimal() +
#'     theme(
#'       plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
#'       plot.subtitle = element_text(size = 12, hjust = 0.5),
#'       axis.title = element_text(size = 12),
#'       axis.text = element_text(size = 10),
#'       legend.title = element_text(size = 11),
#'       legend.text = element_text(size = 10)
#'     ) +
#'     scale_fill_brewer(type = "qual", palette = "Set2")
#'   
#'   return(p)
#' }
#' 
#' #' Create gene-level peptide count comparison plot
#' #' @param tmp1_result Results for TPM threshold 1
#' #' @param file_name Name of the input file
#' #' @return ggplot object
#' create_gene_count_plot <- function(tpm1_result, file_name) {
#'   
#'   # Calculate gene-level counts
#'   gene_counts <- tpm1_result$peptides %>%
#'     group_by(primary_gene) %>%
#'     summarise(
#'       total_peptides = n(),
#'       not_expressed_all = sum(!expressed_in_gtex_all),
#'       not_expressed_filtered = sum(!expressed_in_gtex_filtered),
#'       .groups = 'drop'
#'     ) %>%
#'     arrange(desc(total_peptides)) %>%
#'     head(20)  # Top 20 genes by peptide count
#'   
#'   # Reshape for plotting
#'   gene_counts_long <- gene_counts %>%
#'     pivot_longer(cols = c(total_peptides, not_expressed_all, not_expressed_filtered),
#'                  names_to = "Category", values_to = "Count") %>%
#'     mutate(Category = case_when(
#'       Category == "total_peptides" ~ "Total Peptides",
#'       Category == "not_expressed_all" ~ "Not in GTEx (All Tissues)",
#'       Category == "not_expressed_filtered" ~ "Not in GTEx (Filtered Tissues)",
#'       TRUE ~ Category
#'     ))
#'   
#'   # Create the plot
#'   p <- ggplot(gene_counts_long, aes(x = reorder(primary_gene, Count), y = Count, fill = Category)) +
#'     geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
#'     coord_flip() +
#'     labs(title = "Peptide Counts by Gene (TPM ≥ 1)",
#'          subtitle = paste("Top 20 genes -", gsub("_", " ", file_name)),
#'          x = "Gene", 
#'          y = "Number of Peptides", 
#'          fill = "Category") +
#'     theme_minimal() +
#'     theme(
#'       plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
#'       plot.subtitle = element_text(size = 12, hjust = 0.5),
#'       axis.title = element_text(size = 12),
#'       axis.text.x = element_text(size = 10),
#'       axis.text.y = element_text(size = 9),
#'       legend.title = element_text(size = 11),
#'       legend.text = element_text(size = 10),
#'       legend.position = "bottom"
#'     ) +
#'     scale_fill_brewer(type = "qual", palette = "Set1") +
#'     guides(fill = guide_legend(ncol = 2))
#'   
#'   return(p)
#' }
#' 
#' #' Create expression status pie chart
#' #' @param tpm1_result Results for TPM threshold 1
#' #' @param file_name Name of the input file
#' #' @return ggplot object
#' create_expression_pie_chart <- function(tpm1_result, file_name) {
#'   
#'   # Calculate expression counts
#'   expressed_count <- sum(tpm1_result$peptides$expressed_in_gtex_filtered)
#'   not_expressed_count <- sum(!tmp1_result$peptides$expressed_in_gtex_filtered)
#'   total_count <- nrow(tmp1_result$peptides)
#'   
#'   # Create summary data
#'   venn_summary <- data.frame(
#'     Category = c("Expressed in GTEx", "Not Expressed in GTEx"),
#'     Count = c(expressed_count, not_expressed_count),
#'     Percentage = c(
#'       round(expressed_count / total_count * 100, 1),
#'       round(not_expressed_count / total_count * 100, 1)
#'     )
#'   )
#'   
#'   # Create the pie chart
#'   p <- ggplot(venn_summary, aes(x = "", y = Count, fill = Category)) +
#'     geom_bar(stat = "identity", width = 1, alpha = 0.8) +
#'     coord_polar("y", start = 0) +
#'     labs(title = "GTEx Expression Status (TPM ≥ 1, Filtered Tissues)",
#'          subtitle = paste("Analysis:", gsub("_", " ", file_name)),
#'          fill = "Expression Status") +
#'     theme_void() +
#'     theme(
#'       plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
#'       plot.subtitle = element_text(size = 12, hjust = 0.5),
#'       legend.title = element_text(size = 11),
#'       legend.text = element_text(size = 10),
#'       legend.position = "bottom"
#'     ) +
#'     scale_fill_manual(values = c("Expressed in GTEx" = "#FF6B6B", 
#'                                  "Not Expressed in GTEx" = "#4ECDC4")) +
#'     geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
#'               position = position_stack(vjust = 0.5),
#'               size = 4, fontface = "bold")
#'   
#'   return(p)
#' }
#' 
#' #' Create intensity comparison box plot
#' #' @param tmp1_result Results for TPM threshold 1
#' #' @param file_name Name of the input file
#' #' @return ggplot object or NULL if no intensity data
#' create_intensity_comparison_plot <- function(tmp1_result, file_name) {
#'   
#'   # Check for intensity columns (multiple possible names)
#'   intensity_cols <- c("total_intensity_148T", "intensity", "total_intensity", "final_intensity")
#'   available_intensity <- intersect(intensity_cols, names(tmp1_result$peptides))
#'   
#'   if (length(available_intensity) == 0) {
#'     return(NULL)  # No intensity data available
#'   }
#'   
#'   # Use the first available intensity column
#'   intensity_col <- available_intensity[1]
#'   
#'   # Prepare data for plotting
#'   intensity_comparison <- tmp1_result$peptides %>%
#'     filter(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] > 0) %>%
#'     mutate(
#'       GTEx_Status = ifelse(expressed_in_gtex_filtered, "Expressed in GTEx", "Not Expressed in GTEx"),
#'       log_intensity = log10(.data[[intensity_col]])
#'     )
#'   
#'   if (nrow(intensity_comparison) == 0) {
#'     return(NULL)  # No valid intensity data
#'   }
#'   
#'   # Calculate summary statistics for annotation
#'   summary_stats <- intensity_comparison %>%
#'     group_by(GTEx_Status) %>%
#'     summarise(
#'       count = n(),
#'       median_val = median(log_intensity, na.rm = TRUE),
#'       mean_val = mean(log_intensity, na.rm = TRUE),
#'       .groups = 'drop'
#'     )
#'   
#'   # Create the box plot
#'   p <- ggplot(intensity_comparison, aes(x = GTEx_Status, y = log_intensity, fill = GTEx_Status)) +
#'     geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
#'     geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
#'     stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "black") +
#'     labs(title = "Immunopeptidome Intensity by GTEx Expression Status",
#'          subtitle = paste("Analysis:", gsub("_", " ", file_name), "- TPM ≥ 1, Filtered Tissues"),
#'          x = "GTEx Expression Status", 
#'          y = paste("Log10(", gsub("_", " ", intensity_col), ")"),
#'          fill = "Status") +
#'     theme_minimal() +
#'     theme(
#'       plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
#'       plot.subtitle = element_text(size = 12, hjust = 0.5),
#'       axis.title = element_text(size = 12),
#'       axis.text = element_text(size = 10),
#'       legend.position = "none",  # Remove legend since x-axis is self-explanatory
#'       strip.text = element_text(size = 11)
#'     ) +
#'     scale_fill_manual(values = c("Expressed in GTEx" = "#FF6B6B", 
#'                                  "Not Expressed in GTEx" = "#4ECDC4")) +
#'     # Add sample count annotations
#'     geom_text(data = summary_stats, 
#'               aes(x = GTEx_Status, y = max(intensity_comparison$log_intensity) * 1.05, 
#'                   label = paste("n =", count)),
#'               inherit.aes = FALSE, size = 3.5, fontface = "bold")
#'   
#'   # Add statistical test if both groups have data
#'   if (nrow(summary_stats) == 2 && all(summary_stats$count > 5)) {
#'     # Perform Wilcoxon test
#'     test_result <- wilcox.test(
#'       log_intensity ~ GTEx_Status, 
#'       data = intensity_comparison, 
#'       alternative = "two.sided"
#'     )
#'     
#'     # Add p-value annotation
#'     p_val_text <- ifelse(test_result$p.value < 0.001, "p < 0.001", 
#'                          paste("p =", round(test_result$p.value, 3)))
#'     
#'     p <- p + annotate("text", 
#'                       x = 1.5, 
#'                       y = max(intensity_comparison$log_intensity) * 1.15,
#'                       label = paste("Wilcoxon test:", p_val_text),
#'                       size = 3.5, fontface = "italic")
#'   }
#'   
#'   return(p)
#' }
#' 
#' #' Create tissue expression heatmap for top genes
#' #' @param gtex_annotated GTEx annotated data
#' #' @param tpm1_result Results for TPM threshold 1
#' #' @param file_name Name of the input file
#' #' @param max_genes Maximum number of genes to display (default: 15)
#' #' @param max_tissues Maximum number of tissues to display (default: 25)
#' #' @return ggplot object
#' create_tissue_expression_heatmap <- function(gtex_annotated, tpm1_result, file_name, 
#'                                              max_genes = 15, max_tissues = 25) {
#'   
#'   # Get top genes by peptide count
#'   top_genes <- tmp1_result$peptides %>%
#'     group_by(primary_gene) %>%
#'     summarise(count = n(), .groups = 'drop') %>%
#'     arrange(desc(count)) %>%
#'     head(max_genes) %>%
#'     pull(primary_gene)
#'   
#'   # Get tissue columns
#'   tissue_columns <- setdiff(names(gtex_annotated), c("Name", "Description", "ensembl_id", "hgnc_symbol"))
#'   
#'   # Select top tissues by variance (more informative)
#'   tissue_variance <- gtex_annotated %>%
#'     filter(hgnc_symbol %in% top_genes) %>%
#'     select(all_of(tissue_columns)) %>%
#'     summarise_all(var, na.rm = TRUE) %>%
#'     pivot_longer(everything(), names_to = "Tissue", values_to = "Variance") %>%
#'     arrange(desc(Variance)) %>%
#'     head(max_tissues) %>%
#'     pull(Tissue)
#'   
#'   # Prepare heatmap data
#'   tissue_expression <- gtex_annotated %>%
#'     filter(hgnc_symbol %in% top_genes) %>%
#'     select(hgnc_symbol, all_of(tissue_variance)) %>%
#'     pivot_longer(cols = -hgnc_symbol, names_to = "Tissue", values_to = "TPM") %>%
#'     mutate(
#'       Tissue = str_wrap(Tissue, 25),  # Wrap long tissue names
#'       log_TPM = log10(TPM + 0.1)     # Log transform with pseudocount
#'     )
#'   
#'   # Create the heatmap
#'   p <- ggplot(tissue_expression, aes(x = Tissue, y = hgnc_symbol, fill = log_TPM)) +
#'     geom_tile(color = "white", size = 0.1) +
#'     scale_fill_gradient2(
#'       low = "#313695", 
#'       mid = "#F7F7F7", 
#'       high = "#A50026",
#'       midpoint = 0,
#'       name = "Log10(TPM)",
#'       labels = function(x) ifelse(x <= 0, "≤0.1", as.character(round(10^x, 1)))
#'     ) +
#'     labs(title = "GTEx Expression Heatmap - Top Genes by Peptide Count",
#'          subtitle = paste("Analysis:", gsub("_", " ", file_name), "- Most variable tissues shown"),
#'          x = "Tissue", 
#'          y = "Gene") +
#'     theme_minimal() +
#'     theme(
#'       plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
#'       plot.subtitle = element_text(size = 12, hjust = 0.5),
#'       axis.title = element_text(size = 12),
#'       axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
#'       axis.text.y = element_text(size = 9),
#'       legend.title = element_text(size = 11),
#'       legend.text = element_text(size = 9),
#'       panel.grid = element_blank()
#'     )
#'   
#'   return(p)
#' }
#' 
#' #' Create comprehensive multi-panel summary plot
#' #' @param file_results List of results for a specific input file
#' #' @param gtex_annotated GTEx annotated data (optional, for heatmap)
#' #' @param file_name Name of the input file
#' #' @return ggplot object with multiple panels
#' create_gtex_summary_plot <- function(file_results, file_name, gtex_annotated = NULL) {
#'   
#'   # Get TPM 1 result for detailed plots
#'   tpm1_result <- file_results[[paste(file_name, "tpm", "1", sep = "_")]]
#'   
#'   if (is.null(tmp1_result)) {
#'     return(NULL)
#'   }
#'   
#'   # Create individual plots
#'   p1 <- create_threshold_comparison_plot(file_results, file_name) +
#'     theme(legend.position = "bottom") +
#'     labs(title = "A. Threshold Comparison")
#'   
#'   p2 <- create_expression_pie_chart(tmp1_result, file_name) +
#'     labs(title = "B. Expression Status")
#'   
#'   p3 <- create_intensity_comparison_plot(tmp1_result, file_name)
#'   if (!is.null(p3)) {
#'     p3 <- p3 + labs(title = "C. Intensity Comparison")
#'   }
#'   
#'   # Arrange plots
#'   if (!is.null(p3)) {
#'     combined_plot <- grid.arrange(p1, p2, p3, ncol = 2, nrow = 2,
#'                                   top = textGrob(paste("GTEx Analysis Summary:", gsub("_", " ", file_name)),
#'                                                  gp = gpar(fontsize = 16, fontface = "bold")))
#'   } else {
#'     combined_plot <- grid.arrange(p1, p2, ncol = 2,
#'                                   top = textGrob(paste("GTEx Analysis Summary:", gsub("_", " ", file_name)),
#'                                                  gp = gpar(fontsize = 16, fontface = "bold")))
#'   }
#'   
#'   return(combined_plot)
#' }
#' 
#' #' Create peptide filtering waterfall plot
#' #' @param file_results List of results for a specific input file
#' #' @param file_name Name of the input file
#' #' @return ggplot object
#' create_filtering_waterfall_plot <- function(file_results, file_name) {
#'   
#'   # Use TPM 1 result as example
#'   tpm1_result <- file_results[[paste(file_name, "tpm", "1", sep = "_")]]
#'   
#'   if (is.null(tpm1_result)) {
#'     return(NULL)
#'   }
#'   
#'   # Calculate filtering steps
#'   filtering_steps <- data.frame(
#'     Step = c("Total Peptides", "After GTEx Filter\n(All Tissues)", 
#'              "After GTEx Filter\n(Filtered Tissues)", "Final Candidates"),
#'     Count = c(
#'       tpm1_result$summary$total_peptides,
#'       tpm1_result$summary$peptides_not_expressed_all,
#'       tpm1_result$summary$peptides_not_expressed_filtered,
#'       tmp1_result$summary$peptides_not_expressed_filtered  # Same as previous for now
#'     ),
#'     Type = c("Original", "Filtered", "Filtered", "Final")
#'   )
#'   
#'   # Calculate differences for waterfall effect
#'   filtering_steps$Change <- c(0, diff(-filtering_steps$Count))
#'   filtering_steps$Start <- c(0, head(filtering_steps$Count, -1))
#'   
#'   # Create waterfall plot
#'   p <- ggplot(filtering_steps, aes(x = Step)) +
#'     geom_col(aes(y = Count, fill = Type), alpha = 0.8) +
#'     geom_text(aes(y = Count + max(Count) * 0.02, label = Count), 
#'               size = 4, fontface = "bold") +
#'     labs(title = "Peptide Filtering Waterfall (TPM ≥ 1)",
#'          subtitle = paste("Analysis:", gsub("_", " ", file_name)),
#'          x = "Filtering Step", 
#'          y = "Number of Peptides",
#'          fill = "Category") +
#'     theme_minimal() +
#'     theme(
#'       plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
#'       plot.subtitle = element_text(size = 12, hjust = 0.5),
#'       axis.title = element_text(size = 12),
#'       axis.text.x = element_text(size = 10, angle = 15, hjust = 1),
#'       axis.text.y = element_text(size = 10),
#'       legend.title = element_text(size = 11),
#'       legend.text = element_text(size = 10)
#'     ) +
#'     scale_fill_manual(values = c("Original" = "#3498DB", 
#'                                  "Filtered" = "#E74C3C", 
#'                                  "Final" = "#2ECC71")) +
#'     scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
#'   
#'   return(p)
#'   }
