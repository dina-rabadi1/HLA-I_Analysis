# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# Set file path - modify this to match your actual path
file_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/dina_attempt_fixed_Normalized_Gene_counts_FLCdb_Panel_1.xlsx"

# Try to read the Excel file
tryCatch({
  # Read the Excel file
  cat("Reading Excel file...\n")
  gene_data <- read_excel(file_path, sheet = 1)
  
  # Print the first two rows to verify data
  cat("\n=== First Two Rows of Data ===\n")
  print(head(gene_data, 2))
  
  # Find the index of the first column containing "mean" (case insensitive)
  mean_col_index <- which(str_detect(tolower(names(gene_data)), "mean"))[1]
  if(is.na(mean_col_index)) {
    stop("No 'mean' column found in the dataset")
  }
  
  # Identify normal and tumor columns (before the mean column)
  all_cols <- names(gene_data)[1:mean_col_index-1]
  
  # First column usually contains gene IDs/names, so exclude it from analysis
  analysis_cols <- all_cols[-1]
  
  # Identify normal columns (containing 'n' or 'N')
  normal_cols <- analysis_cols[str_detect(analysis_cols, "[nN]")]
  
  # Identify tumor columns (not containing 'n' or 'N')
  tumor_cols <- setdiff(analysis_cols, normal_cols)
  
  # Sanity check - print column counts
  cat("\n=== Column Counts (Sanity Check) ===\n")
  cat("Total columns before mean:", length(analysis_cols), "\n")
  cat("Normal columns:", length(normal_cols), "\n")
  cat("Tumor columns:", length(tumor_cols), "\n")
  cat("Normal columns:", paste(normal_cols, collapse = " "), "\n")
  cat("Tumor columns:", paste(tumor_cols, collapse = " "), "\n")
  
  # Calculate means for each gene
  cat("\n=== Calculating Means and Fold Change ===\n")
  
  # Create a results data frame - using explicit dplyr namespace to avoid conflicts
  results <- gene_data %>%
    # Select only the necessary columns for calculation
    dplyr::select(1, all_of(c(normal_cols, tumor_cols))) 
  
  # Calculate mean normal
  normal_data <- results %>% dplyr::select(1, all_of(normal_cols))
  results$mean_normal <- rowMeans(as.matrix(normal_data[,-1]), na.rm = TRUE)
  
  # Calculate mean tumor
  tumor_data <- results %>% dplyr::select(1, all_of(tumor_cols))
  results$mean_tumor <- rowMeans(as.matrix(tumor_data[,-1]), na.rm = TRUE)
  
  # Calculate log2 fold change
  results$log2_fold_change <- log2(results$mean_tumor / results$mean_normal)
  
  # Print summary of results
  cat("Results summary:\n")
  print(head(results %>% dplyr::select(1, mean_normal, mean_tumor, log2_fold_change), 10))
  
  # Save results to CSV
  output_file <- "transcriptome_expression_results.csv"
  write.csv(results, output_file, row.names = FALSE)
  cat("\nResults saved to:", output_file, "\n")
  
}, error = function(e) {
  cat("Error:", e$message, "\n")
  
  # If the Excel file couldn't be read, try CSV format
  cat("Trying CSV format instead...\n")
  csv_path <- sub("\\.xlsx$", ".csv", file_path)
  tryCatch({
    gene_data <- read.csv(csv_path)
    cat("CSV file read successfully.\n")
    # Continue with the same analysis as above
    # (You would duplicate the code here)
  }, error = function(e2) {
    cat("Error reading CSV format too:", e2$message, "\n")
    cat("Please check the file path and format.\n")
  })
})