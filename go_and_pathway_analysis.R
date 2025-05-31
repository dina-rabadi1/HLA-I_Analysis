# GO and Pathway Analysis for HLA-I Immunopeptidome
# go_and_pathway_analysis.R

library(readxl)
library(ggplot2)       # For advanced plotting
library(dplyr)         # For data manipulation
library(RColorBrewer)  # For color palettes

# Try installing and loading enrichR which often works with older R versions
if (!require("enrichR")) {
  install.packages("enrichR")
  library(enrichR)
}

setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis")

# 1. Read and preprocess the Excel file
peptide_data <- read_excel("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/shared50_filtered_out_148N_manual.xlsx", sheet = "Sheet3")

# Verify data structure
print("Dimensions of loaded data:")
print(dim(peptide_data))
print("Column names:")
print(colnames(peptide_data))
print("First few rows:")
print(head(peptide_data, 3))

# Check for key columns
required_cols <- c("Peptide", "genes")
missing_cols <- required_cols[!required_cols %in% colnames(peptide_data)]
if (length(missing_cols) > 0) {
  cat("WARNING: Missing required columns:", missing_cols, "\n")
  cat("Available columns are:", colnames(peptide_data), "\n")
}

# Check for NA values in key columns
if ("Peptide" %in% colnames(peptide_data)) {
  cat("NA values in Peptide column:", sum(is.na(peptide_data$Peptide)), "out of", nrow(peptide_data), "\n")
}
if ("genes" %in% colnames(peptide_data)) {
  cat("NA values in genes column:", sum(is.na(peptide_data$genes)), "out of", nrow(peptide_data), "\n")
}

# 2. Data validation and sanity checks
# Check dimensions and structure
dim(peptide_data)
str(peptide_data)
summary(peptide_data)

# Check peptide length distribution (should be ~8-12 for HLA-I)
ggplot(peptide_data, aes(x = `Peptide Length`)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Peptide Length Distribution", 
       x = "Length", 
       y = "Count") +
  scale_x_continuous(breaks = seq(8, 15, by = 1))

# 3. Extract gene symbols
# Process gene column - assuming genes are separated by semicolons
gene_col <- which(colnames(peptide_data) == "genes")
genes_raw <- peptide_data[[gene_col]]

# Initialize an empty vector to store genes
genes <- character(0)

# Process each entry, handling NAs
for (gene_entry in genes_raw) {
  if (!is.na(gene_entry) && gene_entry != "") {
    # Split genes by semicolon
    genes_split <- strsplit(gene_entry, ";")[[1]]
    genes_split <- trimws(genes_split)  # Remove spaces
    
    # Add valid genes to the list
    valid_genes <- genes_split[genes_split != ""]
    genes <- c(genes, valid_genes)
  }
}

# Get unique genes
genes <- unique(genes)

# Remove any empty strings that might remain
genes <- genes[genes != ""]

# Sanity check: How many unique genes do we have?
cat("Number of unique genes:", length(genes), "\n")

# # 3. Extract gene symbols
# # Process gene column - assuming genes are separated by semicolons
# gene_col <- which(colnames(peptide_data) == "genes")
# genes_raw <- peptide_data[[gene_col]]
# genes_list <- strsplit(genes_raw, ";")
# genes <- unique(unlist(genes_list))
# genes <- gsub(" ", "", genes)  # Remove spaces if any
# 
# # Remove any empty strings
# genes <- genes[genes != ""]
# 
# # Sanity check: How many unique genes do we have?
# cat("Number of unique genes:", length(genes), "\n")

# 4. Gene frequency analysis
# Initialize an empty data frame if genes is empty
if (length(genes) == 0) {
  cat("Warning: No valid genes found in the dataset.\n")
  gene_counts_df <- data.frame(Gene = character(0), Frequency = integer(0))
  top_genes <- gene_counts_df
} else {
  # Count gene frequencies
  gene_counts <- table(genes)
  gene_counts_df <- data.frame(
    Gene = names(gene_counts),
    Frequency = as.numeric(gene_counts),
    stringsAsFactors = FALSE
  )
  
  # Sort by frequency
  gene_counts_df <- gene_counts_df[order(-gene_counts_df$Frequency), ]
  
  # Get top genes (limited to what's available)
  top_genes <- if(nrow(gene_counts_df) > 0) {
    head(gene_counts_df, min(20, nrow(gene_counts_df)))
  } else {
    # If no genes, create an empty data frame with the right structure
    data.frame(Gene = character(0), Frequency = integer(0))
  }
}

# Add diagnostic output
cat("Gene counts data frame dimensions:", dim(gene_counts_df), "\n")
cat("Top genes data frame dimensions:", dim(top_genes), "\n")

# Plot top genes only if we have data
if (nrow(top_genes) > 0) {
  ggplot(top_genes, aes(x = reorder(Gene, Frequency), y = Frequency)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top 20 Most Frequent Genes", x = "Gene", y = "Frequency")
} else {
  cat("No genes available for plotting.\n")
}

# 5. Use enrichR for web-based GO enrichment if available
if (exists("enrichR")) {
  # List available databases in enrichR
  dbs <- enrichR::listEnrichrDbs()
  print("Available enrichR databases:")
  print(head(dbs))
  
  # Select relevant databases
  selected_dbs <- c("GO_Biological_Process_2021", 
                    "GO_Molecular_Function_2021", 
                    "GO_Cellular_Component_2021",
                    "KEGG_2021_Human",
                    "Reactome_2022")
  
  # Check which databases are available in your enrichR version
  available_dbs <- intersect(selected_dbs, dbs$libraryName)
  
  if (length(available_dbs) > 0) {
    # Run enrichment analysis
    enrichment_results <- enrichR::enrichr(genes, available_dbs)
    
    # Print top results from each database
    for (db in available_dbs) {
      cat("\nTop enriched terms in", db, ":\n")
      result_df <- enrichment_results[[db]]
      result_df <- result_df[order(result_df$P.value),]
      print(head(result_df[, c("Term", "P.value", "Adjusted.P.value", "Genes")], 10))
      
      # Save full results to CSV
      write.csv(result_df, paste0("results/", gsub(" ", "_", db), "_enrichment.csv"), row.names = FALSE)
      
      # Plot top 10 terms
      if (nrow(result_df) > 0) {
        top_terms <- head(result_df, 10)
        top_terms$Term <- factor(top_terms$Term, levels = rev(top_terms$Term))
        
        p <- ggplot(top_terms, aes(x = Term, y = -log10(P.value))) +
          geom_bar(stat = "identity", fill = "steelblue") +
          coord_flip() +
          theme_minimal() +
          labs(title = paste("Top 10 Enriched Terms -", db),
               x = "Term",
               y = "-log10(P-value)")
        
        print(p)
        
        # Save plot
        ggsave(paste0("results/", gsub(" ", "_", db), "_top10.pdf"), p, width = 10, height = 6)
      }
    }
  } else {
    cat("No matching databases found in enrichR. Please update enrichR or use external web tools.\n")
  }
} else {
  cat("enrichR package not available. Saving gene list for external analysis.\n")
}

# 6. Save gene list for external GO analysis (fallback)
write.csv(gene_counts_df, file = "results/gene_frequency.csv", row.names = FALSE)
write.table(genes, file = "results/gene_list_for_GO.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 7. Sample analysis (which samples contribute the most peptides)
sample_col <- which(colnames(peptide_data) == "sample_list")
sample_lists <- strsplit(peptide_data[[sample_col]], ",")

# Extract all samples
all_samples <- c()
for (sample_set in sample_lists) {
  samples <- trimws(sample_set)  # Remove spaces
  all_samples <- c(all_samples, samples)
}

# Count sample frequencies
sample_counts <- table(all_samples)
sample_counts_df <- data.frame(
  Sample = names(sample_counts),
  PeptideCount = as.numeric(sample_counts)
)
sample_counts_df <- sample_counts_df[order(-sample_counts_df$PeptideCount),]

# Plot sample frequencies
ggplot(sample_counts_df, aes(x = reorder(Sample, PeptideCount), y = PeptideCount)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Peptide Frequency by Sample", 
       x = "Sample ID", 
       y = "Number of Peptides")

# Save sample statistics
write.csv(sample_counts_df, file = "results/sample_peptide_counts.csv", row.names = FALSE)

# 8. If enrichR not available, provide instructions for external GO analysis
if (!exists("enrichR")) {
  cat("===============================================\n")
  cat("GENE LIST SAVED FOR EXTERNAL GO/PATHWAY ANALYSIS\n")
  cat("===============================================\n")
  cat("Since clusterProfiler is not available, please use one of these web tools:\n")
  cat("1. DAVID: https://david.ncifcrf.gov/\n")
  cat("2. g:Profiler: https://biit.cs.ut.ee/gprofiler/gost\n")
  cat("3. Enrichr: https://maayanlab.cloud/Enrichr/\n")
  cat("\nUpload the generated file 'gene_list_for_GO.txt' to any of these tools\n")
  cat("for GO term and pathway enrichment analysis.\n")
}