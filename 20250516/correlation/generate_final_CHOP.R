# Integrated CHOP/MSKCC Peptide Data Processing Script
# This script combines CHOP and MSKCC immunopeptidomics data sets

# Load required libraries
library(dplyr)
library(tidyr)
library(readr)

# Set working directory
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation")

# =====================================================================
# PART 1: INTEGRATE CHOP FILES (51 AND 88) - MODIFIED VERSION
# =====================================================================

cat("=====================================================================\n")
cat("PART 1: INTEGRATING CHOP FILES WITH SAMPLES 51 AND 88\n")
cat("=====================================================================\n\n")

# Define file paths for CHOP integration
original_chop_file <- "CHOP_combined_peptide.tsv"
chop_51_file <- "51_CHOP_peptides.tsv"
chop_88_file <- "88_CHOP_peptides.tsv"
integrated_chop_file <- "final_CHOP_combined_peptide.tsv"

# Read input files
cat("Reading CHOP data files...\n")
chop_data <- read.delim(original_chop_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
chop_51_data <- read.delim(chop_51_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
chop_88_data <- read.delim(chop_88_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Print column names to debug the issue
cat("\nColumn names in chop_data:\n")
print(colnames(chop_data))

cat("\nColumn names in chop_51_data:\n")
print(colnames(chop_51_data))

cat("\nColumn names in chop_88_data:\n")
print(colnames(chop_88_data))

# Ensure CHOP data has "Peptide Sequence" column instead of "Peptide"
if("Peptide" %in% colnames(chop_data) && !"Peptide Sequence" %in% colnames(chop_data)) {
  colnames(chop_data)[colnames(chop_data) == "Peptide"] <- "Peptide Sequence"
}

# Process sample 51 data - inspect column names and adjust accordingly
cat("\nProcessing sample 51 data...\n")
# Find the correct column names
spectral_count_col_51 <- grep("spectral|count", colnames(chop_51_data), ignore.case=TRUE, value=TRUE)[1]
intensity_col_51 <- grep("intensity", colnames(chop_51_data), ignore.case=TRUE, value=TRUE)[1]
peptide_col_51 <- grep("peptide", colnames(chop_51_data), ignore.case=TRUE, value=TRUE)[1]
gene_col_51 <- grep("gene", colnames(chop_51_data), ignore.case=TRUE, value=TRUE)[1]
protein_col_51 <- grep("protein$", colnames(chop_51_data), ignore.case=TRUE, value=TRUE)[1]

cat("Using these columns for sample 51:\n")
cat("Peptide column:", peptide_col_51, "\n")
cat("Spectral count column:", spectral_count_col_51, "\n")
cat("Intensity column:", intensity_col_51, "\n")
cat("Gene column:", gene_col_51, "\n")
cat("Protein column:", protein_col_51, "\n")

# Use the identified column names
chop_51_processed <- chop_51_data %>%
  select(all_of(c(peptide_col_51, spectral_count_col_51, intensity_col_51, gene_col_51, protein_col_51)))

# Rename columns to standardized names with spaces, not periods
names(chop_51_processed) <- c("Peptide Sequence", "FL51 Spectral Count", "FL51 Intensity", "Gene", "Protein")

# Process sample 88 data - similar approach
cat("\nProcessing sample 88 data...\n")
# Find the correct column names
spectral_count_col_88 <- grep("spectral|count", colnames(chop_88_data), ignore.case=TRUE, value=TRUE)[1]
intensity_col_88 <- grep("intensity", colnames(chop_88_data), ignore.case=TRUE, value=TRUE)[1]
peptide_col_88 <- grep("peptide", colnames(chop_88_data), ignore.case=TRUE, value=TRUE)[1]
gene_col_88 <- grep("gene", colnames(chop_88_data), ignore.case=TRUE, value=TRUE)[1]
protein_col_88 <- grep("protein$", colnames(chop_88_data), ignore.case=TRUE, value=TRUE)[1]

cat("Using these columns for sample 88:\n")
cat("Peptide column:", peptide_col_88, "\n")
cat("Spectral count column:", spectral_count_col_88, "\n")
cat("Intensity column:", intensity_col_88, "\n")
cat("Gene column:", gene_col_88, "\n")
cat("Protein column:", protein_col_88, "\n")

# Use the identified column names
chop_88_processed <- chop_88_data %>%
  select(all_of(c(peptide_col_88, spectral_count_col_88, intensity_col_88, gene_col_88, protein_col_88)))

# Rename columns to standardized names with spaces, not periods
names(chop_88_processed) <- c("Peptide Sequence", "FL88 Spectral Count", "FL88 Intensity", "Gene", "Protein")

# Integrate the CHOP datasets
cat("Integrating CHOP datasets...\n")
# First merge with 51 data
chop_integrated <- full_join(
  chop_data, 
  chop_51_processed,
  by = "Peptide Sequence",
  suffix = c("", ".51")
)

# Then merge with 88 data
chop_integrated <- full_join(
  chop_integrated, 
  chop_88_processed,
  by = "Peptide Sequence",
  suffix = c("", ".88")
)

# Handle duplicate columns - consolidate gene and protein info
cat("Handling duplicate columns...\n")
# Function to consolidate data from multiple columns
consolidate_columns <- function(df, base_col) {
  # Check for duplicate columns
  cols <- grep(paste0("^", base_col, "(\\.\\d+)?$"), names(df), value = TRUE)
  
  if (length(cols) > 1) {
    # Create a consolidated column
    df[[base_col]] <- apply(df[, cols, drop = FALSE], 1, function(row) {
      non_na <- row[!is.na(row) & row != ""]
      if (length(non_na) > 0) non_na[1] else NA
    })
    
    # Remove duplicate columns
    for (col in cols) {
      if (col != base_col) {
        df[[col]] <- NULL
      }
    }
  }
  
  return(df)
}

# Consolidate Gene and Protein columns
chop_integrated <- consolidate_columns(chop_integrated, "Gene")
chop_integrated <- consolidate_columns(chop_integrated, "Protein")

# Check for any remaining .x or .y or .51 or .88 columns
duplicate_columns <- grep("\\.x$|\\.y$|\\.51$|\\.88$", names(chop_integrated), value = TRUE)
if (length(duplicate_columns) > 0) {
  cat("Warning: Some duplicate columns remain:", paste(duplicate_columns, collapse=", "), "\n")
}

# Fix column names - replace periods with spaces in all column names
names(chop_integrated) <- gsub("\\.", " ", names(chop_integrated))

# Ensure proper column order - move Gene and Protein to the end
# First identify all columns except Gene and Protein
non_gene_protein_cols <- setdiff(colnames(chop_integrated), c("Gene", "Protein"))
# Then reorder columns
chop_integrated <- chop_integrated[, c(non_gene_protein_cols, "Gene", "Protein")]

# Save the integrated CHOP file
cat("Saving integrated CHOP file...\n")
write.table(chop_integrated, file=integrated_chop_file, sep="\t", quote=FALSE, row.names=FALSE)

# Print CHOP integration summary
cat("\nCHOP integration completed successfully!\n")
cat("Original CHOP file:", nrow(chop_data), "rows,", ncol(chop_data), "columns\n")
cat("CHOP 51 file:", nrow(chop_51_data), "rows\n")
cat("CHOP 88 file:", nrow(chop_88_data), "rows\n")
cat("Integrated CHOP file:", nrow(chop_integrated), "rows,", ncol(chop_integrated), "columns\n")

# Verify 51 and 88 samples
fl51_cols <- grep("FL51", colnames(chop_integrated), value=TRUE)
fl88_cols <- grep("FL88", colnames(chop_integrated), value=TRUE)

cat("\nFL51 columns in integrated file:", paste(fl51_cols, collapse=", "), "\n")
cat("FL88 columns in integrated file:", paste(fl88_cols, collapse=", "), "\n")

# Count peptides by sample
peptides_51 <- sum(chop_integrated$`FL51 Spectral Count` > 0, na.rm=TRUE)
peptides_88 <- sum(chop_integrated$`FL88 Spectral Count` > 0, na.rm=TRUE)

cat("\nPeptides detected in sample 51:", peptides_51, "\n")
cat("Peptides detected in sample 88:", peptides_88, "\n")

# =====================================================================
# PART 2: MERGE MSKCC AND INTEGRATED CHOP DATA - MODIFIED SECTION
# =====================================================================

# Define file paths for merging
mskcc_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation_results_unique_peptides_unmodified_20250519_183234/unique_peptides_unmodified.tsv"
merged_output_file <- "merged_immunopeptidomics_data.tsv"

# Read MSKCC data
cat("Reading MSKCC data...\n")
mskcc_data <- read.delim(mskcc_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Print column names to debug
cat("\nColumn names in mskcc_data:\n")
print(colnames(mskcc_data))

# Create comprehensive peptide annotations table
cat("Creating comprehensive peptide annotations table...\n")

# Function to concatenate unique values with a delimiter
concatenate_unique <- function(x, delimiter = "; ") {
  if (all(is.na(x)) || length(x) == 0) return(NA_character_)
  unique_vals <- unique(x[!is.na(x) & x != ""])
  if (length(unique_vals) == 0) return(NA_character_)
  paste(sort(unique_vals), collapse = delimiter)
}

# 1. Collect all MSKCC peptide annotations with sample information
mskcc_annotations <- mskcc_data %>%
  select(Peptide, Gene, Protein, SampleID) %>%
  filter(!is.na(Gene) | !is.na(Protein)) %>%
  distinct() %>%
  mutate(Source = "MSKCC")

# 2. Collect all CHOP peptide annotations
# We need to handle the integrated CHOP data format which might have different column names
if ("Peptide Sequence" %in% colnames(chop_integrated)) {
  peptide_col_chop <- "Peptide Sequence"
} else {
  peptide_col_chop <- grep("Peptide", colnames(chop_integrated), ignore.case=TRUE, value=TRUE)[1]
}

# Get sample IDs from CHOP data
chop_sample_ids <- c()
for (col in colnames(chop_integrated)) {
  # Extract sample IDs from spectral count column names
  if (grepl("Spectral Count", col, ignore.case=TRUE)) {
    sample_match <- regexpr("FL[0-9]+|H[0-9]_[0-9]+", col)
    if (sample_match > 0) {
      sample_id <- regmatches(col, sample_match)
      chop_sample_ids <- c(chop_sample_ids, sample_id)
    }
  }
}
chop_sample_ids <- unique(chop_sample_ids)

# Create a lookup table for CHOP samples
chop_annotations <- chop_integrated %>%
  select(all_of(c(peptide_col_chop, "Gene", "Protein"))) %>%
  filter(!is.na(Gene) | !is.na(Protein)) %>%
  rename(Peptide = all_of(peptide_col_chop)) %>%
  distinct() %>%
  mutate(Source = "CHOP", SampleID = NA)

# 3. Combine all annotations
all_annotations <- bind_rows(mskcc_annotations, chop_annotations)

# 4. Create a comprehensive peptide-to-annotation mapping table
peptide_annotation_mapping <- all_annotations %>%
  group_by(Peptide) %>%
  summarize(
    All_Genes = concatenate_unique(Gene),
    All_Proteins = concatenate_unique(Protein),
    MSKCC_SampleIDs = concatenate_unique(SampleID[Source == "MSKCC"]),
    Sources = concatenate_unique(Source),
    Annotation_Count = n_distinct(paste(Gene, Protein, sep = "_"), na.rm = TRUE),
    .groups = 'drop'
  )

# 5. Save the comprehensive annotations table
annotations_output_file <- "peptide_annotations_all.tsv"
write.table(peptide_annotation_mapping, file = annotations_output_file, 
            sep="\t", quote=FALSE, row.names=FALSE)

cat("Saved comprehensive peptide annotations to:", annotations_output_file, "\n")
cat("This file contains all possible gene and protein annotations for each peptide.\n")

# Print some stats about the annotations
multi_gene_count <- sum(grepl(";", peptide_annotation_mapping$All_Genes, fixed = TRUE), na.rm = TRUE)
multi_protein_count <- sum(grepl(";", peptide_annotation_mapping$All_Proteins, fixed = TRUE), na.rm = TRUE)
cat("Peptides with multiple gene annotations:", multi_gene_count, "\n")
cat("Peptides with multiple protein annotations:", multi_protein_count, "\n")

# Prepare MSKCC data for merging
cat("Preparing MSKCC data for merging...\n")
# Ensure numeric values
spectral_count_col_mskcc <- grep("Spectral", colnames(mskcc_data), ignore.case=TRUE, value=TRUE)[1]
intensity_col_mskcc <- grep("Intensity", colnames(mskcc_data), ignore.case=TRUE, value=TRUE)[1]

cat("Using these columns for MSKCC:\n")
cat("Spectral count column:", spectral_count_col_mskcc, "\n")
cat("Intensity column:", intensity_col_mskcc, "\n")

mskcc_data[[spectral_count_col_mskcc]] <- as.numeric(mskcc_data[[spectral_count_col_mskcc]])
mskcc_data[[intensity_col_mskcc]] <- as.numeric(mskcc_data[[intensity_col_mskcc]])

# Aggregate data (in case of duplicates)
mskcc_agg <- mskcc_data %>%
  group_by(Peptide, SampleID) %>%
  summarize(
    Spectral_Count = sum(as.numeric(!!sym(spectral_count_col_mskcc)), na.rm = TRUE),
    Intensity_Value = sum(as.numeric(!!sym(intensity_col_mskcc)), na.rm = TRUE),
    Gene = first(Gene),
    Protein = first(Protein),
    .groups = 'drop'
  )

# Create new columns for the wide format names - USING PREVIOUS APPROACH
cat("Creating column names for MSKCC wide format...\n")
mskcc_data_prep <- mskcc_agg %>%
  mutate(
    SC_ColName = paste0("MSKCC_", SampleID, "_Spectral_Count"),
    Int_ColName = paste0("MSKCC_", SampleID, "_Intensity")
  )

# Create spectral count wide format
cat("Creating spectral count wide format...\n")
mskcc_sc_wide <- mskcc_data_prep %>%
  select(Peptide, Gene, Protein, SC_ColName, Spectral_Count) %>%
  pivot_wider(
    id_cols = c(Peptide, Gene, Protein),
    names_from = SC_ColName,
    values_from = Spectral_Count,
    values_fill = list(Spectral_Count = 0)
  )

# Create intensity wide format
cat("Creating intensity wide format...\n")
mskcc_int_wide <- mskcc_data_prep %>%
  select(Peptide, Int_ColName, Intensity_Value) %>%
  pivot_wider(
    id_cols = Peptide,
    names_from = Int_ColName,
    values_from = Intensity_Value,
    values_fill = list(Intensity_Value = 0)
  )

# Combine all MSKCC data
mskcc_wide <- left_join(mskcc_sc_wide, mskcc_int_wide, by = "Peptide")

# Now that mskcc_wide exists, we can check for duplicates
cat("Checking for duplicate peptides in MSKCC_wide...\n")
if (any(duplicated(mskcc_wide$Peptide))) {
  cat("WARNING: Found duplicate peptides in MSKCC_wide, deduplicating...\n")
  
  # Deduplicate the mskcc_wide data
  mskcc_wide <- mskcc_wide %>%
    group_by(Peptide) %>%
    summarize(
      # For Gene and Protein, keep one value (first non-NA)
      Gene = first(na.omit(Gene)),
      Protein = first(na.omit(Protein)),
      # For all other columns, keep maximum values
      across(where(is.numeric), ~max(., na.rm = TRUE)),
      .groups = 'drop'
    )
  
  cat("After deduplication, mskcc_wide has", nrow(mskcc_wide), "rows\n")
}

# Prepare integrated CHOP data for merging
cat("Preparing integrated CHOP data for merging...\n")
# Identify peptide sequence column in CHOP data
if ("Peptide" %in% colnames(chop_integrated)) {
  chop_prep <- chop_integrated
} else {
  peptide_col_chop <- grep("Peptide", colnames(chop_integrated), ignore.case=TRUE, value=TRUE)[1]
  chop_prep <- chop_integrated
  colnames(chop_prep)[colnames(chop_prep) == peptide_col_chop] <- "Peptide"
}

# Identify spectral count and intensity columns
chop_sc_cols <- grep("Spectral Count", colnames(chop_prep), value = TRUE)
chop_int_cols <- grep("Intensity$", colnames(chop_prep), value = TRUE)

# Select only relevant columns
chop_ready <- chop_prep %>%
  select(Peptide, Gene, Protein, all_of(chop_sc_cols), all_of(chop_int_cols))

# Check for duplicates in chop_ready
cat("Checking for duplicate peptides in chop_ready...\n")
if (any(duplicated(chop_ready$Peptide))) {
  cat("WARNING: Found duplicate peptides in chop_ready, deduplicating...\n")
  
  # Deduplicate the chop_ready data
  chop_ready <- chop_ready %>%
    group_by(Peptide) %>%
    summarize(
      # For Gene and Protein, keep one value (first non-NA)
      Gene = first(na.omit(Gene)),
      Protein = first(na.omit(Protein)),
      # For all other columns, keep maximum values (conservative approach)
      across(where(is.numeric), ~max(., na.rm = TRUE)),
      .groups = 'drop'
    )
  
  cat("After deduplication, chop_ready has", nrow(chop_ready), "rows\n")
}

# Merge the datasets
cat("Merging MSKCC and CHOP datasets...\n")
merged_data <- full_join(mskcc_wide, chop_ready, by = "Peptide", suffix = c("_MSKCC", "_CHOP"))

# Replace NA values with zeros for numeric columns
cat("Handling missing values...\n")
numeric_cols <- c(
  grep("Spectral_Count", names(merged_data)),
  grep("Intensity", names(merged_data))
)

for (col in numeric_cols) {
  merged_data[[col]] <- ifelse(is.na(merged_data[[col]]), 0, merged_data[[col]])
}

# After merging, add a final safety check
if (any(duplicated(merged_data$Peptide))) {
  cat("WARNING: Found", sum(duplicated(merged_data$Peptide)), "duplicate peptides in merged data\n")
  cat("Performing final deduplication...\n")
  
  # Deduplicate the merged data
  merged_data <- merged_data %>%
    group_by(Peptide) %>%
    summarize(
      # For Gene and Protein columns, keep first non-NA value
      Gene_MSKCC = first(na.omit(Gene_MSKCC)),
      Protein_MSKCC = first(na.omit(Protein_MSKCC)),
      Gene_CHOP = first(na.omit(Gene_CHOP)),
      Protein_CHOP = first(na.omit(Protein_CHOP)),
      # For all other columns, take the maximum value
      across(where(is.numeric), ~max(., na.rm = TRUE)),
      .groups = 'drop'
    )
  
  cat("After final deduplication, merged data has", nrow(merged_data), "rows\n")
}

# Save the merged file
cat("Saving merged file...\n")
write.table(merged_data, file = merged_output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Print merging summary
cat("\nMerging completed successfully!\n")
cat("MSKCC dataset:", nrow(mskcc_data), "rows,", length(unique(mskcc_data$Peptide)), "unique peptides\n")
cat("CHOP dataset:", nrow(chop_integrated), "rows,", length(unique(chop_prep$Peptide)), "unique peptides\n")
cat("Merged dataset:", nrow(merged_data), "rows (unique peptides)\n")

# Calculate dataset overlap
mskcc_peptides <- unique(mskcc_data$Peptide)
chop_peptides <- unique(chop_prep$Peptide)
common_peptides <- intersect(mskcc_peptides, chop_peptides)

cat("\nPeptide overlap:\n")
cat("Peptides only in MSKCC:", length(setdiff(mskcc_peptides, chop_peptides)), "\n")
cat("Peptides only in CHOP:", length(setdiff(chop_peptides, mskcc_peptides)), "\n")
cat("Peptides in both datasets:", length(common_peptides), "\n")
cat("Overlap percentage:",
    round(100 * length(common_peptides) / length(union(mskcc_peptides, chop_peptides)), 1), "%\n")

# =====================================================================
# FINAL SUMMARY
# =====================================================================

cat("\n=====================================================================\n")
cat("INTEGRATION COMPLETE\n")
cat("=====================================================================\n\n")

cat("Files created:\n")
cat("1. final_CHOP_combined_peptide.tsv - Integrated CHOP dataset with samples 51 and 88\n")
cat("2. merged_immunopeptidomics_data.tsv - Merged MSKCC and CHOP data for visualization\n")
cat("3. peptide_annotations_all.tsv - Comprehensive table with all gene/protein annotations\n\n")

cat("Next steps:\n")
cat("- Run your visualization script to create correlation plots and overlap statistics\n")
cat("- These files contain all the data needed for your analysis\n")
cat("- For neoantigen analysis, use the peptide_annotations_all.tsv file to see all possible annotations\n")

cat("\nAll data preparation tasks completed successfully!\n")

# 
# 
# setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation")
# 
# # =====================================================================
# # PART 1: INTEGRATE CHOP FILES (51 AND 88) - MODIFIED VERSION
# # =====================================================================
# 
# cat("=====================================================================\n")
# cat("PART 1: INTEGRATING CHOP FILES WITH SAMPLES 51 AND 88\n")
# cat("=====================================================================\n\n")
# 
# # Define file paths for CHOP integration
# original_chop_file <- "CHOP_combined_peptide.tsv"
# chop_51_file <- "51_CHOP_peptides.tsv"
# chop_88_file <- "88_CHOP_peptides.tsv"
# integrated_chop_file <- "final_CHOP_combined_peptide.tsv"
# 
# # Read input files
# cat("Reading CHOP data files...\n")
# chop_data <- read.delim(original_chop_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# chop_51_data <- read.delim(chop_51_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# chop_88_data <- read.delim(chop_88_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# 
# # Print column names to debug the issue
# cat("\nColumn names in chop_data:\n")
# print(colnames(chop_data))
# cat("\nColumn names in chop_51_data:\n")
# 
# print(colnames(chop_51_data))
# 
# cat("\nColumn names in chop_88_data:\n")
# print(colnames(chop_88_data))
# 
# # Ensure CHOP data has "Peptide Sequence" column instead of "Peptide"
# if("Peptide" %in% colnames(chop_data) && !"Peptide Sequence" %in% colnames(chop_data)) {
#   colnames(chop_data)[colnames(chop_data) == "Peptide"] <- "Peptide Sequence"
# }
# 
# # Process sample 51 data - inspect column names and adjust accordingly
# cat("\nProcessing sample 51 data...\n")
# # Find the correct column names
# spectral_count_col_51 <- grep("spectral|count", colnames(chop_51_data), ignore.case=TRUE, value=TRUE)[1]
# intensity_col_51 <- grep("intensity", colnames(chop_51_data), ignore.case=TRUE, value=TRUE)[1]
# peptide_col_51 <- grep("peptide", colnames(chop_51_data), ignore.case=TRUE, value=TRUE)[1]
# gene_col_51 <- grep("gene", colnames(chop_51_data), ignore.case=TRUE, value=TRUE)[1]
# protein_col_51 <- grep("protein$", colnames(chop_51_data), ignore.case=TRUE, value=TRUE)[1]
# 
# cat("Using these columns for sample 51:\n")
# cat("Peptide column:", peptide_col_51, "\n")
# cat("Spectral count column:", spectral_count_col_51, "\n")
# cat("Intensity column:", intensity_col_51, "\n")
# cat("Gene column:", gene_col_51, "\n")
# cat("Protein column:", protein_col_51, "\n")
# 
# # Use the identified column names
# chop_51_processed <- chop_51_data %>%
#   select(all_of(c(peptide_col_51, spectral_count_col_51, intensity_col_51, gene_col_51, protein_col_51)))
# 
# # Rename columns to standardized names with spaces, not periods
# names(chop_51_processed) <- c("Peptide Sequence", "FL51 Spectral Count", "FL51 Intensity", "Gene", "Protein")
# 
# # Process sample 88 data - similar approach
# cat("\nProcessing sample 88 data...\n")
# # Find the correct column names
# spectral_count_col_88 <- grep("spectral|count", colnames(chop_88_data), ignore.case=TRUE, value=TRUE)[1]
# intensity_col_88 <- grep("intensity", colnames(chop_88_data), ignore.case=TRUE, value=TRUE)[1]
# peptide_col_88 <- grep("peptide", colnames(chop_88_data), ignore.case=TRUE, value=TRUE)[1]
# gene_col_88 <- grep("gene", colnames(chop_88_data), ignore.case=TRUE, value=TRUE)[1]
# protein_col_88 <- grep("protein$", colnames(chop_88_data), ignore.case=TRUE, value=TRUE)[1]
# 
# cat("Using these columns for sample 88:\n")
# cat("Peptide column:", peptide_col_88, "\n")
# cat("Spectral count column:", spectral_count_col_88, "\n")
# cat("Intensity column:", intensity_col_88, "\n")
# cat("Gene column:", gene_col_88, "\n")
# cat("Protein column:", protein_col_88, "\n")
# 
# # Use the identified column names
# chop_88_processed <- chop_88_data %>%
#   select(all_of(c(peptide_col_88, spectral_count_col_88, intensity_col_88, gene_col_88, protein_col_88)))
# 
# # Rename columns to standardized names with spaces, not periods
# names(chop_88_processed) <- c("Peptide Sequence", "FL88 Spectral Count", "FL88 Intensity", "Gene", "Protein")
# 
# # Integrate the CHOP datasets
# cat("Integrating CHOP datasets...\n")
# # First merge with 51 data
# chop_integrated <- full_join(
#   chop_data, 
#   chop_51_processed,
#   by = "Peptide Sequence",
#   suffix = c("", ".51")
# )
# 
# # Then merge with 88 data
# chop_integrated <- full_join(
#   chop_integrated, 
#   chop_88_processed,
#   by = "Peptide Sequence",
#   suffix = c("", ".88")
# )
# 
# # Handle duplicate columns - consolidate gene and protein info
# cat("Handling duplicate columns...\n")
# # Function to consolidate data from multiple columns
# consolidate_columns <- function(df, base_col) {
#   # Check for duplicate columns
#   cols <- grep(paste0("^", base_col, "(\\.\\d+)?$"), names(df), value = TRUE)
#   
#   if (length(cols) > 1) {
#     # Create a consolidated column
#     df[[base_col]] <- apply(df[, cols, drop = FALSE], 1, function(row) {
#       non_na <- row[!is.na(row) & row != ""]
#       if (length(non_na) > 0) non_na[1] else NA
#     })
#     
#     # Remove duplicate columns
#     for (col in cols) {
#       if (col != base_col) {
#         df[[col]] <- NULL
#       }
#     }
#   }
#   
#   return(df)
# }
# 
# # Consolidate Gene and Protein columns
# chop_integrated <- consolidate_columns(chop_integrated, "Gene")
# chop_integrated <- consolidate_columns(chop_integrated, "Protein")
# 
# # Check for any remaining .x or .y or .51 or .88 columns
# duplicate_columns <- grep("\\.x$|\\.y$|\\.51$|\\.88$", names(chop_integrated), value = TRUE)
# if (length(duplicate_columns) > 0) {
#   cat("Warning: Some duplicate columns remain:", paste(duplicate_columns, collapse=", "), "\n")
# }
# 
# # Fix column names - replace periods with spaces in all column names
# names(chop_integrated) <- gsub("\\.", " ", names(chop_integrated))
# 
# # Ensure proper column order - move Gene and Protein to the end
# # First identify all columns except Gene and Protein
# non_gene_protein_cols <- setdiff(colnames(chop_integrated), c("Gene", "Protein"))
# # Then reorder columns
# chop_integrated <- chop_integrated[, c(non_gene_protein_cols, "Gene", "Protein")]
# 
# # Save the integrated CHOP file
# cat("Saving integrated CHOP file...\n")
# write.table(chop_integrated, file=integrated_chop_file, sep="\t", quote=FALSE, row.names=FALSE)
# 
# # Print CHOP integration summary
# cat("\nCHOP integration completed successfully!\n")
# cat("Original CHOP file:", nrow(chop_data), "rows,", ncol(chop_data), "columns\n")
# cat("CHOP 51 file:", nrow(chop_51_data), "rows\n")
# cat("CHOP 88 file:", nrow(chop_88_data), "rows\n")
# cat("Integrated CHOP file:", nrow(chop_integrated), "rows,", ncol(chop_integrated), "columns\n")
# 
# # Verify 51 and 88 samples
# fl51_cols <- grep("FL51", colnames(chop_integrated), value=TRUE)
# fl88_cols <- grep("FL88", colnames(chop_integrated), value=TRUE)
# 
# cat("\nFL51 columns in integrated file:", paste(fl51_cols, collapse=", "), "\n")
# cat("FL88 columns in integrated file:", paste(fl88_cols, collapse=", "), "\n")
# 
# # Count peptides by sample
# peptides_51 <- sum(chop_integrated$`FL51 Spectral Count` > 0, na.rm=TRUE)
# peptides_88 <- sum(chop_integrated$`FL88 Spectral Count` > 0, na.rm=TRUE)
# 
# cat("\nPeptides detected in sample 51:", peptides_51, "\n")
# cat("Peptides detected in sample 88:", peptides_88, "\n")
# 
# # =====================================================================
# # PART 2: MERGE MSKCC AND INTEGRATED CHOP DATA - MODIFIED SECTION
# # =====================================================================
# 
# # Add this right before the final merge operation
# cat("Preparing datasets for clean merging...\n")
# 
# # First, check if mskcc_wide has duplicate peptides
# if (any(duplicated(mskcc_wide$Peptide))) {
#   cat("WARNING: Found duplicate peptides in MSKCC_wide, deduplicating...\n")
#   
#   # Deduplicate the mskcc_wide data
#   mskcc_wide <- mskcc_wide %>%
#     group_by(Peptide) %>%
#     summarize(
#       # For Gene and Protein, keep one value (first non-NA)
#       Gene = first(na.omit(Gene)),
#       Protein = first(na.omit(Protein)),
#       # For all other columns, keep them (they already have unique names)
#       across(where(is.numeric), ~max(., na.rm = TRUE)),
#       .groups = 'drop'
#     )
#   
#   cat("After deduplication, mskcc_wide has", nrow(mskcc_wide), "rows\n")
# }
# 
# # Also check if chop_ready has duplicate peptides
# if (any(duplicated(chop_ready$Peptide))) {
#   cat("WARNING: Found duplicate peptides in chop_ready, deduplicating...\n")
#   
#   # Deduplicate the chop_ready data
#   chop_ready <- chop_ready %>%
#     group_by(Peptide) %>%
#     summarize(
#       # For Gene and Protein, keep one value (first non-NA)
#       Gene = first(na.omit(Gene)),
#       Protein = first(na.omit(Protein)),
#       # For all other columns, keep maximum values (conservative approach)
#       across(where(is.numeric), ~max(., na.rm = TRUE)),
#       .groups = 'drop'
#     )
#   
#   cat("After deduplication, chop_ready has", nrow(chop_ready), "rows\n")
# }
# 
# # Prepare integrated CHOP data for merging
# cat("Preparing integrated CHOP data for merging...\n")
# # Identify peptide sequence column in CHOP data
# if ("Peptide Sequence" %in% colnames(chop_integrated)) {
#   chop_prep <- chop_integrated
# } else {
#   peptide_col_chop <- grep("Peptide", colnames(chop_integrated), ignore.case=TRUE, value=TRUE)[1]
#   chop_prep <- chop_integrated
#   colnames(chop_prep)[colnames(chop_prep) == peptide_col_chop] <- "Peptide Sequence"
# }
# 
# # Now when working with MSKCC data, ensure it uses "Peptide Sequence" instead of "Peptide"
# if ("Peptide" %in% colnames(mskcc_data) && !"Peptide Sequence" %in% colnames(mskcc_data)) {
#   colnames(mskcc_data)[colnames(mskcc_data) == "Peptide"] <- "Peptide Sequence"
# }
# 
# # Define file paths for merging
# mskcc_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation_results_unique_peptides_unmodified_20250519_183234/unique_peptides_unmodified.tsv"
# merged_output_file <- "merged_immunopeptidomics_data.tsv"
# 
# # Read MSKCC data
# cat("Reading MSKCC data...\n")
# mskcc_data <- read.delim(mskcc_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# 
# # Print column names to debug
# cat("\nColumn names in mskcc_data:\n")
# print(colnames(mskcc_data))
# 
# # Add this after reading the data files but before the final merge
# # This section creates a separate table with all peptide annotations
# cat("Creating comprehensive peptide annotations table...\n")
# 
# # Function to concatenate unique values with a delimiter
# concatenate_unique <- function(x, delimiter = "; ") {
#   if (all(is.na(x)) || length(x) == 0) return(NA_character_)
#   unique_vals <- unique(x[!is.na(x) & x != ""])
#   if (length(unique_vals) == 0) return(NA_character_)
#   paste(sort(unique_vals), collapse = delimiter)
# }
# 
# # 1. Collect all MSKCC peptide annotations with sample information
# mskcc_annotations <- mskcc_data %>%
#   select(Peptide, Gene, Protein, SampleID) %>%
#   distinct() %>%
#   mutate(Source = "MSKCC")
# 
# # 2. Collect all CHOP peptide annotations
# # We need to handle the integrated CHOP data format which might have different column names
# if ("Peptide Sequence" %in% colnames(chop_integrated)) {
#   peptide_col_chop <- "Peptide Sequence"
# } else {
#   peptide_col_chop <- grep("Peptide", colnames(chop_integrated), ignore.case=TRUE, value=TRUE)[1]
# }
# 
# # Get sample IDs from CHOP data
# chop_sample_ids <- c()
# for (col in colnames(chop_integrated)) {
#   # Extract sample IDs from spectral count column names
#   if (grepl("Spectral Count", col, ignore.case=TRUE)) {
#     sample_match <- regexpr("FL[0-9]+|H[0-9]_[0-9]+", col)
#     if (sample_match > 0) {
#       sample_id <- regmatches(col, sample_match)
#       chop_sample_ids <- c(chop_sample_ids, sample_id)
#     }
#   }
# }
# chop_sample_ids <- unique(chop_sample_ids)
# 
# # Create a lookup table for CHOP samples
# chop_annotations <- chop_integrated %>%
#   select(all_of(c(peptide_col_chop, "Gene", "Protein"))) %>%
#   rename(Peptide = all_of(peptide_col_chop)) %>%
#   distinct() %>%
#   mutate(Source = "CHOP", SampleID = NA)
# 
# # 3. Combine all annotations
# all_annotations <- bind_rows(mskcc_annotations, chop_annotations)
# 
# # 4. Create a comprehensive peptide-to-annotation mapping table
# peptide_annotation_mapping <- all_annotations %>%
#   group_by(Peptide) %>%
#   summarize(
#     All_Genes = concatenate_unique(Gene),
#     All_Proteins = concatenate_unique(Protein),
#     MSKCC_SampleIDs = concatenate_unique(SampleID[Source == "MSKCC"]),
#     Sources = concatenate_unique(Source),
#     .groups = 'drop'
#   )
# 
# # 5. Save the comprehensive annotations table
# annotations_output_file <- "peptide_annotations_all.tsv"
# write.table(peptide_annotation_mapping, file = annotations_output_file, 
#             sep="\t", quote=FALSE, row.names=FALSE)
# 
# cat("Saved comprehensive peptide annotations to:", annotations_output_file, "\n")
# cat("This file contains all possible gene and protein annotations for each peptide.\n")
# 
# # Prepare MSKCC data for merging
# cat("Preparing MSKCC data for merging...\n")
# # Ensure numeric values
# spectral_count_col_mskcc <- grep("Spectral", colnames(mskcc_data), ignore.case=TRUE, value=TRUE)[1]
# intensity_col_mskcc <- grep("Intensity", colnames(mskcc_data), ignore.case=TRUE, value=TRUE)[1]
# 
# cat("Using these columns for MSKCC:\n")
# cat("Spectral count column:", spectral_count_col_mskcc, "\n")
# cat("Intensity column:", intensity_col_mskcc, "\n")
# 
# mskcc_data[[spectral_count_col_mskcc]] <- as.numeric(mskcc_data[[spectral_count_col_mskcc]])
# mskcc_data[[intensity_col_mskcc]] <- as.numeric(mskcc_data[[intensity_col_mskcc]])
# 
# # Aggregate data (in case of duplicates)
# mskcc_agg <- mskcc_data %>%
#   group_by(Peptide, SampleID) %>%
#   summarize(
#     Spectral_Count = sum(as.numeric(!!sym(spectral_count_col_mskcc)), na.rm = TRUE),
#     Intensity_Value = sum(as.numeric(!!sym(intensity_col_mskcc)), na.rm = TRUE),
#     Gene = first(Gene),
#     Protein = first(Protein),
#     .groups = 'drop'
#   )
# 
# # Create new columns for the wide format names - USING PREVIOUS APPROACH
# cat("Creating column names for MSKCC wide format...\n")
# mskcc_data_prep <- mskcc_agg %>%
#   mutate(
#     SC_ColName = paste0("MSKCC_", SampleID, "_Spectral_Count"),
#     Int_ColName = paste0("MSKCC_", SampleID, "_Intensity")
#   )
# 
# # Create spectral count wide format
# cat("Creating spectral count wide format...\n")
# mskcc_sc_wide <- mskcc_data_prep %>%
#   select(Peptide, Gene, Protein, SC_ColName, Spectral_Count) %>%
#   pivot_wider(
#     id_cols = c(Peptide, Gene, Protein),
#     names_from = SC_ColName,
#     values_from = Spectral_Count,
#     values_fill = list(Spectral_Count = 0)
#   )
# 
# # Create intensity wide format
# cat("Creating intensity wide format...\n")
# mskcc_int_wide <- mskcc_data_prep %>%
#   select(Peptide, Int_ColName, Intensity_Value) %>%
#   pivot_wider(
#     id_cols = Peptide,
#     names_from = Int_ColName,
#     values_from = Intensity_Value,
#     values_fill = list(Intensity_Value = 0)
#   )
# 
# # Combine all MSKCC data
# mskcc_wide <- left_join(mskcc_sc_wide, mskcc_int_wide, by = "Peptide")
# 
# # Prepare integrated CHOP data for merging
# cat("Preparing integrated CHOP data for merging...\n")
# # Identify peptide sequence column in CHOP data
# if ("Peptide" %in% colnames(chop_integrated)) {
#   chop_prep <- chop_integrated
# } else {
#   peptide_col_chop <- grep("Peptide", colnames(chop_integrated), ignore.case=TRUE, value=TRUE)[1]
#   chop_prep <- chop_integrated
#   colnames(chop_prep)[colnames(chop_prep) == peptide_col_chop] <- "Peptide"
# }
# 
# # Identify spectral count and intensity columns
# chop_sc_cols <- grep("Spectral Count", colnames(chop_prep), value = TRUE)
# chop_int_cols <- grep("Intensity$", colnames(chop_prep), value = TRUE)
# 
# # Select only relevant columns
# chop_ready <- chop_prep %>%
#   select(Peptide, Gene, Protein, all_of(chop_sc_cols), all_of(chop_int_cols))
# 
# # Merge the datasets
# cat("Merging MSKCC and CHOP datasets...\n")
# merged_data <- full_join(mskcc_wide, chop_ready, by = "Peptide", suffix = c("_MSKCC", "_CHOP"))
# 
# # Replace NA values with zeros for numeric columns
# cat("Handling missing values...\n")
# numeric_cols <- c(
#   grep("Spectral_Count", names(merged_data)),
#   grep("Intensity", names(merged_data))
# )
# 
# for (col in numeric_cols) {
#   merged_data[[col]] <- ifelse(is.na(merged_data[[col]]), 0, merged_data[[col]])
# }
# 
# # After merging, add a final safety check
# if (any(duplicated(merged_data$Peptide))) {
#   cat("WARNING: Found", sum(duplicated(merged_data$Peptide)), "duplicate peptides in merged data\n")
#   cat("Performing final deduplication...\n")
#   
#   # Deduplicate the merged data
#   merged_data <- merged_data %>%
#     group_by(Peptide) %>%
#     summarize(
#       # For Gene and Protein columns, keep first non-NA value
#       Gene_MSKCC = first(na.omit(Gene_MSKCC)),
#       Protein_MSKCC = first(na.omit(Protein_MSKCC)),
#       Gene_CHOP = first(na.omit(Gene_CHOP)),
#       Protein_CHOP = first(na.omit(Protein_CHOP)),
#       # For all other columns, take the maximum value
#       across(where(is.numeric), ~max(., na.rm = TRUE)),
#       .groups = 'drop'
#     )
#   
#   cat("After final deduplication, merged data has", nrow(merged_data), "rows\n")
# }
# 
# # Save the merged file
# cat("Saving merged file...\n")
# write.table(merged_data, file = merged_output_file, sep = "\t", quote = FALSE, row.names = FALSE)
# 
# # Print merging summary
# cat("\nMerging completed successfully!\n")
# cat("MSKCC dataset:", nrow(mskcc_data), "rows,", length(unique(mskcc_data$Peptide)), "unique peptides\n")
# cat("CHOP dataset:", nrow(chop_integrated), "rows,", length(unique(chop_prep$Peptide)), "unique peptides\n")
# cat("Merged dataset:", nrow(merged_data), "rows (unique peptides)\n")
# 
# # Calculate dataset overlap
# mskcc_peptides <- unique(mskcc_data$Peptide)
# chop_peptides <- unique(chop_prep$Peptide)
# common_peptides <- intersect(mskcc_peptides, chop_peptides)
# 
# cat("\nPeptide overlap:\n")
# cat("Peptides only in MSKCC:", length(setdiff(mskcc_peptides, chop_peptides)), "\n")
# cat("Peptides only in CHOP:", length(setdiff(chop_peptides, mskcc_peptides)), "\n")
# cat("Peptides in both datasets:", length(common_peptides), "\n")
# cat("Overlap percentage:",
#     round(100 * length(common_peptides) / length(union(mskcc_peptides, chop_peptides)), 1), "%\n")
# 
# # =====================================================================
# # FINAL SUMMARY
# # =====================================================================
# 
# cat("\n=====================================================================\n")
# cat("INTEGRATION COMPLETE\n")
# cat("=====================================================================\n\n")
# 
# cat("Files created:\n")
# cat("1. final_CHOP_combined_peptide.tsv - Integrated CHOP dataset with samples 51 and 88\n")
# cat("2. merged_immunopeptidomics_data.tsv - Merged MSKCC and CHOP data for visualization\n\n")
# 
# cat("Next steps:\n")
# cat("- Run your visualization script to create correlation plots and overlap statistics\n")
# cat("- These files contain all the data needed for your analysis\n")
# 
# cat("\nAll data preparation tasks completed successfully!\n")
