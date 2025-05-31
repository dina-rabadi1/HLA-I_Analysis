# Multi-Sample Immunopeptidome Analysis Script
# Adapted from RU148 comparison script for analyzing multiple samples
# This script analyzes immunopeptidome data and compares with transcriptome, LFQ proteome, and TMT proteome

# Setting directory
setwd("~/Documents/Github/HLA-I_Analysis/")

# Create output directory
output_dir <- "MultiSample_Analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Created output directory:", output_dir, "\n")
}

# Create a directory for visualizations
viz_dir <- file.path(output_dir, "MultiSample_visualizations")
if (!dir.exists(viz_dir)) {
  dir.create(viz_dir)
  cat("Created visualizations directory:", viz_dir, "\n")
}

# Load required packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(writexl)   # For Excel output
library(openxlsx)  # For better Excel formatting
library(readxl)    # For reading Excel files
library(VennDiagram) # For Venn diagram visualization
library(gridExtra)  # For arranging multiple plots
library(RColorBrewer) # For enhanced color palettes

# Define the path to your data files
immunopeptidome_path <- "combined_peptides.tsv"
transcriptome_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/Normalized_Gene_counts_FLCdb_Panel_1.xlsx"
lfq_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/Levin2023/adg7038_Table_S2_LFQ.xlsx"
tmt_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/Levin2023/adg7038_Table_S1_TMT.xlsx"

#--------------------------------------------------
# PART 1: Process the immunopeptidome data
#--------------------------------------------------

# Read the combined peptides TSV file
cat("Reading immunopeptidome data from:", immunopeptidome_path, "\n")
immunopeptidome_data <- read.delim(immunopeptidome_path, stringsAsFactors = FALSE)

# Print column names to ensure correct reading
cat("Immunopeptidome data columns:", paste(colnames(immunopeptidome_data), collapse=", "), "\n")

# Extract unique sample IDs
sample_ids <- unique(immunopeptidome_data$SampleID)
cat("Found", length(sample_ids), "unique sample IDs:", paste(sample_ids, collapse=", "), "\n")

# Clean up column names if necessary (remove any "X." prefixes or similar)
colnames(immunopeptidome_data) <- gsub("^X\\.", "", colnames(immunopeptidome_data))

# Create a summary of peptide detection for each sample
peptide_summary <- immunopeptidome_data %>%
  # Add column for peptide length
  mutate(peptide_length = nchar(Peptide)) %>%
  # FILTER: Keep only peptides with lengths 8-12 amino acids
  filter(peptide_length >= 8 & peptide_length <= 12) %>%
  group_by(SampleID, Peptide) %>%
  summarize(
    peptide_length = first(peptide_length),
    spectral_count = sum(Spectral.Count),
    intensity = sum(Intensity),
    protein_ids = paste(unique(Protein.ID), collapse = "; "),
    genes = paste(unique(Gene), collapse = "; "),
    source_filenames = paste(unique(SourceFile), collapse = "; "),
    .groups = "drop"
  )

# Print summary statistics
cat("\nSummary of filtered peptides (length 8-12):\n")
cat("Total unique peptides:", length(unique(peptide_summary$Peptide)), "\n")
cat("Peptides per sample:\n")
peptide_counts <- peptide_summary %>%
  group_by(SampleID) %>%
  summarise(peptide_count = n_distinct(Peptide), .groups = "drop")
print(peptide_counts)

#--------------------------------------------------
# PART 1.2: Create presence/absence matrix for peptides across samples
#--------------------------------------------------

# Create a presence/absence matrix
presence_matrix <- peptide_summary %>%
  # Create a binary indicator for presence
  mutate(present = 1) %>%
  # Spread to wide format
  pivot_wider(
    id_cols = Peptide,
    names_from = SampleID,
    values_from = present,
    values_fill = 0
  )

# Count the number of samples each peptide appears in
presence_matrix <- presence_matrix %>%
  mutate(samples_present = rowSums(select(., -Peptide)))

# Create filters for peptides present in at least X samples
peptides_in_5plus <- presence_matrix %>% filter(samples_present >= 5)
peptides_in_6plus <- presence_matrix %>% filter(samples_present >= 6)
peptides_in_7plus <- presence_matrix %>% filter(samples_present >= 7)
peptides_in_8plus <- presence_matrix %>% filter(samples_present >= 8)
peptides_in_9plus <- presence_matrix %>% filter(samples_present >= 9)

# Print summary of multi-sample peptides
cat("\nPeptides present in multiple samples:\n")
cat("In 5+ samples:", nrow(peptides_in_5plus), "\n")
cat("In 6+ samples:", nrow(peptides_in_6plus), "\n")
cat("In 7+ samples:", nrow(peptides_in_7plus), "\n")
cat("In 8+ samples:", nrow(peptides_in_8plus), "\n")
cat("In 9+ samples:", nrow(peptides_in_9plus), "\n")
cat("In all 10 samples:", sum(presence_matrix$samples_present == length(sample_ids)), "\n")

# Create a column to indicate presence in multiple samples
presence_categories <- presence_matrix %>%
  select(Peptide, samples_present) %>%
  mutate(
    in_5plus_samples = samples_present >= 5,
    in_6plus_samples = samples_present >= 6,
    in_7plus_samples = samples_present >= 7,
    in_8plus_samples = samples_present >= 8,
    in_9plus_samples = samples_present >= 9,
    in_all_samples = samples_present == length(sample_ids),
    presence_category = case_when(
      samples_present == length(sample_ids) ~ "All samples",
      samples_present >= 9 ~ "9+ samples",
      samples_present >= 8 ~ "8+ samples",
      samples_present >= 7 ~ "7+ samples",
      samples_present >= 6 ~ "6+ samples",
      samples_present >= 5 ~ "5+ samples",
      samples_present >= 3 ~ "3-4 samples",
      samples_present == 2 ~ "2 samples",
      TRUE ~ "1 sample"
    )
  )

# Create intensity matrix for heatmap visualization
intensity_matrix <- peptide_summary %>%
  select(Peptide, SampleID, intensity) %>%
  pivot_wider(
    id_cols = Peptide,
    names_from = SampleID,
    values_from = intensity,
    values_fill = 0
  )

#--------------------------------------------------
# PART 2: Define the fusion protein sequence and detection function
#--------------------------------------------------

# Define the fusion protein sequence
fusion_protein <- "RKREIFDRYGEEVKEFLAKAKEDF"
# Identify the fusion junction - between E (DNAJB1) and V (PRKACA)
dnajb1_part <- "RKREIFDRYGEE"  # Ends with YGEE (full DNAJB1 part)
prkaca_part <- "VKEFLAKAKEDF"  # Starts with VKEF (full PRKACA part)
junction_position <- nchar(dnajb1_part)  # Position after the last E in YGEE

# Print info about the fusion protein for verification
cat("\nFusion protein information:\n")
cat("DNAJB1 part:", dnajb1_part, "\n")
cat("PRKACA part:", prkaca_part, "\n")
cat("Junction position:", junction_position, "\n")
cat("Fusion protein:", fusion_protein, "\n")
cat("Checking junction: character at junction position is", substr(fusion_protein, junction_position, junction_position), "\n")
cat("Checking junction: next character is", substr(fusion_protein, junction_position+1, junction_position+1), "\n\n")

# Function to check if a peptide spans the fusion junction
is_fusion_junction_peptide <- function(peptide_seq) {
  # Check if the peptide spans the fusion junction
  spans_junction <- FALSE
  
  if (nchar(peptide_seq) >= 6) { # Require at least 6 AA to be meaningful
    for (i in 3:(nchar(peptide_seq) - 3)) { # At least 3 AA from each side
      left_part <- substr(peptide_seq, 1, i)
      right_part <- substr(peptide_seq, i + 1, nchar(peptide_seq))
      
      # Check if left part is in DNAJB1 and right part in PRKACA
      if (grepl(left_part, dnajb1_part, fixed = TRUE) && 
          grepl(right_part, prkaca_part, fixed = TRUE)) {
        
        # Additional check to ensure left part aligns with end of DNAJB1
        left_pos <- gregexpr(left_part, dnajb1_part, fixed = TRUE)[[1]]
        if (length(left_pos) > 0 && any(left_pos + nchar(left_part) - 1 >= nchar(dnajb1_part) - 5)) {
          
          # Additional check to ensure right part aligns with start of PRKACA
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
  from_dnajb1 <- grepl(peptide_seq, dnajb1_part, fixed = TRUE)
  from_prkaca <- grepl(peptide_seq, prkaca_part, fixed = TRUE)
  from_fusion <- grepl(peptide_seq, fusion_protein, fixed = TRUE)
  
  return(list(
    spans_junction = spans_junction,
    from_dnajb1 = from_dnajb1,
    from_prkaca = from_prkaca,
    from_fusion = from_fusion | spans_junction
  ))
}

# Add fusion protein information to all peptides
peptide_presence_fusion <- presence_categories %>%
  rowwise() %>%
  mutate(
    fusion_info = list(is_fusion_junction_peptide(Peptide)),
    from_fusion = fusion_info$from_fusion,
    spans_junction = fusion_info$spans_junction,
    from_dnajb1 = fusion_info$from_dnajb1,
    from_prkaca = fusion_info$from_prkaca,
    fusion_peptide_type = case_when(
      spans_junction ~ "Junction-spanning",
      from_dnajb1 ~ "DNAJB1 part",
      from_prkaca ~ "PRKACA part",
      TRUE ~ "Not from fusion"
    )
  ) %>%
  select(-fusion_info)

# Create a separate analysis specifically for fusion peptides
fusion_peptides_analysis <- peptide_presence_fusion %>%
  filter(from_fusion) %>%
  arrange(desc(spans_junction), desc(samples_present))

# Print summary of fusion peptides
cat("\nFusion peptide analysis:\n")
cat("Total fusion-related peptides:", nrow(fusion_peptides_analysis), "\n")
if(nrow(fusion_peptides_analysis) > 0) {
  cat("Junction-spanning peptides:", sum(fusion_peptides_analysis$spans_junction), "\n")
  cat("DNAJB1-part peptides:", sum(fusion_peptides_analysis$from_dnajb1), "\n")
  cat("PRKACA-part peptides:", sum(fusion_peptides_analysis$from_prkaca), "\n")
  
  # Summary by sample count
  fusion_sample_summary <- fusion_peptides_analysis %>%
    group_by(samples_present, fusion_peptide_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(samples_present), fusion_peptide_type)
  
  cat("\nFusion peptides by sample count:\n")
  print(fusion_sample_summary)
}

#--------------------------------------------------
# PART 3: Process the transcriptome data
#--------------------------------------------------

# Read the transcriptome data
cat("\nReading transcriptome data from:", transcriptome_path, "\n")
transcriptome_data <- read_excel(transcriptome_path)

# Print column names to determine available samples
cat("Transcriptome data contains", ncol(transcriptome_data), "columns\n")

# Map immunopeptidome sample IDs to transcriptome column names
# Create mapping for each sample ID
sample_mapping <- c(
  "51" = "RU51_M49",    # Using M49 for sample 51
  "57" = "RU59_M6",     # Using nearest available match for 57 (not exact)
  "59" = "RU59_M6",     # Using M6 for sample 59
  "62" = "RU63_R2",     # Using nearest available match for 62 (not exact)
  "63" = "RU63_R2",     # Using R2 for sample 63
  "88" = "RU88_M15",    # Using M15 for sample 88
  "117" = "RU117_M10",  # Using M10 for sample 117
  "123" = "RU123_MT",   # Using MT for sample 123
  "148T" = "RU148_T8",  # Using T8 for sample 148T
  "148N" = "RU148_N"    # Using N for sample 148N
)

# Verify all mapped columns exist
missing_columns <- setdiff(unname(sample_mapping), colnames(transcriptome_data))
if(length(missing_columns) > 0) {
  warning("Missing transcriptome columns: ", paste(missing_columns, collapse=", "))
} else {
  cat("All mapped transcriptome columns found\n")
}

# Extract the relevant transcriptome columns
relevant_cols <- c("symbol", unname(sample_mapping))
relevant_cols <- relevant_cols[relevant_cols %in% colnames(transcriptome_data)]

transcriptome_filtered <- transcriptome_data %>%
  select(all_of(relevant_cols))

# Calculate mean expression across all samples for each gene
transcriptome_filtered <- transcriptome_filtered %>%
  rowwise() %>%
  mutate(
    mean_expression = mean(c_across(all_of(unname(sample_mapping)[unname(sample_mapping) %in% colnames(transcriptome_data)])), na.rm = TRUE)
  ) %>%
  ungroup()

# Print summary of processed transcriptome data
cat("Processed", nrow(transcriptome_filtered), "genes from transcriptome data\n")
cat("Sample columns kept:", paste(colnames(transcriptome_filtered)[2:(length(sample_mapping)+1)], collapse=", "), "\n")

#--------------------------------------------------
# PART 3.1: Process the LFQ proteome data
#--------------------------------------------------

# Read the LFQ data
cat("\nReading LFQ proteome data from:", lfq_path, "\n")
lfq_data <- tryCatch({
  read_excel(lfq_path)
}, error = function(e) {
  # Try reading a specific sheet if the file has multiple sheets
  sheets <- excel_sheets(lfq_path)
  if(length(sheets) > 0) {
    cat("Attempting to read first sheet:", sheets[1], "\n")
    read_excel(lfq_path, sheet = sheets[1])
  } else {
    stop("Error reading LFQ data: ", e$message)
  }
})

# Print column names to confirm structure
cat("LFQ proteome data columns:", paste(colnames(lfq_data), collapse=", "), "\n")

# Process LFQ data
# If the column structures match what was described:
if("Gene Name" %in% colnames(lfq_data) && "Log2 Difference" %in% colnames(lfq_data)) {
  lfq_processed <- lfq_data %>%
    # Ensure consistent column names
    rename_with(~ gsub(" ", "_", .), everything()) %>%
    # Filter for significant entries if p-value available
    filter(if("P.value" %in% colnames(.)) P.value <= 0.05 else TRUE) %>%
    # Clean up and standardize
    mutate(
      Gene_Name = trimws(Gene_Name),
      log2_fold_change_lfq = Log2_Difference,
      protein_category_lfq = case_when(
        Log2_Difference > 1 ~ "Up (FC > 2)",
        Log2_Difference < -1 ~ "Down (FC < 0.5)",
        TRUE ~ "Similar (-1 < log2FC < 1)"
      )
    ) %>%
    # Select relevant columns
    select(Gene_Name, log2_fold_change_lfq, Protein_Name, protein_category_lfq)
} else {
  # If column names don't match exactly, try to adapt
  # You may need to adjust this based on your actual column names
  cat("Column names don't match expected structure, attempting to adapt\n")
  lfq_processed <- lfq_data %>%
    rename_with(~ gsub(" ", "_", .), everything())
  
  # Try to find gene name column
  gene_col <- grep("Gene|gene", colnames(lfq_processed), ignore.case = TRUE)[1]
  fc_col <- grep("Difference|FC|log2|fold", colnames(lfq_processed), ignore.case = TRUE)[1]
  
  if(!is.na(gene_col) && !is.na(fc_col)) {
    gene_col_name <- colnames(lfq_processed)[gene_col]
    fc_col_name <- colnames(lfq_processed)[fc_col]
    
    lfq_processed <- lfq_processed %>%
      rename(Gene_Name = !!gene_col_name,
             log2_fold_change_lfq = !!fc_col_name) %>%
      mutate(
        Gene_Name = trimws(Gene_Name),
        protein_category_lfq = case_when(
          log2_fold_change_lfq > 1 ~ "Up (FC > 2)",
          log2_fold_change_lfq < -1 ~ "Down (FC < 0.5)",
          TRUE ~ "Similar (-1 < log2FC < 1)"
        )
      ) %>%
      select(Gene_Name, log2_fold_change_lfq, protein_category_lfq)
  } else {
    warning("Could not identify proper columns in LFQ data, creating placeholder")
    lfq_processed <- data.frame(
      Gene_Name = character(),
      log2_fold_change_lfq = numeric(),
      protein_category_lfq = character(),
      stringsAsFactors = FALSE
    )
  }
}

cat("Processed", nrow(lfq_processed), "entries from LFQ proteome data\n")

#--------------------------------------------------
# PART 3.2: Process the TMT proteome data
#--------------------------------------------------

# Read the TMT data
cat("\nReading TMT proteome data from:", tmt_path, "\n")
tmt_data <- tryCatch({
  read_excel(tmt_path)
}, error = function(e) {
  # Try reading a specific sheet if the file has multiple sheets
  sheets <- excel_sheets(tmt_path)
  if(length(sheets) > 0) {
    cat("Attempting to read first sheet:", sheets[1], "\n")
    read_excel(tmt_path, sheet = sheets[1])
  } else {
    stop("Error reading TMT data: ", e$message)
  }
})

# Print column names to confirm structure
cat("TMT proteome data columns:", paste(colnames(tmt_data), collapse=", "), "\n")

# Process TMT data
# If the column structures match what was described:
if("Gene Name" %in% colnames(tmt_data) && "Log2 Difference" %in% colnames(tmt_data)) {
  tmt_processed <- tmt_data %>%
    # Ensure consistent column names
    rename_with(~ gsub(" ", "_", .), everything()) %>%
    # Filter for significant entries if p-value available
    filter(if("P.value" %in% colnames(.)) P.value <= 0.05 else TRUE) %>%
    # Clean up and standardize
    mutate(
      Gene_Name = trimws(Gene_Name),
      log2_fold_change_tmt = Log2_Difference,
      protein_category_tmt = case_when(
        Log2_Difference > 1 ~ "Up (FC > 2)",
        Log2_Difference < -1 ~ "Down (FC < 0.5)",
        TRUE ~ "Similar (-1 < log2FC < 1)"
      )
    ) %>%
    # Select relevant columns
    select(Gene_Name, log2_fold_change_tmt, Protein_Name, protein_category_tmt)
} else {
  # If column names don't match exactly, try to adapt
  cat("Column names don't match expected structure, attempting to adapt\n")
  tmt_processed <- tmt_data %>%
    rename_with(~ gsub(" ", "_", .), everything())
  
  # Try to find gene name column
  gene_col <- grep("Gene|gene", colnames(tmt_processed), ignore.case = TRUE)[1]
  fc_col <- grep("Difference|FC|log2|fold", colnames(tmt_processed), ignore.case = TRUE)[1]
  
  if(!is.na(gene_col) && !is.na(fc_col)) {
    gene_col_name <- colnames(tmt_processed)[gene_col]
    fc_col_name <- colnames(tmt_processed)[fc_col]
    
    tmt_processed <- tmt_processed %>%
      rename(Gene_Name = !!gene_col_name,
             log2_fold_change_tmt = !!fc_col_name) %>%
      mutate(
        Gene_Name = trimws(Gene_Name),
        protein_category_tmt = case_when(
          log2_fold_change_tmt > 1 ~ "Up (FC > 2)",
          log2_fold_change_tmt < -1 ~ "Down (FC < 0.5)",
          TRUE ~ "Similar (-1 < log2FC < 1)"
        )
      ) %>%
      select(Gene_Name, log2_fold_change_tmt, protein_category_tmt)
  } else {
    warning("Could not identify proper columns in TMT data, creating placeholder")
    tmt_processed <- data.frame(
      Gene_Name = character(),
      log2_fold_change_tmt = numeric(),
      protein_category_tmt = character(),
      stringsAsFactors = FALSE
    )
  }
}

cat("Processed", nrow(tmt_processed), "entries from TMT proteome data\n")

#--------------------------------------------------
# PART 4: Extract gene information for peptides and merge with transcriptome/proteome data
#--------------------------------------------------

# Extract a list of unique peptides with gene information
peptide_gene_mapping <- peptide_summary %>%
  select(Peptide, genes) %>%
  distinct() %>%
  # Extract primary gene for each peptide
  mutate(
    primary_gene = sapply(strsplit(genes, ";\\s*"), function(x) trimws(x[1]))
  )

# Add this primary gene information to the presence matrix
presence_with_genes <- peptide_presence_fusion %>%
  left_join(peptide_gene_mapping, by = "Peptide")

# Add intensity information from each sample for visualization
# First, ensure presence_with_genes has peptides in the same order as intensity_matrix
# Then add mean intensity across samples
peptide_intensity_summary <- peptide_summary %>%
  group_by(Peptide) %>%
  summarise(
    mean_intensity = mean(intensity, na.rm = TRUE),
    max_intensity = max(intensity, na.rm = TRUE),
    .groups = "drop"
  )

presence_with_genes <- presence_with_genes %>%
  left_join(peptide_intensity_summary, by = "Peptide")

#--------------------------------------------------
# PART 4.1: Merge with transcriptome data
#--------------------------------------------------

# Merge with transcriptome data based on primary gene
multi_omics_data <- presence_with_genes %>%
  left_join(
    transcriptome_filtered,
    by = c("primary_gene" = "symbol"),
    relationship = "many-to-many"  # Explicitly acknowledge many-to-many relationship
  )

# Add transcriptome expression categorization
multi_omics_data <- multi_omics_data %>%
  mutate(
    transcriptome_expression_category = case_when(
      is.na(mean_expression) ~ "No transcriptome data",
      mean_expression > median(transcriptome_filtered$mean_expression, na.rm = TRUE) * 2 ~ "High expression",
      mean_expression > median(transcriptome_filtered$mean_expression, na.rm = TRUE) ~ "Medium-high expression",
      mean_expression > median(transcriptome_filtered$mean_expression, na.rm = TRUE) / 2 ~ "Medium-low expression",
      TRUE ~ "Low expression"
    )
  )

#--------------------------------------------------
# PART 4.2: Merge with LFQ proteome data
#--------------------------------------------------

multi_omics_data <- multi_omics_data %>%
  left_join(
    lfq_processed,
    by = c("primary_gene" = "Gene_Name"),
    relationship = "many-to-many"
  )

#--------------------------------------------------
# PART 4.3: Merge with TMT proteome data
#--------------------------------------------------

multi_omics_data <- multi_omics_data %>%
  left_join(
    tmt_processed,
    by = c("primary_gene" = "Gene_Name"),
    relationship = "many-to-many"
  )

#--------------------------------------------------
# PART 5: Define criteria for potential public neoantigens
#--------------------------------------------------

multi_omics_data <- multi_omics_data %>%
  mutate(
    # Criteria 1: Present in at least 5 samples and high expression in transcriptome
    high_presence_high_transcriptome = case_when(
      samples_present >= 5 & 
        !is.na(mean_expression) & 
        transcriptome_expression_category %in% c("High expression", "Medium-high expression") ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Criteria 2: Present in at least 5 samples and upregulated in LFQ
    high_presence_up_lfq = case_when(
      samples_present >= 5 & 
        !is.na(log2_fold_change_lfq) & 
        log2_fold_change_lfq > 1 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Criteria 3: Present in at least 5 samples and upregulated in TMT
    high_presence_up_tmt = case_when(
      samples_present >= 5 & 
        !is.na(log2_fold_change_tmt) & 
        log2_fold_change_tmt > 1 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Criteria 4: Present in most samples (8+) and detected in either proteome dataset
    very_high_presence_and_proteome = case_when(
      samples_present >= 8 & 
        (!is.na(log2_fold_change_lfq) | !is.na(log2_fold_change_tmt)) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Criteria 5: FUSION PEPTIDE present in 3+ samples - lower threshold for fusion due to importance
    fusion_in_multiple_samples = case_when(
      from_fusion & samples_present >= 3 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Combined neoantigen score (weighted sum of criteria)
    public_neoantigen_score = as.integer(high_presence_high_transcriptome) * 2 + 
      as.integer(high_presence_up_lfq) * 2 + 
      as.integer(high_presence_up_tmt) * 2 + 
      as.integer(very_high_presence_and_proteome) * 3 +
      as.integer(fusion_in_multiple_samples) * 5,
    
    # Classify as potential public neoantigen if score above threshold
    potential_public_neoantigen = public_neoantigen_score >= 3,
    
    # Classification label
    public_neoantigen_classification = case_when(
      fusion_in_multiple_samples & spans_junction ~ "Tier 1A (Junction-spanning fusion peptide in 3+ samples)",
      fusion_in_multiple_samples ~ "Tier 1B (Fusion peptide in 3+ samples)",
      public_neoantigen_score >= 7 ~ "Tier 2 (Multiple strong criteria)",
      public_neoantigen_score >= 5 ~ "Tier 3 (Several criteria met)",
      public_neoantigen_score >= 3 ~ "Tier 4 (Some criteria met)",
      TRUE ~ "Not a potential public neoantigen"
    )
  )

# Make sure peptide_length is available in multi_omics_data
# It's likely that this column was not properly carried over during the joining process
# Let's add it from the presence_with_genes dataset which should have it

# First, check if peptide_length exists in multi_omics_data
if(!"peptide_length" %in% colnames(multi_omics_data)) {
  # If not, join it from peptide_summary
  # Create a mapping of peptide to peptide_length
  peptide_length_mapping <- peptide_summary %>%
    distinct(Peptide, peptide_length)
  
  # Add this to multi_omics_data
  multi_omics_data <- multi_omics_data %>%
    left_join(peptide_length_mapping, by = "Peptide")
}

# Now create the filtered dataset for potential public neoantigens
public_neoantigens <- multi_omics_data %>%
  filter(potential_public_neoantigen) %>%
  arrange(desc(public_neoantigen_score), desc(samples_present)) %>%
  select(
    Peptide, 
    primary_gene, 
    genes,
    peptide_length,
    samples_present,
    presence_category,
    public_neoantigen_classification,
    public_neoantigen_score,
    mean_expression,
    transcriptome_expression_category,
    log2_fold_change_lfq,
    protein_category_lfq,
    log2_fold_change_tmt,
    protein_category_tmt,
    from_fusion,
    spans_junction,
    mean_intensity
  )

# Count the number of potential public neoantigens
cat("\nPotential public neoantigens found:", nrow(public_neoantigens), "\n")
if(nrow(public_neoantigens) > 0) {
  tier_counts <- public_neoantigens %>%
    group_by(public_neoantigen_classification) %>%
    summarise(count = n(), .groups = "drop")
  
  print(tier_counts)
}

# Check for fusion-derived public neoantigens
fusion_public_neoantigens <- public_neoantigens %>%
  filter(from_fusion)

if(nrow(fusion_public_neoantigens) > 0) {
  cat("\nFound", nrow(fusion_public_neoantigens), "fusion-derived potential public neoantigens\n")
  cat("Junction-spanning fusion public neoantigens:", sum(fusion_public_neoantigens$spans_junction), "\n")
  
  # Print the fusion-derived public neoantigens
  cat("\nFusion-derived public neoantigens:\n")
  print(fusion_public_neoantigens %>% 
          select(Peptide, peptide_length, samples_present, spans_junction, public_neoantigen_classification))
}

#--------------------------------------------------
# PART 6: Create final multi-omics dataset with key metrics
#--------------------------------------------------

# Create a comprehensive dataset with all peptides and their multi-omics information
multi_omics_final <- multi_omics_data %>%
  select(
    # Peptide information
    Peptide,
    peptide_length,
    primary_gene,
    genes,
    
    # Sample presence information
    samples_present,
    presence_category,
    in_5plus_samples,
    in_6plus_samples,
    in_7plus_samples,
    in_8plus_samples,
    in_9plus_samples,
    
    # Intensity information
    mean_intensity,
    max_intensity,
    
    # Transcriptome data
    mean_expression,
    transcriptome_expression_category,
    
    # LFQ proteome data
    log2_fold_change_lfq,
    protein_category_lfq,
    
    # TMT proteome data
    log2_fold_change_tmt,
    protein_category_tmt,
    
    # Fusion information
    from_fusion,
    spans_junction,
    from_dnajb1,
    from_prkaca,
    fusion_peptide_type,
    
    # Public neoantigen info
    potential_public_neoantigen,
    public_neoantigen_classification,
    public_neoantigen_score,
    high_presence_high_transcriptome,
    high_presence_up_lfq,
    high_presence_up_tmt,
    very_high_presence_and_proteome,
    fusion_in_multiple_samples
  )

# Now we can safely arrange by the newly created variables
multi_omics_final <- multi_omics_final %>%
  arrange(desc(potential_public_neoantigen), desc(samples_present), desc(mean_intensity))

#--------------------------------------------------
# PART 7: Create visualizations
#--------------------------------------------------

#--------------------------------------------------
# PART 7.1: Sample presence heatmap
#--------------------------------------------------

# Create a heatmap of peptide presen#--------------------------------------------------
# Utility function for safely setting row names in heatmaps
#--------------------------------------------------

#--------------------------------------------------
# Utility function for safely setting row names in heatmaps
#--------------------------------------------------

# Function to safely set row names for heatmaps
safe_set_rownames <- function(data, id_column) {
  if(any(duplicated(data[[id_column]]))) {
    warning("Duplicate values found in ", id_column, ". Making row names unique.")
    data$unique_id <- make.unique(as.character(data[[id_column]]))
    rownames(data) <- data$unique_id
    data <- data %>% select(-all_of(c(id_column, "unique_id")))
  } else {
    rownames(data) <- data[[id_column]]
    data <- data %>% select(-all_of(id_column))
  }
  return(data)
}

# First, identify high-value peptides for the heatmap (public neoantigens and/or fusion peptides)
highlight_peptides <- multi_omics_final %>%
  filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
  pull(Peptide)

# If there are too many peptides, limit to top ones by score
if(length(highlight_peptides) > 100) {
  highlight_peptides <- multi_omics_final %>%
    filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
    arrange(desc(public_neoantigen_score), desc(samples_present)) %>%
    head(100) %>%
    pull(Peptide)
}

# Create a presence matrix just for these highlighted peptides
# Function to safely set row names for heatmaps
safe_set_rownames <- function(data, id_column) {
  if(any(duplicated(data[[id_column]]))) {
    warning("Duplicate values found in ", id_column, ". Making row names unique.")
    data$unique_id <- make.unique(as.character(data[[id_column]]))
    rownames(data) <- data$unique_id
    data <- data %>% select(-all_of(c(id_column, "unique_id")))
  } else {
    rownames(data) <- data[[id_column]]
    data <- data %>% select(-all_of(id_column))
  }
  return(data)
}

# First, identify high-value peptides for the heatmap (public neoantigens and/or fusion peptides)
highlight_peptides <- multi_omics_final %>%
  filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
  pull(Peptide)

# If there are too many peptides, limit to top ones by score
if(length(highlight_peptides) > 100) {
  highlight_peptides <- multi_omics_final %>%
    filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
    arrange(desc(public_neoantigen_score), desc(samples_present)) %>%
    head(100) %>%
    pull(Peptide)
}

# Create a presence matrix just for these highlighted peptides
# Function to safely set row names for heatmaps
safe_set_rownames <- function(data, id_column) {
  if(any(duplicated(data[[id_column]]))) {
    warning("Duplicate values found in ", id_column, ". Making row names unique.")
    data$unique_id <- make.unique(as.character(data[[id_column]]))
    rownames(data) <- data$unique_id
    data <- data %>% select(-all_of(c(id_column, "unique_id")))
  } else {
    rownames(data) <- data[[id_column]]
    data <- data %>% select(-all_of(id_column))
  }
  return(data)
}

# First, identify high-value peptides for the heatmap (public neoantigens and/or fusion peptides)
highlight_peptides <- multi_omics_final %>%
  filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
  pull(Peptide)

# If there are too many peptides, limit to top ones by score
if(length(highlight_peptides) > 100) {
  highlight_peptides <- multi_omics_final %>%
    filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
    arrange(desc(public_neoantigen_score), desc(samples_present)) %>%
    head(100) %>%
    pull(Peptide)
}

# Create a presence matrix just for these highlighted peptides
# Function to safely set row names for heatmaps
safe_set_rownames <- function(data, id_column) {
  if(any(duplicated(data[[id_column]]))) {
    warning("Duplicate values found in ", id_column, ". Making row names unique.")
    data$unique_id <- make.unique(as.character(data[[id_column]]))
    rownames(data) <- data$unique_id
    data <- data %>% select(-all_of(c(id_column, "unique_id")))
  } else {
    rownames(data) <- data[[id_column]]
    data <- data %>% select(-all_of(id_column))
  }
  return(data)
}

# First, identify high-value peptides for the heatmap (public neoantigens and/or fusion peptides)
highlight_peptides <- multi_omics_final %>%
  filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
  pull(Peptide)

# If there are too many peptides, limit to top ones by score
if(length(highlight_peptides) > 100) {
  highlight_peptides <- multi_omics_final %>%
    filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
    arrange(desc(public_neoantigen_score), desc(samples_present)) %>%
    head(100) %>%
    pull(Peptide)
}

# Create a presence matrix just for these highlighted peptides
# Function to safely set row names for heatmaps
safe_set_rownames <- function(data, id_column) {
  if(any(duplicated(data[[id_column]]))) {
    warning("Duplicate values found in ", id_column, ". Making row names unique.")
    data$unique_id <- make.unique(as.character(data[[id_column]]))
    rownames(data) <- data$unique_id
    data <- data %>% select(-all_of(c(id_column, "unique_id")))
  } else {
    rownames(data) <- data[[id_column]]
    data <- data %>% select(-all_of(id_column))
  }
  return(data)
}

# First, identify high-value peptides for the heatmap (public neoantigens and/or fusion peptides)
highlight_peptides <- multi_omics_final %>%
  filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
  pull(Peptide)

# If there are too many peptides, limit to top ones by score
if(length(highlight_peptides) > 100) {
  highlight_peptides <- multi_omics_final %>%
    filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
    arrange(desc(public_neoantigen_score), desc(samples_present)) %>%
    head(100) %>%
    pull(Peptide)
}

# Create a presence matrix just for these highlighted peptides
# Function to safely set row names for heatmaps
safe_set_rownames <- function(data, id_column) {
  if(any(duplicated(data[[id_column]]))) {
    warning("Duplicate values found in ", id_column, ". Making row names unique.")
    data$unique_id <- make.unique(as.character(data[[id_column]]))
    rownames(data) <- data$unique_id
    data <- data %>% select(-all_of(c(id_column, "unique_id")))
  } else {
    rownames(data) <- data[[id_column]]
    data <- data %>% select(-all_of(id_column))
  }
  return(data)
}

# First, identify high-value peptides for the heatmap (public neoantigens and/or fusion peptides)
highlight_peptides <- multi_omics_final %>%
  filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
  pull(Peptide)

# If there are too many peptides, limit to top ones by score
if(length(highlight_peptides) > 100) {
  highlight_peptides <- multi_omics_final %>%
    filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
    arrange(desc(public_neoantigen_score), desc(samples_present)) %>%
    head(100) %>%
    pull(Peptide)
}

# Create a presence matrix just for these highlighted peptides
if(length(highlight_peptides) > 0) {
  # Get unique peptides to prevent duplicates
  highlight_peptides <- unique(highlight_peptides)
  
  # Get the annotation data first with distinct peptides
  annotation_data <- multi_omics_final %>% 
    filter(Peptide %in% highlight_peptides) %>%
    select(Peptide, public_neoantigen_classification, from_fusion, spans_junction) %>%
    distinct(Peptide, .keep_all = TRUE)
  
  # Make sure annotation data has no NA values
  annotation_data <- annotation_data %>%
    mutate(
      public_neoantigen_classification = if_else(
        is.na(public_neoantigen_classification), 
        "Not classified", 
        public_neoantigen_classification
      ),
      from_fusion = if_else(is.na(from_fusion), FALSE, from_fusion),
      spans_junction = if_else(is.na(spans_junction), FALSE, spans_junction)
    )
  
  # Now join with the presence matrix
  highlight_presence <- presence_matrix %>%
    filter(Peptide %in% highlight_peptides) %>%
    left_join(annotation_data, by = "Peptide")
  
  # Prepare data for heatmap
  heatmap_data <- highlight_presence %>%
    select(-samples_present, -public_neoantigen_classification, -from_fusion, -spans_junction)
  
  # Use our safe function to set row names
  heatmap_data <- safe_set_rownames(heatmap_data, "Peptide")
  
  # Create annotation for the rows with matching rownames
  row_annotation <- data.frame(
    Classification = annotation_data$public_neoantigen_classification,
    From_Fusion = annotation_data$from_fusion,
    Spans_Junction = annotation_data$spans_junction,
    row.names = annotation_data$Peptide
  )
  
  # Make sure the row names in both data frames match
  common_peptides <- intersect(rownames(heatmap_data), rownames(row_annotation))
  heatmap_data <- heatmap_data[common_peptides, ]
  row_annotation <- row_annotation[common_peptides, ]
  
  # Define annotation colors explicitly
  annotation_colors <- list(
    Classification = c(
      "Tier 1A (Junction-spanning fusion peptide in 3+ samples)" = "red",
      "Tier 1B (Fusion peptide in 3+ samples)" = "orange",
      "Tier 2 (Multiple strong criteria)" = "yellow",
      "Tier 3 (Several criteria met)" = "green",
      "Tier 4 (Some criteria met)" = "blue",
      "Not a potential public neoantigen" = "gray",
      "Not classified" = "lightgray"
    ),
    From_Fusion = c("TRUE" = "purple", "FALSE" = "lightblue"),
    Spans_Junction = c("TRUE" = "brown", "FALSE" = "beige")
  )
  
  # Create the heatmap PDF
  pdf(file.path(viz_dir, "Peptide_Presence_Heatmap.pdf"), 
      width = 12, height = max(8, nrow(heatmap_data)/5))
  
  # Heatmap with annotations and explicit colors
  heatmap_plot <- pheatmap(
    heatmap_data,
    main = "Peptide Presence Across Samples",
    color = colorRampPalette(c("white", "steelblue"))(2),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    annotation_colors = annotation_colors,
    fontsize_row = 8,
    fontsize_col = 10
  )
  
  print(heatmap_plot)
  dev.off()
  
  # Create the heatmap PNG
  png(file.path(viz_dir, "Peptide_Presence_Heatmap.png"), 
      width = 1200, height = max(800, nrow(heatmap_data)*15), res = 120)
  
  print(heatmap_plot)
  dev.off()
}

#--------------------------------------------------
# PART 7.2: Intensity heatmap
#--------------------------------------------------

# Create intensity matrix for the highlighted peptides
# Function to safely set row names for heatmaps
safe_set_rownames <- function(data, id_column) {
  if(any(duplicated(data[[id_column]]))) {
    warning("Duplicate values found in ", id_column, ". Making row names unique.")
    data$unique_id <- make.unique(as.character(data[[id_column]]))
    rownames(data) <- data$unique_id
    data <- data %>% select(-all_of(c(id_column, "unique_id")))
  } else {
    rownames(data) <- data[[id_column]]
    data <- data %>% select(-all_of(id_column))
  }
  return(data)
}

# First, identify high-value peptides for the heatmap (public neoantigens and/or fusion peptides)
highlight_peptides <- multi_omics_final %>%
  filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
  pull(Peptide)

# If there are too many peptides, limit to top ones by score
if(length(highlight_peptides) > 100) {
  highlight_peptides <- multi_omics_final %>%
    filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
    arrange(desc(public_neoantigen_score), desc(samples_present)) %>%
    head(100) %>%
    pull(Peptide)
}

# Create a presence matrix just for these highlighted peptides
# Function to safely set row names for heatmaps
safe_set_rownames <- function(data, id_column) {
  if(any(duplicated(data[[id_column]]))) {
    warning("Duplicate values found in ", id_column, ". Making row names unique.")
    data$unique_id <- make.unique(as.character(data[[id_column]]))
    rownames(data) <- data$unique_id
    data <- data %>% select(-all_of(c(id_column, "unique_id")))
  } else {
    rownames(data) <- data[[id_column]]
    data <- data %>% select(-all_of(id_column))
  }
  return(data)
}

# First, identify high-value peptides for the heatmap (public neoantigens and/or fusion peptides)
highlight_peptides <- multi_omics_final %>%
  filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
  pull(Peptide)

# If there are too many peptides, limit to top ones by score
if(length(highlight_peptides) > 100) {
  highlight_peptides <- multi_omics_final %>%
    filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
    arrange(desc(public_neoantigen_score), desc(samples_present)) %>%
    head(100) %>%
    pull(Peptide)
}

# Create a presence matrix just for these highlighted peptides
# Function to safely set row names for heatmaps
safe_set_rownames <- function(data, id_column) {
  if(any(duplicated(data[[id_column]]))) {
    warning("Duplicate values found in ", id_column, ". Making row names unique.")
    data$unique_id <- make.unique(as.character(data[[id_column]]))
    rownames(data) <- data$unique_id
    data <- data %>% select(-all_of(c(id_column, "unique_id")))
  } else {
    rownames(data) <- data[[id_column]]
    data <- data %>% select(-all_of(id_column))
  }
  return(data)
}

# First, identify high-value peptides for the heatmap (public neoantigens and/or fusion peptides)
highlight_peptides <- multi_omics_final %>%
  filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
  pull(Peptide)

# If there are too many peptides, limit to top ones by score
if(length(highlight_peptides) > 100) {
  highlight_peptides <- multi_omics_final %>%
    filter(potential_public_neoantigen | fusion_in_multiple_samples) %>%
    arrange(desc(public_neoantigen_score), desc(samples_present)) %>%
    head(100) %>%
    pull(Peptide)
}

# Create a presence matrix just for these highlighted peptides
if(length(highlight_peptides) > 0) {
  # Get unique peptides to prevent duplicates
  highlight_peptides <- unique(highlight_peptides)
  
  # Get the annotation data first with distinct peptides
  annotation_data <- multi_omics_final %>% 
    filter(Peptide %in% highlight_peptides) %>%
    select(Peptide, public_neoantigen_classification, from_fusion, spans_junction) %>%
    distinct(Peptide, .keep_all = TRUE)
  
  # Make sure annotation data has no NA values
  annotation_data <- annotation_data %>%
    mutate(
      public_neoantigen_classification = if_else(
        is.na(public_neoantigen_classification), 
        "Not classified", 
        public_neoantigen_classification
      ),
      from_fusion = if_else(is.na(from_fusion), FALSE, from_fusion),
      spans_junction = if_else(is.na(spans_junction), FALSE, spans_junction)
    )
  
  # Now join with the presence matrix
  highlight_presence <- presence_matrix %>%
    filter(Peptide %in% highlight_peptides) %>%
    left_join(annotation_data, by = "Peptide")
  
  # Prepare data for heatmap
  heatmap_data <- highlight_presence %>%
    select(-samples_present, -public_neoantigen_classification, -from_fusion, -spans_junction)
  
  # Use our safe function to set row names
  heatmap_data <- safe_set_rownames(heatmap_data, "Peptide")
  
  # Create annotation for the rows with matching rownames
  row_annotation <- data.frame(
    Classification = annotation_data$public_neoantigen_classification,
    From_Fusion = annotation_data$from_fusion,
    Spans_Junction = annotation_data$spans_junction,
    row.names = annotation_data$Peptide
  )
  
  # Make sure the row names in both data frames match
  common_peptides <- intersect(rownames(heatmap_data), rownames(row_annotation))
  heatmap_data <- heatmap_data[common_peptides, ]
  row_annotation <- row_annotation[common_peptides, ]
  
  # Define annotation colors explicitly
  annotation_colors <- list(
    Classification = c(
      "Tier 1A (Junction-spanning fusion peptide in 3+ samples)" = "red",
      "Tier 1B (Fusion peptide in 3+ samples)" = "orange",
      "Tier 2 (Multiple strong criteria)" = "yellow",
      "Tier 3 (Several criteria met)" = "green",
      "Tier 4 (Some criteria met)" = "blue",
      "Not a potential public neoantigen" = "gray",
      "Not classified" = "lightgray"
    ),
    From_Fusion = c("TRUE" = "purple", "FALSE" = "lightblue"),
    Spans_Junction = c("TRUE" = "brown", "FALSE" = "beige")
  )
  
  # Create the heatmap PDF
  pdf(file.path(viz_dir, "Peptide_Presence_Heatmap.pdf"), 
      width = 12, height = max(8, nrow(heatmap_data)/5))
  
  # Heatmap with annotations and explicit colors
  heatmap_plot <- pheatmap(
    heatmap_data,
    main = "Peptide Presence Across Samples",
    color = colorRampPalette(c("white", "steelblue"))(2),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    annotation_colors = annotation_colors,
    fontsize_row = 8,
    fontsize_col = 10
  )
  
  print(heatmap_plot)
  dev.off()
  
  # Create the heatmap PNG
  png(file.path(viz_dir, "Peptide_Presence_Heatmap.png"), 
      width = 1200, height = max(800, nrow(heatmap_data)*15), res = 120)
  
  print(heatmap_plot)
  dev.off()
}

#--------------------------------------------------
# PART 7.3: Fusion peptide visualizations
#--------------------------------------------------

# Visualizations for fusion peptides (if any)
if(nrow(fusion_peptides_analysis) > 0) {
  # Create barplot of fusion peptides by sample presence
  pdf(file.path(viz_dir, "Fusion_Peptides_By_Samples.pdf"), width = 10, height = 7)
  
  fusion_sample_plot <- ggplot(fusion_peptides_analysis, 
                               aes(x = samples_present, fill = fusion_peptide_type)) +
    geom_histogram(binwidth = 1, position = "stack", color = "black") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    labs(
      title = "DNAJB1-PRKACA Fusion Peptides by Sample Count",
      x = "Number of Samples Present",
      y = "Number of Peptides",
      fill = "Fusion Peptide Type"
    )
  
  print(fusion_sample_plot)
  dev.off()
  
  png(file.path(viz_dir, "Fusion_Peptides_By_Samples.png"), width = 800, height = 600, res = 100)
  print(fusion_sample_plot)
  dev.off()
  
  # Create a heatmap specifically for fusion peptides
  fusion_presence <- presence_matrix %>%
    filter(Peptide %in% fusion_peptides_analysis$Peptide) %>%
    left_join(
      fusion_peptides_analysis %>% 
        select(Peptide, fusion_peptide_type, spans_junction) %>%
        distinct(Peptide, .keep_all = TRUE),  # Ensure uniqueness
      by = "Peptide"
    )
  
  if(nrow(fusion_presence) > 0) {
    # Prepare data for heatmap
    fusion_heatmap_data <- fusion_presence %>%
      select(-samples_present, -fusion_peptide_type, -spans_junction)
    
    # Safely set row names
    fusion_heatmap_data <- safe_set_rownames(fusion_heatmap_data, "Peptide")
    
    # Create annotation for the rows with matching rownames
    fusion_row_annotation <- data.frame(
      Peptide_Type = fusion_presence$fusion_peptide_type,
      Spans_Junction = fusion_presence$spans_junction,
      row.names = rownames(fusion_heatmap_data)
    )
    
    # Create the heatmap PDF
    pdf(file.path(viz_dir, "Fusion_Peptides_Presence_Heatmap.pdf"), 
        width = 12, height = max(8, nrow(fusion_heatmap_data)/3))
    
    fusion_heatmap_plot <- pheatmap(
      fusion_heatmap_data,
      main = "DNAJB1-PRKACA Fusion Peptides Presence Across Samples",
      color = colorRampPalette(c("white", "darkred"))(2),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = fusion_row_annotation,
      fontsize_row = 10,
      fontsize_col = 10
    )
    
    print(fusion_heatmap_plot)
    dev.off()
    
    # Create the heatmap PNG
    png(file.path(viz_dir, "Fusion_Peptides_Presence_Heatmap.png"), 
        width = 1200, height = max(600, nrow(fusion_heatmap_data)*40), res = 120)
    
    print(fusion_heatmap_plot)
    dev.off()
  }
  
  # If there are fusion peptides with transcriptome data, create a plot
  fusion_with_transcriptome <- multi_omics_final %>%
    filter(from_fusion & !is.na(mean_expression))
  
  if(nrow(fusion_with_transcriptome) > 0) {
    pdf(file.path(viz_dir, "Fusion_Peptides_Transcriptome.pdf"), width = 10, height = 8)
    
    fusion_transcriptome_plot <- ggplot(fusion_with_transcriptome, 
                                        aes(x = mean_expression, 
                                            y = samples_present, 
                                            color = fusion_peptide_type,
                                            shape = spans_junction)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_text(aes(label = Peptide), hjust = -0.1, vjust = 0.2, size = 3, check_overlap = TRUE) +
      scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) +
      scale_color_brewer(palette = "Set1") +
      theme_minimal() +
      labs(
        title = "Fusion Peptides: Transcriptome Expression vs. Sample Presence",
        subtitle = "Triangles indicate junction-spanning peptides",
        x = "Mean Expression (Transcriptome)",
        y = "Number of Samples Present",
        color = "Fusion Peptide Type",
        shape = "Spans Junction"
      )
    
    print(fusion_transcriptome_plot)
    dev.off()
    
    png(file.path(viz_dir, "Fusion_Peptides_Transcriptome.png"), width = 800, height = 600, res = 100)
    print(fusion_transcriptome_plot)
    dev.off()
  }
}

#--------------------------------------------------
# PART 7.4: Public neoantigen visualizations
#--------------------------------------------------

# Create visualizations for potential public neoantigens
if(nrow(public_neoantigens) > 0) {
  # Create a summary of neoantigen tiers
  tier_summary <- public_neoantigens %>%
    group_by(public_neoantigen_classification) %>%
    summarise(
      count = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(count))
  
  # Barplot of public neoantigen tiers
  pdf(file.path(viz_dir, "Public_Neoantigen_Tiers.pdf"), width = 10, height = 6)
  
  tier_plot <- ggplot(tier_summary, 
                      aes(x = reorder(public_neoantigen_classification, -count), 
                          y = count, 
                          fill = public_neoantigen_classification)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), vjust = -0.5) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Distribution of Potential Public Neoantigens",
      x = "Classification Tier",
      y = "Count",
      fill = "Classification"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  print(tier_plot)
  dev.off()
  
  png(file.path(viz_dir, "Public_Neoantigen_Tiers.png"), width = 800, height = 600, res = 100)
  print(tier_plot)
  dev.off()
  
  # Create a scatter plot of sample presence vs transcriptome expression
  public_neo_scatter <- public_neoantigens %>%
    filter(!is.na(mean_expression))
  
  if(nrow(public_neo_scatter) > 0) {
    pdf(file.path(viz_dir, "Public_Neoantigens_Scatter.pdf"), width = 10, height = 8)
    
    neo_scatter_plot <- ggplot(public_neo_scatter, 
                               aes(x = mean_expression, 
                                   y = samples_present, 
                                   color = public_neoantigen_classification,
                                   shape = from_fusion)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_text(aes(label = Peptide), hjust = -0.1, vjust = 0.2, size = 3, check_overlap = TRUE) +
      scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) +
      theme_minimal() +
      labs(
        title = "Potential Public Neoantigens: Transcriptome Expression vs. Sample Presence",
        subtitle = "Triangles indicate fusion-derived peptides",
        x = "Mean Expression (Transcriptome)",
        y = "Number of Samples Present",
        color = "Neoantigen Classification",
        shape = "From Fusion"
      )
    
    print(neo_scatter_plot)
    dev.off()
    
    png(file.path(viz_dir, "Public_Neoantigens_Scatter.png"), width = 800, height = 600, res = 100)
    print(neo_scatter_plot)
    dev.off()
  }
  
  # Create a heatmap of public neoantigens across different omics
  # First, create a matrix of key metrics
  neo_heatmap_data <- public_neoantigens %>%
    select(
      Peptide,
      samples_present,
      mean_expression,
      log2_fold_change_lfq,
      log2_fold_change_tmt
    )
  
  # Handle potentially duplicate peptides
  neo_heatmap_data <- safe_set_rownames(neo_heatmap_data, "Peptide")
  neo_heatmap_matrix <- as.matrix(neo_heatmap_data)
  
  # Replace NA with 0 for visualization purposes
  neo_heatmap_matrix[is.na(neo_heatmap_matrix)] <- 0
  
  # Scale the columns to be comparable
  neo_heatmap_matrix_scaled <- scale(neo_heatmap_matrix)
  
  # Create annotation for the rows - ensuring matching rownames
  row_info <- public_neoantigens %>%
    select(Peptide, public_neoantigen_classification, from_fusion) %>%
    distinct(Peptide, .keep_all = TRUE)
  
  # Create a mapping from original peptides to unique rownames if needed
  if(exists("unique_id", where = neo_heatmap_data)) {
    peptide_to_rowname <- setNames(rownames(neo_heatmap_data), neo_heatmap_data$Peptide)
    row_info$unique_id <- peptide_to_rowname[row_info$Peptide]
    neo_row_annotation <- data.frame(
      Classification = row_info$public_neoantigen_classification,
      From_Fusion = row_info$from_fusion,
      row.names = row_info$unique_id
    )
  } else {
    neo_row_annotation <- data.frame(
      Classification = row_info$public_neoantigen_classification,
      From_Fusion = row_info$from_fusion,
      row.names = row_info$Peptide
    )
  }
  
  # Create the heatmap PDF
  pdf(file.path(viz_dir, "Public_Neoantigens_Multi_Omics_Heatmap.pdf"), 
      width = 12, height = max(8, nrow(neo_heatmap_matrix)/5))
  
  neo_heatmap_plot <- pheatmap(
    neo_heatmap_matrix_scaled,
    main = "Potential Public Neoantigens: Multi-Omics Profile",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-3, 3, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_row = neo_row_annotation,
    display_numbers = FALSE,
    fontsize_row = 8,
    fontsize_col = 10,
    labels_col = c("Samples Present", "Transcriptome Expression", "LFQ log2FC", "TMT log2FC")
  )
  
  print(neo_heatmap_plot)
  dev.off()
  
  # Create the heatmap PNG
  png(file.path(viz_dir, "Public_Neoantigens_Multi_Omics_Heatmap.png"), 
      width = 1200, height = max(800, nrow(neo_heatmap_matrix)*15), res = 120)
  
  print(neo_heatmap_plot)
  dev.off()

#--------------------------------------------------
# PART 7.5: General summary visualizations
#--------------------------------------------------

# Create a histogram of sample counts
pdf(file.path(viz_dir, "Peptide_Sample_Count_Histogram.pdf"), width = 10, height = 6)

sample_hist <- ggplot(presence_categories, aes(x = samples_present, fill = presence_category)) +
  geom_histogram(binwidth = 1, color = "black") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(
    title = "Distribution of Peptides by Sample Count",
    x = "Number of Samples Present",
    y = "Number of Peptides",
    fill = "Presence Category"
  )

print(sample_hist)
dev.off()

png(file.path(viz_dir, "Peptide_Sample_Count_Histogram.png"), width = 800, height = 600, res = 100)
print(sample_hist)
dev.off()

# Create a bar chart of peptide presence categories
pdf(file.path(viz_dir, "Peptide_Presence_Categories.pdf"), width = 10, height = 6)

presence_cat_counts <- presence_categories %>%
  group_by(presence_category) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

presence_cat_plot <- ggplot(presence_cat_counts, 
                            aes(x = reorder(presence_category, -count), 
                                y = count, 
                                fill = presence_category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = -0.5) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Distribution of Peptides by Presence Category",
    x = "Presence Category",
    y = "Number of Peptides",
    fill = "Category"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(presence_cat_plot)
dev.off()

png(file.path(viz_dir, "Peptide_Presence_Categories.png"), width = 800, height = 600, res = 100)
print(presence_cat_plot)
dev.off()

# Create a Venn diagram of genes in different omics datasets
# First, extract gene lists
genes_immuno <- multi_omics_final %>%
  filter(samples_present >= 5) %>%
  pull(primary_gene) %>%
  unique()

genes_transcriptome <- transcriptome_filtered %>%
  filter(mean_expression > median(transcriptome_filtered$mean_expression, na.rm = TRUE)) %>%
  pull(symbol) %>%
  unique()

genes_lfq <- lfq_processed %>%
  filter(log2_fold_change_lfq > 1) %>%
  pull(Gene_Name) %>%
  unique()

genes_tmt <- tmt_processed %>%
  filter(log2_fold_change_tmt > 1) %>%
  pull(Gene_Name) %>%
  unique()

# Create Venn diagram if all gene lists have entries
if(length(genes_immuno) > 0 && length(genes_transcriptome) > 0 && 
   length(genes_lfq) > 0 && length(genes_tmt) > 0) {
  
  venn_colors <- brewer.pal(4, "Set1")
  venn_output <- file.path(viz_dir, "Multi_Omics_Genes_Venn.png")
  
  png(venn_output, width = 800, height = 800, res = 120)
  
  venn_plot <- draw.quad.venn(
    area1 = length(genes_immuno),
    area2 = length(genes_transcriptome),
    area3 = length(genes_lfq),
    area4 = length(genes_tmt),
    n12 = length(intersect(genes_immuno, genes_transcriptome)),
    n13 = length(intersect(genes_immuno, genes_lfq)),
    n14 = length(intersect(genes_immuno, genes_tmt)),
    n23 = length(intersect(genes_transcriptome, genes_lfq)),
    n24 = length(intersect(genes_transcriptome, genes_tmt)),
    n34 = length(intersect(genes_lfq, genes_tmt)),
    n123 = length(intersect(intersect(genes_immuno, genes_transcriptome), genes_lfq)),
    n124 = length(intersect(intersect(genes_immuno, genes_transcriptome), genes_tmt)),
    n134 = length(intersect(intersect(genes_immuno, genes_lfq), genes_tmt)),
    n234 = length(intersect(intersect(genes_transcriptome, genes_lfq), genes_tmt)),
    n1234 = length(intersect(intersect(intersect(genes_immuno, genes_transcriptome), genes_lfq), genes_tmt)),
    category = c("Immunopeptidome (5+ samples)", "Transcriptome (High Expr)", "LFQ Proteome (Upregulated)", "TMT Proteome (Upregulated)"),
    fill = venn_colors,
    alpha = 0.5,
    lty = "blank",
    cex = 1.5,
    cat.cex = 1.2,
    cat.col = venn_colors,
    cat.dist = 0.09,
    cat.pos = c(0, 0, 0, 0),
    euler.d = TRUE,
    scaled = TRUE
  )
  
  dev.off()
}

#--------------------------------------------------
# PART 8: Save results to Excel file
#--------------------------------------------------

# Create Excel output with multiple sheets
excel_sheets <- list(
  "Multi_Omics_All_Peptides" = multi_omics_final,
  "Peptides_in_5plus_Samples" = multi_omics_final %>% filter(in_5plus_samples),
  "Peptides_in_6plus_Samples" = multi_omics_final %>% filter(in_6plus_samples),
  "Peptides_in_7plus_Samples" = multi_omics_final %>% filter(in_7plus_samples),
  "Peptides_in_8plus_Samples" = multi_omics_final %>% filter(in_8plus_samples),
  "Peptides_in_9plus_Samples" = multi_omics_final %>% filter(in_9plus_samples),
  "Peptides_in_All_Samples" = multi_omics_final %>% filter(samples_present == length(sample_ids))
)

# Add public neoantigens sheet
if(nrow(public_neoantigens) > 0) {
  excel_sheets[["Public_Neoantigens"]] = public_neoantigens
  
  # Add tier-specific sheets
  tier_groups <- split(public_neoantigens, public_neoantigens$public_neoantigen_classification)
  for(tier in names(tier_groups)) {
    sheet_name <- paste0("Neoantigen_", gsub("[^a-zA-Z0-9]", "_", substr(tier, 1, 10)))
    excel_sheets[[sheet_name]] <- tier_groups[[tier]]
  }
}

# Add fusion peptide sheets
if(nrow(fusion_peptides_analysis) > 0) {
  excel_sheets[["Fusion_Peptides"]] <- fusion_peptides_analysis
  
  # Add junction-spanning fusion peptides
  junction_peptides <- fusion_peptides_analysis %>% filter(spans_junction)
  if(nrow(junction_peptides) > 0) {
    excel_sheets[["Junction_Spanning_Peptides"]] <- junction_peptides
  }
  
  # Add fusion public neoantigens
  if(nrow(fusion_public_neoantigens) > 0) {
    excel_sheets[["Fusion_Public_Neoantigens"]] <- fusion_public_neoantigens
  }
}

# Create a more formatted Excel workbook
wb <- createWorkbook()

# Add sheets with formatting
for (sheet_name in names(excel_sheets)) {
  # Add a worksheet
  addWorksheet(wb, sheet_name)
  
  # Write data
  writeData(wb, sheet_name, excel_sheets[[sheet_name]], headerStyle = createStyle(textDecoration = "bold"))
  
  # Auto-adjust column widths
  setColWidths(wb, sheet_name, cols = 1:ncol(excel_sheets[[sheet_name]]), widths = "auto")
  
  # Freeze the header row
  freezePane(wb, sheet_name, firstRow = TRUE)
  
  # Add conditional formatting for public neoantigens sheets
  if(grepl("Neoantigen", sheet_name)) {
    # Highlight by public neoantigen score
    conditionalFormatting(wb, sheet_name, 
                          cols = which(colnames(excel_sheets[[sheet_name]]) == "public_neoantigen_score"), 
                          rows = 2:(nrow(excel_sheets[[sheet_name]]) + 1), 
                          rule = ">5", 
                          style = createStyle(bgFill = "#E2EFDA")) # Light green
    
    # Highlight by sample presence
    conditionalFormatting(wb, sheet_name, 
                          cols = which(colnames(excel_sheets[[sheet_name]]) == "samples_present"), 
                          rows = 2:(nrow(excel_sheets[[sheet_name]]) + 1), 
                          rule = ">7", 
                          style = createStyle(bgFill = "#FCE4D6")) # Light orange
  }
  
  # Add conditional formatting for fusion peptides
  if(grepl("Fusion|Junction", sheet_name)) {
    # Highlight spanning junction
    conditionalFormatting(wb, sheet_name, 
                          cols = which(colnames(excel_sheets[[sheet_name]]) == "spans_junction"), 
                          rows = 2:(nrow(excel_sheets[[sheet_name]]) + 1), 
                          rule = "==TRUE", 
                          style = createStyle(bgFill = "#FFCCCC")) # Light red
    
    # Highlight by sample presence
    conditionalFormatting(wb, sheet_name, 
                          cols = which(colnames(excel_sheets[[sheet_name]]) == "samples_present"), 
                          rows = 2:(nrow(excel_sheets[[sheet_name]]) + 1), 
                          rule = ">5", 
                          style = createStyle(bgFill = "#FCE4D6")) # Light orange
  }
}

# Save the Excel file
output_file <- file.path(output_dir, "MultiSample_Omics_Comparison.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

#--------------------------------------------------
# PART 9: Print summary information
#--------------------------------------------------

cat("\n----------------------------------------\n")
cat("SUMMARY OF MULTI-SAMPLE OMICS ANALYSIS\n")
cat("----------------------------------------\n\n")

cat("Total peptides analyzed:", nrow(presence_matrix), "\n")
cat("Peptides present in 5+ samples:", nrow(peptides_in_5plus), "\n")
cat("Peptides present in 6+ samples:", nrow(peptides_in_6plus), "\n")
cat("Peptides present in 7+ samples:", nrow(peptides_in_7plus), "\n")
cat("Peptides present in 8+ samples:", nrow(peptides_in_8plus), "\n")
cat("Peptides present in 9+ samples:", nrow(peptides_in_9plus), "\n")
cat("Peptides present in all samples:", sum(presence_matrix$samples_present == length(sample_ids)), "\n")

if(exists("transcriptome_filtered")) {
  cat("\nTranscriptome data summary:\n")
  cat("Total genes in transcriptome:", nrow(transcriptome_filtered), "\n")
}

if(exists("lfq_processed")) {
  cat("\nLFQ proteome data summary:\n")
  cat("Total genes in LFQ proteome:", nrow(lfq_processed), "\n")
  cat("Upregulated genes (log2FC > 1):", sum(lfq_processed$log2_fold_change_lfq > 1, na.rm = TRUE), "\n")
}

if(exists("tmt_processed")) {
  cat("\nTMT proteome data summary:\n")
  cat("Total genes in TMT proteome:", nrow(tmt_processed), "\n")
  cat("Upregulated genes (log2FC > 1):", sum(tmt_processed$log2_fold_change_tmt > 1, na.rm = TRUE), "\n")
}

# Print public neoantigen summary
cat("\nPUBLIC NEOANTIGEN SUMMARY:\n")
cat("Total potential public neoantigens found:", nrow(public_neoantigens), "\n")

if(nrow(public_neoantigens) > 0) {
  tier_counts <- public_neoantigens %>%
    group_by(public_neoantigen_classification) %>%
    summarise(count = n(), .groups = "drop")
  
  print(tier_counts)
  
  cat("\nTop public neoantigens (by score):\n")
  top_neoantigens <- public_neoantigens %>%
    arrange(desc(public_neoantigen_score), desc(samples_present)) %>%
    head(10) %>%
    select(Peptide, primary_gene, public_neoantigen_classification, 
           samples_present, mean_expression, 
           log2_fold_change_lfq, log2_fold_change_tmt)
  
  print(top_neoantigens)
}

# Print fusion peptide information if any found
if(nrow(fusion_peptides_analysis) > 0) {
  cat("\nFUSION PEPTIDE SUMMARY:\n")
  cat("Fusion peptides found:", nrow(fusion_peptides_analysis), "\n")
  cat("Junction-spanning peptides:", sum(fusion_peptides_analysis$spans_junction), "\n")
  
  cat("\nFusion peptide presence summary:\n")
  fusion_presence <- fusion_peptides_analysis %>%
    group_by(samples_present, fusion_peptide_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(samples_present), fusion_peptide_type)
  
  print(fusion_presence)
  
  if(nrow(fusion_public_neoantigens) > 0) {
    cat("\nFusion-derived public neoantigens:", nrow(fusion_public_neoantigens), "\n")
    print(fusion_public_neoantigens %>% 
            select(Peptide, fusion_peptide_type, spans_junction, 
                   public_neoantigen_classification, samples_present))
  } else {
    cat("\nNo fusion-derived public neoantigens identified.\n")
  }
} else {
  cat("\nNo fusion peptides found in the analysis.\n")
}

# Print multi-omics overlap summary
if(length(genes_immuno) > 0 && length(genes_transcriptome) > 0 && 
   length(genes_lfq) > 0 && length(genes_tmt) > 0) {
  
  cat("\nMULTI-OMICS OVERLAP SUMMARY:\n")
  cat("Genes from immunopeptidome (5+ samples):", length(genes_immuno), "\n")
  cat("Genes with high expression in transcriptome:", length(genes_transcriptome), "\n")
  cat("Upregulated genes in LFQ proteome:", length(genes_lfq), "\n")
  cat("Upregulated genes in TMT proteome:", length(genes_tmt), "\n")
  
  cat("\nOverlap statistics:\n")
  cat("Immunopeptidome + Transcriptome:", length(intersect(genes_immuno, genes_transcriptome)), "\n")
  cat("Immunopeptidome + LFQ proteome:", length(intersect(genes_immuno, genes_lfq)), "\n")
  cat("Immunopeptidome + TMT proteome:", length(intersect(genes_immuno, genes_tmt)), "\n")
  cat("Common in all 4 datasets:", 
      length(intersect(intersect(intersect(genes_immuno, genes_transcriptome), genes_lfq), genes_tmt)), "\n")
}

cat("\nAnalysis complete! Results saved to:", output_dir, "\n")
cat("Main output file:", output_file, "\n")
cat("Visualizations saved to:", viz_dir, "\n")