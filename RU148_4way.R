# Enhanced script to analyze immunopeptidome data and compare with 
# transcriptome, LFQ proteome, and TMT proteome for identifying public neoantigens

# Setting directory
setwd("~/Documents/Github/HLA-I_Analysis/")

# Create output directory in HLA-I_Analysis
output_dir <- "RU148_4wayanalysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Created output directory:", output_dir, "\n")
}

# Create a directory for visualizations
viz_dir <- file.path(output_dir, "RU148_4way_visualizations")
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
data_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/imp_014_rawdata"
transcriptome_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/Normalized_Gene_counts_FLCdb_Panel_1.xlsx"
lfq_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/Levin2023/adg7038_Table_S2_LFQ.xlsx"
tmt_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/Levin2023/adg7038_Table_S1_TMT.xlsx"

#--------------------------------------------------
# PART 1: Process the immunopeptidome data
#--------------------------------------------------

# Read all TSV files in the directory and extract sample IDs from filenames
files <- list.files(path = data_path, pattern = ".*_peptides\\.tsv$|.*_untargeted_peptide\\.tsv$", full.names = TRUE)

# Filter for only 148T and 148N files
tumor_normal_files <- files[grepl("148[TN]", files)]

if (length(tumor_normal_files) == 0) {
  stop("No 148T or 148N peptide files found in ", data_path)
}

# Split files into 2CV and 3CV categories
files_2cv <- tumor_normal_files[grepl("2CV", tumor_normal_files)]
files_3cv <- tumor_normal_files[grepl("3CV", tumor_normal_files)]

cat("Found", length(files_2cv), "2CV files and", length(files_3cv), "3CV files for analysis\n")

# Function to read and process immunopeptidome files
process_immunopeptidome_files <- function(file_list) {
  all_data <- list()
  
  for (file in file_list) {
    filename <- basename(file)
    
    # Extract sample ID from filename (148T or 148N)
    sample_id <- ifelse(grepl("148T", filename), "148T", "148N")
    
    cat("Reading file:", filename, "- Sample ID:", sample_id, "\n")
    
    # Read the file
    data <- read.delim(file, stringsAsFactors = FALSE)
    
    # Add a column for sample ID
    data$Sample_ID <- sample_id
    data$Filename <- filename
    
    # Add to our list
    all_data[[filename]] <- data
  }
  
  return(bind_rows(all_data))
}

# Process 2CV and 3CV files separately
data_2cv <- process_immunopeptidome_files(files_2cv)
data_3cv <- process_immunopeptidome_files(files_3cv)

# Combine 2CV and 3CV data
combined_data <- bind_rows(
  data_2cv %>% mutate(CV_type = "2CV"),
  data_3cv %>% mutate(CV_type = "3CV")
)

# Create a summary of peptide detection for tumor and normal
peptide_summary <- combined_data %>%
  # Add column for peptide length
  mutate(peptide_length = nchar(Peptide)) %>%
  # FILTER 2: Keep only peptides with lengths 8-12 amino acids
  filter(peptide_length >= 8 & peptide_length <= 12) %>%
  group_by(Sample_ID, Peptide) %>%
  summarize(
    peptide_length = first(peptide_length),
    spectral_count = sum(Spectral.Count),
    total_intensity = sum(Intensity),
    protein_ids = paste(unique(Protein.ID), collapse = "; "),
    genes = paste(unique(Gene), collapse = "; "),
    source_filenames = paste(unique(Filename), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Peptide, Sample_ID)

# Create a wide format table with tumor and normal side by side
peptide_comparison <- peptide_summary %>%
  select(Sample_ID, Peptide, peptide_length, spectral_count, total_intensity, protein_ids, genes) %>%
  pivot_wider(
    names_from = Sample_ID,
    values_from = c(spectral_count, total_intensity, protein_ids, genes),
    values_fill = list(spectral_count = 0, total_intensity = 0, protein_ids = NA, genes = NA)
  )

# Calculate fold changes and identify tumor-specific and normal-specific peptides
immunopeptidome_analysis <- peptide_comparison %>%
  mutate(
    # Replace zero with small value to prevent division by zero or Inf
    total_intensity_148N_adj = ifelse(total_intensity_148N == 0, 0.1, total_intensity_148N),
    total_intensity_148T_adj = ifelse(total_intensity_148T == 0, 0.1, total_intensity_148T),
    
    # Calculate fold changes (log2)
    log2_fold_change_immuno = log2(total_intensity_148T_adj / total_intensity_148N_adj),
    
    # Determine if peptide is specific to tumor or normal
    detection_status = case_when(
      total_intensity_148T > 0 & total_intensity_148N == 0 ~ "Tumor-specific",
      total_intensity_148N > 0 & total_intensity_148T == 0 ~ "Normal-specific",
      total_intensity_148T > 0 & total_intensity_148N > 0 ~ "Detected in both",
      TRUE ~ "Not detected"
    ),
    
    # Add peptide length
    peptide_length = nchar(Peptide),
    
    # Simplified category for plotting
    peptide_category = case_when(
      log2_fold_change_immuno > 1 ~ "Up in Tumor (FC > 2)",
      log2_fold_change_immuno < -1 ~ "Down in Tumor (FC < 0.5)",
      TRUE ~ "Similar (-1 < log2FC < 1)"
    )
  ) %>%
  # Clean up protein and gene info
  mutate(
    genes_combined = coalesce(genes_148T, genes_148N),
    proteins_combined = coalesce(protein_ids_148T, protein_ids_148N)
  ) %>%
  # Extract primary gene for later comparison with transcriptome
  mutate(
    primary_gene = sapply(strsplit(genes_combined, ";\\s*"), function(x) trimws(x[1]))
  ) %>%
  # Sort by fold change for easier viewing
  arrange(desc(log2_fold_change_immuno))

#--------------------------------------------------
# PART 2: Define the fusion protein sequence and detection function
#--------------------------------------------------

# Define the fusion protein sequence - UPDATED based on provided sequences
fusion_protein <- "RKREIFDRYGEEVKEFLAKAKEDF"
# Identify the fusion junction - between E (DNAJB1) and V (PRKACA)
dnajb1_part <- "RKREIFDRYGEE"  # Ends with YGEE (full DNAJB1 part)
prkaca_part <- "VKEFLAKAKEDF"  # Starts with VKEF (full PRKACA part)
junction_position <- nchar(dnajb1_part)  # Position after the last E in YGEE

# Print info about the fusion protein for verification
cat("DNAJB1 part:", dnajb1_part, "\n")
cat("PRKACA part:", prkaca_part, "\n")
cat("Junction position:", junction_position, "\n")
cat("Fusion protein:", fusion_protein, "\n")
cat("Checking junction: character at junction position is", substr(fusion_protein, junction_position, junction_position), "\n")
cat("Checking junction: next character is", substr(fusion_protein, junction_position+1, junction_position+1), "\n\n")

# Modified function to check if a peptide spans the fusion junction with the updated sequences
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

# Add fusion protein information to the immunopeptidome analysis
immunopeptidome_analysis <- immunopeptidome_analysis %>%
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
fusion_peptides_analysis <- immunopeptidome_analysis %>%
  filter(from_fusion) %>%
  arrange(desc(spans_junction), desc(log2_fold_change_immuno))

#--------------------------------------------------
# PART 2.5: Identify peptides exclusive to tumor or normal
#--------------------------------------------------

# Find peptides that are exclusively in tumor (not detected in normal)
tumor_exclusive_peptides <- immunopeptidome_analysis %>%
  filter(total_intensity_148T > 0 & total_intensity_148N == 0) %>%
  arrange(desc(total_intensity_148T))

# Find peptides that are exclusively in normal (not detected in tumor)
normal_exclusive_peptides <- immunopeptidome_analysis %>%
  filter(total_intensity_148N > 0 & total_intensity_148T == 0) %>%
  arrange(desc(total_intensity_148N))

# Count of exclusive peptides 
cat("\nExclusive peptide counts:\n")
cat("Peptides found only in tumor:", nrow(tumor_exclusive_peptides), "\n")
cat("Peptides found only in normal:", nrow(normal_exclusive_peptides), "\n")

# Add fusion protein analysis for exclusive peptides
tumor_exclusive_fusion <- tumor_exclusive_peptides %>%
  filter(from_fusion) %>%
  arrange(desc(spans_junction), desc(total_intensity_148T))

normal_exclusive_fusion <- normal_exclusive_peptides %>%
  filter(from_fusion) %>%
  arrange(desc(spans_junction), desc(total_intensity_148N))

#--------------------------------------------------
# PART 3: Process the transcriptome data
#--------------------------------------------------

# First, examine the structure of the transcriptome data
cat("Reading transcriptome data from:", transcriptome_path, "\n")
transcriptome_data <- read_excel(transcriptome_path)

# Print column names to determine what's available
cat("Transcriptome data columns:", paste(colnames(transcriptome_data), collapse=", "), "\n")

# Process transcriptome data based on actual column names
# We'll be more flexible with column naming and existence
transcriptome_processed <- transcriptome_data %>%
  # Create default placeholder columns if they don't exist
  mutate(
    Mean.Normal = NA_real_,
    Mean.Tumor = NA_real_
  )

# Check if RU148 columns exist and compute averages if they do
if(all(c("RU148_T8", "RU148_T11") %in% colnames(transcriptome_data))) {
  transcriptome_processed <- transcriptome_processed %>%
    mutate(
      RU148_T_Average = (RU148_T8 + RU148_T11) / 2
    )
} else {
  # If columns don't exist, create placeholder
  transcriptome_processed$RU148_T_Average <- NA_real_
  cat("Warning: RU148_T8 and/or RU148_T11 columns not found in transcriptome data\n")
}

# Calculate log2 fold change if possible
if(all(c("RU148_N", "RU148_T_Average") %in% colnames(transcriptome_processed)) && 
   !all(is.na(transcriptome_processed$RU148_N)) && 
   !all(is.na(transcriptome_processed$RU148_T_Average))) {
  
  transcriptome_processed <- transcriptome_processed %>%
    mutate(
      log2_fold_change_transcriptome = log2(
        ifelse(RU148_T_Average == 0, 0.1, RU148_T_Average) / 
          ifelse(RU148_N == 0, 0.1, RU148_N)
      )
    )
} else {
  transcriptome_processed$log2_fold_change_transcriptome <- NA_real_
  cat("Warning: Unable to calculate transcriptome log2 fold change\n")
}

# Ensure we have a symbol column for joining
if("symbol" %in% colnames(transcriptome_processed)) {
  cat("Using 'symbol' column for joining\n")
} else if("gene_symbol" %in% colnames(transcriptome_processed)) {
  transcriptome_processed <- transcriptome_processed %>%
    rename(symbol = gene_symbol)
  cat("Renamed 'gene_symbol' to 'symbol' for joining\n")
} else if("Symbol" %in% colnames(transcriptome_processed)) {
  transcriptome_processed <- transcriptome_processed %>%
    rename(symbol = Symbol)
  cat("Renamed 'Symbol' to 'symbol' for joining\n")
} else {
  cat("Warning: No suitable symbol column found for joining\n")
  # Create an empty symbol column to avoid join errors
  transcriptome_processed$symbol <- NA_character_
}

#--------------------------------------------------
# PART 3.1: Process the LFQ proteome data
#--------------------------------------------------

# Read the LFQ data
cat("Reading LFQ proteome data from:", lfq_path, "\n")
lfq_data <- read_excel(lfq_path, sheet = "Significant and 1.5x_2")

# Print column names to confirm structure
cat("LFQ proteome data columns:", paste(colnames(lfq_data), collapse=", "), "\n")

# Process LFQ data
lfq_processed <- lfq_data %>%
  # Ensure consistent column names
  rename_with(~ gsub(" ", "_", .), everything()) %>%
  # Filter for significant entries (p-value <= 0.05)
  filter(P.value <= 0.05) %>%
  # Clean up and standardize
  mutate(
    Gene_Name = trimws(Gene_Name),
    log2_fold_change_lfq = Log2_Difference,
    p_value_lfq = P.value,
    protein_category_lfq = case_when(
      Log2_Difference > 1 ~ "Up in Tumor (FC > 2)",
      Log2_Difference < -1 ~ "Down in Tumor (FC < 0.5)",
      TRUE ~ "Similar (-1 < log2FC < 1)"
    )
  ) %>%
  # Select relevant columns
  select(Gene_Name, log2_fold_change_lfq, p_value_lfq, Protein_Name, protein_category_lfq)

cat("Processed", nrow(lfq_processed), "significant entries from LFQ proteome data\n")

#--------------------------------------------------
# PART 3.2: Process the TMT proteome data
#--------------------------------------------------

# Read the TMT data
cat("Reading TMT proteome data from:", tmt_path, "\n")
tmt_data <- read_excel(tmt_path, sheet = "Significant and 1.5x_2")

# Print column names to confirm structure
cat("TMT proteome data columns:", paste(colnames(tmt_data), collapse=", "), "\n")

# Process TMT data
tmt_processed <- tmt_data %>%
  # Ensure consistent column names
  rename_with(~ gsub(" ", "_", .), everything()) %>%
  # Filter for significant entries (p-value <= 0.05)
  filter(P.value <= 0.05) %>%
  # Clean up and standardize
  mutate(
    Gene_Name = trimws(Gene_Name),
    log2_fold_change_tmt = Log2_Difference,
    p_value_tmt = P.value,
    protein_category_tmt = case_when(
      Log2_Difference > 1 ~ "Up in Tumor (FC > 2)",
      Log2_Difference < -1 ~ "Down in Tumor (FC < 0.5)",
      TRUE ~ "Similar (-1 < log2FC < 1)"
    )
  ) %>%
  # Select relevant columns
  select(Gene_Name, log2_fold_change_tmt, p_value_tmt, Protein_Name, protein_category_tmt)

cat("Processed", nrow(tmt_processed), "significant entries from TMT proteome data\n")

#--------------------------------------------------
# PART 4: Combine immunopeptidome and transcriptome data
#--------------------------------------------------

# Join the immunopeptidome and transcriptome data based on gene symbol
# Address the many-to-many relationship warning by explicitly setting the relationship
combined_analysis <- immunopeptidome_analysis %>%
  left_join(
    transcriptome_processed,
    by = c("primary_gene" = "symbol"),
    relationship = "many-to-many"  # Explicitly acknowledge many-to-many relationship
  )

# Add comparison metrics if log2 fold changes are available
if("log2_fold_change_transcriptome" %in% colnames(combined_analysis) && 
   !all(is.na(combined_analysis$log2_fold_change_transcriptome))) {
  
  combined_analysis <- combined_analysis %>%
    mutate(
      immuno_trans_correlation = log2_fold_change_immuno * log2_fold_change_transcriptome,
      expression_category_transcriptome = case_when(
        is.na(log2_fold_change_transcriptome) ~ "No transcriptome data",
        log2_fold_change_transcriptome > 1 & log2_fold_change_immuno > 1 ~ "Up in both",
        log2_fold_change_transcriptome < -1 & log2_fold_change_immuno < -1 ~ "Down in both",
        log2_fold_change_transcriptome > 1 & log2_fold_change_immuno < -1 ~ "Up in transcriptome, down in immunopeptidome",
        log2_fold_change_transcriptome < -1 & log2_fold_change_immuno > 1 ~ "Down in transcriptome, up in immunopeptidome",
        TRUE ~ "No significant change"
      )
    )
} else {
  combined_analysis$immuno_trans_correlation <- NA_real_
  combined_analysis$expression_category_transcriptome <- "No transcriptome data"
}

#--------------------------------------------------
# PART 4.1: Combine with LFQ proteome data
#--------------------------------------------------

# Join the combined analysis with LFQ data
combined_analysis <- combined_analysis %>%
  left_join(
    lfq_processed,
    by = c("primary_gene" = "Gene_Name"),
    relationship = "many-to-many"
  )

# Add comparison metrics for LFQ
combined_analysis <- combined_analysis %>%
  mutate(
    immuno_lfq_correlation = case_when(
      !is.na(log2_fold_change_lfq) ~ log2_fold_change_immuno * log2_fold_change_lfq,
      TRUE ~ NA_real_
    ),
    expression_category_lfq = case_when(
      is.na(log2_fold_change_lfq) ~ "No LFQ proteome data",
      log2_fold_change_lfq > 1 & log2_fold_change_immuno > 1 ~ "Up in both",
      log2_fold_change_lfq < -1 & log2_fold_change_immuno < -1 ~ "Down in both",
      log2_fold_change_lfq > 1 & log2_fold_change_immuno < -1 ~ "Up in LFQ, down in immunopeptidome",
      log2_fold_change_lfq < -1 & log2_fold_change_immuno > 1 ~ "Down in LFQ, up in immunopeptidome",
      TRUE ~ "No significant change"
    )
  )

#--------------------------------------------------
# PART 4.2: Combine with TMT proteome data
#--------------------------------------------------

# Join the combined analysis with TMT data
combined_analysis <- combined_analysis %>%
  left_join(
    tmt_processed,
    by = c("primary_gene" = "Gene_Name"),
    relationship = "many-to-many"
  )

# Add comparison metrics for TMT
combined_analysis <- combined_analysis %>%
  mutate(
    immuno_tmt_correlation = case_when(
      !is.na(log2_fold_change_tmt) ~ log2_fold_change_immuno * log2_fold_change_tmt,
      TRUE ~ NA_real_
    ),
    expression_category_tmt = case_when(
      is.na(log2_fold_change_tmt) ~ "No TMT proteome data",
      log2_fold_change_tmt > 1 & log2_fold_change_immuno > 1 ~ "Up in both",
      log2_fold_change_tmt < -1 & log2_fold_change_immuno < -1 ~ "Down in both",
      log2_fold_change_tmt > 1 & log2_fold_change_immuno < -1 ~ "Up in TMT, down in immunopeptidome",
      log2_fold_change_tmt < -1 & log2_fold_change_immuno > 1 ~ "Down in TMT, up in immunopeptidome",
      TRUE ~ "No significant change"
    )
  )

#--------------------------------------------------
# PART 4.3: Create a 4-way comparison for potential public neoantigens
#--------------------------------------------------

# Define public neoantigen criteria
combined_analysis <- combined_analysis %>%
  mutate(
    # Criteria 1: Upregulated in all 4 datasets (log2FC > 1)
    upregulated_in_all = case_when(
      log2_fold_change_immuno > 1 & 
        !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > 1 &
        !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > 1 &
        !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > 1 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Criteria 2: Upregulated in immunopeptidome and at least 2 other datasets
    upregulated_in_3_datasets = case_when(
      log2_fold_change_immuno > 1 & 
        sum(
          !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > 1,
          !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > 1,
          !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > 1
        ) >= 2 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Criteria 3: Tumor-specific in immunopeptidome and upregulated in at least 2 other datasets
    tumor_specific_upregulated = case_when(
      detection_status == "Tumor-specific" & 
        sum(
          !is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > 1,
          !is.na(log2_fold_change_lfq) & log2_fold_change_lfq > 1,
          !is.na(log2_fold_change_tmt) & log2_fold_change_tmt > 1
        ) >= 2 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Combined public neoantigen score (sum of criteria)
    public_neoantigen_score = as.integer(upregulated_in_all) * 3 + 
      as.integer(upregulated_in_3_datasets) * 2 + 
      as.integer(tumor_specific_upregulated) * 1,
    
    # Classify as potential public neoantigen if any criteria are met
    potential_public_neoantigen = case_when(
      public_neoantigen_score > 0 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Classification label
    public_neoantigen_classification = case_when(
      upregulated_in_all ~ "Tier 1 (Upregulated in all 4 datasets)",
      upregulated_in_3_datasets ~ "Tier 2 (Upregulated in immunopeptidome + 2 others)",
      tumor_specific_upregulated ~ "Tier 3 (Tumor-specific + upregulated in 2 others)",
      TRUE ~ "Not a potential public neoantigen"
    )
  )

# Create a filtered dataset for potential public neoantigens
public_neoantigens <- combined_analysis %>%
  filter(potential_public_neoantigen) %>%
  arrange(desc(public_neoantigen_score), desc(log2_fold_change_immuno)) %>%
  select(
    Peptide, 
    primary_gene, 
    genes_combined,
    peptide_length,
    public_neoantigen_classification,
    public_neoantigen_score,
    log2_fold_change_immuno,
    log2_fold_change_transcriptome,
    log2_fold_change_lfq,
    log2_fold_change_tmt,
    detection_status,
    from_fusion,
    spans_junction
  )

# Count the number of potential public neoantigens
cat("\nPotential public neoantigens found:", nrow(public_neoantigens), "\n")
cat("Tier 1 (Upregulated in all 4 datasets):", sum(public_neoantigens$public_neoantigen_classification == "Tier 1 (Upregulated in all 4 datasets)"), "\n")
cat("Tier 2 (Upregulated in immunopeptidome + 2 others):", sum(public_neoantigens$public_neoantigen_classification == "Tier 2 (Upregulated in immunopeptidome + 2 others)"), "\n")
cat("Tier 3 (Tumor-specific + upregulated in 2 others):", sum(public_neoantigens$public_neoantigen_classification == "Tier 3 (Tumor-specific + upregulated in 2 others)"), "\n")

# Check for fusion-derived public neoantigens
fusion_public_neoantigens <- public_neoantigens %>%
  filter(from_fusion)

if(nrow(fusion_public_neoantigens) > 0) {
  cat("\nFound", nrow(fusion_public_neoantigens), "fusion-derived potential public neoantigens\n")
  cat("Junction-spanning fusion public neoantigens:", sum(fusion_public_neoantigens$spans_junction), "\n")
}

# Create a subset of combined_analysis with key columns for the final output
combined_analysis_final <- combined_analysis %>%
  select(
    # Peptide information
    Peptide, 
    peptide_length,
    primary_gene,
    genes_combined,
    proteins_combined,
    
    # Immunopeptidome data
    spectral_count_148N,
    spectral_count_148T,
    total_intensity_148N,
    total_intensity_148T,
    log2_fold_change_immuno,
    detection_status,
    
    # Transcriptome data
    log2_fold_change_transcriptome,
    
    # LFQ proteome data
    log2_fold_change_lfq,
    p_value_lfq,
    
    # TMT proteome data
    log2_fold_change_tmt,
    p_value_tmt,
    
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
    upregulated_in_all,
    upregulated_in_3_datasets,
    tumor_specific_upregulated
  )

# Regenerate expression categories for plots
combined_analysis_final <- combined_analysis_final %>%
  mutate(
    # Recreate transcriptome expression category
    expression_category_transcriptome = case_when(
      is.na(log2_fold_change_transcriptome) ~ "No transcriptome data",
      log2_fold_change_transcriptome > 1 & log2_fold_change_immuno > 1 ~ "Up in both",
      log2_fold_change_transcriptome < -1 & log2_fold_change_immuno < -1 ~ "Down in both",
      log2_fold_change_transcriptome > 1 & log2_fold_change_immuno < -1 ~ "Up in transcriptome, down in immunopeptidome",
      log2_fold_change_transcriptome < -1 & log2_fold_change_immuno > 1 ~ "Down in transcriptome, up in immunopeptidome",
      TRUE ~ "No significant change"
    ),
    
    # Recreate LFQ expression category
    expression_category_lfq = case_when(
      is.na(log2_fold_change_lfq) ~ "No LFQ proteome data",
      log2_fold_change_lfq > 1 & log2_fold_change_immuno > 1 ~ "Up in both",
      log2_fold_change_lfq < -1 & log2_fold_change_immuno < -1 ~ "Down in both",
      log2_fold_change_lfq > 1 & log2_fold_change_immuno < -1 ~ "Up in LFQ, down in immunopeptidome",
      log2_fold_change_lfq < -1 & log2_fold_change_immuno > 1 ~ "Down in LFQ, up in immunopeptidome",
      TRUE ~ "No significant change"
    ),
    
    # Recreate TMT expression category
    expression_category_tmt = case_when(
      is.na(log2_fold_change_tmt) ~ "No TMT proteome data",
      log2_fold_change_tmt > 1 & log2_fold_change_immuno > 1 ~ "Up in both",
      log2_fold_change_tmt < -1 & log2_fold_change_immuno < -1 ~ "Down in both",
      log2_fold_change_tmt > 1 & log2_fold_change_immuno < -1 ~ "Up in TMT, down in immunopeptidome",
      log2_fold_change_tmt < -1 & log2_fold_change_immuno > 1 ~ "Down in TMT, up in immunopeptidome",
      TRUE ~ "No significant change"
    )
  )

# Fix for the public neoantigens heatmap
if(nrow(public_neoantigens) > 0) {
  # Create a matrix for the heatmap
  heatmap_data <- public_neoantigens %>%
    select(
      Peptide,
      primary_gene,
      log2_fold_change_immuno,
      log2_fold_change_transcriptome,
      log2_fold_change_lfq,
      log2_fold_change_tmt,
      public_neoantigen_classification
    )
  
  # Create the matrix for visualization (replace NA with 0 for visualization)
  heatmap_matrix <- heatmap_data %>%
    select(log2_fold_change_immuno, log2_fold_change_transcriptome, 
           log2_fold_change_lfq, log2_fold_change_tmt) %>%
    as.matrix()
  
  # Replace NA with 0 for visualization purposes
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  
  # Make sure row names are unique by adding a sequence number for duplicates
  row_labels <- paste0(heatmap_data$primary_gene, " (", heatmap_data$Peptide, ")")
  make_unique_rownames <- function(labels) {
    result <- character(length(labels))
    counts <- table(labels)
    for (name in names(counts)) {
      if (counts[name] == 1) {
        # If there's only one occurrence, use the original name
        result[labels == name] <- name
      } else {
        # If there are multiple occurrences, add a counter
        counter <- 1
        for (i in which(labels == name)) {
          result[i] <- paste0(name, "_", counter)
          counter <- counter + 1
        }
      }
    }
    return(result)
  }
  
  # Create unique row names
  unique_rownames <- make_unique_rownames(row_labels)
  rownames(heatmap_matrix) <- unique_rownames
  
  # Create row annotations - use the unique row names for the annotation data frame as well
  row_annotation <- data.frame(
    Classification = heatmap_data$public_neoantigen_classification,
    row.names = unique_rownames
  )
  
  # Create the heatmap PDF
  pdf(file.path(viz_dir, "Public_Neoantigens_Heatmap.pdf"), 
      width = 12, height = max(8, nrow(heatmap_matrix)/3))
  
  pheatmap(
    heatmap_matrix,
    main = "Potential Public Neoantigens: Log2 Fold Changes Across Omics Datasets",
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
  
  # Create the heatmap PNG
  png(file.path(viz_dir, "Public_Neoantigens_Heatmap.png"), 
      width = 1200, height = max(800, nrow(heatmap_matrix)*40), res = 120)
  
  pheatmap(
    heatmap_matrix,
    main = "Potential Public Neoantigens: Log2 Fold Changes Across Omics Datasets",
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
}

#--------------------------------------------------
# PART 5: Save outputs to the new directory
#--------------------------------------------------

# Create Excel output with multiple sheets
excel_sheets <- list(
  "Combined_Analysis" = combined_analysis_final,
  "Immunopeptidome_Only" = immunopeptidome_analysis %>% 
    select(-total_intensity_148N_adj, -total_intensity_148T_adj)
)

# Add the 4 individual data types sheets
if(exists("transcriptome_processed") && nrow(transcriptome_processed) > 0) {
  excel_sheets[["Transcriptome_Only"]] <- transcriptome_processed
}

if(exists("lfq_processed") && nrow(lfq_processed) > 0) {
  excel_sheets[["LFQ_Proteome_Only"]] <- lfq_processed
}

if(exists("tmt_processed") && nrow(tmt_processed) > 0) {
  excel_sheets[["TMT_Proteome_Only"]] <- tmt_processed
}

# Add public neoantigens sheet
if(nrow(public_neoantigens) > 0) {
  excel_sheets[["Public_Neoantigens"]] <- public_neoantigens
}

# Add fusion peptide sheets
if(nrow(fusion_peptides_analysis) > 0) {
  excel_sheets[["Fusion_Peptides"]] <- fusion_peptides_analysis %>% 
    select(-total_intensity_148N_adj, -total_intensity_148T_adj)
}

if(exists("fusion_peptides_with_transcriptome") && nrow(fusion_peptides_with_transcriptome) > 0) {
  excel_sheets[["Fusion_Peptides_With_Transcriptome"]] <- fusion_peptides_with_transcriptome
}

if(exists("fusion_upregulated_both") && nrow(fusion_upregulated_both) > 0) {
  excel_sheets[["Fusion_Upregulated_Both"]] <- fusion_upregulated_both
}

if(exists("fusion_public_neoantigens") && nrow(fusion_public_neoantigens) > 0) {
  excel_sheets[["Fusion_Public_Neoantigens"]] <- fusion_public_neoantigens
}

# Add tumor/normal-specific peptides sheets
excel_sheets[["Tumor_Specific_Peptides"]] <- combined_analysis_final %>% 
  filter(detection_status == "Tumor-specific")

excel_sheets[["Normal_Specific_Peptides"]] <- combined_analysis_final %>% 
  filter(detection_status == "Normal-specific")

# Add tier-specific public neoantigen sheets
if(sum(public_neoantigens$public_neoantigen_classification == "Tier 1 (Upregulated in all 4 datasets)") > 0) {
  excel_sheets[["Tier1_Public_Neoantigens"]] <- public_neoantigens %>%
    filter(public_neoantigen_classification == "Tier 1 (Upregulated in all 4 datasets)")
}

if(sum(public_neoantigens$public_neoantigen_classification == "Tier 2 (Upregulated in immunopeptidome + 2 others)") > 0) {
  excel_sheets[["Tier2_Public_Neoantigens"]] <- public_neoantigens %>%
    filter(public_neoantigen_classification == "Tier 2 (Upregulated in immunopeptidome + 2 others)")
}

if(sum(public_neoantigens$public_neoantigen_classification == "Tier 3 (Tumor-specific + upregulated in 2 others)") > 0) {
  excel_sheets[["Tier3_Public_Neoantigens"]] <- public_neoantigens %>%
    filter(public_neoantigen_classification == "Tier 3 (Tumor-specific + upregulated in 2 others)")
}

# Add pairwise comparison sheets
# Immunopeptidome vs Transcriptome - with column existence check
excel_sheets[["Immuno_vs_Transcriptome"]] <- combined_analysis_final %>%
  filter(!is.na(log2_fold_change_transcriptome)) %>%
  # Check if expression_category_transcriptome exists, otherwise create it
  {
    df <- .
    if (!"expression_category_transcriptome" %in% colnames(df)) {
      # Create the column if it doesn't exist
      df <- df %>%
        mutate(
          expression_category_transcriptome = case_when(
            is.na(log2_fold_change_transcriptome) ~ "No transcriptome data",
            log2_fold_change_transcriptome > 1 & log2_fold_change_immuno > 1 ~ "Up in both",
            log2_fold_change_transcriptome < -1 & log2_fold_change_immuno < -1 ~ "Down in both",
            log2_fold_change_transcriptome > 1 & log2_fold_change_immuno < -1 ~ "Up in transcriptome, down in immunopeptidome",
            log2_fold_change_transcriptome < -1 & log2_fold_change_immuno > 1 ~ "Down in transcriptome, up in immunopeptidome",
            TRUE ~ "No significant change"
          )
        )
    }
    df
  } %>%
  select(
    Peptide, primary_gene, peptide_length,
    log2_fold_change_immuno, log2_fold_change_transcriptome,
    expression_category_transcriptome, detection_status
  ) %>%
  arrange(desc(log2_fold_change_immuno))

# Immunopeptidome vs LFQ - with column existence check
excel_sheets[["Immuno_vs_LFQ"]] <- combined_analysis_final %>%
  filter(!is.na(log2_fold_change_lfq)) %>%
  {
    df <- .
    if (!"expression_category_lfq" %in% colnames(df)) {
      # Create the column if it doesn't exist
      df <- df %>%
        mutate(
          expression_category_lfq = case_when(
            is.na(log2_fold_change_lfq) ~ "No LFQ proteome data",
            log2_fold_change_lfq > 1 & log2_fold_change_immuno > 1 ~ "Up in both",
            log2_fold_change_lfq < -1 & log2_fold_change_immuno < -1 ~ "Down in both",
            log2_fold_change_lfq > 1 & log2_fold_change_immuno < -1 ~ "Up in LFQ, down in immunopeptidome",
            log2_fold_change_lfq < -1 & log2_fold_change_immuno > 1 ~ "Down in LFQ, up in immunopeptidome",
            TRUE ~ "No significant change"
          )
        )
    }
    df
  } %>%
  select(
    Peptide, primary_gene, peptide_length,
    log2_fold_change_immuno, log2_fold_change_lfq,
    p_value_lfq, expression_category_lfq, detection_status
  ) %>%
  arrange(desc(log2_fold_change_immuno))

# Immunopeptidome vs TMT - with column existence check
excel_sheets[["Immuno_vs_TMT"]] <- combined_analysis_final %>%
  filter(!is.na(log2_fold_change_tmt)) %>%
  {
    df <- .
    if (!"expression_category_tmt" %in% colnames(df)) {
      # Create the column if it doesn't exist
      df <- df %>%
        mutate(
          expression_category_tmt = case_when(
            is.na(log2_fold_change_tmt) ~ "No TMT proteome data",
            log2_fold_change_tmt > 1 & log2_fold_change_immuno > 1 ~ "Up in both",
            log2_fold_change_tmt < -1 & log2_fold_change_immuno < -1 ~ "Down in both",
            log2_fold_change_tmt > 1 & log2_fold_change_immuno < -1 ~ "Up in TMT, down in immunopeptidome",
            log2_fold_change_tmt < -1 & log2_fold_change_immuno > 1 ~ "Down in TMT, up in immunopeptidome",
            TRUE ~ "No significant change"
          )
        )
    }
    df
  } %>%
  select(
    Peptide, primary_gene, peptide_length,
    log2_fold_change_immuno, log2_fold_change_tmt,
    p_value_tmt, expression_category_tmt, detection_status
  ) %>%
  arrange(desc(log2_fold_change_immuno))

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
  if(grepl("Public_Neoantigens|Tier", sheet_name)) {
    # Highlight upregulated genes
    conditionalFormatting(wb, sheet_name, 
                          cols = which(colnames(excel_sheets[[sheet_name]]) %in% 
                                         c("log2_fold_change_immuno", "log2_fold_change_transcriptome", 
                                           "log2_fold_change_lfq", "log2_fold_change_tmt")), 
                          rows = 2:(nrow(excel_sheets[[sheet_name]]) + 1), 
                          rule = ">1", 
                          style = createStyle(bgFill = "#E2EFDA")) # Light green
  }
  
  # Add conditional formatting for fold change columns
  if(!grepl("Only", sheet_name)) {
    # Highlight upregulated (red)
    if("log2_fold_change_immuno" %in% colnames(excel_sheets[[sheet_name]])) {
      conditionalFormatting(wb, sheet_name, 
                            cols = which(colnames(excel_sheets[[sheet_name]]) == "log2_fold_change_immuno"), 
                            rows = 2:(nrow(excel_sheets[[sheet_name]]) + 1), 
                            rule = ">1", 
                            style = createStyle(bgFill = "#FFCCCC")) # Light red
    }
    
    # Highlight downregulated (blue)
    if("log2_fold_change_immuno" %in% colnames(excel_sheets[[sheet_name]])) {
      conditionalFormatting(wb, sheet_name, 
                            cols = which(colnames(excel_sheets[[sheet_name]]) == "log2_fold_change_immuno"), 
                            rows = 2:(nrow(excel_sheets[[sheet_name]]) + 1), 
                            rule = "<-1", 
                            style = createStyle(bgFill = "#CCCCFF")) # Light blue
    }
  }
}

# Save the Excel file
output_file <- file.path(output_dir, "RU148_4Way_Omics_Comparison.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

#--------------------------------------------------
# PART 6: Create visualizations
#--------------------------------------------------

#--------------------------------------------------
# PART 6.1: Pairwise scatter plots for all omics comparisons
#--------------------------------------------------

# Function to create standardized scatter plots
create_scatter_plot <- function(data, x_col, y_col, x_label, y_label, title, subtitle, color_col = NULL, color_label = NULL, highlight_points = NULL) {
  # Filter out NAs
  plot_data <- data %>%
    filter(!is.na(!!sym(x_col))) %>%
    filter(!is.na(!!sym(y_col)))
  
  if(nrow(plot_data) == 0) {
    cat("Warning: No data available for scatter plot:", title, "\n")
    return(NULL)
  }
  
  # Base plot
  p <- ggplot(plot_data, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    # Add lines at +/- 1 log2FC
    geom_hline(yintercept = 1, linetype = "dotted", color = "darkgray") +
    geom_hline(yintercept = -1, linetype = "dotted", color = "darkgray") +
    geom_vline(xintercept = 1, linetype = "dotted", color = "darkgray") +
    geom_vline(xintercept = -1, linetype = "dotted", color = "darkgray") +
    theme_minimal() +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  # Add color if specified
  if(!is.null(color_col)) {
    p <- p + geom_point(aes(color = !!sym(color_col)), alpha = 0.7) +
      labs(color = color_label)
  } else {
    p <- p + geom_point(alpha = 0.7, color = "steelblue")
  }
  
  # Highlight specific points if specified
  if(!is.null(highlight_points)) {
    # Filter only highlighted points
    highlighted_data <- plot_data %>%
      filter(!!sym(highlight_points))
    
    if(nrow(highlighted_data) > 0) {
      p <- p + 
        geom_point(data = highlighted_data, 
                   aes(x = !!sym(x_col), y = !!sym(y_col)), 
                   color = "red", size = 3, shape = 17) +
        geom_text(data = highlighted_data,
                  aes(x = !!sym(x_col), y = !!sym(y_col), label = primary_gene),
                  hjust = -0.1, vjust = -0.1, size = 3)
    }
  }
  
  return(p)
}

# 1. Immunopeptidome vs Transcriptome
immuno_trans_scatter <- create_scatter_plot(
  combined_analysis_final, 
  "log2_fold_change_transcriptome", 
  "log2_fold_change_immuno",
  "Log2 Fold Change Transcriptome (Tumor/Normal)",
  "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
  "Comparison of RU148 Tumor/Normal Fold Changes",
  "Immunopeptidome vs Transcriptome",
  "expression_category_transcriptome",
  "Expression Category",
  "potential_public_neoantigen"
)

if(!is.null(immuno_trans_scatter)) {
  pdf(file.path(viz_dir, "Immuno_vs_Transcriptome_Scatter.pdf"), width = 10, height = 8)
  print(immuno_trans_scatter)
  dev.off()
  
  png(file.path(viz_dir, "Immuno_vs_Transcriptome_Scatter.png"), width = 800, height = 600, res = 100)
  print(immuno_trans_scatter)
  dev.off()
}

# 2. Immunopeptidome vs LFQ Proteome
immuno_lfq_scatter <- create_scatter_plot(
  combined_analysis_final, 
  "log2_fold_change_lfq", 
  "log2_fold_change_immuno",
  "Log2 Fold Change LFQ Proteome (Tumor/Normal)",
  "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
  "Comparison of RU148 Tumor/Normal Fold Changes",
  "Immunopeptidome vs LFQ Proteome",
  "expression_category_lfq",
  "Expression Category",
  "potential_public_neoantigen"
)

if(!is.null(immuno_lfq_scatter)) {
  pdf(file.path(viz_dir, "Immuno_vs_LFQ_Scatter.pdf"), width = 10, height = 8)
  print(immuno_lfq_scatter)
  dev.off()
  
  png(file.path(viz_dir, "Immuno_vs_LFQ_Scatter.png"), width = 800, height = 600, res = 100)
  print(immuno_lfq_scatter)
  dev.off()
}

# 3. Immunopeptidome vs TMT Proteome
immuno_tmt_scatter <- create_scatter_plot(
  combined_analysis_final, 
  "log2_fold_change_tmt", 
  "log2_fold_change_immuno",
  "Log2 Fold Change TMT Proteome (Tumor/Normal)",
  "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
  "Comparison of RU148 Tumor/Normal Fold Changes",
  "Immunopeptidome vs TMT Proteome",
  "expression_category_tmt",
  "Expression Category",
  "potential_public_neoantigen"
)

if(!is.null(immuno_tmt_scatter)) {
  pdf(file.path(viz_dir, "Immuno_vs_TMT_Scatter.pdf"), width = 10, height = 8)
  print(immuno_tmt_scatter)
  dev.off()
  
  png(file.path(viz_dir, "Immuno_vs_TMT_Scatter.png"), width = 800, height = 600, res = 100)
  print(immuno_tmt_scatter)
  dev.off()
}

# 4. LFQ vs TMT Proteome
lfq_tmt_scatter <- create_scatter_plot(
  combined_analysis_final, 
  "log2_fold_change_tmt", 
  "log2_fold_change_lfq",
  "Log2 Fold Change TMT Proteome (Tumor/Normal)",
  "Log2 Fold Change LFQ Proteome (Tumor/Normal)",
  "Comparison of RU148 Tumor/Normal Fold Changes",
  "LFQ Proteome vs TMT Proteome",
  NULL,
  NULL,
  "potential_public_neoantigen"
)

if(!is.null(lfq_tmt_scatter)) {
  pdf(file.path(viz_dir, "LFQ_vs_TMT_Scatter.pdf"), width = 10, height = 8)
  print(lfq_tmt_scatter)
  dev.off()
  
  png(file.path(viz_dir, "LFQ_vs_TMT_Scatter.png"), width = 800, height = 600, res = 100)
  print(lfq_tmt_scatter)
  dev.off()
}

# 5. Transcriptome vs LFQ Proteome
trans_lfq_scatter <- create_scatter_plot(
  combined_analysis_final, 
  "log2_fold_change_lfq", 
  "log2_fold_change_transcriptome",
  "Log2 Fold Change LFQ Proteome (Tumor/Normal)",
  "Log2 Fold Change Transcriptome (Tumor/Normal)",
  "Comparison of RU148 Tumor/Normal Fold Changes",
  "Transcriptome vs LFQ Proteome",
  NULL,
  NULL,
  "potential_public_neoantigen"
)

if(!is.null(trans_lfq_scatter)) {
  pdf(file.path(viz_dir, "Transcriptome_vs_LFQ_Scatter.pdf"), width = 10, height = 8)
  print(trans_lfq_scatter)
  dev.off()
  
  png(file.path(viz_dir, "Transcriptome_vs_LFQ_Scatter.png"), width = 800, height = 600, res = 100)
  print(trans_lfq_scatter)
  dev.off()
}

# 6. Transcriptome vs TMT Proteome
trans_tmt_scatter <- create_scatter_plot(
  combined_analysis_final, 
  "log2_fold_change_tmt", 
  "log2_fold_change_transcriptome",
  "Log2 Fold Change TMT Proteome (Tumor/Normal)",
  "Log2 Fold Change Transcriptome (Tumor/Normal)",
  "Comparison of RU148 Tumor/Normal Fold Changes",
  "Transcriptome vs TMT Proteome",
  NULL,
  NULL,
  "potential_public_neoantigen"
)

if(!is.null(trans_tmt_scatter)) {
  pdf(file.path(viz_dir, "Transcriptome_vs_TMT_Scatter.pdf"), width = 10, height = 8)
  print(trans_tmt_scatter)
  dev.off()
  
  png(file.path(viz_dir, "Transcriptome_vs_TMT_Scatter.png"), width = 800, height = 600, res = 100)
  print(trans_tmt_scatter)
  dev.off()
}

#--------------------------------------------------
# PART 6.2: Public neoantigen visualizations
#--------------------------------------------------

# Create summary of public neoantigen tiers
if(nrow(public_neoantigens) > 0) {
  # Create summary of tier counts
  tier_summary <- public_neoantigens %>%
    group_by(public_neoantigen_classification) %>%
    summarise(
      count = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(count))
  
  # Barplot of public neoantigen tiers
  if(nrow(tier_summary) > 0) {
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
  }
}

#--------------------------------------------------
# PART 6.3: Fusion peptide visualizations
#--------------------------------------------------

# 1. Heatmap of fusion peptides (if any found)
if (nrow(fusion_peptides_analysis) > 0) {
  # Create a matrix for the heatmap
  fusion_intensity_data <- fusion_peptides_analysis %>%
    select(Peptide, total_intensity_148T, total_intensity_148N, fusion_peptide_type, spans_junction) %>%
    pivot_longer(
      cols = c(total_intensity_148T, total_intensity_148N),
      names_to = "Sample",
      values_to = "Intensity"
    ) %>%
    mutate(Sample = gsub("total_intensity_", "", Sample)) %>%
    pivot_wider(
      names_from = Sample,
      values_from = Intensity
    ) %>%
    arrange(desc(spans_junction), Peptide)
  
  # Log transform the values
  fusion_matrix <- as.matrix(fusion_intensity_data[, c("148T", "148N")])
  rownames(fusion_matrix) <- fusion_intensity_data$Peptide
  log_fusion_matrix <- log10(fusion_matrix + 1)
  
  # Create annotation for the rows
  row_annotation <- data.frame(
    Peptide_Type = fusion_intensity_data$fusion_peptide_type,
    Spans_Junction = fusion_intensity_data$spans_junction,
    row.names = fusion_intensity_data$Peptide
  )
  
  # Create a PDF of the heatmap
  pdf(file.path(viz_dir, "Fusion_Peptides_Heatmap.pdf"), 
      width = 10, height = max(8, nrow(fusion_matrix)/3))
  pheatmap(
    log_fusion_matrix,
    main = "Intensity of DNAJB1_PRKACA Fusion Peptides in 148T vs 148N (log10)",
    color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    display_numbers = TRUE,
    number_format = "%.1f",
    fontsize_row = 10,
    fontsize_col = 10
  )
  dev.off()
  
  # Create a PNG of the heatmap
  png(file.path(viz_dir, "Fusion_Peptides_Heatmap.png"), 
      width = 800, height = max(600, nrow(fusion_matrix)*40), res = 100)
  pheatmap(
    log_fusion_matrix,
    main = "Intensity of DNAJB1_PRKACA Fusion Peptides in 148T vs 148N (log10)",
    color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    display_numbers = TRUE,
    number_format = "%.1f",
    fontsize_row = 10,
    fontsize_col = 10
  )
  dev.off()
  
  # Create a heatmap of fusion peptides across all omics datasets
  fusion_with_omics <- combined_analysis_final %>%
    filter(from_fusion) %>%
    filter(!is.na(log2_fold_change_transcriptome) | !is.na(log2_fold_change_lfq) | !is.na(log2_fold_change_tmt))
  
  if(nrow(fusion_with_omics) > 0) {
    # Prepare data for heatmap
    fusion_omics_data <- fusion_with_omics %>%
      select(
        Peptide,
        primary_gene,
        log2_fold_change_immuno,
        log2_fold_change_transcriptome,
        log2_fold_change_lfq,
        log2_fold_change_tmt,
        fusion_peptide_type,
        spans_junction
      )
    
    # Replace NA with 0 for visualization purposes
    fusion_omics_data[is.na(fusion_omics_data)] <- 0
    
    # Create matrix for heatmap
    fusion_omics_matrix <- as.matrix(fusion_omics_data[, c("log2_fold_change_immuno", 
                                                           "log2_fold_change_transcriptome",
                                                           "log2_fold_change_lfq", 
                                                           "log2_fold_change_tmt")])
    rownames(fusion_omics_matrix) <- fusion_omics_data$Peptide
    
    # Create annotation for the rows
    fusion_row_annotation <- data.frame(
      Peptide_Type = fusion_omics_data$fusion_peptide_type,
      Spans_Junction = fusion_omics_data$spans_junction,
      row.names = rownames(fusion_omics_matrix)
    )
    
    # Create a PDF of the heatmap
    pdf(file.path(viz_dir, "Fusion_Peptides_Omics_Heatmap.pdf"), 
        width = 12, height = max(8, nrow(fusion_omics_matrix)/3))
    
    pheatmap(
      fusion_omics_matrix,
      main = "DNAJB1-PRKACA Fusion Peptides: Log2 Fold Changes Across Omics Datasets",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-3, 3, length.out = 101),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = fusion_row_annotation,
      display_numbers = TRUE,
      number_format = "%.1f",
      fontsize_row = 10,
      fontsize_col = 10,
      labels_col = c("Immunopeptidome", "Transcriptome", "LFQ Proteome", "TMT Proteome")
    )
    dev.off()
    
    # Create a PNG of the heatmap
    png(file.path(viz_dir, "Fusion_Peptides_Omics_Heatmap.png"), 
        width = 1200, height = max(800, nrow(fusion_omics_matrix)*40), res = 120)
    
    pheatmap(
      fusion_omics_matrix,
      main = "DNAJB1-PRKACA Fusion Peptides: Log2 Fold Changes Across Omics Datasets",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-3, 3, length.out = 101),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = fusion_row_annotation,
      display_numbers = TRUE,
      number_format = "%.1f",
      fontsize_row = 10,
      fontsize_col = 10,
      labels_col = c("Immunopeptidome", "Transcriptome", "LFQ Proteome", "TMT Proteome")
    )
    dev.off()
  }
  
  # Generate scatter plots for fusion peptides
  if(any(!is.na(fusion_with_omics$log2_fold_change_transcriptome))) {
    # Immunopeptidome vs Transcriptome for fusion peptides
    fusion_scatter_data <- fusion_with_omics %>%
      filter(!is.na(log2_fold_change_transcriptome))
    
    if(nrow(fusion_scatter_data) > 0) {
      pdf(file.path(viz_dir, "Fusion_Immuno_vs_Transcriptome_Scatter.pdf"), width = 10, height = 8)
      fusion_scatter_plot <- ggplot(fusion_scatter_data, 
                                    aes(x = log2_fold_change_transcriptome, 
                                        y = log2_fold_change_immuno, 
                                        color = fusion_peptide_type,
                                        shape = spans_junction)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_text(aes(label = Peptide), hjust = -0.1, vjust = 0.2, size = 3) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
        scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) +
        theme_minimal() +
        labs(
          title = "Fusion Peptides: Immunopeptidome vs Transcriptome Fold Changes",
          subtitle = "Triangles indicate junction-spanning peptides",
          x = "Log2 Fold Change Transcriptome (Tumor/Normal)",
          y = "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
          color = "Fusion Peptide Type",
          shape = "Spans Junction"
        ) +
        theme(
          legend.position = "right",
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12)
        )
      print(fusion_scatter_plot)
      dev.off()
      
      png(file.path(viz_dir, "Fusion_Immuno_vs_Transcriptome_Scatter.png"), 
          width = 800, height = 600, res = 100)
      print(fusion_scatter_plot)
      dev.off()
    }
  }
}

#--------------------------------------------------
# PART 6.4: Additional visualizations
#--------------------------------------------------

# 1. Volcano plot of immunopeptidome data
volcano_data <- immunopeptidome_analysis %>%
  filter(!is.na(log2_fold_change_immuno))

pdf(file.path(viz_dir, "Immunopeptidome_Volcano_Plot.pdf"), width = 10, height = 8)
volcano_plot <- ggplot(volcano_data, aes(x = log2_fold_change_immuno, y = -log10(0.05), 
                                         color = peptide_category)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +
  scale_color_manual(values = c("Up in Tumor (FC > 2)" = "red", 
                                "Down in Tumor (FC < 0.5)" = "blue", 
                                "Similar (-1 < log2FC < 1)" = "gray")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Peptides in 148T vs 148N",
    subtitle = "Red: Upregulated in Tumor, Blue: Downregulated in Tumor",
    x = "Log2 Fold Change (Tumor/Normal)",
    y = "-log10(p-value) [placeholder]",
    color = "Peptide Category"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )
print(volcano_plot)
dev.off()

png(file.path(viz_dir, "Immunopeptidome_Volcano_Plot.png"), width = 800, height = 600, res = 100)
print(volcano_plot)
dev.off()

# 2. Create Venn diagram of significant genes across datasets
# Extract gene lists that are upregulated in each dataset
genes_immuno_up <- combined_analysis_final %>%
  filter(log2_fold_change_immuno > 1) %>%
  pull(primary_gene) %>%
  unique()

genes_trans_up <- combined_analysis_final %>%
  filter(!is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > 1) %>%
  pull(primary_gene) %>%
  unique()

genes_lfq_up <- combined_analysis_final %>%
  filter(!is.na(log2_fold_change_lfq) & log2_fold_change_lfq > 1) %>%
  pull(primary_gene) %>%
  unique()

genes_tmt_up <- combined_analysis_final %>%
  filter(!is.na(log2_fold_change_tmt) & log2_fold_change_tmt > 1) %>%
  pull(primary_gene) %>%
  unique()

# Create a Venn diagram PDF
venn_colors <- brewer.pal(4, "Set1")
venn_output <- file.path(viz_dir, "Upregulated_Genes_Venn.png")

# Create Venn diagram
png(venn_output, width = 800, height = 800, res = 120)
draw.quad.venn(
  area1 = length(genes_immuno_up),
  area2 = length(genes_trans_up),
  area3 = length(genes_lfq_up),
  area4 = length(genes_tmt_up),
  n12 = length(intersect(genes_immuno_up, genes_trans_up)),
  n13 = length(intersect(genes_immuno_up, genes_lfq_up)),
  n14 = length(intersect(genes_immuno_up, genes_tmt_up)),
  n23 = length(intersect(genes_trans_up, genes_lfq_up)),
  n24 = length(intersect(genes_trans_up, genes_tmt_up)),
  n34 = length(intersect(genes_lfq_up, genes_tmt_up)),
  n123 = length(intersect(intersect(genes_immuno_up, genes_trans_up), genes_lfq_up)),
  n124 = length(intersect(intersect(genes_immuno_up, genes_trans_up), genes_tmt_up)),
  n134 = length(intersect(intersect(genes_immuno_up, genes_lfq_up), genes_tmt_up)),
  n234 = length(intersect(intersect(genes_trans_up, genes_lfq_up), genes_tmt_up)),
  n1234 = length(intersect(intersect(intersect(genes_immuno_up, genes_trans_up), genes_lfq_up), genes_tmt_up)),
  category = c("Immunopeptidome", "Transcriptome", "LFQ Proteome", "TMT Proteome"),
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

# 3. Create bar chart of overlapping genes
overlap_summary <- data.frame(
  Category = c(
    "Immuno only", 
    "Transcriptome only", 
    "LFQ only", 
    "TMT only",
    "Immuno + Trans", 
    "Immuno + LFQ", 
    "Immuno + TMT",
    "Trans + LFQ", 
    "Trans + TMT", 
    "LFQ + TMT",
    "Immuno + Trans + LFQ",
    "Immuno + Trans + TMT",
    "Immuno + LFQ + TMT",
    "Trans + LFQ + TMT",
    "All 4 datasets"
  ),
  Count = c(
    length(setdiff(genes_immuno_up, union(union(genes_trans_up, genes_lfq_up), genes_tmt_up))),
    length(setdiff(genes_trans_up, union(union(genes_immuno_up, genes_lfq_up), genes_tmt_up))),
    length(setdiff(genes_lfq_up, union(union(genes_immuno_up, genes_trans_up), genes_tmt_up))),
    length(setdiff(genes_tmt_up, union(union(genes_immuno_up, genes_trans_up), genes_lfq_up))),
    
    length(setdiff(intersect(genes_immuno_up, genes_trans_up), 
                   union(genes_lfq_up, genes_tmt_up))),
    length(setdiff(intersect(genes_immuno_up, genes_lfq_up), 
                   union(genes_trans_up, genes_tmt_up))),
    length(setdiff(intersect(genes_immuno_up, genes_tmt_up), 
                   union(genes_trans_up, genes_lfq_up))),
    length(setdiff(intersect(genes_trans_up, genes_lfq_up), 
                   union(genes_immuno_up, genes_tmt_up))),
    length(setdiff(intersect(genes_trans_up, genes_tmt_up), 
                   union(genes_immuno_up, genes_lfq_up))),
    length(setdiff(intersect(genes_lfq_up, genes_tmt_up), 
                   union(genes_immuno_up, genes_trans_up))),
    
    length(setdiff(intersect(intersect(genes_immuno_up, genes_trans_up), genes_lfq_up), 
                   genes_tmt_up)),
    length(setdiff(intersect(intersect(genes_immuno_up, genes_trans_up), genes_tmt_up), 
                   genes_lfq_up)),
    length(setdiff(intersect(intersect(genes_immuno_up, genes_lfq_up), genes_tmt_up), 
                   genes_trans_up)),
    length(setdiff(intersect(intersect(genes_trans_up, genes_lfq_up), genes_tmt_up), 
                   genes_immuno_up)),
    
    length(intersect(intersect(intersect(genes_immuno_up, genes_trans_up), genes_lfq_up), genes_tmt_up))
  ),
  Type = c(
    rep("Single dataset", 4),
    rep("Two datasets", 6),
    rep("Three datasets", 4),
    "All datasets"
  )
)

# Order the categories by count
overlap_summary$Category <- factor(overlap_summary$Category, 
                                   levels = overlap_summary$Category[order(overlap_summary$Count, decreasing = TRUE)])

# Create bar chart
pdf(file.path(viz_dir, "Upregulated_Genes_Overlap.pdf"), width = 12, height = 8)
overlap_plot <- ggplot(overlap_summary, aes(x = Category, y = Count, fill = Type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Overlap of Upregulated Genes Across Omics Datasets",
    x = "",
    y = "Number of Genes",
    fill = "Overlap Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold")
  )
print(overlap_plot)
dev.off()

png(file.path(viz_dir, "Upregulated_Genes_Overlap.png"), width = 1200, height = 800, res = 120)
print(overlap_plot)
dev.off()

# 4. Histogram of fold changes
pdf(file.path(viz_dir, "Immunopeptidome_Fold_Change_Histogram.pdf"), width = 10, height = 6)
hist_plot <- ggplot(volcano_data, aes(x = log2_fold_change_immuno, fill = peptide_category)) +
  geom_histogram(bins = 50, color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("Up in Tumor (FC > 2)" = "red", 
                               "Down in Tumor (FC < 0.5)" = "blue", 
                               "Similar (-1 < log2FC < 1)" = "gray")) +
  theme_minimal() +
  labs(
    title = "Distribution of Peptide Fold Changes in 148T vs 148N",
    x = "Log2 Fold Change (Tumor/Normal)",
    y = "Count",
    fill = "Peptide Category"
  )
print(hist_plot)
dev.off()

png(file.path(viz_dir, "Immunopeptidome_Fold_Change_Histogram.png"), width = 800, height = 600, res = 100)
print(hist_plot)
dev.off()

# 5. Barplot of detection status
detection_summary <- immunopeptidome_analysis %>%
  group_by(detection_status) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(count))

pdf(file.path(viz_dir, "Immunopeptidome_Detection_Status.pdf"), width = 8, height = 6)
detection_plot <- ggplot(detection_summary, aes(x = detection_status, y = count, fill = detection_status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = -0.5) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Peptide Detection Status in 148T vs 148N",
    x = "Detection Status",
    y = "Count",
    fill = "Status"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(detection_plot)
dev.off()

png(file.path(viz_dir, "Immunopeptidome_Detection_Status.png"), width = 800, height = 600, res = 100)
print(detection_plot)
dev.off()

#--------------------------------------------------
# PART 7: Print summary information
#--------------------------------------------------

# Print summary information
cat("\n----------------------------------------\n")
cat("SUMMARY OF 4-WAY OMICS COMPARISON ANALYSIS\n")
cat("----------------------------------------\n\n")

cat("Total peptides analyzed:", nrow(immunopeptidome_analysis), "\n")

if(exists("transcriptome_processed")) {
  cat("Total genes in transcriptome:", nrow(transcriptome_processed), "\n")
}

if(exists("lfq_processed")) {
  cat("Total significant genes in LFQ proteome:", nrow(lfq_processed), "\n")
}

if(exists("tmt_processed")) {
  cat("Total significant genes in TMT proteome:", nrow(tmt_processed), "\n")
}

cat("\nMatched data counts:\n")
cat("Peptides with matching transcriptome data:", 
    sum(!is.na(combined_analysis_final$log2_fold_change_transcriptome)), "\n")
cat("Peptides with matching LFQ proteome data:", 
    sum(!is.na(combined_analysis_final$log2_fold_change_lfq)), "\n")
cat("Peptides with matching TMT proteome data:", 
    sum(!is.na(combined_analysis_final$log2_fold_change_tmt)), "\n")

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
    arrange(desc(public_neoantigen_score), desc(log2_fold_change_immuno)) %>%
    head(10) %>%
    select(primary_gene, Peptide, public_neoantigen_classification, 
           log2_fold_change_immuno, log2_fold_change_transcriptome, 
           log2_fold_change_lfq, log2_fold_change_tmt)
  
  print(top_neoantigens)
}

# Print fusion peptide information if any found
if(nrow(fusion_peptides_analysis) > 0) {
  cat("\nFUSION PEPTIDE SUMMARY:\n")
  cat("Fusion peptides found:", nrow(fusion_peptides_analysis), "\n")
  cat("Junction-spanning peptides:", sum(fusion_peptides_analysis$spans_junction), "\n")
  
  cat("\nFusion peptide detection summary:\n")
  fusion_detection <- fusion_peptides_analysis %>%
    group_by(fusion_peptide_type, detection_status) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(fusion_peptide_type, detection_status)
  
  print(fusion_detection)
  
  if(nrow(fusion_public_neoantigens) > 0) {
    cat("\nFusion-derived public neoantigens:", nrow(fusion_public_neoantigens), "\n")
    print(fusion_public_neoantigens %>% 
            select(Peptide, fusion_peptide_type, spans_junction, 
                   public_neoantigen_classification, log2_fold_change_immuno))
  } else {
    cat("\nNo fusion-derived public neoantigens identified.\n")
  }
} else {
  cat("\nNo fusion peptides found in the analysis.\n")
}

# Print multi-omics overlap summary
cat("\nMULTI-OMICS OVERLAP SUMMARY:\n")
cat("Genes upregulated in all 4 datasets:", 
    length(intersect(intersect(intersect(genes_immuno_up, genes_trans_up), genes_lfq_up), genes_tmt_up)), "\n")
cat("Genes upregulated in immunopeptidome + at least 2 other datasets:", 
    length(filter(combined_analysis_final, upregulated_in_3_datasets)$primary_gene %>% unique()), "\n")
cat("Genes upregulated in immunopeptidome + transcriptome:", 
    length(intersect(genes_immuno_up, genes_trans_up)), "\n")
cat("Genes upregulated in immunopeptidome + LFQ proteome:", 
    length(intersect(genes_immuno_up, genes_lfq_up)), "\n")
cat("Genes upregulated in immunopeptidome + TMT proteome:", 
    length(intersect(genes_immuno_up, genes_tmt_up)), "\n")

cat("\nAnalysis complete! Results saved to:", output_dir, "\n")
cat("\nOUTPUT FILES:\n")
cat("1. RU148_4Way_Omics_Comparison.xlsx - Excel file with comprehensive analysis results\n")
cat("2. Multiple visualizations in the", viz_dir, "directory\n")

