# Modified script to analyze peptide differences and compare with transcriptome data

# Setting directory
setwd("~/Documents/Github/HLA-I_Analysis/")

# Create output directory in HLA-I_Analysis
output_dir <- "RU148_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Created output directory:", output_dir, "\n")
}

# Load required packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(writexl)   # For Excel output
library(openxlsx)  # For better Excel formatting
library(readxl)    # For reading Excel files

# Define the path to your data files
data_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/imp_014_rawdata"
transcriptome_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/Normalized_Gene_counts_FLCdb_Panel_1.xlsx"

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
  group_by(Sample_ID, Peptide) %>%
  summarize(
    peptide_length = first(nchar(Peptide)),
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

# After creating the immunopeptidome_analysis dataframe, add this section:

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

# Create new sheets in the Excel workbook for exclusive peptides
excel_sheets[["Tumor_Exclusive_Peptides"]] <- tumor_exclusive_peptides %>%
  select(-total_intensity_148N_adj, -total_intensity_148T_adj) %>%
  arrange(desc(total_intensity_148T))

excel_sheets[["Normal_Exclusive_Peptides"]] <- normal_exclusive_peptides %>%
  select(-total_intensity_148N_adj, -total_intensity_148T_adj) %>%
  arrange(desc(total_intensity_148N))

# Add fusion protein analysis for exclusive peptides
tumor_exclusive_fusion <- tumor_exclusive_peptides %>%
  filter(from_fusion) %>%
  arrange(desc(spans_junction), desc(total_intensity_148T))

normal_exclusive_fusion <- normal_exclusive_peptides %>%
  filter(from_fusion) %>%
  arrange(desc(spans_junction), desc(total_intensity_148N))

if(nrow(tumor_exclusive_fusion) > 0) {
  excel_sheets[["Tumor_Exclusive_Fusion"]] <- tumor_exclusive_fusion %>%
    select(-total_intensity_148N_adj, -total_intensity_148T_adj)
}

if(nrow(normal_exclusive_fusion) > 0) {
  excel_sheets[["Normal_Exclusive_Fusion"]] <- normal_exclusive_fusion %>%
    select(-total_intensity_148N_adj, -total_intensity_148T_adj)
}

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
      expression_category = case_when(
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
  combined_analysis$expression_category <- "No transcriptome data"
}

# Create a filtered dataset of fusion peptides with transcriptome information
fusion_peptides_with_transcriptome <- combined_analysis %>%
  filter(from_fusion) %>%
  arrange(desc(spans_junction), desc(log2_fold_change_immuno))

# Find fusion peptides upregulated in both immunopeptidome and transcriptome
fusion_upregulated_both <- fusion_peptides_with_transcriptome %>%
  filter(
    !is.na(log2_fold_change_transcriptome),
    log2_fold_change_immuno > 1,
    log2_fold_change_transcriptome > 1
  ) %>%
  arrange(desc(log2_fold_change_immuno))

# Select columns that definitely exist for the final combined analysis
combined_analysis_final <- combined_analysis %>%
  select(
    # Peptide information (these should always exist)
    Peptide, 
    peptide_length,
    primary_gene,
    genes_combined,
    proteins_combined,
    
    # Immunopeptidome data (these should always exist)
    spectral_count_148N,
    spectral_count_148T,
    total_intensity_148N,
    total_intensity_148T,
    log2_fold_change_immuno,
    detection_status,
    
    # Fusion information 
    from_fusion,
    spans_junction,
    from_dnajb1,
    from_prkaca,
    fusion_peptide_type
  )

# Add transcriptome columns if they exist
for(col in c("Mean.Normal", "Mean.Tumor", "RU148_N", "RU148_T8", "RU148_T11", 
             "RU148_T_Average", "log2_fold_change_transcriptome")) {
  if(col %in% colnames(combined_analysis)) {
    combined_analysis_final[[col]] <- combined_analysis[[col]]
  }
}

# Add gene information columns if they exist
for(col in c("geneID", "biotype", "chromosome", "gene_start", "gene_end", 
             "gene_length", "description")) {
  if(col %in% colnames(combined_analysis)) {
    combined_analysis_final[[col]] <- combined_analysis[[col]]
  }
}

# Add comparison metrics if they exist
for(col in c("expression_category", "immuno_trans_correlation")) {
  if(col %in% colnames(combined_analysis)) {
    combined_analysis_final[[col]] <- combined_analysis[[col]]
  }
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

# Add transcriptome sheet if we have data
if(exists("transcriptome_processed") && nrow(transcriptome_processed) > 0) {
  excel_sheets[["Transcriptome_Only"]] <- transcriptome_processed
}

# Add fusion peptide sheets
if(nrow(fusion_peptides_analysis) > 0) {
  excel_sheets[["Fusion_Peptides"]] <- fusion_peptides_analysis %>% 
    select(-total_intensity_148N_adj, -total_intensity_148T_adj)
}

if(nrow(fusion_peptides_with_transcriptome) > 0) {
  excel_sheets[["Fusion_Peptides_With_Transcriptome"]] <- fusion_peptides_with_transcriptome
}

if(nrow(fusion_upregulated_both) > 0) {
  excel_sheets[["Fusion_Upregulated_Both"]] <- fusion_upregulated_both
}

# Add additional analysis sheets if expression category exists
if("expression_category" %in% colnames(combined_analysis_final)) {
  # Tumor-specific peptides
  excel_sheets[["Tumor_Specific_Peptides"]] <- combined_analysis_final %>% 
    filter(detection_status == "Tumor-specific")
  
  # Normal-specific peptides
  excel_sheets[["Normal_Specific_Peptides"]] <- combined_analysis_final %>% 
    filter(detection_status == "Normal-specific")
  
  # Expression category-based sheets
  if("Up in both" %in% unique(combined_analysis_final$expression_category)) {
    excel_sheets[["Up_In_Both"]] <- combined_analysis_final %>% 
      filter(expression_category == "Up in both") %>%
      arrange(desc(log2_fold_change_immuno))
  }
  
  if("Down in both" %in% unique(combined_analysis_final$expression_category)) {
    excel_sheets[["Down_In_Both"]] <- combined_analysis_final %>% 
      filter(expression_category == "Down in both") %>%
      arrange(log2_fold_change_immuno)
  }
  
  # Check if we have any discordant expression
  discordant_categories <- c("Up in transcriptome, down in immunopeptidome", 
                             "Down in transcriptome, up in immunopeptidome")
  
  if(any(discordant_categories %in% unique(combined_analysis_final$expression_category))) {
    excel_sheets[["Discordant_Expression"]] <- combined_analysis_final %>% 
      filter(expression_category %in% discordant_categories)
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
}

# Save the Excel file
output_file <- file.path(output_dir, "RU148_Immuno_Transcriptome_Comparison.xlsx")
saveWorkbook(wb, output_file, overwrite = TRUE)

#--------------------------------------------------
# PART 6: Create visualizations
#--------------------------------------------------

# 1. Scatter plot comparing immunopeptidome and transcriptome fold changes
if("log2_fold_change_transcriptome" %in% colnames(combined_analysis_final) && 
   !all(is.na(combined_analysis_final$log2_fold_change_transcriptome))) {
  
  scatter_data <- combined_analysis_final %>%
    filter(!is.na(log2_fold_change_transcriptome)) %>%
    filter(!is.na(log2_fold_change_immuno))
  
  if(nrow(scatter_data) > 0) {
    pdf(file.path(output_dir, "Immuno_vs_Transcriptome_Scatter.pdf"), width = 10, height = 8)
    scatter_plot <- ggplot(scatter_data, 
                           aes(x = log2_fold_change_transcriptome, 
                               y = log2_fold_change_immuno, 
                               color = expression_category)) +
      geom_point(alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
      theme_minimal() +
      labs(
        title = "Comparison of RU148 Tumor/Normal Fold Changes",
        subtitle = "Immunopeptidome vs Transcriptome",
        x = "Log2 Fold Change Transcriptome (Tumor/Normal)",
        y = "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
        color = "Expression Category"
      ) +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    print(scatter_plot)
    dev.off()
    
    png(file.path(output_dir, "Immuno_vs_Transcriptome_Scatter.png"), width = 800, height = 600, res = 100)
    print(scatter_plot)
    dev.off()
  }
  
  # 2. Barplot of expression categories
  if("expression_category" %in% colnames(combined_analysis_final)) {
    expression_summary <- combined_analysis_final %>%
      group_by(expression_category) %>%
      summarise(
        count = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(count))
    
    if(nrow(expression_summary) > 0) {
      pdf(file.path(output_dir, "Expression_Category_Counts.pdf"), width = 10, height = 6)
      barplot <- ggplot(expression_summary, 
                        aes(x = reorder(expression_category, -count), 
                            y = count, 
                            fill = expression_category)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = count), vjust = -0.5) +
        theme_minimal() +
        labs(
          title = "Distribution of Expression Categories",
          x = "Category",
          y = "Count",
          fill = "Expression Category"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(barplot)
      dev.off()
      
      png(file.path(output_dir, "Expression_Category_Counts.png"), width = 800, height = 600, res = 100)
      print(barplot)
      dev.off()
    }
  }
}

# 3. Heatmap of fusion peptides (if any found)
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
    row.names = fusion_intensity_data$Peptide
  )
  
  # Create a PDF of the heatmap
  pdf(file.path(output_dir, "148T_vs_148N_fusion_peptides_heatmap.pdf"), 
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
  png(file.path(output_dir, "148T_vs_148N_fusion_peptides_heatmap.png"), 
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
  
  # 4. If we have transcriptome data for fusion peptides, create a scatter plot
  if(nrow(fusion_peptides_with_transcriptome) > 0 && 
     !all(is.na(fusion_peptides_with_transcriptome$log2_fold_change_transcriptome))) {
    
    fusion_scatter_data <- fusion_peptides_with_transcriptome %>%
      filter(!is.na(log2_fold_change_transcriptome))
    
    if(nrow(fusion_scatter_data) > 0) {
      pdf(file.path(output_dir, "Fusion_Immuno_vs_Transcriptome_Scatter.pdf"), width = 10, height = 8)
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
          plot.subtitle = element_text(size = 12)) +
        theme(
          legend.position = "right",
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12)
        )
      print(fusion_scatter_plot)
      dev.off()
      
      png(file.path(output_dir, "Fusion_Immuno_vs_Transcriptome_Scatter.png"), 
          width = 800, height = 600, res = 100)
      print(fusion_scatter_plot)
      dev.off()
    }
  }
}

# 5. Volcano plot of immunopeptidome data
volcano_data <- immunopeptidome_analysis %>%
  filter(!is.na(log2_fold_change_immuno))

pdf(file.path(output_dir, "148T_vs_148N_volcano_plot.pdf"), width = 10, height = 8)
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

png(file.path(output_dir, "148T_vs_148N_volcano_plot.png"), width = 800, height = 600, res = 100)
print(volcano_plot)
dev.off()

# 6. Histogram of fold changes
pdf(file.path(output_dir, "148T_vs_148N_fold_change_histogram.pdf"), width = 10, height = 6)
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

png(file.path(output_dir, "148T_vs_148N_fold_change_histogram.png"), width = 800, height = 600, res = 100)
print(hist_plot)
dev.off()

# 7. Barplot of detection status
detection_summary <- immunopeptidome_analysis %>%
  group_by(detection_status) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(count))

pdf(file.path(output_dir, "148T_vs_148N_detection_status.pdf"), width = 8, height = 6)
detection_plot <- ggplot(detection_summary, aes(x = detection_status, y = count, fill = detection_status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = -0.5) +
  theme_minimal() +
  labs(
    title = "Peptide Detection Status in 148T vs 148N",
    x = "Detection Status",
    y = "Count",
    fill = "Status"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(detection_plot)
dev.off()

png(file.path(output_dir, "148T_vs_148N_detection_status.png"), width = 800, height = 600, res = 100)
print(detection_plot)
dev.off()

#--------------------------------------------------
# PART 7: Print summary information
#--------------------------------------------------

# Print summary information
cat("\nSummary of RU148 Immunopeptidome and Transcriptome Analysis:\n")
cat("Total peptides analyzed:", nrow(immunopeptidome_analysis), "\n")

if(exists("transcriptome_processed")) {
  cat("Total genes in transcriptome:", nrow(transcriptome_processed), "\n")
}

if("log2_fold_change_transcriptome" %in% colnames(combined_analysis_final)) {
  cat("Peptides with matching transcriptome data:", 
      sum(!is.na(combined_analysis_final$log2_fold_change_transcriptome)), "\n")
}

# Print detection status summary
cat("\nDetection status summary:\n")
print(detection_summary)

# Print fusion peptide information if any found
if(nrow(fusion_peptides_analysis) > 0) {
  cat("\nFusion peptides found:", nrow(fusion_peptides_analysis), "\n")
  cat("Junction-spanning peptides:", sum(fusion_peptides_analysis$spans_junction), "\n")
  
  cat("\nFusion peptide detection summary:\n")
  fusion_peptides_analysis %>%
    group_by(fusion_peptide_type, detection_status) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(fusion_peptide_type, detection_status) %>%
    print(n = Inf)
  
  if(nrow(fusion_upregulated_both) > 0) {
    cat("\nFusion peptides upregulated in both immunopeptidome and transcriptome:", 
        nrow(fusion_upregulated_both), "\n")
    
    cat("\nList of fusion peptides upregulated in both datasets:\n")
    fusion_upregulated_both %>% 
      select(Peptide, fusion_peptide_type, log2_fold_change_immuno, log2_fold_change_transcriptome) %>%
      print(n = Inf)
  } else {
    cat("\nNo fusion peptides found to be upregulated in both immunopeptidome and transcriptome.\n")
  }
} else {
  cat("\nNo fusion peptides found in the analysis.\n")
}

# Print expression categories if available
if("expression_category" %in% colnames(combined_analysis_final)) {
  expression_summary <- combined_analysis_final %>%
    group_by(expression_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(count))
  
  cat("\nExpression categories:\n")
  print(expression_summary)
}

cat("\nAnalysis complete! Results saved to:", output_dir, "\n")
cat("\nThe following files were generated:\n")
cat("1. RU148_Immuno_Transcriptome_Comparison.xlsx - Excel file with comprehensive analysis results\n")

# List visualizations that were created
visualization_count = 2
if(exists("scatter_plot")) {
  cat(visualization_count, ". Immuno_vs_Transcriptome_Scatter.pdf/png - Scatter plot comparing fold changes\n", sep="")
  visualization_count = visualization_count + 1
}

if(exists("barplot")) {
  cat(visualization_count, ". Expression_Category_Counts.pdf/png - Barplot of expression categories\n", sep="")
  visualization_count = visualization_count + 1
}

if(nrow(fusion_peptides_analysis) > 0) {
  cat(visualization_count, ". 148T_vs_148N_fusion_peptides_heatmap.pdf/png - Heatmap of fusion peptides\n", sep="")
  visualization_count = visualization_count + 1
  
  if(exists("fusion_scatter_plot")) {
    cat(visualization_count, ". Fusion_Immuno_vs_Transcriptome_Scatter.pdf/png - Scatter plot of fusion peptides fold changes\n", sep="")
    visualization_count = visualization_count + 1
  }
}

cat(visualization_count, ". 148T_vs_148N_volcano_plot.pdf/png - Volcano plot showing peptide fold changes\n", sep="")
visualization_count = visualization_count + 1

cat(visualization_count, ". 148T_vs_148N_fold_change_histogram.pdf/png - Histogram of fold changes\n", sep="")
visualization_count = visualization_count + 1

cat(visualization_count, ". 148T_vs_148N_detection_status.pdf/png - Barplot of peptide detection status\n", sep="")

#--------------------------------------------------
# PART 8: Add sanity checks to validate the analysis
#--------------------------------------------------

# Add at the end of the script:

#--------------------------------------------------
# Sanity Checks
#--------------------------------------------------
cat("\n----- RUNNING SANITY CHECKS -----\n")

# 1. Check consistency in detection status vs intensity values
detection_sanity_check <- immunopeptidome_analysis %>%
  mutate(
    status_check = case_when(
      detection_status == "Tumor-specific" & total_intensity_148N > 0 ~ "FAIL",
      detection_status == "Normal-specific" & total_intensity_148T > 0 ~ "FAIL",
      detection_status == "Detected in both" & (total_intensity_148T == 0 | total_intensity_148N == 0) ~ "FAIL",
      TRUE ~ "PASS"
    )
  )

failed_detection <- detection_sanity_check %>% filter(status_check == "FAIL")
if(nrow(failed_detection) > 0) {
  cat("WARNING: Found", nrow(failed_detection), "peptides with inconsistent detection status\n")
  print(failed_detection %>% select(Peptide, detection_status, total_intensity_148T, total_intensity_148N))
} else {
  cat("Sanity check PASSED: Detection status is consistent with intensity values\n")
}

# 2. Check fold change calculations
fold_change_sanity_check <- immunopeptidome_analysis %>%
  mutate(
    expected_fold_change = log2(total_intensity_148T_adj / total_intensity_148N_adj),
    fold_change_diff = abs(expected_fold_change - log2_fold_change_immuno)
  ) %>%
  filter(fold_change_diff > 0.001)  # Allow for small floating point differences

if(nrow(fold_change_sanity_check) > 0) {
  cat("WARNING: Found", nrow(fold_change_sanity_check), "peptides with inconsistent fold change calculations\n")
  print(fold_change_sanity_check %>% select(Peptide, expected_fold_change, log2_fold_change_immuno, fold_change_diff))
} else {
  cat("Sanity check PASSED: Fold change calculations are consistent\n")
}

# 3. Check peptide category assignment
category_sanity_check <- immunopeptidome_analysis %>%
  mutate(
    expected_category = case_when(
      log2_fold_change_immuno > 1 ~ "Up in Tumor (FC > 2)",
      log2_fold_change_immuno < -1 ~ "Down in Tumor (FC < 0.5)",
      TRUE ~ "Similar (-1 < log2FC < 1)"
    ),
    category_check = ifelse(expected_category == peptide_category, "PASS", "FAIL")
  ) %>%
  filter(category_check == "FAIL")

if(nrow(category_sanity_check) > 0) {
  cat("WARNING: Found", nrow(category_sanity_check), "peptides with inconsistent category assignment\n")
  print(category_sanity_check %>% select(Peptide, log2_fold_change_immuno, expected_category, peptide_category))
} else {
  cat("Sanity check PASSED: Peptide category assignments are consistent\n")
}

# 4. Check fusion peptide analysis
if(nrow(fusion_peptides_analysis) > 0) {
  fusion_sanity_check <- fusion_peptides_analysis %>%
    filter(
      (from_dnajb1 & !grepl(dnajb1_part, Peptide, fixed = TRUE)) |
        (from_prkaca & !grepl(prkaca_part, Peptide, fixed = TRUE)) |
        (spans_junction & !(grepl(substr(dnajb1_part, nchar(dnajb1_part)-3, nchar(dnajb1_part)), Peptide, fixed = TRUE) &
                              grepl(substr(prkaca_part, 1, 4), Peptide, fixed = TRUE)))
    )
  
  if(nrow(fusion_sanity_check) > 0) {
    cat("WARNING: Found", nrow(fusion_sanity_check), "fusion peptides with inconsistent classification\n")
    print(fusion_sanity_check %>% select(Peptide, from_dnajb1, from_prkaca, spans_junction, fusion_peptide_type))
  } else {
    cat("Sanity check PASSED: Fusion peptide classifications are consistent\n")
  }
}

# 5. Check for NAs in key analytical columns
na_check <- immunopeptidome_analysis %>%
  summarise(
    NA_in_fold_change = sum(is.na(log2_fold_change_immuno)),
    NA_in_detection = sum(is.na(detection_status)),
    NA_in_category = sum(is.na(peptide_category))
  )

if(any(unlist(na_check) > 0)) {
  cat("WARNING: Found missing values in key analytical columns:\n")
  print(na_check)
} else {
  cat("Sanity check PASSED: No missing values in key analytical columns\n")
}

# 6. Check transcriptome analysis (if applicable)
if("log2_fold_change_transcriptome" %in% colnames(combined_analysis_final)) {
  trans_category_check <- combined_analysis_final %>%
    filter(!is.na(log2_fold_change_transcriptome)) %>%
    mutate(
      expected_category = case_when(
        is.na(log2_fold_change_transcriptome) ~ "No transcriptome data",
        log2_fold_change_transcriptome > 1 & log2_fold_change_immuno > 1 ~ "Up in both",
        log2_fold_change_transcriptome < -1 & log2_fold_change_immuno < -1 ~ "Down in both",
        log2_fold_change_transcriptome > 1 & log2_fold_change_immuno < -1 ~ "Up in transcriptome, down in immunopeptidome",
        log2_fold_change_transcriptome < -1 & log2_fold_change_immuno > 1 ~ "Down in transcriptome, up in immunopeptidome",
        TRUE ~ "No significant change"
      ),
      category_check = ifelse(expected_category == expression_category, "PASS", "FAIL")
    ) %>%
    filter(category_check == "FAIL")
  
  if(nrow(trans_category_check) > 0) {
    cat("WARNING: Found", nrow(trans_category_check), "peptides with inconsistent expression category\n")
    print(trans_category_check %>% select(Peptide, log2_fold_change_immuno, log2_fold_change_transcriptome, 
                                          expected_category, expression_category))
  } else {
    cat("Sanity check PASSED: Expression categories are consistent\n")
  }
}

# 7. Basic statistical sanity check
basic_stats <- immunopeptidome_analysis %>%
  summarise(
    min_fold_change = min(log2_fold_change_immuno, na.rm = TRUE),
    max_fold_change = max(log2_fold_change_immuno, na.rm = TRUE),
    mean_fold_change = mean(log2_fold_change_immuno, na.rm = TRUE),
    median_fold_change = median(log2_fold_change_immuno, na.rm = TRUE),
    min_intensity_tumor = min(total_intensity_148T, na.rm = TRUE),
    max_intensity_tumor = max(total_intensity_148T, na.rm = TRUE),
    min_intensity_normal = min(total_intensity_148N, na.rm = TRUE),
    max_intensity_normal = max(total_intensity_148N, na.rm = TRUE)
  )

cat("\nBasic statistics for sanity check:\n")
print(basic_stats)

# 8. Summary counts for validation
count_summary <- immunopeptidome_analysis %>%
  summarise(
    total_peptides = n(),
    tumor_specific = sum(detection_status == "Tumor-specific"),
    normal_specific = sum(detection_status == "Normal-specific"),
    in_both = sum(detection_status == "Detected in both"),
    not_detected = sum(detection_status == "Not detected"),
    up_in_tumor = sum(peptide_category == "Up in Tumor (FC > 2)"),
    down_in_tumor = sum(peptide_category == "Down in Tumor (FC < 0.5)"),
    similar = sum(peptide_category == "Similar (-1 < log2FC < 1)")
  )

cat("\nCount summary for validation:\n")
print(count_summary)

# Fix: Use a more direct approach that avoids the dataframe issue
total_count <- as.numeric(count_summary[1, "total_peptides"])
detection_sum <- as.numeric(count_summary[1, "tumor_specific"]) + 
  as.numeric(count_summary[1, "normal_specific"]) + 
  as.numeric(count_summary[1, "in_both"]) + 
  as.numeric(count_summary[1, "not_detected"])

category_sum <- as.numeric(count_summary[1, "up_in_tumor"]) + 
  as.numeric(count_summary[1, "down_in_tumor"]) + 
  as.numeric(count_summary[1, "similar"])

# Check that counts add up correctly
if(total_count != detection_sum) {
  cat("WARNING: Detection status counts don't add up to total peptides\n")
  cat("Total peptides:", total_count, "Sum of detection categories:", detection_sum, "\n")
} else {
  cat("Sanity check PASSED: Detection status counts add up correctly\n")
}

if(total_count != category_sum) {
  cat("WARNING: Peptide category counts don't add up to total peptides\n")
  cat("Total peptides:", total_count, "Sum of fold change categories:", category_sum, "\n")
} else {
  cat("Sanity check PASSED: Peptide category counts add up correctly\n")
}

# # Modified script that integrates both fusion protein analysis and transcriptome comparison
# 
# # [Keep all the initial setup and data loading code the same as the current script]
# 
# #--------------------------------------------------
# # Add fusion protein analysis from original script
# #--------------------------------------------------
# 
# # Identify DNAJB1-PRKACA fusion peptides
# # Define the sequences of DNAJB1 and PRKACA parts of the fusion protein
# dnajb1_seq <- "GKDYYQTLGLARGASDEEIKRAYRRQALRYHPDKNKEPGAEEKFKEIAEAYDVLSDPRKREIFDRYGEE"
# prkaca_seq <- "VKEFLAKAKEDFLKKWESPAQNTAHLDQFERIKTLGTGSFGRVMLVKHKETGNHYAMKILDKQKVVKLKQIEHTLNEKRILQAVNFPFLVKLEFSFKDNSNLYMVMEYVPGGEMFSHLRRIGRFSEPHARFYAAQIVLTFEYLHSLDLIYRDLKPENLLIDQQGYIQVTDFGFAKRVKGRTWTLCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFATTDWIAIYQRKVEAPFIPKFKGPGDTSNFDDYEEEEIRVSINEKCGKEFSEF"
# fusion_protein <- paste0(dnajb1_seq, prkaca_seq)
# 
# # Function to check if a peptide spans the fusion junction
# is_fusion_junction_peptide <- function(peptide_seq) {
#   # Define the end of DNAJB1 and start of PRKACA for the fusion
#   dnajb1_end <- "IFDRYGEE"
#   prkaca_start <- "VKEFLAK"
#   
#   # Check if the peptide spans the fusion junction
#   spans_junction <- FALSE
#   
#   if (nchar(peptide_seq) >= 6) { # Require at least 6 AA to be meaningful
#     for (i in 3:(nchar(peptide_seq) - 3)) { # At least 3 AA from each side
#       left_part <- substr(peptide_seq, 1, i)
#       right_part <- substr(peptide_seq, i + 1, nchar(peptide_seq))
#       
#       # Check if left part is in DNAJB1 (at the end) and right part in PRKACA (at the beginning)
#       if (grepl(left_part, dnajb1_seq, fixed = TRUE) && 
#           grepl(right_part, prkaca_seq, fixed = TRUE)) {
#         
#         # Additional check to ensure left part aligns with end of DNAJB1
#         left_pos <- gregexpr(left_part, dnajb1_seq, fixed = TRUE)[[1]]
#         if (length(left_pos) > 0 && any(left_pos + nchar(left_part) - 1 >= nchar(dnajb1_seq) - 10)) {
#           
#           # Additional check to ensure right part aligns with start of PRKACA
#           right_pos <- gregexpr(right_part, prkaca_seq, fixed = TRUE)[[1]]
#           if (length(right_pos) > 0 && any(right_pos <= 10)) {
#             spans_junction <- TRUE
#             break
#           }
#         }
#       }
#     }
#   }
#   
#   # Check if peptide is from either part of the fusion protein
#   from_dnajb1 <- grepl(peptide_seq, dnajb1_seq, fixed = TRUE)
#   from_prkaca <- grepl(peptide_seq, prkaca_seq, fixed = TRUE)
#   
#   return(list(
#     spans_junction = spans_junction,
#     from_dnajb1 = from_dnajb1,
#     from_prkaca = from_prkaca,
#     from_fusion = from_dnajb1 | from_prkaca | spans_junction
#   ))
# }
# 
# # Add fusion protein information to the immunopeptidome analysis
# immunopeptidome_analysis <- immunopeptidome_analysis %>%
#   rowwise() %>%
#   mutate(
#     fusion_info = list(is_fusion_junction_peptide(Peptide)),
#     from_fusion = fusion_info$from_fusion,
#     spans_junction = fusion_info$spans_junction,
#     from_dnajb1 = fusion_info$from_dnajb1,
#     from_prkaca = fusion_info$from_prkaca,
#     fusion_peptide_type = case_when(
#       spans_junction ~ "Junction-spanning",
#       from_dnajb1 ~ "DNAJB1 part",
#       from_prkaca ~ "PRKACA part",
#       TRUE ~ "Not from fusion"
#     )
#   ) %>%
#   select(-fusion_info)
# 
# # Create a separate analysis specifically for fusion peptides
# fusion_peptides_analysis <- immunopeptidome_analysis %>%
#   filter(from_fusion) %>%
#   arrange(desc(spans_junction), desc(log2_fold_change_immuno))
# 
# # Add fusion information to the combined analysis
# combined_analysis <- combined_analysis %>%
#   rowwise() %>%
#   mutate(
#     fusion_info = list(is_fusion_junction_peptide(Peptide)),
#     from_fusion = fusion_info$from_fusion,
#     spans_junction = fusion_info$spans_junction,
#     from_dnajb1 = fusion_info$from_dnajb1,
#     from_prkaca = fusion_info$from_prkaca,
#     fusion_peptide_type = case_when(
#       spans_junction ~ "Junction-spanning",
#       from_dnajb1 ~ "DNAJB1 part",
#       from_prkaca ~ "PRKACA part",
#       TRUE ~ "Not from fusion"
#     )
#   ) %>%
#   select(-fusion_info)
# 
# # Also add this information to the combined_analysis_final
# combined_analysis_final <- combined_analysis_final %>%
#   left_join(
#     combined_analysis %>% 
#       select(Peptide, from_fusion, spans_junction, from_dnajb1, from_prkaca, fusion_peptide_type),
#     by = "Peptide"
#   )
# 
# # Create fusion peptides with transcriptome analysis
# fusion_peptides_with_transcriptome <- combined_analysis_final %>%
#   filter(from_fusion) %>%
#   arrange(desc(spans_junction), desc(log2_fold_change_immuno))
# 
# # Find fusion peptides upregulated in both immunopeptidome and transcriptome
# fusion_upregulated_both <- fusion_peptides_with_transcriptome %>%
#   filter(
#     !is.na(log2_fold_change_transcriptome),
#     log2_fold_change_immuno > 1,
#     log2_fold_change_transcriptome > 1
#   ) %>%
#   arrange(desc(log2_fold_change_immuno))
# 
# #--------------------------------------------------
# # Update Excel output to include fusion peptide sheets
# #--------------------------------------------------
# 
# # Add fusion peptide sheets
# excel_sheets[["Fusion_Peptides"]] <- fusion_peptides_analysis %>% 
#   select(-total_intensity_148N_adj, -total_intensity_148T_adj)
# 
# if(nrow(fusion_peptides_with_transcriptome) > 0) {
#   excel_sheets[["Fusion_Peptides_With_Transcriptome"]] <- fusion_peptides_with_transcriptome
# }
# 
# if(nrow(fusion_upregulated_both) > 0) {
#   excel_sheets[["Fusion_Upregulated_Both"]] <- fusion_upregulated_both
# }
# 
# #--------------------------------------------------
# # Add fusion peptide visualizations
# #--------------------------------------------------
# 
# # Heatmap of fusion peptides (if any found)
# if (nrow(fusion_peptides_analysis) > 0) {
#   # Create a matrix for the heatmap
#   fusion_intensity_data <- fusion_peptides_analysis %>%
#     select(Peptide, total_intensity_148T, total_intensity_148N, fusion_peptide_type, spans_junction) %>%
#     pivot_longer(
#       cols = c(total_intensity_148T, total_intensity_148N),
#       names_to = "Sample",
#       values_to = "Intensity"
#     ) %>%
#     mutate(Sample = gsub("total_intensity_", "", Sample)) %>%
#     pivot_wider(
#       names_from = Sample,
#       values_from = Intensity
#     ) %>%
#     arrange(desc(spans_junction), Peptide)
#   
#   # Log transform the values
#   fusion_matrix <- as.matrix(fusion_intensity_data[, c("148T", "148N")])
#   rownames(fusion_matrix) <- fusion_intensity_data$Peptide
#   log_fusion_matrix <- log10(fusion_matrix + 1)
#   
#   # Create annotation for the rows
#   row_annotation <- data.frame(
#     Peptide_Type = fusion_intensity_data$fusion_peptide_type,
#     row.names = fusion_intensity_data$Peptide
#   )
#   
#   # Create a PDF of the heatmap
#   pdf(file.path(output_dir, "148T_vs_148N_fusion_peptides_heatmap.pdf"), 
#       width = 10, height = max(8, nrow(fusion_matrix)/3))
#   pheatmap(
#     log_fusion_matrix,
#     main = "Intensity of DNAJB1_PRKACA Fusion Peptides in 148T vs 148N (log10)",
#     color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
#     cluster_rows = FALSE,
#     cluster_cols = FALSE,
#     annotation_row = row_annotation,
#     display_numbers = TRUE,
#     number_format = "%.1f",
#     fontsize_row = 10,
#     fontsize_col = 10
#   )
#   dev.off()
#   
#   # Create a PNG of the heatmap
#   png(file.path(output_dir, "148T_vs_148N_fusion_peptides_heatmap.png"), 
#       width = 800, height = max(600, nrow(fusion_matrix)*40), res = 100)
#   pheatmap(
#     log_fusion_matrix,
#     main = "Intensity of DNAJB1_PRKACA Fusion Peptides in 148T vs 148N (log10)",
#     color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
#     cluster_rows = FALSE,
#     cluster_cols = FALSE,
#     annotation_row = row_annotation,
#     display_numbers = TRUE,
#     number_format = "%.1f",
#     fontsize_row = 10,
#     fontsize_col = 10
#   )
#   dev.off()
#   
#   # If we have transcriptome data for fusion peptides, create a scatter plot
#   if(nrow(fusion_peptides_with_transcriptome) > 0 && 
#      !all(is.na(fusion_peptides_with_transcriptome$log2_fold_change_transcriptome))) {
#     
#     fusion_scatter_data <- fusion_peptides_with_transcriptome %>%
#       filter(!is.na(log2_fold_change_transcriptome))
#     
#     if(nrow(fusion_scatter_data) > 0) {
#       pdf(file.path(output_dir, "Fusion_Immuno_vs_Transcriptome_Scatter.pdf"), width = 10, height = 8)
#       fusion_scatter_plot <- ggplot(fusion_scatter_data, 
#                                     aes(x = log2_fold_change_transcriptome, 
#                                         y = log2_fold_change_immuno, 
#                                         color = fusion_peptide_type,
#                                         shape = spans_junction)) +
#         geom_point(size = 3, alpha = 0.8) +
#         geom_text(aes(label = Peptide), hjust = -0.1, vjust = 0.2, size = 3) +
#         geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
#         geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
#         scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) +
#         theme_minimal() +
#         labs(
#           title = "Fusion Peptides: Immunopeptidome vs Transcriptome Fold Changes",
#           subtitle = "Triangles indicate junction-spanning peptides",
#           x = "Log2 Fold Change Transcriptome (Tumor/Normal)",
#           y = "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
#           color = "Fusion Peptide Type",
#           shape = "Spans Junction"
#         ) +
#         theme(
#           legend.position = "right",
#           plot.title = element_text(size = 14, face = "bold"),
#           plot.subtitle = element_text(size = 12)
#         )
#       print(fusion_scatter_plot)
#       dev.off()
#       
#       png(file.path(output_dir, "Fusion_Immuno_vs_Transcriptome_Scatter.png"), 
#           width = 800, height = 600, res = 100)
#       print(fusion_scatter_plot)
#       dev.off()
#     }
#   }
# }
# 
# #--------------------------------------------------
# # Update summary information to include fusion peptides
# #--------------------------------------------------
# 
# # Add fusion peptide information to summary
# if(nrow(fusion_peptides_analysis) > 0) {
#   cat("\nFusion peptides found:", nrow(fusion_peptides_analysis), "\n")
#   cat("Junction-spanning peptides:", sum(fusion_peptides_analysis$spans_junction), "\n")
#   
#   cat("\nFusion peptide detection summary:\n")
#   fusion_peptides_analysis %>%
#     group_by(fusion_peptide_type, detection_status) %>%
#     summarise(count = n(), .groups = "drop") %>%
#     arrange(fusion_peptide_type, detection_status) %>%
#     print(n = Inf)
#   
#   if(nrow(fusion_upregulated_both) > 0) {
#     cat("\nFusion peptides upregulated in both immunopeptidome and transcriptome:", 
#         nrow(fusion_upregulated_both), "\n")
#     
#     cat("\nList of fusion peptides upregulated in both datasets:\n")
#     print(fusion_upregulated_both %>% 
#             select(Peptide, fusion_peptide_type, log2_fold_change_immuno, log2_fold_change_transcriptome))
#   } else {
#     cat("\nNo fusion peptides found to be upregulated in both immunopeptidome and transcriptome.\n")
#   }
# } else {
#   cat("\nNo fusion peptides found in the analysis.\n")
# }
# 
# # Add fusion peptide visualizations to the output files list
# if(nrow(fusion_peptides_analysis) > 0) {
#   cat("5. 148T_vs_148N_fusion_peptides_heatmap.pdf/png - Heatmap of fusion peptides\n")
#   
#   if(exists("fusion_scatter_plot")) {
#     cat("6. Fusion_Immuno_vs_Transcriptome_Scatter.pdf/png - Scatter plot of fusion peptides fold changes\n")
#   }
# }
# 
# # # Modified script to analyze peptide differences and compare with transcriptome data
# # 
# # # Setting directory
# # setwd("~/Documents/Github/HLA-I_Analysis/")
# # 
# # # Create output directory in HLA-I_Analysis
# # output_dir <- "RU148_analysis"
# # if (!dir.exists(output_dir)) {
# #   dir.create(output_dir)
# #   cat("Created output directory:", output_dir, "\n")
# # }
# # 
# # # Load required packages
# # library(tidyverse)
# # library(ggplot2)
# # library(pheatmap)
# # library(writexl)   # For Excel output
# # library(openxlsx)  # For better Excel formatting
# # library(readxl)    # For reading Excel files
# # 
# # # Define the path to your data files
# # data_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/imp_014_rawdata"
# # transcriptome_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/Normalized_Gene_counts_FLCdb_Panel_1.xlsx"
# # 
# # #--------------------------------------------------
# # # PART 1: Process the immunopeptidome data
# # #--------------------------------------------------
# # 
# # # Read all TSV files in the directory and extract sample IDs from filenames
# # files <- list.files(path = data_path, pattern = ".*_peptides\\.tsv$|.*_untargeted_peptide\\.tsv$", full.names = TRUE)
# # 
# # # Filter for only 148T and 148N files
# # tumor_normal_files <- files[grepl("148[TN]", files)]
# # 
# # if (length(tumor_normal_files) == 0) {
# #   stop("No 148T or 148N peptide files found in ", data_path)
# # }
# # 
# # # Split files into 2CV and 3CV categories
# # files_2cv <- tumor_normal_files[grepl("2CV", tumor_normal_files)]
# # files_3cv <- tumor_normal_files[grepl("3CV", tumor_normal_files)]
# # 
# # cat("Found", length(files_2cv), "2CV files and", length(files_3cv), "3CV files for analysis\n")
# # 
# # # Function to read and process immunopeptidome files
# # process_immunopeptidome_files <- function(file_list) {
# #   all_data <- list()
# #   
# #   for (file in file_list) {
# #     filename <- basename(file)
# #     
# #     # Extract sample ID from filename (148T or 148N)
# #     sample_id <- ifelse(grepl("148T", filename), "148T", "148N")
# #     
# #     cat("Reading file:", filename, "- Sample ID:", sample_id, "\n")
# #     
# #     # Read the file
# #     data <- read.delim(file, stringsAsFactors = FALSE)
# #     
# #     # Add a column for sample ID
# #     data$Sample_ID <- sample_id
# #     data$Filename <- filename
# #     
# #     # Add to our list
# #     all_data[[filename]] <- data
# #   }
# #   
# #   return(bind_rows(all_data))
# # }
# # 
# # # Process 2CV and 3CV files separately
# # data_2cv <- process_immunopeptidome_files(files_2cv)
# # data_3cv <- process_immunopeptidome_files(files_3cv)
# # 
# # # Combine 2CV and 3CV data
# # combined_data <- bind_rows(
# #   data_2cv %>% mutate(CV_type = "2CV"),
# #   data_3cv %>% mutate(CV_type = "3CV")
# # )
# # 
# # # Create a summary of peptide detection for tumor and normal
# # peptide_summary <- combined_data %>%
# #   group_by(Sample_ID, Peptide) %>%
# #   summarize(
# #     peptide_length = first(nchar(Peptide)),
# #     spectral_count = sum(Spectral.Count),
# #     total_intensity = sum(Intensity),
# #     protein_ids = paste(unique(Protein.ID), collapse = "; "),
# #     genes = paste(unique(Gene), collapse = "; "),
# #     source_filenames = paste(unique(Filename), collapse = "; "),
# #     .groups = "drop"
# #   ) %>%
# #   arrange(Peptide, Sample_ID)
# # 
# # # Create a wide format table with tumor and normal side by side
# # peptide_comparison <- peptide_summary %>%
# #   select(Sample_ID, Peptide, peptide_length, spectral_count, total_intensity, protein_ids, genes) %>%
# #   pivot_wider(
# #     names_from = Sample_ID,
# #     values_from = c(spectral_count, total_intensity, protein_ids, genes),
# #     values_fill = list(spectral_count = 0, total_intensity = 0, protein_ids = NA, genes = NA)
# #   )
# # 
# # # Calculate fold changes and identify tumor-specific and normal-specific peptides
# # immunopeptidome_analysis <- peptide_comparison %>%
# #   mutate(
# #     # Replace zero with small value to prevent division by zero or Inf
# #     total_intensity_148N_adj = ifelse(total_intensity_148N == 0, 0.1, total_intensity_148N),
# #     total_intensity_148T_adj = ifelse(total_intensity_148T == 0, 0.1, total_intensity_148T),
# #     
# #     # Calculate fold changes (log2)
# #     log2_fold_change_immuno = log2(total_intensity_148T_adj / total_intensity_148N_adj),
# #     
# #     # Determine if peptide is specific to tumor or normal
# #     detection_status = case_when(
# #       total_intensity_148T > 0 & total_intensity_148N == 0 ~ "Tumor-specific",
# #       total_intensity_148N > 0 & total_intensity_148T == 0 ~ "Normal-specific",
# #       total_intensity_148T > 0 & total_intensity_148N > 0 ~ "Detected in both",
# #       TRUE ~ "Not detected"
# #     ),
# #     
# #     # Add peptide length
# #     peptide_length = nchar(Peptide),
# #     
# #     # Simplified category for plotting
# #     peptide_category = case_when(
# #       log2_fold_change_immuno > 1 ~ "Up in Tumor (FC > 2)",
# #       log2_fold_change_immuno < -1 ~ "Down in Tumor (FC < 0.5)",
# #       TRUE ~ "Similar (-1 < log2FC < 1)"
# #     )
# #   ) %>%
# #   # Clean up protein and gene info
# #   mutate(
# #     genes_combined = coalesce(genes_148T, genes_148N),
# #     proteins_combined = coalesce(protein_ids_148T, protein_ids_148N)
# #   ) %>%
# #   # Extract primary gene for later comparison with transcriptome
# #   mutate(
# #     primary_gene = sapply(strsplit(genes_combined, ";\\s*"), function(x) trimws(x[1]))
# #   ) %>%
# #   # Sort by fold change for easier viewing
# #   arrange(desc(log2_fold_change_immuno))
# # 
# # #--------------------------------------------------
# # # PART 2: Process the transcriptome data
# # #--------------------------------------------------
# # 
# # # First, examine the structure of the transcriptome data
# # cat("Reading transcriptome data from:", transcriptome_path, "\n")
# # transcriptome_data <- read_excel(transcriptome_path)
# # 
# # # Print column names to determine what's available
# # cat("Transcriptome data columns:", paste(colnames(transcriptome_data), collapse=", "), "\n")
# # 
# # # Process transcriptome data based on actual column names
# # # We'll be more flexible with column naming and existence
# # transcriptome_processed <- transcriptome_data %>%
# #   # Create default placeholder columns if they don't exist
# #   mutate(
# #     Mean.Normal = NA_real_,
# #     Mean.Tumor = NA_real_
# #   )
# # 
# # # Check if RU148 columns exist and compute averages if they do
# # if(all(c("RU148_T8", "RU148_T11") %in% colnames(transcriptome_data))) {
# #   transcriptome_processed <- transcriptome_processed %>%
# #     mutate(
# #       RU148_T_Average = (RU148_T8 + RU148_T11) / 2
# #     )
# # } else {
# #   # If columns don't exist, create placeholder
# #   transcriptome_processed$RU148_T_Average <- NA_real_
# #   cat("Warning: RU148_T8 and/or RU148_T11 columns not found in transcriptome data\n")
# # }
# # 
# # # Calculate log2 fold change if possible
# # if(all(c("RU148_N", "RU148_T_Average") %in% colnames(transcriptome_processed)) && 
# #    !all(is.na(transcriptome_processed$RU148_N)) && 
# #    !all(is.na(transcriptome_processed$RU148_T_Average))) {
# #   
# #   transcriptome_processed <- transcriptome_processed %>%
# #     mutate(
# #       log2_fold_change_transcriptome = log2(
# #         ifelse(RU148_T_Average == 0, 0.1, RU148_T_Average) / 
# #           ifelse(RU148_N == 0, 0.1, RU148_N)
# #       )
# #     )
# # } else {
# #   transcriptome_processed$log2_fold_change_transcriptome <- NA_real_
# #   cat("Warning: Unable to calculate transcriptome log2 fold change\n")
# # }
# # 
# # # Ensure we have a symbol column for joining
# # if("symbol" %in% colnames(transcriptome_processed)) {
# #   cat("Using 'symbol' column for joining\n")
# # } else if("gene_symbol" %in% colnames(transcriptome_processed)) {
# #   transcriptome_processed <- transcriptome_processed %>%
# #     rename(symbol = gene_symbol)
# #   cat("Renamed 'gene_symbol' to 'symbol' for joining\n")
# # } else if("Symbol" %in% colnames(transcriptome_processed)) {
# #   transcriptome_processed <- transcriptome_processed %>%
# #     rename(symbol = Symbol)
# #   cat("Renamed 'Symbol' to 'symbol' for joining\n")
# # } else {
# #   cat("Warning: No suitable symbol column found for joining\n")
# #   # Create an empty symbol column to avoid join errors
# #   transcriptome_processed$symbol <- NA_character_
# # }
# # 
# # #--------------------------------------------------
# # # PART 3: Combine immunopeptidome and transcriptome data
# # #--------------------------------------------------
# # 
# # # Join the immunopeptidome and transcriptome data based on gene symbol
# # # Address the many-to-many relationship warning by explicitly setting the relationship
# # combined_analysis <- immunopeptidome_analysis %>%
# #   left_join(
# #     transcriptome_processed,
# #     by = c("primary_gene" = "symbol"),
# #     relationship = "many-to-many"  # Explicitly acknowledge many-to-many relationship
# #   )
# # 
# # # Add comparison metrics if log2 fold changes are available
# # if("log2_fold_change_transcriptome" %in% colnames(combined_analysis) && 
# #    !all(is.na(combined_analysis$log2_fold_change_transcriptome))) {
# #   
# #   combined_analysis <- combined_analysis %>%
# #     mutate(
# #       immuno_trans_correlation = log2_fold_change_immuno * log2_fold_change_transcriptome,
# #       expression_category = case_when(
# #         is.na(log2_fold_change_transcriptome) ~ "No transcriptome data",
# #         log2_fold_change_transcriptome > 1 & log2_fold_change_immuno > 1 ~ "Up in both",
# #         log2_fold_change_transcriptome < -1 & log2_fold_change_immuno < -1 ~ "Down in both",
# #         log2_fold_change_transcriptome > 1 & log2_fold_change_immuno < -1 ~ "Up in transcriptome, down in immunopeptidome",
# #         log2_fold_change_transcriptome < -1 & log2_fold_change_immuno > 1 ~ "Down in transcriptome, up in immunopeptidome",
# #         TRUE ~ "No significant change"
# #       )
# #     )
# # } else {
# #   combined_analysis$immuno_trans_correlation <- NA_real_
# #   combined_analysis$expression_category <- "No transcriptome data"
# # }
# # 
# # # Select columns that definitely exist
# # combined_analysis_final <- combined_analysis %>%
# #   select(
# #     # Peptide information (these should always exist)
# #     Peptide, 
# #     peptide_length,
# #     primary_gene,
# #     genes_combined,
# #     proteins_combined,
# #     
# #     # Immunopeptidome data (these should always exist)
# #     spectral_count_148N,
# #     spectral_count_148T,
# #     total_intensity_148N,
# #     total_intensity_148T,
# #     log2_fold_change_immuno,
# #     detection_status
# #   )
# # 
# # # Add transcriptome columns if they exist
# # for(col in c("Mean.Normal", "Mean.Tumor", "RU148_N", "RU148_T8", "RU148_T11", 
# #              "RU148_T_Average", "log2_fold_change_transcriptome")) {
# #   if(col %in% colnames(combined_analysis)) {
# #     combined_analysis_final[[col]] <- combined_analysis[[col]]
# #   }
# # }
# # 
# # # Add gene information columns if they exist
# # for(col in c("geneID", "biotype", "chromosome", "gene_start", "gene_end", 
# #              "gene_length", "description")) {
# #   if(col %in% colnames(combined_analysis)) {
# #     combined_analysis_final[[col]] <- combined_analysis[[col]]
# #   }
# # }
# # 
# # # Add comparison metrics if they exist
# # for(col in c("expression_category", "immuno_trans_correlation")) {
# #   if(col %in% colnames(combined_analysis)) {
# #     combined_analysis_final[[col]] <- combined_analysis[[col]]
# #   }
# # }
# # 
# # #--------------------------------------------------
# # # PART 4: Save outputs to the new directory
# # #--------------------------------------------------
# # 
# # # Create Excel output with multiple sheets
# # excel_sheets <- list(
# #   "Combined_Analysis" = combined_analysis_final,
# #   "Immunopeptidome_Only" = immunopeptidome_analysis %>% 
# #     select(-total_intensity_148N_adj, -total_intensity_148T_adj)
# # )
# # 
# # # Add transcriptome sheet if we have data
# # if(exists("transcriptome_processed") && nrow(transcriptome_processed) > 0) {
# #   excel_sheets[["Transcriptome_Only"]] <- transcriptome_processed
# # }
# # 
# # # Add additional analysis sheets if expression category exists
# # if("expression_category" %in% colnames(combined_analysis_final)) {
# #   # Tumor-specific peptides
# #   excel_sheets[["Tumor_Specific_Peptides"]] <- combined_analysis_final %>% 
# #     filter(detection_status == "Tumor-specific")
# #   
# #   # Normal-specific peptides
# #   excel_sheets[["Normal_Specific_Peptides"]] <- combined_analysis_final %>% 
# #     filter(detection_status == "Normal-specific")
# #   
# #   # Expression category-based sheets
# #   if("Up in both" %in% unique(combined_analysis_final$expression_category)) {
# #     excel_sheets[["Up_In_Both"]] <- combined_analysis_final %>% 
# #       filter(expression_category == "Up in both") %>%
# #       arrange(desc(log2_fold_change_immuno))
# #   }
# #   
# #   if("Down in both" %in% unique(combined_analysis_final$expression_category)) {
# #     excel_sheets[["Down_In_Both"]] <- combined_analysis_final %>% 
# #       filter(expression_category == "Down in both") %>%
# #       arrange(log2_fold_change_immuno)
# #   }
# #   
# #   # Check if we have any discordant expression
# #   discordant_categories <- c("Up in transcriptome, down in immunopeptidome", 
# #                              "Down in transcriptome, up in immunopeptidome")
# #   
# #   if(any(discordant_categories %in% unique(combined_analysis_final$expression_category))) {
# #     excel_sheets[["Discordant_Expression"]] <- combined_analysis_final %>% 
# #       filter(expression_category %in% discordant_categories)
# #   }
# # }
# # 
# # # Create a more formatted Excel workbook
# # wb <- createWorkbook()
# # 
# # # Add sheets with formatting
# # for (sheet_name in names(excel_sheets)) {
# #   # Add a worksheet
# #   addWorksheet(wb, sheet_name)
# #   
# #   # Write data
# #   writeData(wb, sheet_name, excel_sheets[[sheet_name]], headerStyle = createStyle(textDecoration = "bold"))
# #   
# #   # Auto-adjust column widths
# #   setColWidths(wb, sheet_name, cols = 1:ncol(excel_sheets[[sheet_name]]), widths = "auto")
# #   
# #   # Freeze the header row
# #   freezePane(wb, sheet_name, firstRow = TRUE)
# # }
# # 
# # # Save the Excel file
# # output_file <- file.path(output_dir, "RU148_Immuno_Transcriptome_Comparison.xlsx")
# # saveWorkbook(wb, output_file, overwrite = TRUE)
# # 
# # #--------------------------------------------------
# # # PART 5: Create visualizations
# # #--------------------------------------------------
# # 
# # # Only create visualizations if we have transcriptome data
# # if("log2_fold_change_transcriptome" %in% colnames(combined_analysis_final) && 
# #    !all(is.na(combined_analysis_final$log2_fold_change_transcriptome))) {
# #   
# #   # 1. Scatter plot comparing immunopeptidome and transcriptome fold changes
# #   scatter_data <- combined_analysis_final %>%
# #     filter(!is.na(log2_fold_change_transcriptome)) %>%
# #     filter(!is.na(log2_fold_change_immuno))
# #   
# #   if(nrow(scatter_data) > 0) {
# #     pdf(file.path(output_dir, "Immuno_vs_Transcriptome_Scatter.pdf"), width = 10, height = 8)
# #     scatter_plot <- ggplot(scatter_data, 
# #                            aes(x = log2_fold_change_transcriptome, 
# #                                y = log2_fold_change_immuno, 
# #                                color = expression_category)) +
# #       geom_point(alpha = 0.7) +
# #       geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
# #       geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
# #       theme_minimal() +
# #       labs(
# #         title = "Comparison of RU148 Tumor/Normal Fold Changes",
# #         subtitle = "Immunopeptidome vs Transcriptome",
# #         x = "Log2 Fold Change Transcriptome (Tumor/Normal)",
# #         y = "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
# #         color = "Expression Category"
# #       ) +
# #       theme(
# #         legend.position = "right",
# #         plot.title = element_text(size = 14, face = "bold"),
# #         plot.subtitle = element_text(size = 12)
# #       )
# #     print(scatter_plot)
# #     dev.off()
# #     
# #     png(file.path(output_dir, "Immuno_vs_Transcriptome_Scatter.png"), width = 800, height = 600, res = 100)
# #     print(scatter_plot)
# #     dev.off()
# #   }
# #   
# #   # 2. Barplot of expression categories
# #   if("expression_category" %in% colnames(combined_analysis_final)) {
# #     expression_summary <- combined_analysis_final %>%
# #       group_by(expression_category) %>%
# #       summarise(
# #         count = n(),
# #         .groups = "drop"
# #       ) %>%
# #       arrange(desc(count))
# #     
# #     if(nrow(expression_summary) > 0) {
# #       pdf(file.path(output_dir, "Expression_Category_Counts.pdf"), width = 10, height = 6)
# #       barplot <- ggplot(expression_summary, 
# #                         aes(x = reorder(expression_category, -count), 
# #                             y = count, 
# #                             fill = expression_category)) +
# #         geom_bar(stat = "identity") +
# #         geom_text(aes(label = count), vjust = -0.5) +
# #         theme_minimal() +
# #         labs(
# #           title = "Distribution of Expression Categories",
# #           x = "Category",
# #           y = "Count",
# #           fill = "Expression Category"
# #         ) +
# #         theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #       print(barplot)
# #       dev.off()
# #       
# #       png(file.path(output_dir, "Expression_Category_Counts.png"), width = 800, height = 600, res = 100)
# #       print(barplot)
# #       dev.off()
# #     }
# #   }
# # }
# # 
# # # Print summary information
# # cat("\nSummary of RU148 Immunopeptidome and Transcriptome Analysis:\n")
# # cat("Total peptides analyzed:", nrow(immunopeptidome_analysis), "\n")
# # 
# # if(exists("transcriptome_processed")) {
# #   cat("Total genes in transcriptome:", nrow(transcriptome_processed), "\n")
# # }
# # 
# # if("log2_fold_change_transcriptome" %in% colnames(combined_analysis_final)) {
# #   cat("Peptides with matching transcriptome data:", 
# #       sum(!is.na(combined_analysis_final$log2_fold_change_transcriptome)), "\n")
# # }
# # 
# # if("expression_category" %in% colnames(combined_analysis_final)) {
# #   expression_summary <- combined_analysis_final %>%
# #     group_by(expression_category) %>%
# #     summarise(count = n(), .groups = "drop") %>%
# #     arrange(desc(count))
# #   
# #   cat("\nExpression categories:\n")
# #   print(expression_summary)
# # }
# # 
# # cat("\nAnalysis complete! Results saved to:", output_dir, "\n")
# # cat("\nThe following files were generated:\n")
# # cat("1. RU148_Immuno_Transcriptome_Comparison.xlsx - Excel file with comprehensive analysis results\n")
# # 
# # if(exists("scatter_plot")) {
# #   cat("2. Immuno_vs_Transcriptome_Scatter.pdf/png - Scatter plot comparing fold changes\n")
# # }
# # 
# # if(exists("barplot")) {
# #   cat("3. Expression_Category_Counts.pdf/png - Barplot of expression categories\n")
# # }
# # # # Script to analyze peptide differences between 148T (tumor) and 148N (normal) samples
# # # # focusing on fold change analysis
# # # 
# # # # Setting directory
# # # setwd("~/Documents/Github/HLA-I_Analysis/")
# # # #make data file in HLA-I_Analysis then save it
# # # 
# # # # Load required packages
# # # library(tidyverse)
# # # library(ggplot2)
# # # library(pheatmap)
# # # library(writexl)  # For Excel output
# # # library(openxlsx) # For better Excel formatting
# # # 
# # # # Define the path to your data files
# # # data_path <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/imp_014_rawdata"
# # # 
# # # # Read all TSV files in the directory and extract sample IDs from filenames
# # # files <- list.files(path = data_path, pattern = ".*_peptides\\.tsv$|.*_untargeted_peptide\\.tsv$", full.names = TRUE)
# # # 
# # # # Filter for only 148T and 148N files
# # # tumor_normal_files <- files[grepl("148[TN]", files)]
# # # 
# # # if (length(tumor_normal_files) == 0) {
# # #   stop("No 148T or 148N peptide files found in ", data_path)
# # # }
# # # 
# # # cat("Found", length(tumor_normal_files), "tumor/normal files for analysis\n")
# # # 
# # # # Read the 148T and 148N files
# # # all_data <- list()
# # # sample_ids <- c()
# # # 
# # # for (file in tumor_normal_files) {
# # #   filename <- basename(file)
# # #   
# # #   # Extract sample ID from filename (148T or 148N)
# # #   sample_id <- ifelse(grepl("148T", filename), "148T", "148N")
# # #   
# # #   cat("Reading file:", filename, "- Sample ID:", sample_id, "\n")
# # #   
# # #   # Read the file
# # #   data <- read.delim(file, stringsAsFactors = FALSE)
# # #   
# # #   # Add a column for sample ID
# # #   data$Sample_ID <- sample_id
# # #   data$Filename <- filename
# # #   
# # #   # Add to our list
# # #   all_data[[sample_id]] <- data
# # #   sample_ids <- c(sample_ids, sample_id)
# # # }
# # # 
# # # # Get unique sample IDs (should be 148T and 148N)
# # # sample_ids <- unique(sample_ids)
# # # cat("Found data for", length(sample_ids), "unique sample IDs:", paste(sample_ids, collapse = ", "), "\n")
# # # 
# # # # Ensure we have both tumor and normal samples
# # # if (!all(c("148T", "148N") %in% sample_ids)) {
# # #   stop("Missing either tumor (148T) or normal (148N) sample data")
# # # }
# # # 
# # # # Combine all data frames
# # # combined_data <- bind_rows(all_data)
# # # 
# # # # Create a summary of peptide detection for tumor and normal
# # # peptide_summary <- combined_data %>%
# # #   group_by(Sample_ID, Peptide) %>%
# # #   summarize(
# # #     peptide_length = first(nchar(Peptide)),
# # #     spectral_count = sum(Spectral.Count),
# # #     total_intensity = sum(Intensity),
# # #     protein_ids = paste(unique(Protein.ID), collapse = "; "),
# # #     genes = paste(unique(Gene), collapse = "; "),
# # #     source_filename = first(Filename),
# # #     .groups = "drop"
# # #   ) %>%
# # #   arrange(Peptide, Sample_ID)
# # # 
# # # # Create a wide format table with tumor and normal side by side
# # # peptide_comparison <- peptide_summary %>%
# # #   select(Sample_ID, Peptide, peptide_length, spectral_count, total_intensity, protein_ids, genes) %>%
# # #   pivot_wider(
# # #     names_from = Sample_ID,
# # #     values_from = c(spectral_count, total_intensity, protein_ids, genes),
# # #     values_fill = list(spectral_count = 0, total_intensity = 0, protein_ids = NA, genes = NA)
# # #   )
# # # 
# # # # Calculate fold changes and identify tumor-specific and normal-specific peptides
# # # peptide_analysis <- peptide_comparison %>%
# # #   mutate(
# # #     # Replace zero with small value to prevent division by zero or Inf
# # #     total_intensity_148N_adj = ifelse(total_intensity_148N == 0, 0.1, total_intensity_148N),
# # #     total_intensity_148T_adj = ifelse(total_intensity_148T == 0, 0.1, total_intensity_148T),
# # #     
# # #     # Calculate fold changes (log2)
# # #     log2_fold_change = log2(total_intensity_148T_adj / total_intensity_148N_adj),
# # #     
# # #     # Determine if peptide is specific to tumor or normal
# # #     detection_status = case_when(
# # #       total_intensity_148T > 0 & total_intensity_148N == 0 ~ "Tumor-specific",
# # #       total_intensity_148N > 0 & total_intensity_148T == 0 ~ "Normal-specific",
# # #       total_intensity_148T > 0 & total_intensity_148N > 0 ~ "Detected in both",
# # #       TRUE ~ "Not detected"
# # #     ),
# # #     
# # #     # Add peptide length
# # #     peptide_length = nchar(Peptide),
# # #     
# # #     # Simplified category for plotting
# # #     peptide_category = case_when(
# # #       log2_fold_change > 1 ~ "Up in Tumor (FC > 2)",
# # #       log2_fold_change < -1 ~ "Down in Tumor (FC < 0.5)",
# # #       TRUE ~ "Similar (-1 < log2FC < 1)"
# # #     )
# # #   ) %>%
# # #   # Clean up protein and gene info
# # #   mutate(
# # #     genes_combined = coalesce(genes_148T, genes_148N),
# # #     proteins_combined = coalesce(protein_ids_148T, protein_ids_148N)
# # #   ) %>%
# # #   # Sort by fold change for easier viewing
# # #   arrange(desc(log2_fold_change))
# # # 
# # # # Identify DNAJB1-PRKACA fusion peptides
# # # # Define the sequences of DNAJB1 and PRKACA parts of the fusion protein
# # # dnajb1_seq <- "GKDYYQTLGLARGASDEEIKRAYRRQALRYHPDKNKEPGAEEKFKEIAEAYDVLSDPRKREIFDRYGEE"
# # # prkaca_seq <- "VKEFLAKAKEDFLKKWESPAQNTAHLDQFERIKTLGTGSFGRVMLVKHKETGNHYAMKILDKQKVVKLKQIEHTLNEKRILQAVNFPFLVKLEFSFKDNSNLYMVMEYVPGGEMFSHLRRIGRFSEPHARFYAAQIVLTFEYLHSLDLIYRDLKPENLLIDQQGYIQVTDFGFAKRVKGRTWTLCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFATTDWIAIYQRKVEAPFIPKFKGPGDTSNFDDYEEEEIRVSINEKCGKEFSEF"
# # # fusion_protein <- paste0(dnajb1_seq, prkaca_seq)
# # # 
# # # # Function to check if a peptide spans the fusion junction
# # # is_fusion_junction_peptide <- function(peptide_seq) {
# # #   # Define the end of DNAJB1 and start of PRKACA for the fusion
# # #   dnajb1_end <- "IFDRYGEE"
# # #   prkaca_start <- "VKEFLAK"
# # #   
# # #   # Check if the peptide spans the fusion junction
# # #   spans_junction <- FALSE
# # #   
# # #   if (nchar(peptide_seq) >= 6) { # Require at least 6 AA to be meaningful
# # #     for (i in 3:(nchar(peptide_seq) - 3)) { # At least 3 AA from each side
# # #       left_part <- substr(peptide_seq, 1, i)
# # #       right_part <- substr(peptide_seq, i + 1, nchar(peptide_seq))
# # #       
# # #       # Check if left part is in DNAJB1 (at the end) and right part in PRKACA (at the beginning)
# # #       if (grepl(left_part, dnajb1_seq, fixed = TRUE) && 
# # #           grepl(right_part, prkaca_seq, fixed = TRUE)) {
# # #         
# # #         # Additional check to ensure left part aligns with end of DNAJB1
# # #         left_pos <- gregexpr(left_part, dnajb1_seq, fixed = TRUE)[[1]]
# # #         if (length(left_pos) > 0 && any(left_pos + nchar(left_part) - 1 >= nchar(dnajb1_seq) - 10)) {
# # #           
# # #           # Additional check to ensure right part aligns with start of PRKACA
# # #           right_pos <- gregexpr(right_part, prkaca_seq, fixed = TRUE)[[1]]
# # #           if (length(right_pos) > 0 && any(right_pos <= 10)) {
# # #             spans_junction <- TRUE
# # #             break
# # #           }
# # #         }
# # #       }
# # #     }
# # #   }
# # #   
# # #   # Check if peptide is from either part of the fusion protein
# # #   from_dnajb1 <- grepl(peptide_seq, dnajb1_seq, fixed = TRUE)
# # #   from_prkaca <- grepl(peptide_seq, prkaca_seq, fixed = TRUE)
# # #   
# # #   return(list(
# # #     spans_junction = spans_junction,
# # #     from_dnajb1 = from_dnajb1,
# # #     from_prkaca = from_prkaca,
# # #     from_fusion = from_dnajb1 | from_prkaca | spans_junction
# # #   ))
# # # }
# # # 
# # # # Add fusion protein information to the peptide analysis
# # # peptide_analysis <- peptide_analysis %>%
# # #   rowwise() %>%
# # #   mutate(
# # #     fusion_info = list(is_fusion_junction_peptide(Peptide)),
# # #     from_fusion = fusion_info$from_fusion,
# # #     spans_junction = fusion_info$spans_junction,
# # #     from_dnajb1 = fusion_info$from_dnajb1,
# # #     from_prkaca = fusion_info$from_prkaca,
# # #     fusion_peptide_type = case_when(
# # #       spans_junction ~ "Junction-spanning",
# # #       from_dnajb1 ~ "DNAJB1 part",
# # #       from_prkaca ~ "PRKACA part",
# # #       TRUE ~ "Not from fusion"
# # #     )
# # #   ) %>%
# # #   select(-fusion_info)
# # # 
# # # # Create a separate analysis specifically for fusion peptides
# # # fusion_peptides_analysis <- peptide_analysis %>%
# # #   filter(from_fusion) %>%
# # #   arrange(desc(spans_junction), desc(log2_fold_change))
# # # 
# # # # Create Excel output with multiple sheets
# # # excel_sheets <- list(
# # #   "All_Peptides_Analysis" = peptide_analysis %>% 
# # #     select(-total_intensity_148N_adj, -total_intensity_148T_adj),
# # #   "Fusion_Peptides" = fusion_peptides_analysis %>% 
# # #     select(-total_intensity_148N_adj, -total_intensity_148T_adj),
# # #   "Tumor_Specific_Peptides" = peptide_analysis %>% 
# # #     filter(detection_status == "Tumor-specific") %>% 
# # #     select(-total_intensity_148N_adj, -total_intensity_148T_adj),
# # #   "Normal_Specific_Peptides" = peptide_analysis %>% 
# # #     filter(detection_status == "Normal-specific") %>% 
# # #     select(-total_intensity_148N_adj, -total_intensity_148T_adj),
# # #   "Up_In_Tumor" = peptide_analysis %>% 
# # #     filter(log2_fold_change > 1) %>% 
# # #     arrange(desc(log2_fold_change)) %>% 
# # #     select(-total_intensity_148N_adj, -total_intensity_148T_adj),
# # #   "Down_In_Tumor" = peptide_analysis %>% 
# # #     filter(log2_fold_change < -1) %>% 
# # #     arrange(log2_fold_change) %>% 
# # #     select(-total_intensity_148N_adj, -total_intensity_148T_adj)
# # # )
# # # 
# # # # Create a more formatted Excel workbook
# # # wb <- createWorkbook()
# # # 
# # # # Add sheets with formatting
# # # for (sheet_name in names(excel_sheets)) {
# # #   # Add a worksheet
# # #   addWorksheet(wb, sheet_name)
# # #   
# # #   # Write data
# # #   writeData(wb, sheet_name, excel_sheets[[sheet_name]], headerStyle = createStyle(textDecoration = "bold"))
# # #   
# # #   # Auto-adjust column widths
# # #   setColWidths(wb, sheet_name, cols = 1:ncol(excel_sheets[[sheet_name]]), widths = "auto")
# # #   
# # #   # Freeze the header row
# # #   freezePane(wb, sheet_name, firstRow = TRUE)
# # # }
# # # 
# # # # Save the Excel file
# # # saveWorkbook(wb, file.path(data_path, "148T_vs_148N_Peptide_Analysis.xlsx"), overwrite = TRUE)
# # # 
# # # # Write simple Excel file as backup
# # # write_xlsx(excel_sheets, path = file.path(data_path, "148T_vs_148N_Peptide_Analysis_simple.xlsx"))
# # # 
# # # # Create visualizations
# # # 
# # # # 1. Volcano plot of all peptides
# # # volcano_data <- peptide_analysis %>%
# # #   filter(!is.na(log2_fold_change))
# # # 
# # # pdf(file.path(data_path, "148T_vs_148N_volcano_plot.pdf"), width = 10, height = 8)
# # # volcano_plot <- ggplot(volcano_data, aes(x = log2_fold_change, y = -log10(0.05), color = peptide_category)) +
# # #   geom_point(alpha = 0.7) +
# # #   geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +
# # #   scale_color_manual(values = c("Up in Tumor (FC > 2)" = "red", 
# # #                                 "Down in Tumor (FC < 0.5)" = "blue", 
# # #                                 "Similar (-1 < log2FC < 1)" = "gray")) +
# # #   theme_minimal() +
# # #   labs(
# # #     title = "Volcano Plot of Peptides in 148T vs 148N",
# # #     subtitle = "Red: Upregulated in Tumor, Blue: Downregulated in Tumor",
# # #     x = "Log2 Fold Change (Tumor/Normal)",
# # #     y = "-log10(p-value) [placeholder]",
# # #     color = "Peptide Category"
# # #   ) +
# # #   theme(
# # #     legend.position = "right",
# # #     plot.title = element_text(size = 14, face = "bold"),
# # #     plot.subtitle = element_text(size = 12)
# # #   )
# # # print(volcano_plot)
# # # dev.off()
# # # 
# # # png(file.path(data_path, "148T_vs_148N_volcano_plot.png"), width = 800, height = 600, res = 100)
# # # print(volcano_plot)
# # # dev.off()
# # # 
# # # # 2. Histogram of fold changes
# # # pdf(file.path(data_path, "148T_vs_148N_fold_change_histogram.pdf"), width = 10, height = 6)
# # # hist_plot <- ggplot(volcano_data, aes(x = log2_fold_change, fill = peptide_category)) +
# # #   geom_histogram(bins = 50, color = "black", alpha = 0.7) +
# # #   scale_fill_manual(values = c("Up in Tumor (FC > 2)" = "red", 
# # #                                "Down in Tumor (FC < 0.5)" = "blue", 
# # #                                "Similar (-1 < log2FC < 1)" = "gray")) +
# # #   theme_minimal() +
# # #   labs(
# # #     title = "Distribution of Peptide Fold Changes in 148T vs 148N",
# # #     x = "Log2 Fold Change (Tumor/Normal)",
# # #     y = "Count",
# # #     fill = "Peptide Category"
# # #   )
# # # print(hist_plot)
# # # dev.off()
# # # 
# # # png(file.path(data_path, "148T_vs_148N_fold_change_histogram.png"), width = 800, height = 600, res = 100)
# # # print(hist_plot)
# # # dev.off()
# # # 
# # # # 3. Barplot of detection status
# # # detection_summary <- peptide_analysis %>%
# # #   group_by(detection_status) %>%
# # #   summarise(
# # #     count = n(),
# # #     .groups = "drop"
# # #   ) %>%
# # #   arrange(desc(count))
# # # 
# # # pdf(file.path(data_path, "148T_vs_148N_detection_status.pdf"), width = 8, height = 6)
# # # detection_plot <- ggplot(detection_summary, aes(x = detection_status, y = count, fill = detection_status)) +
# # #   geom_bar(stat = "identity") +
# # #   geom_text(aes(label = count), vjust = -0.5) +
# # #   theme_minimal() +
# # #   labs(
# # #     title = "Peptide Detection Status in 148T vs 148N",
# # #     x = "Detection Status",
# # #     y = "Count",
# # #     fill = "Status"
# # #   ) +
# # #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# # # print(detection_plot)
# # # dev.off()
# # # 
# # # png(file.path(data_path, "148T_vs_148N_detection_status.png"), width = 800, height = 600, res = 100)
# # # print(detection_plot)
# # # dev.off()
# # # 
# # # # 4. Heatmap of fusion peptides (if any found)
# # # if (nrow(fusion_peptides_analysis) > 0) {
# # #   # Create a matrix for the heatmap
# # #   fusion_intensity_data <- fusion_peptides_analysis %>%
# # #     select(Peptide, total_intensity_148T, total_intensity_148N, fusion_peptide_type, spans_junction) %>%
# # #     pivot_longer(
# # #       cols = c(total_intensity_148T, total_intensity_148N),
# # #       names_to = "Sample",
# # #       values_to = "Intensity"
# # #     ) %>%
# # #     mutate(Sample = gsub("total_intensity_", "", Sample)) %>%
# # #     pivot_wider(
# # #       names_from = Sample,
# # #       values_from = Intensity
# # #     ) %>%
# # #     arrange(desc(spans_junction), Peptide)
# # #   
# # #   # Log transform the values
# # #   fusion_matrix <- as.matrix(fusion_intensity_data[, c("148T", "148N")])
# # #   rownames(fusion_matrix) <- fusion_intensity_data$Peptide
# # #   log_fusion_matrix <- log10(fusion_matrix + 1)
# # #   
# # #   # Create annotation for the rows
# # #   row_annotation <- data.frame(
# # #     Peptide_Type = fusion_intensity_data$fusion_peptide_type,
# # #     row.names = fusion_intensity_data$Peptide
# # #   )
# # #   
# # #   # Create a PDF of the heatmap
# # #   pdf(file.path(data_path, "148T_vs_148N_fusion_peptides_heatmap.pdf"), width = 10, height = max(8, nrow(fusion_matrix)/3))
# # #   pheatmap(
# # #     log_fusion_matrix,
# # #     main = "Intensity of DNAJB1_PRKACA Fusion Peptides in 148T vs 148N (log10)",
# # #     color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
# # #     cluster_rows = FALSE,
# # #     cluster_cols = FALSE,
# # #     annotation_row = row_annotation,
# # #     display_numbers = TRUE,
# # #     number_format = "%.1f",
# # #     fontsize_row = 10,
# # #     fontsize_col = 10
# # #   )
# # #   dev.off()
# # #   
# # #   # Create a PNG of the heatmap
# # #   png(file.path(data_path, "148T_vs_148N_fusion_peptides_heatmap.png"), width = 800, height = max(600, nrow(fusion_matrix)*40), res = 100)
# # #   pheatmap(
# # #     log_fusion_matrix,
# # #     main = "Intensity of DNAJB1_PRKACA Fusion Peptides in 148T vs 148N (log10)",
# # #     color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
# # #     cluster_rows = FALSE,
# # #     cluster_cols = FALSE,
# # #     annotation_row = row_annotation,
# # #     display_numbers = TRUE,
# # #     number_format = "%.1f",
# # #     fontsize_row = 10,
# # #     fontsize_col = 10
# # #   )
# # #   dev.off()
# # # }
# # # 
# # # # Print summary information
# # # cat("\nSummary of 148T vs 148N peptide analysis:\n")
# # # cat("Total peptides analyzed:", nrow(peptide_analysis), "\n")
# # # cat("\nDetection status:\n")
# # # print(detection_summary)
# # # 
# # # cat("\nFold change categories:\n")
# # # peptide_analysis %>%
# # #   group_by(peptide_category) %>%
# # #   summarise(count = n(), .groups = "drop") %>%
# # #   arrange(desc(count)) %>%
# # #   print(n = Inf)
# # # 
# # # if (nrow(fusion_peptides_analysis) > 0) {
# # #   cat("\nFusion peptides found:", nrow(fusion_peptides_analysis), "\n")
# # #   cat("Junction-spanning peptides:", sum(fusion_peptides_analysis$spans_junction), "\n")
# # #   
# # #   cat("\nFusion peptide detection summary:\n")
# # #   fusion_peptides_analysis %>%
# # #     group_by(fusion_peptide_type, detection_status) %>%
# # #     summarise(count = n(), .groups = "drop") %>%
# # #     arrange(fusion_peptide_type, detection_status) %>%
# # #     print(n = Inf)
# # # } else {
# # #   cat("\nNo fusion peptides found in the analysis.\n")
# # # }
# # # 
# # # cat("\nAnalysis complete! Results saved to:", data_path, "\n")
# # # cat("\nThe following files were generated:\n")
# # # cat("1. 148T_vs_148N_Peptide_Analysis.xlsx - Excel file with comprehensive analysis results\n")
# # # cat("2. 148T_vs_148N_volcano_plot.pdf/png - Volcano plot showing peptide fold changes\n")
# # # cat("3. 148T_vs_148N_fold_change_histogram.pdf/png - Histogram of fold changes\n")
# # # cat("4. 148T_vs_148N_detection_status.pdf/png - Barplot of peptide detection status\n")
# # # if (nrow(fusion_peptides_analysis) > 0) {
# # #   cat("5. 148T_vs_148N_fusion_peptides_heatmap.pdf/png - Heatmap of fusion peptides\n")
# # # }