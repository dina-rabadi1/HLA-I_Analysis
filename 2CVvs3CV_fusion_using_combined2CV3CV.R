# Streamlined Script to analyze DNAJB1-PRKACA fusion protein peptides in 2CV vs 3CV samples
# Looking for 8-12mers that span the fusion junction
# Uses pre-processed combined dataset

# Setting directory
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/")

# Load required packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(writexl)
library(stringr)
# Try to load ggrepel - install first if not available
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggrepel)

# Define the fusion protein sequence
fusion_protein <- "RKREIFDRYGEEVKEFLAKAKEDF"

# Identify the fusion junction - between E (DNAJB1) and V (PRKACA)
dnajb1_part <- "RKREIFDRYGEE"  # Ends with YGEE (full DNAJB1 part)
prkaca_part <- "VKEFLAKAKEDF"  # Starts with VKEF (full PRKACA part)
junction_position <- nchar(dnajb1_part)  # Position after the last E in YGEE

cat("DNAJB1 part:", dnajb1_part, "\n")
cat("PRKACA part:", prkaca_part, "\n")
cat("Junction position:", junction_position, "\n")
cat("Fusion protein:", fusion_protein, "\n")
cat("Checking junction: character at junction position is", substr(fusion_protein, junction_position, junction_position), "\n")
cat("Checking junction: next character is", substr(fusion_protein, junction_position+1, junction_position+1), "\n\n")

# Generate theoretical junction-spanning peptides (8-12mers)
theoretical_peptides <- list()
peptide_lengths <- 8:12  # Looking for 8-12mers

for (length in peptide_lengths) {
  for (start_pos in 1:(nchar(fusion_protein) - length + 1)) {
    peptide <- substr(fusion_protein, start_pos, start_pos + length - 1)
    
    # Check if this peptide spans the junction
    # It spans if it includes at least one residue from both proteins
    peptide_end_pos <- start_pos + length - 1
    spans_junction <- start_pos <= junction_position && peptide_end_pos > junction_position
    
    if (spans_junction) {
      theoretical_peptides[[length(theoretical_peptides) + 1]] <- list(
        Peptide = peptide,
        Length = length,
        Start_Position = start_pos,
        End_Position = peptide_end_pos,
        DNAJB1_Part = paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
        PRKACA_Part = paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = ""),
        Visualization = paste0(
          paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
          paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = "")
        )
      )
    }
  }
}

# Convert to dataframe
theoretical_junction_peptides <- do.call(rbind, lapply(theoretical_peptides, function(x) {
  data.frame(
    Peptide = x$Peptide,
    Length = x$Length,
    Start_Position = x$Start_Position,
    End_Position = x$End_Position,
    DNAJB1_Part = x$DNAJB1_Part,
    PRKACA_Part = x$PRKACA_Part,
    Visualization = x$Visualization,
    stringsAsFactors = FALSE
  )
}))

cat("Generated", nrow(theoretical_junction_peptides), "theoretical junction-spanning peptides\n")

# Count by length
for (length in peptide_lengths) {
  count <- sum(theoretical_junction_peptides$Length == length)
  cat("  Length", length, ":", count, "peptides\n")
}

# Print all theoretical peptides in a neat table
cat("\nAll theoretical junction-spanning peptides:\n")
for (i in 1:nrow(theoretical_junction_peptides)) {
  peptide <- theoretical_junction_peptides[i,]
  cat(sprintf("%2d. %s (Length: %d, Start: %d, End: %d, Vis: %s)\n", 
              i, peptide$Peptide, peptide$Length, peptide$Start_Position, 
              peptide$End_Position, peptide$Visualization))
}

# Read the combined dataset
cat("Reading combined dataset...\n")
data_file <- "unique_peptides_unmodified.tsv"
combined_data <- read.delim(data_file, stringsAsFactors = FALSE)

cat("Dataset loaded with", nrow(combined_data), "peptide records\n")
cat("Columns available:", paste(colnames(combined_data), collapse = ", "), "\n\n")

# DIAGNOSTIC: Check sample IDs in the dataset
cat("=== DIAGNOSTIC INFORMATION ===\n")
cat("Unique Sample IDs in dataset:\n")
unique_samples <- unique(combined_data$SampleID)
cat(paste(unique_samples, collapse = ", "), "\n\n")

# Check for any patterns that might cause 51S vs 51 issue
samples_with_51 <- unique_samples[grepl("51", unique_samples)]
cat("Sample IDs containing '51':", paste(samples_with_51, collapse = ", "), "\n")

# Check if there are issues with sample ID extraction
cat("\nSample ID patterns:\n")
sample_id_table <- table(combined_data$SampleID)
print(sample_id_table)

# Check source files to understand the naming
if("SourceFile" %in% colnames(combined_data)) {
  cat("\nSource file patterns:\n")
  source_files <- unique(combined_data$SourceFile)
  files_with_51 <- source_files[grepl("51", source_files)]
  cat("Files containing '51':\n")
  for(file in files_with_51) {
    cat("  ", file, "\n")
  }
}

cat("\n=== END DIAGNOSTIC ===\n\n")

# Create output directory for results
results_dir <- "Fusion_Junction_Analysis"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

#===============================#
# Function to check if a peptide spans the fusion junction
#===============================#

is_junction_spanning <- function(peptide) {
  # First check if it's entirely within the fusion protein
  if (!grepl(peptide, fusion_protein, fixed = TRUE)) {
    return(FALSE)
  }
  
  # Find position in fusion protein
  start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
  if (start_pos == -1) {
    return(FALSE)
  }
  
  end_pos <- start_pos + nchar(peptide) - 1
  
  # Check if it spans the junction
  spans_junction <- start_pos <= junction_position && end_pos > junction_position
  
  return(spans_junction)
}

# Filter data for 8-12mers and identify junction-spanning peptides
combined_data <- combined_data %>%
  mutate(
    Peptide_Length = nchar(Peptide),
    Spans_Junction = sapply(Peptide, is_junction_spanning)
  ) %>%
  filter(Peptide_Length >= 8 & Peptide_Length <= 12)  # Only include 8-12mers

# Extract detected fusion-spanning peptides
detected_junction_peptides <- combined_data %>%
  filter(Spans_Junction == TRUE) %>%
  group_by(Peptide, Peptide_Length) %>%
  summarize(
    Total_Samples = n(),
    Detected_2CV_Count = sum(!is.na(detected_2cv) & detected_2cv == TRUE, na.rm = TRUE),
    Detected_3CV_Count = sum(!is.na(detected_3cv) & detected_3cv == TRUE, na.rm = TRUE),
    Detected_Both_Count = sum(!is.na(detected_both) & detected_both == TRUE, na.rm = TRUE),
    Avg_Final_Intensity = mean(final_intensity, na.rm = TRUE),
    Avg_Intensity_2CV = mean(Intensity_2cv, na.rm = TRUE),
    Avg_Intensity_3CV = mean(Intensity_3cv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Peptide_Length, Peptide)

cat("Found", nrow(detected_junction_peptides), "detected peptides that span the fusion junction\n")

# Get summary statistics about the dataset
total_samples <- length(unique(combined_data$SampleID))
samples_with_2cv <- sum(!is.na(combined_data$detected_2cv) & combined_data$detected_2cv == TRUE)
samples_with_3cv <- sum(!is.na(combined_data$detected_3cv) & combined_data$detected_3cv == TRUE)

cat("Dataset summary:\n")
cat("  Total unique samples:", total_samples, "\n")
cat("  Detections in 2CV:", samples_with_2cv, "\n")
cat("  Detections in 3CV:", samples_with_3cv, "\n\n")

if (nrow(detected_junction_peptides) > 0) {
  # Create detailed analysis of each detected junction peptide
  detected_details <- data.frame()
  
  for (i in 1:nrow(detected_junction_peptides)) {
    peptide <- detected_junction_peptides$Peptide[i]
    start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
    end_pos <- start_pos + nchar(peptide) - 1
    
    # Calculate how many residues come from each protein
    dnajb1_residues <- max(0, min(junction_position - start_pos + 1, nchar(peptide)))
    prkaca_residues <- max(0, min(end_pos - junction_position, nchar(peptide)))
    
    # Visualization string (D for DNAJB1, P for PRKACA)
    vis_string <- paste0(
      paste(rep("D", dnajb1_residues), collapse = ""),
      paste(rep("P", prkaca_residues), collapse = "")
    )
    
    # Add to details dataframe
    detected_details <- rbind(detected_details, data.frame(
      Peptide = peptide,
      Length = nchar(peptide),
      Start_Position = start_pos,
      End_Position = end_pos,
      DNAJB1_Residues = dnajb1_residues,
      PRKACA_Residues = prkaca_residues,
      Visualization = vis_string,
      Total_Samples = detected_junction_peptides$Total_Samples[i],
      Detected_2CV_Count = detected_junction_peptides$Detected_2CV_Count[i],
      Detected_3CV_Count = detected_junction_peptides$Detected_3CV_Count[i],
      Detected_Both_Count = detected_junction_peptides$Detected_Both_Count[i],
      Avg_Final_Intensity = detected_junction_peptides$Avg_Final_Intensity[i],
      Avg_Intensity_2CV = detected_junction_peptides$Avg_Intensity_2CV[i],
      Avg_Intensity_3CV = detected_junction_peptides$Avg_Intensity_3CV[i],
      stringsAsFactors = FALSE
    ))
  }
  
  # Create heatmaps showing detection patterns
  
  # Prepare data for heatmaps - focusing on detected junction peptides
  junction_data <- combined_data %>%
    filter(Spans_Junction == TRUE) %>%
    select(Peptide, SampleID, detected_2cv, detected_3cv, detected_both, 
           final_intensity, Intensity_2cv, Intensity_3cv)
  
  # Debug script to find why EIFDRYGEEV is missing from heatmap
  # Add this section right after you create the detected_junction_peptides dataframe
  
  cat("\n=== DEBUGGING EIFDRYGEEV PEPTIDE ===\n")
  
  # Define the missing peptide
  missing_peptide <- "EIFDRYGEEV"
  
  # 1. Check if the peptide exists in the raw dataset
  cat("1. Checking if", missing_peptide, "exists in raw dataset...\n")
  raw_matches <- combined_data[combined_data$Peptide == missing_peptide, ]
  cat("   Found", nrow(raw_matches), "raw matches for", missing_peptide, "\n")
  
  if(nrow(raw_matches) > 0) {
    cat("   Sample details:\n")
    print(raw_matches[, c("SampleID", "Peptide", "detected_2cv", "detected_3cv", "detected_both")])
  }
  
  # 2. Check if it spans the junction according to our function
  cat("\n2. Checking if", missing_peptide, "spans junction...\n")
  spans_result <- is_junction_spanning(missing_peptide)
  cat("   is_junction_spanning(", missing_peptide, ") =", spans_result, "\n")
  
  # 3. Manual junction check
  cat("\n3. Manual junction check for", missing_peptide, ":\n")
  fusion_protein <- "RKREIFDRYGEEVKEFLAKAKEDF"
  junction_position <- 12  # Position after the last E in YGEE
  
  # Find position in fusion protein
  start_pos <- regexpr(missing_peptide, fusion_protein, fixed = TRUE)[1]
  if (start_pos != -1) {
    end_pos <- start_pos + nchar(missing_peptide) - 1
    cat("   Found at position", start_pos, "to", end_pos, "\n")
    cat("   Junction position is:", junction_position, "\n")
    cat("   Spans junction?", start_pos <= junction_position && end_pos > junction_position, "\n")
    
    # Show which part is from which protein
    dnajb1_residues <- max(0, min(junction_position - start_pos + 1, nchar(missing_peptide)))
    prkaca_residues <- max(0, min(end_pos - junction_position, nchar(missing_peptide)))
    cat("   DNAJB1 residues:", dnajb1_residues, "\n")
    cat("   PRKACA residues:", prkaca_residues, "\n")
    
    # Show visual breakdown
    peptide_chars <- strsplit(missing_peptide, "")[[1]]
    protein_source <- c(rep("D", dnajb1_residues), rep("P", prkaca_residues))
    cat("   Peptide breakdown:\n")
    cat("   ", paste(peptide_chars, collapse = " "), "\n")
    cat("   ", paste(protein_source, collapse = " "), "\n")
  } else {
    cat("   ERROR: Peptide not found in fusion protein!\n")
  }
  
  # 4. Check the filtering steps
  cat("\n4. Checking filtering steps...\n")
  
  # Check length filter
  peptide_length <- nchar(missing_peptide)
  cat("   Peptide length:", peptide_length, "(should be 8-12) - PASS:", peptide_length >= 8 && peptide_length <= 12, "\n")
  
  # Check if it's in the filtered dataset
  filtered_matches <- combined_data[combined_data$Peptide == missing_peptide & 
                                      combined_data$Peptide_Length >= 8 & 
                                      combined_data$Peptide_Length <= 12, ]
  cat("   Matches after length filter:", nrow(filtered_matches), "\n")
  
  # Check if it's marked as junction spanning
  junction_matches <- combined_data[combined_data$Peptide == missing_peptide & 
                                      combined_data$Spans_Junction == TRUE, ]
  cat("   Matches marked as junction-spanning:", nrow(junction_matches), "\n")
  
  # 5. Check the grouping/summarization step
  cat("\n5. Checking grouping step...\n")
  if(nrow(junction_matches) > 0) {
    cat("   Junction matches details:\n")
    print(junction_matches[, c("SampleID", "Peptide", "detected_2cv", "detected_3cv", "detected_both", "Spans_Junction")])
    
    # Manual grouping calculation
    manual_summary <- junction_matches %>%
      group_by(Peptide, Peptide_Length) %>%
      summarize(
        Total_Samples = n(),
        Detected_2CV_Count = sum(!is.na(detected_2cv) & detected_2cv == TRUE, na.rm = TRUE),
        Detected_3CV_Count = sum(!is.na(detected_3cv) & detected_3cv == TRUE, na.rm = TRUE),
        Detected_Both_Count = sum(!is.na(detected_both) & detected_both == TRUE, na.rm = TRUE),
        .groups = "drop"
      )
    
    cat("   Manual summary:\n")
    print(manual_summary)
  }
  
  # 6. Check if it's in the detected_junction_peptides dataframe
  cat("\n6. Checking detected_junction_peptides dataframe...\n")
  detected_match <- detected_junction_peptides[detected_junction_peptides$Peptide == missing_peptide, ]
  cat("   Found in detected_junction_peptides:", nrow(detected_match), "rows\n")
  if(nrow(detected_match) > 0) {
    print(detected_match)
  }
  
  # 7. Check the heatmap data preparation
  cat("\n7. Checking heatmap data preparation...\n")
  
  # Check junction_data
  junction_data_matches <- junction_data[junction_data$Peptide == missing_peptide, ]
  cat("   Found in junction_data:", nrow(junction_data_matches), "rows\n")
  if(nrow(junction_data_matches) > 0) {
    print(junction_data_matches)
  }
  
  # Check old_style_combined
  if(exists("old_style_combined")) {
    old_style_matches <- old_style_combined[old_style_combined$Peptide == missing_peptide, ]
    cat("   Found in old_style_combined:", nrow(old_style_matches), "rows\n")
    if(nrow(old_style_matches) > 0) {
      print(old_style_matches)
    }
  }
  
  # 8. Check if the peptide has any detections at all
  cat("\n8. Summary check - does", missing_peptide, "have any TRUE detections?\n")
  if(nrow(raw_matches) > 0) {
    cat("   Any detected_2cv = TRUE?", any(raw_matches$detected_2cv == TRUE, na.rm = TRUE), "\n")
    cat("   Any detected_3cv = TRUE?", any(raw_matches$detected_3cv == TRUE, na.rm = TRUE), "\n")
    cat("   Any detected_both = TRUE?", any(raw_matches$detected_both == TRUE, na.rm = TRUE), "\n")
    
    # Count non-NA detections
    cat("   Non-NA detected_2cv count:", sum(!is.na(raw_matches$detected_2cv)), "\n")
    cat("   Non-NA detected_3cv count:", sum(!is.na(raw_matches$detected_3cv)), "\n")
    cat("   TRUE detected_2cv count:", sum(raw_matches$detected_2cv == TRUE, na.rm = TRUE), "\n")
    cat("   TRUE detected_3cv count:", sum(raw_matches$detected_3cv == TRUE, na.rm = TRUE), "\n")
  }
  
  cat("\n=== END DEBUGGING ===\n\n")
  
  # COMPREHENSIVE DEBUG: Find ALL instances of EIFDRYGEEV
  cat("\n=== COMPREHENSIVE EIFDRYGEEV DEBUG ===\n")
  
  missing_peptide <- "EIFDRYGEEV"
  
  # 1. Find ALL raw instances in the original dataset (before any filtering)
  cat("1. ALL instances of", missing_peptide, "in the ORIGINAL dataset:\n")
  all_raw_matches <- read.delim("unique_peptides_unmodified.tsv", stringsAsFactors = FALSE)
  all_eifdrygeev_raw <- all_raw_matches[all_raw_matches$Peptide == missing_peptide, ]
  cat("   Found", nrow(all_eifdrygeev_raw), "total raw instances\n")
  
  if(nrow(all_eifdrygeev_raw) > 0) {
    cat("   Raw instances details:\n")
    for(i in 1:nrow(all_eifdrygeev_raw)) {
      row <- all_eifdrygeev_raw[i, ]
      cat(sprintf("     Row %d: SampleID=%s, detected_2cv=%s, detected_3cv=%s, detected_both=%s, SourceFile=%s\n",
                  i, row$SampleID, row$detected_2cv, row$detected_3cv, row$detected_both, 
                  substr(row$SourceFile, 1, 50))) # Truncate long filenames
    }
  }
  
  # 2. Check what happens after initial data loading and filtering
  cat("\n2. After initial data loading (combined_data):\n")
  combined_eifdrygeev <- combined_data[combined_data$Peptide == missing_peptide, ]
  cat("   Found", nrow(combined_eifdrygeev), "instances after initial processing\n")
  
  if(nrow(combined_eifdrygeev) > 0) {
    cat("   After initial processing:\n")
    print(combined_eifdrygeev[, c("SampleID", "Peptide", "detected_2cv", "detected_3cv", "detected_both", "Peptide_Length", "Spans_Junction")])
  }
  
  # 3. Check what's in junction_data
  cat("\n3. In junction_data (junction-spanning peptides only):\n")
  if(exists("junction_data")) {
    junction_eifdrygeev <- junction_data[junction_data$Peptide == missing_peptide, ]
    cat("   Found", nrow(junction_eifdrygeev), "instances in junction_data\n")
    if(nrow(junction_eifdrygeev) > 0) {
      print(junction_eifdrygeev)
    }
  } else {
    cat("   junction_data does not exist yet\n")
  }
  
  # 4. Check the sample ID mapping issue
  cat("\n4. Investigating sample ID discrepancy:\n")
  cat("   You mentioned row 61230 has SampleID = 51S\n")
  cat("   But debug showed SampleID = 62\n")
  cat("   Let's check if there are multiple EIFDRYGEEV entries with different sample IDs:\n")
  
  if(nrow(all_eifdrygeev_raw) > 0) {
    sample_ids <- unique(all_eifdrygeev_raw$SampleID)
    cat("   Unique Sample IDs for EIFDRYGEEV:", paste(sample_ids, collapse = ", "), "\n")
    
    for(sid in sample_ids) {
      subset_data <- all_eifdrygeev_raw[all_eifdrygeev_raw$SampleID == sid, ]
      cat(sprintf("   Sample %s: %d rows, detected_2cv=%s, detected_3cv=%s\n", 
                  sid, nrow(subset_data), 
                  paste(unique(subset_data$detected_2cv), collapse = ","),
                  paste(unique(subset_data$detected_3cv), collapse = ",")))
    }
  }
  
  # 5. Check if there are any data processing steps that might be changing sample IDs
  cat("\n5. Checking for data processing issues:\n")
  cat("   Are there any duplicate rows being removed?\n")
  cat("   Are sample IDs being modified during processing?\n")
  
  # Check for exact row that was mentioned
  if(nrow(all_raw_matches) >= 61230) {
    cat("   Checking row 61230 specifically:\n")
    row_61230 <- all_raw_matches[61230, ]
    cat("   Row 61230 Peptide:", row_61230$Peptide, "\n")
    cat("   Row 61230 SampleID:", row_61230$SampleID, "\n")
    cat("   Row 61230 detected_2cv:", row_61230$detected_2cv, "\n")
    cat("   Row 61230 detected_3cv:", row_61230$detected_3cv, "\n")
  } else {
    cat("   Dataset has fewer than 61230 rows\n")
  }
  
  # 6. Check if the issue is in the heatmap matrix creation
  cat("\n6. Checking matrix creation logic:\n")
  if(exists("old_style_matrix_data")) {
    matrix_eifdrygeev <- old_style_matrix_data[old_style_matrix_data$Peptide == missing_peptide, ]
    cat("   EIFDRYGEEV in matrix data:", nrow(matrix_eifdrygeev), "rows\n")
    if(nrow(matrix_eifdrygeev) > 0) {
      # Check which samples have TRUE detections
      true_detections <- matrix_eifdrygeev[matrix_eifdrygeev$detected == TRUE, ]
      cat("   TRUE detections in matrix:\n")
      if(nrow(true_detections) > 0) {
        print(true_detections[, c("SampleID", "CV_Type", "Sample_CV", "detected")])
      } else {
        cat("   No TRUE detections found in matrix data!\n")
      }
    }
  }
  
  # 7. Final verification - check the heatmap row
  cat("\n7. Final heatmap verification:\n")
  if(exists("old_style_presence_numeric")) {
    if("EIFDRYGEEV" %in% rownames(old_style_presence_numeric)) {
      eifdrygeev_row <- old_style_presence_numeric["EIFDRYGEEV", ]
      cat("   EIFDRYGEEV row in heatmap:\n")
      print(eifdrygeev_row)
      
      # Check which columns have 1 (TRUE detection)
      true_cols <- names(eifdrygeev_row)[eifdrygeev_row == 1]
      cat("   Columns with TRUE detection:", paste(true_cols, collapse = ", "), "\n")
    } else {
      cat("   EIFDRYGEEV NOT FOUND in final heatmap matrix!\n")
    }
  }
  
  cat("\n=== END COMPREHENSIVE DEBUG ===\n")
  
  # TARGETED DEBUG: Find the missing EIFDRYGEEV in sample 51S
  cat("\n=== TARGETED DEBUG: MISSING 51S INSTANCE ===\n")
  
  # 1. Read the raw file again and search more comprehensively
  cat("1. Comprehensive search for EIFDRYGEEV in raw data:\n")
  raw_data <- read.delim("unique_peptides_unmodified.tsv", stringsAsFactors = FALSE)
  
  # Search for ALL instances of EIFDRYGEEV (case insensitive)
  all_eifdrygeev <- raw_data[grepl("EIFDRYGEEV", raw_data$Peptide, ignore.case = TRUE, fixed = TRUE), ]
  cat("   Total EIFDRYGEEV instances found:", nrow(all_eifdrygeev), "\n")
  
  if(nrow(all_eifdrygeev) > 0) {
    cat("   Details of ALL EIFDRYGEEV instances:\n")
    for(i in 1:nrow(all_eifdrygeev)) {
      row <- all_eifdrygeev[i, ]
      cat(sprintf("   Instance %d: Row %d, SampleID='%s', detected_2cv=%s, detected_3cv=%s, SourceFile=%s\n", 
                  i, which(raw_data$Peptide == "EIFDRYGEEV")[i], 
                  row$SampleID, row$detected_2cv, row$detected_3cv, 
                  basename(row$SourceFile)))
    }
  }
  
  # 2. Specifically search for 51S samples
  cat("\n2. Searching for ALL 51S samples:\n")
  samples_51s <- raw_data[raw_data$SampleID == "51S", ]
  cat("   Total rows with SampleID = '51S':", nrow(samples_51s), "\n")
  
  if(nrow(samples_51s) > 0) {
    # Check if any of the 51S samples contain EIFDRYGEEV
    eifdrygeev_51s <- samples_51s[samples_51s$Peptide == "EIFDRYGEEV", ]
    cat("   EIFDRYGEEV instances in 51S:", nrow(eifdrygeev_51s), "\n")
    
    if(nrow(eifdrygeev_51s) > 0) {
      cat("   FOUND 51S INSTANCE:\n")
      print(eifdrygeev_51s[, c("SampleID", "Peptide", "detected_2cv", "detected_3cv", "SourceFile")])
    } else {
      cat("   NO EIFDRYGEEV found in 51S samples\n")
      
      # Check what peptides ARE in 51S around row 61230 area
      cat("   Sample of peptides in 51S (first 10):\n")
      print(head(samples_51s$Peptide, 10))
    }
  }
  
  # 3. Check row 61230 specifically and surrounding rows
  cat("\n3. Checking row 61230 and surrounding rows:\n")
  if(nrow(raw_data) >= 61235) {
    for(row_num in 61225:61235) {
      row_data <- raw_data[row_num, ]
      cat(sprintf("   Row %d: SampleID='%s', Peptide='%s'\n", 
                  row_num, row_data$SampleID, row_data$Peptide))
    }
  } else {
    cat("   Dataset has fewer rows than 61235\n")
  }
  
  # 4. Check if there are different variations of the peptide sequence
  cat("\n4. Checking for peptide sequence variations:\n")
  # Look for similar sequences (in case there are modifications or case differences)
  similar_peptides <- raw_data[grepl("IFDRYGEEV", raw_data$Peptide, ignore.case = TRUE), ]
  cat("   Peptides containing 'IFDRYGEEV':", nrow(similar_peptides), "\n")
  if(nrow(similar_peptides) > 0) {
    unique_similar <- unique(similar_peptides$Peptide)
    cat("   Unique similar sequences:", paste(unique_similar, collapse = ", "), "\n")
  }
  
  # 5. Check if there are issues with sample ID parsing
  cat("\n5. Checking sample ID patterns:\n")
  unique_samples <- unique(raw_data$SampleID)
  samples_with_51 <- unique_samples[grepl("51", unique_samples)]
  cat("   All sample IDs containing '51':", paste(samples_with_51, collapse = ", "), "\n")
  
  # Check if "51S" might be stored differently (with spaces, different case, etc.)
  samples_similar_51s <- unique_samples[grepl("51.*S|S.*51", unique_samples, ignore.case = TRUE)]
  cat("   Sample IDs similar to '51S':", paste(samples_similar_51s, collapse = ", "), "\n")
  
  # 6. Check the source files for 51S
  cat("\n6. Checking source files:\n")
  source_files_51s <- unique(samples_51s$SourceFile)
  cat("   Source files for 51S samples:\n")
  for(file in source_files_51s) {
    cat("   ", basename(file), "\n")
  }
  
  # 7. Manual verification - read the specific area around your mentioned row
  cat("\n7. Manual verification around your specified row:\n")
  cat("   You mentioned row 61230 has EIFDRYGEEV in sample 51S\n")
  cat("   Let's check if there might be a row numbering difference:\n")
  
  # Search for the exact combination you mentioned
  eifdrygeev_in_51s_files <- raw_data[raw_data$Peptide == "EIFDRYGEEV" & 
                                        grepl("51S", raw_data$SourceFile, ignore.case = TRUE), ]
  cat("   EIFDRYGEEV in files containing '51S':", nrow(eifdrygeev_in_51s_files), "\n")
  if(nrow(eifdrygeev_in_51s_files) > 0) {
    print(eifdrygeev_in_51s_files)
  }
  
  cat("\n=== END TARGETED DEBUG ===\n")
  
  # Create presence matrices for 2CV and 3CV
  presence_2cv_data <- junction_data %>%
    filter(!is.na(detected_2cv)) %>%
    select(Peptide, SampleID, detected_2cv) %>%
    pivot_wider(names_from = SampleID, values_from = detected_2cv, values_fill = FALSE) %>%
    column_to_rownames("Peptide")
  
  presence_3cv_data <- junction_data %>%
    filter(!is.na(detected_3cv)) %>%
    select(Peptide, SampleID, detected_3cv) %>%
    pivot_wider(names_from = SampleID, values_from = detected_3cv, values_fill = FALSE) %>%
    column_to_rownames("Peptide")
  
  # Convert to numeric matrices
  if(nrow(presence_2cv_data) > 0 && ncol(presence_2cv_data) > 0) {
    presence_2cv_matrix <- as.matrix(presence_2cv_data)
    presence_2cv_numeric <- matrix(as.numeric(presence_2cv_matrix), 
                                   nrow = nrow(presence_2cv_matrix),
                                   dimnames = dimnames(presence_2cv_matrix))
    
    # Create 2CV presence heatmap
    pdf(file.path(results_dir, "junction_peptides_presence_2CV.pdf"), width = 12, height = 8)
    pheatmap(
      presence_2cv_numeric,
      main = "Presence of Junction Peptides in 2CV Samples",
      color = c("white", "steelblue"),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      number_format = "%.0f"
    )
    dev.off()
    
    png(file.path(results_dir, "junction_peptides_presence_2CV.png"), width = 1000, height = 600, res = 100)
    pheatmap(
      presence_2cv_numeric,
      main = "Presence of Junction Peptides in 2CV Samples",
      color = c("white", "steelblue"),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      number_format = "%.0f"
    )
    dev.off()
  }
  
  if(nrow(presence_3cv_data) > 0 && ncol(presence_3cv_data) > 0) {
    presence_3cv_matrix <- as.matrix(presence_3cv_data)
    presence_3cv_numeric <- matrix(as.numeric(presence_3cv_matrix), 
                                   nrow = nrow(presence_3cv_matrix),
                                   dimnames = dimnames(presence_3cv_matrix))
    
    # Create 3CV presence heatmap
    pdf(file.path(results_dir, "junction_peptides_presence_3CV.pdf"), width = 12, height = 8)
    pheatmap(
      presence_3cv_numeric,
      main = "Presence of Junction Peptides in 3CV Samples",
      color = c("white", "steelblue"),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      number_format = "%.0f"
    )
    dev.off()
    
    png(file.path(results_dir, "junction_peptides_presence_3CV.png"), width = 1000, height = 600, res = 100)
    pheatmap(
      presence_3cv_numeric,
      main = "Presence of Junction Peptides in 3CV Samples",
      color = c("white", "steelblue"),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      number_format = "%.0f"
    )
    dev.off()
  }
  
  # Create intensity heatmaps
  intensity_2cv_data <- junction_data %>%
    filter(!is.na(Intensity_2cv) & Intensity_2cv > 0) %>%
    select(Peptide, SampleID, Intensity_2cv) %>%
    pivot_wider(names_from = SampleID, values_from = Intensity_2cv, values_fill = 0) %>%
    column_to_rownames("Peptide")
  
  intensity_3cv_data <- junction_data %>%
    filter(!is.na(Intensity_3cv) & Intensity_3cv > 0) %>%
    select(Peptide, SampleID, Intensity_3cv) %>%
    pivot_wider(names_from = SampleID, values_from = Intensity_3cv, values_fill = 0) %>%
    column_to_rownames("Peptide")
  
  if(nrow(intensity_2cv_data) > 0 && ncol(intensity_2cv_data) > 0) {
    log_intensity_2cv <- log10(as.matrix(intensity_2cv_data) + 1)
    
    pdf(file.path(results_dir, "junction_peptides_intensity_2CV.pdf"), width = 12, height = 8)
    pheatmap(
      log_intensity_2cv,
      main = "Intensity of Junction Peptides in 2CV Samples (log10)",
      color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      number_format = "%.1f"
    )
    dev.off()
    
    png(file.path(results_dir, "junction_peptides_intensity_2CV.png"), width = 1000, height = 600, res = 100)
    pheatmap(
      log_intensity_2cv,
      main = "Intensity of Junction Peptides in 2CV Samples (log10)",
      color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      number_format = "%.1f"
    )
    dev.off()
  }
  
  if(nrow(intensity_3cv_data) > 0 && ncol(intensity_3cv_data) > 0) {
    log_intensity_3cv <- log10(as.matrix(intensity_3cv_data) + 1)
    
    pdf(file.path(results_dir, "junction_peptides_intensity_3CV.pdf"), width = 12, height = 8)
    pheatmap(
      log_intensity_3cv,
      main = "Intensity of Junction Peptides in 3CV Samples (log10)",
      color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      number_format = "%.1f"
    )
    dev.off()
    
    png(file.path(results_dir, "junction_peptides_intensity_3CV.png"), width = 1000, height = 600, res = 100)
    pheatmap(
      log_intensity_3cv,
      main = "Intensity of Junction Peptides in 3CV Samples (log10)",
      color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      fontsize_col = 8,
      display_numbers = TRUE,
      number_format = "%.1f"
    )
    dev.off()
  }
  
  # Compare 2CV vs 3CV detection efficiency for junction peptides
  cv_comparison_plot <- ggplot(detected_details, aes(x = Peptide, y = Length)) +
    geom_point(aes(size = Total_Samples, color = Detected_3CV_Count / (Detected_2CV_Count + 0.001))) +
    scale_color_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red", 
      midpoint = 1,
      name = "3CV/2CV\nDetection Ratio"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Fusion Junction Peptide Detection: 2CV vs 3CV",
      x = "Peptide Sequence",
      y = "Peptide Length",
      size = "Total\nSamples"
    )
  
  ggsave(file.path(results_dir, "junction_peptide_cv_comparison.pdf"), cv_comparison_plot, width = 12, height = 8)
  ggsave(file.path(results_dir, "junction_peptide_cv_comparison.png"), cv_comparison_plot, width = 12, height = 8)
  
  # Create a visualization showing peptide coverage across the fusion protein
  peptide_coverage_plot <- ggplot(detected_details, 
                                  aes(x = Start_Position, xend = End_Position, 
                                      y = reorder(Peptide, Start_Position), yend = reorder(Peptide, Start_Position))) +
    geom_segment(aes(color = Total_Samples), linewidth = 5) +
    geom_vline(xintercept = junction_position, linetype = "dashed", color = "red") +
    annotate("text", x = junction_position - 2, y = nrow(detected_details) + 1, 
             label = "DNAJB1", color = "darkgreen", fontface = "bold") +
    annotate("text", x = junction_position + 2, y = nrow(detected_details) + 1, 
             label = "PRKACA", color = "purple", fontface = "bold") +
    scale_color_gradient(low = "lightsteelblue", high = "steelblue") +
    theme_minimal() +
    labs(
      title = "Coverage of Fusion Junction by Detected Peptides",
      x = "Position in Fusion Protein",
      y = "Peptide",
      color = "Total\nSamples"
    )
  
  ggsave(file.path(results_dir, "junction_peptide_coverage.pdf"), peptide_coverage_plot, width = 12, height = 8)
  ggsave(file.path(results_dir, "junction_peptide_coverage.png"), peptide_coverage_plot, width = 12, height = 8)
  
  # FIXED VERSION: Create OLD SCRIPT STYLE combined heatmap for direct comparison
  # The key fix: include ALL junction peptides regardless of detection status
  
  # Get all unique samples that have junction peptides
  samples_with_junction_peptides <- junction_data %>%
    filter(!is.na(SampleID)) %>%
    pull(SampleID) %>%
    unique()
  
  cat("Samples with junction peptides:", length(samples_with_junction_peptides), "\n")
  cat("Sample IDs:", paste(samples_with_junction_peptides, collapse = ", "), "\n")
  
  if(length(samples_with_junction_peptides) > 0 && nrow(detected_details) > 0) {
    
    # FIXED: Include ALL junction peptides, not just those with TRUE detections
    old_style_data_2cv <- junction_data %>%
      filter(!is.na(detected_2cv)) %>%  # Just need non-NA, not necessarily TRUE
      mutate(Sample_CV = paste0(SampleID, "_2CV")) %>%
      select(Peptide, Sample_CV, detected = detected_2cv, Intensity = Intensity_2cv)
    
    old_style_data_3cv <- junction_data %>%
      filter(!is.na(detected_3cv)) %>%  # Just need non-NA, not necessarily TRUE
      mutate(Sample_CV = paste0(SampleID, "_3CV")) %>%
      select(Peptide, Sample_CV, detected = detected_3cv, Intensity = Intensity_3cv)
    
    # Combine both 2CV and 3CV data
    old_style_combined <- bind_rows(old_style_data_2cv, old_style_data_3cv)
    
    # Create all possible combinations (like the old script does)
    all_sample_cv_combinations <- expand.grid(
      Peptide = detected_details$Peptide,
      SampleID = samples_with_junction_peptides,
      CV_Type = c("2CV", "3CV"),
      stringsAsFactors = FALSE
    ) %>%
      mutate(Sample_CV = paste0(SampleID, "_", CV_Type))
    
    # Merge with actual detection data (mimicking old script's left_join)
    old_style_matrix_data <- all_sample_cv_combinations %>%
      left_join(
        old_style_combined,
        by = c("Peptide", "Sample_CV")
      ) %>%
      mutate(
        detected = ifelse(is.na(detected), FALSE, detected),
        Intensity = ifelse(is.na(Intensity), 0, Intensity)
      )
    
    # Debug: Check if EIFDRYGEEV is now in the matrix data
    cat("\n=== CHECKING EIFDRYGEEV IN MATRIX DATA ===\n")
    eifdrygeev_matrix_data <- old_style_matrix_data[old_style_matrix_data$Peptide == "EIFDRYGEEV", ]
    cat("EIFDRYGEEV rows in matrix data:", nrow(eifdrygeev_matrix_data), "\n")
    if(nrow(eifdrygeev_matrix_data) > 0) {
      print(eifdrygeev_matrix_data)
    }
    
    # Create presence matrix (matching old script exactly)
    old_style_presence <- old_style_matrix_data %>%
      select(Peptide, Sample_CV, present = detected) %>%
      pivot_wider(
        names_from = Sample_CV,
        values_from = present,
        values_fill = FALSE
      ) %>%
      column_to_rownames("Peptide")
    
    # Create intensity matrix (matching old script exactly)
    old_style_intensity <- old_style_matrix_data %>%
      select(Peptide, Sample_CV, Intensity) %>%
      pivot_wider(
        names_from = Sample_CV,
        values_from = Intensity,
        values_fill = 0
      ) %>%
      column_to_rownames("Peptide")
    
    # Debug: Check if EIFDRYGEEV is in the final matrices
    cat("\n=== CHECKING EIFDRYGEEV IN FINAL MATRICES ===\n")
    cat("EIFDRYGEEV in presence matrix:", "EIFDRYGEEV" %in% rownames(old_style_presence), "\n")
    cat("EIFDRYGEEV in intensity matrix:", "EIFDRYGEEV" %in% rownames(old_style_intensity), "\n")
    
    if("EIFDRYGEEV" %in% rownames(old_style_presence)) {
      cat("EIFDRYGEEV presence row:\n")
      print(old_style_presence["EIFDRYGEEV", ])
    }
    
    if(nrow(old_style_presence) > 0 && ncol(old_style_presence) > 0) {
      
      # Convert to numeric (exactly like old script)
      old_style_presence_mat <- as.matrix(old_style_presence)
      old_style_presence_numeric <- matrix(as.numeric(old_style_presence_mat), 
                                           nrow = nrow(old_style_presence_mat),
                                           dimnames = dimnames(old_style_presence_mat))
      
      # Log transform intensities (exactly like old script)
      log_old_style_intensity <- log10(as.matrix(old_style_intensity) + 1)
      
      # Create column annotations for CV type (simplified - no colors on top)
      old_style_column_ann <- data.frame(
        CV_Type = str_extract(colnames(old_style_presence), "2CV|3CV")
      )
      rownames(old_style_column_ann) <- colnames(old_style_presence)
      
      # Define minimal colors for annotation (just CV type)
      old_style_ann_colors <- list(
        CV_Type = c("2CV" = "steelblue", "3CV" = "tomato")
      )
      
      # Create FIXED OLD SCRIPT STYLE heatmaps
      pdf(file.path(results_dir, "FIXED_OLD_STYLE_junction_peptides_presence_comparison.pdf"), width = 16, height = 8)
      pheatmap(
        old_style_presence_numeric,
        main = "FIXED OLD SCRIPT STYLE: Presence of Junction Peptides (2CV vs 3CV)",
        color = c("white", "steelblue"),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        fontsize_row = 10,
        fontsize_col = 8,
        display_numbers = TRUE,
        number_format = "%.0f",
        annotation_col = old_style_column_ann,
        annotation_colors = old_style_ann_colors
      )
      dev.off()
      
      png(file.path(results_dir, "FIXED_OLD_STYLE_junction_peptides_presence_comparison.png"), width = 1400, height = 600, res = 100)
      pheatmap(
        old_style_presence_numeric,
        main = "FIXED OLD SCRIPT STYLE: Presence of Junction Peptides (2CV vs 3CV)",
        color = c("white", "steelblue"),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        fontsize_row = 10,
        fontsize_col = 8,
        display_numbers = TRUE,
        number_format = "%.0f",
        annotation_col = old_style_column_ann,
        annotation_colors = old_style_ann_colors
      )
      dev.off()
      
      # Also create the version WITH clustering as requested earlier
      pdf(file.path(results_dir, "FIXED_OLD_STYLE_junction_peptides_presence_CLUSTERED.pdf"), width = 16, height = 8)
      pheatmap(
        old_style_presence_numeric,
        main = "FIXED OLD SCRIPT STYLE: Presence of Junction Peptides (2CV vs 3CV) - WITH CLUSTERING", 
        color = c("white", "steelblue"),
        cluster_rows = TRUE,    # Enable row clustering
        cluster_cols = TRUE,    # Enable column clustering - THIS ADDS THE TREE
        clustering_distance_cols = "euclidean",
        clustering_method = "complete",
        fontsize_row = 10,
        fontsize_col = 8,
        display_numbers = TRUE,
        number_format = "%.0f",
        annotation_col = old_style_column_ann,
        annotation_colors = old_style_ann_colors,
        cutree_cols = 3,       # Cut dendrogram into groups
        show_colnames = TRUE,
        angle_col = 45
      )
      dev.off()
      
      png(file.path(results_dir, "FIXED_OLD_STYLE_junction_peptides_presence_CLUSTERED.png"), width = 1400, height = 600, res = 100)
      pheatmap(
        old_style_presence_numeric,
        main = "FIXED OLD SCRIPT STYLE: Presence of Junction Peptides (2CV vs 3CV) - WITH CLUSTERING",
        color = c("white", "steelblue"),
        cluster_rows = TRUE,
        cluster_cols = TRUE,    # This creates the relationship tree at the top
        clustering_distance_cols = "euclidean",
        clustering_method = "complete",
        fontsize_row = 10,
        fontsize_col = 8,
        display_numbers = TRUE,
        number_format = "%.0f",
        annotation_col = old_style_column_ann,
        annotation_colors = old_style_ann_colors,
        cutree_cols = 3,
        show_colnames = TRUE,
        angle_col = 45
      )
      dev.off()
      
      # Print final comparison information  
      cat("\n=== FIXED COMPARISON INFORMATION ===\n")
      cat("FIXED OLD SCRIPT STYLE Results:\n")
      cat("  Peptides in heatmap:", nrow(old_style_presence_numeric), "\n")
      cat("  Samples in heatmap:", ncol(old_style_presence_numeric), "\n")  
      cat("  Sample-CV combinations:", paste(colnames(old_style_presence_numeric), collapse = ", "), "\n")
      cat("  Peptides detected:", paste(rownames(old_style_presence_numeric), collapse = ", "), "\n")
      cat("  EIFDRYGEEV included:", "EIFDRYGEEV" %in% rownames(old_style_presence_numeric), "\n\n")
    }
  }
  
  # #===============================#
  # # Create OLD SCRIPT STYLE combined heatmap for direct comparison
  # #===============================#
  # 
  # # Mimic the old script's approach: create matrices and fill them manually
  # # This matches exactly what the old script does for paired samples
  # 
  # # First, simulate the old script's paired sample logic using our data
  # # Get all unique samples that have junction peptides
  # samples_with_junction_peptides <- junction_data %>%
  #   filter(!is.na(SampleID)) %>%
  #   pull(SampleID) %>%
  #   unique()
  # 
  # cat("Samples with junction peptides:", length(samples_with_junction_peptides), "\n")
  # cat("Sample IDs:", paste(samples_with_junction_peptides, collapse = ", "), "\n")
  # 
  # if(length(samples_with_junction_peptides) > 0 && nrow(detected_details) > 0) {
  #   
  #   # Create the exact same structure as the old script's paired heatmap
  #   # This simulates: paired_data <- combined_data %>% filter(Sample_ID %in% paired_samples)
  #   old_style_data_2cv <- junction_data %>%
  #     filter(!is.na(detected_2cv) & detected_2cv == TRUE) %>%
  #     mutate(Sample_CV = paste0(SampleID, "_2CV")) %>%
  #     select(Peptide, Sample_CV, detected = detected_2cv, Intensity = Intensity_2cv)
  #   
  #   old_style_data_3cv <- junction_data %>%
  #     filter(!is.na(detected_3cv) & detected_3cv == TRUE) %>%
  #     mutate(Sample_CV = paste0(SampleID, "_3CV")) %>%
  #     select(Peptide, Sample_CV, detected = detected_3cv, Intensity = Intensity_3cv)
  #   
  #   # Combine both 2CV and 3CV data
  #   old_style_combined <- bind_rows(old_style_data_2cv, old_style_data_3cv)
  #   
  #   # Create all possible combinations (like the old script does)
  #   all_sample_cv_combinations <- expand.grid(
  #     Peptide = detected_details$Peptide,
  #     SampleID = samples_with_junction_peptides,
  #     CV_Type = c("2CV", "3CV"),
  #     stringsAsFactors = FALSE
  #   ) %>%
  #     mutate(Sample_CV = paste0(SampleID, "_", CV_Type))
  #   
  #   # Merge with actual detection data (mimicking old script's left_join)
  #   old_style_matrix_data <- all_sample_cv_combinations %>%
  #     left_join(
  #       old_style_combined,
  #       by = c("Peptide", "Sample_CV")
  #     ) %>%
  #     mutate(
  #       detected = ifelse(is.na(detected), FALSE, detected),
  #       Intensity = ifelse(is.na(Intensity), 0, Intensity)
  #     )
  #   
  #   # Create presence matrix (matching old script exactly)
  #   old_style_presence <- old_style_matrix_data %>%
  #     select(Peptide, Sample_CV, present = detected) %>%
  #     pivot_wider(
  #       names_from = Sample_CV,
  #       values_from = present,
  #       values_fill = FALSE
  #     ) %>%
  #     column_to_rownames("Peptide")
  #   
  #   # Create intensity matrix (matching old script exactly)
  #   old_style_intensity <- old_style_matrix_data %>%
  #     select(Peptide, Sample_CV, Intensity) %>%
  #     pivot_wider(
  #       names_from = Sample_CV,
  #       values_from = Intensity,
  #       values_fill = 0
  #     ) %>%
  #     column_to_rownames("Peptide")
  #   
  #   if(nrow(old_style_presence) > 0 && ncol(old_style_presence) > 0) {
  #     
  #     # Convert to numeric (exactly like old script)
  #     old_style_presence_mat <- as.matrix(old_style_presence)
  #     old_style_presence_numeric <- matrix(as.numeric(old_style_presence_mat), 
  #                                          nrow = nrow(old_style_presence_mat),
  #                                          dimnames = dimnames(old_style_presence_mat))
  #     
  #     # Log transform intensities (exactly like old script)
  #     log_old_style_intensity <- log10(as.matrix(old_style_intensity) + 1)
  #     
  #     # Create column annotations for CV type (simplified - no colors on top)
  #     old_style_column_ann <- data.frame(
  #       CV_Type = str_extract(colnames(old_style_presence), "2CV|3CV")
  #     )
  #     rownames(old_style_column_ann) <- colnames(old_style_presence)
  #     
  #     # Define minimal colors for annotation (just CV type)
  #     old_style_ann_colors <- list(
  #       CV_Type = c("2CV" = "steelblue", "3CV" = "tomato")
  #     )
  #     
  #     # Create OLD SCRIPT STYLE heatmaps for direct comparison
  #     pdf(file.path(results_dir, "OLD_STYLE_junction_peptides_presence_comparison.pdf"), width = 16, height = 8)
  #     pheatmap(
  #       old_style_presence_numeric,
  #       main = "OLD SCRIPT STYLE: Presence of Junction Peptides (2CV vs 3CV)",
  #       color = c("white", "steelblue"),
  #       cluster_rows = FALSE,
  #       cluster_cols = FALSE,  # No clustering like old script
  #       fontsize_row = 10,
  #       fontsize_col = 8,
  #       display_numbers = TRUE,
  #       number_format = "%.0f",
  #       annotation_col = old_style_column_ann,
  #       annotation_colors = old_style_ann_colors
  #     )
  #     dev.off()
  #     
  #     png(file.path(results_dir, "OLD_STYLE_junction_peptides_presence_comparison.png"), width = 1400, height = 600, res = 100)
  #     pheatmap(
  #       old_style_presence_numeric,
  #       main = "OLD SCRIPT STYLE: Presence of Junction Peptides (2CV vs 3CV)",
  #       color = c("white", "steelblue"),
  #       cluster_rows = FALSE,
  #       cluster_cols = FALSE,
  #       fontsize_row = 10,
  #       fontsize_col = 8,
  #       display_numbers = TRUE,
  #       number_format = "%.0f",
  #       annotation_col = old_style_column_ann,
  #       annotation_colors = old_style_ann_colors
  #     )
  #     dev.off()
  #     
  #     pdf(file.path(results_dir, "OLD_STYLE_junction_peptides_intensity_comparison.pdf"), width = 16, height = 8)
  #     pheatmap(
  #       log_old_style_intensity,
  #       main = "OLD SCRIPT STYLE: Intensity of Junction Peptides (2CV vs 3CV, log10)",
  #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
  #       cluster_rows = FALSE,
  #       cluster_cols = FALSE,
  #       fontsize_row = 10,
  #       fontsize_col = 8,
  #       display_numbers = TRUE,
  #       number_format = "%.1f",
  #       annotation_col = old_style_column_ann,
  #       annotation_colors = old_style_ann_colors
  #     )
  #     dev.off()
  #     
  #     png(file.path(results_dir, "OLD_STYLE_junction_peptides_intensity_comparison.png"), width = 1400, height = 600, res = 100)
  #     pheatmap(
  #       log_old_style_intensity,
  #       main = "OLD SCRIPT STYLE: Intensity of Junction Peptides (2CV vs 3CV, log10)",
  #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
  #       cluster_rows = FALSE,
  #       cluster_cols = FALSE,
  #       fontsize_row = 10,
  #       fontsize_col = 8,
  #       display_numbers = TRUE,
  #       number_format = "%.1f",
  #       annotation_col = old_style_column_ann,
  #       annotation_colors = old_style_ann_colors
  #     )
  #     dev.off()
  #     
  #     # Print comparison information
  #     cat("\n=== COMPARISON INFORMATION ===\n")
  #     cat("OLD SCRIPT STYLE Results:\n")
  #     cat("  Peptides in heatmap:", nrow(old_style_presence_numeric), "\n")
  #     cat("  Samples in heatmap:", ncol(old_style_presence_numeric), "\n")  
  #     cat("  Sample-CV combinations:", paste(colnames(old_style_presence_numeric), collapse = ", "), "\n")
  #     cat("  Peptides detected:", paste(rownames(old_style_presence_numeric), collapse = ", "), "\n\n")
  #   }
  # }
  
  # Create intensity comparison plot
  intensity_data <- detected_details %>%
    filter(!is.na(Avg_Intensity_2CV) & !is.na(Avg_Intensity_3CV) & 
             Avg_Intensity_2CV > 0 & Avg_Intensity_3CV > 0)
  
  if(nrow(intensity_data) > 0) {
    intensity_comparison_plot <- ggplot(intensity_data, 
                                        aes(x = Avg_Intensity_2CV + 1, y = Avg_Intensity_3CV + 1)) +
      geom_point(aes(size = Total_Samples, color = Length)) +
      geom_text_repel(aes(label = Peptide), size = 3) +
      scale_x_log10() +
      scale_y_log10() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      theme_minimal() +
      labs(
        title = "Average Intensity of Junction Peptides: 2CV vs 3CV",
        x = "Average Intensity in 2CV (log scale)",
        y = "Average Intensity in 3CV (log scale)",
        color = "Peptide\nLength",
        size = "Total\nSamples"
      )
    
    ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.pdf"), intensity_comparison_plot, width = 10, height = 8)
    ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.png"), intensity_comparison_plot, width = 10, height = 8)
  }
  
  #===============================#
  # Create Excel Output
  #===============================#
  
  # Sheet with theoretical junction-spanning peptides
  excel_theoretical <- theoretical_junction_peptides %>%
    arrange(Length, Start_Position)
  
  # Sheet with detected junction peptides details
  excel_detected <- detected_details %>%
    arrange(Length, Start_Position)
  
  # Sheet with sample-level detection of junction peptides
  excel_sample_detection <- combined_data %>%
    filter(Spans_Junction == TRUE) %>%
    select(SampleID, Peptide, Peptide_Length, detected_2cv, detected_3cv, detected_both, 
           final_intensity, Intensity_2cv, Intensity_3cv, Spans_Junction) %>%
    arrange(SampleID, Peptide)
  
  # Sheet with 2CV vs 3CV detection statistics
  excel_cv_comparison <- detected_details %>%
    select(
      Peptide,
      Length,
      DNAJB1_Residues,
      PRKACA_Residues,
      Visualization,
      Total_Samples,
      Detected_2CV_Count,
      Detected_3CV_Count,
      Detected_Both_Count,
      Avg_Final_Intensity,
      Avg_Intensity_2CV,
      Avg_Intensity_3CV
    ) %>%
    mutate(
      `2CV Detection %` = (Detected_2CV_Count / Total_Samples) * 100,
      `3CV Detection %` = (Detected_3CV_Count / Total_Samples) * 100,
      `Both Detection %` = (Detected_Both_Count / Total_Samples) * 100,
      `3CV/2CV Detection Ratio` = (Detected_3CV_Count + 0.001) / (Detected_2CV_Count + 0.001),
      `3CV/2CV Intensity Ratio` = (Avg_Intensity_3CV + 1) / (Avg_Intensity_2CV + 1)
    ) %>%
    arrange(Length, Peptide)
  
  # Create a list of sheets for the Excel file
  excel_sheets <- list(
    "Theoretical_Junction_Peptides" = excel_theoretical,
    "Detected_Junction_Peptides" = excel_detected,
    "Sample_Level_Detection" = excel_sample_detection,
    "2CV_vs_3CV_Comparison" = excel_cv_comparison
  )
  
  # Write Excel file with multiple sheets
  write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
  
} else {
  cat("No junction-spanning peptides were detected in the dataset.\n")
  
  # Create Excel with just theoretical peptides
  excel_sheets <- list(
    "Theoretical_Junction_Peptides" = theoretical_junction_peptides
  )
  
  write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
}

# Print summary information
cat("\nAnalysis complete! Results saved to:", results_dir, "\n")
cat("\nThe following files were generated:\n")
cat("1. Fusion_Junction_Peptides_Analysis.xlsx - Excel file with detailed analysis\n")

if (nrow(detected_junction_peptides) > 0) {
  cat("2. junction_peptides_presence_2CV.pdf/png - Heatmap showing presence of junction peptides in 2CV samples\n")
  cat("3. junction_peptides_presence_3CV.pdf/png - Heatmap showing presence of junction peptides in 3CV samples\n")
  cat("4. junction_peptides_intensity_2CV.pdf/png - Heatmap showing intensity of junction peptides in 2CV samples\n")
  cat("5. junction_peptides_intensity_3CV.pdf/png - Heatmap showing intensity of junction peptides in 3CV samples\n")
  cat("6. OLD_STYLE_junction_peptides_presence_comparison.pdf/png - OLD SCRIPT STYLE heatmap for direct comparison\n")
  cat("7. OLD_STYLE_junction_peptides_intensity_comparison.pdf/png - OLD SCRIPT STYLE intensity heatmap for comparison\n")
}

cat("\nSummary of detected junction-spanning peptides:\n")
if (nrow(detected_junction_peptides) > 0) {
  for (i in 1:nrow(detected_junction_peptides)) {
    peptide_info <- detected_junction_peptides[i,]
    cat(sprintf("  %s (Length: %d) - 2CV: %d, 3CV: %d, Both: %d\n", 
                peptide_info$Peptide, peptide_info$Peptide_Length,
                peptide_info$Detected_2CV_Count, peptide_info$Detected_3CV_Count,
                peptide_info$Detected_Both_Count))
  }
} else {
  cat("  No junction-spanning peptides detected.\n")
}

# # Streamlined Script to analyze DNAJB1-PRKACA fusion protein peptides in 2CV vs 3CV samples
# # Looking for 8-12mers that span the fusion junction
# # Uses pre-processed combined dataset
# 
# # Setting directory
# setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/")
# 
# # Load required packages
# library(tidyverse)
# library(ggplot2)
# library(pheatmap)
# library(writexl)
# library(stringr)
# # Try to load ggrepel - install first if not available
# if (!requireNamespace("ggrepel", quietly = TRUE)) {
#   install.packages("ggrepel")
# }
# library(ggrepel)
# 
# # Define the fusion protein sequence
# fusion_protein <- "RKREIFDRYGEEVKEFLAKAKEDF"
# 
# # Identify the fusion junction - between E (DNAJB1) and V (PRKACA)
# dnajb1_part <- "RKREIFDRYGEE"  # Ends with YGEE (full DNAJB1 part)
# prkaca_part <- "VKEFLAKAKEDF"  # Starts with VKEF (full PRKACA part)
# junction_position <- nchar(dnajb1_part)  # Position after the last E in YGEE
# 
# cat("DNAJB1 part:", dnajb1_part, "\n")
# cat("PRKACA part:", prkaca_part, "\n")
# cat("Junction position:", junction_position, "\n")
# cat("Fusion protein:", fusion_protein, "\n")
# cat("Checking junction: character at junction position is", substr(fusion_protein, junction_position, junction_position), "\n")
# cat("Checking junction: next character is", substr(fusion_protein, junction_position+1, junction_position+1), "\n\n")
# 
# # Generate theoretical junction-spanning peptides (8-12mers)
# theoretical_peptides <- list()
# peptide_lengths <- 8:12  # Looking for 8-12mers
# 
# for (length in peptide_lengths) {
#   for (start_pos in 1:(nchar(fusion_protein) - length + 1)) {
#     peptide <- substr(fusion_protein, start_pos, start_pos + length - 1)
#     
#     # Check if this peptide spans the junction
#     # It spans if it includes at least one residue from both proteins
#     peptide_end_pos <- start_pos + length - 1
#     spans_junction <- start_pos <= junction_position && peptide_end_pos > junction_position
#     
#     if (spans_junction) {
#       theoretical_peptides[[length(theoretical_peptides) + 1]] <- list(
#         Peptide = peptide,
#         Length = length,
#         Start_Position = start_pos,
#         End_Position = peptide_end_pos,
#         DNAJB1_Part = paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
#         PRKACA_Part = paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = ""),
#         Visualization = paste0(
#           paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
#           paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = "")
#         )
#       )
#     }
#   }
# }
# 
# # Convert to dataframe
# theoretical_junction_peptides <- do.call(rbind, lapply(theoretical_peptides, function(x) {
#   data.frame(
#     Peptide = x$Peptide,
#     Length = x$Length,
#     Start_Position = x$Start_Position,
#     End_Position = x$End_Position,
#     DNAJB1_Part = x$DNAJB1_Part,
#     PRKACA_Part = x$PRKACA_Part,
#     Visualization = x$Visualization,
#     stringsAsFactors = FALSE
#   )
# }))
# 
# cat("Generated", nrow(theoretical_junction_peptides), "theoretical junction-spanning peptides\n")
# 
# # Count by length
# for (length in peptide_lengths) {
#   count <- sum(theoretical_junction_peptides$Length == length)
#   cat("  Length", length, ":", count, "peptides\n")
# }
# 
# # Print all theoretical peptides in a neat table
# cat("\nAll theoretical junction-spanning peptides:\n")
# for (i in 1:nrow(theoretical_junction_peptides)) {
#   peptide <- theoretical_junction_peptides[i,]
#   cat(sprintf("%2d. %s (Length: %d, Start: %d, End: %d, Vis: %s)\n", 
#               i, peptide$Peptide, peptide$Length, peptide$Start_Position, 
#               peptide$End_Position, peptide$Visualization))
# }
# 
# # Read the combined dataset
# cat("Reading combined dataset...\n")
# data_file <- "unique_peptides_all.tsv"
# combined_data <- read.delim(data_file, stringsAsFactors = FALSE)
# 
# cat("Dataset loaded with", nrow(combined_data), "peptide records\n")
# cat("Columns available:", paste(colnames(combined_data), collapse = ", "), "\n\n")
# 
# # DIAGNOSTIC: Check sample IDs in the dataset
# cat("=== DIAGNOSTIC INFORMATION ===\n")
# cat("Unique Sample IDs in dataset:\n")
# unique_samples <- unique(combined_data$SampleID)
# cat(paste(unique_samples, collapse = ", "), "\n\n")
# 
# # Check for any patterns that might cause 51S vs 51 issue
# samples_with_51 <- unique_samples[grepl("51", unique_samples)]
# cat("Sample IDs containing '51':", paste(samples_with_51, collapse = ", "), "\n")
# 
# # Check if there are issues with sample ID extraction
# cat("\nSample ID patterns:\n")
# sample_id_table <- table(combined_data$SampleID)
# print(sample_id_table)
# 
# # Check source files to understand the naming
# if("SourceFile" %in% colnames(combined_data)) {
#   cat("\nSource file patterns:\n")
#   source_files <- unique(combined_data$SourceFile)
#   files_with_51 <- source_files[grepl("51", source_files)]
#   cat("Files containing '51':\n")
#   for(file in files_with_51) {
#     cat("  ", file, "\n")
#   }
# }
# 
# cat("\n=== END DIAGNOSTIC ===\n\n")
# 
# # Create output directory for results
# results_dir <- "Fusion_Junction_Analysis"
# dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
# 
# #===============================#
# # Function to check if a peptide spans the fusion junction
# #===============================#
# 
# is_junction_spanning <- function(peptide) {
#   # First check if it's entirely within the fusion protein
#   if (!grepl(peptide, fusion_protein, fixed = TRUE)) {
#     return(FALSE)
#   }
#   
#   # Find position in fusion protein
#   start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
#   if (start_pos == -1) {
#     return(FALSE)
#   }
#   
#   end_pos <- start_pos + nchar(peptide) - 1
#   
#   # Check if it spans the junction
#   spans_junction <- start_pos <= junction_position && end_pos > junction_position
#   
#   return(spans_junction)
# }
# 
# # Filter data for 8-12mers and identify junction-spanning peptides
# combined_data <- combined_data %>%
#   mutate(
#     Peptide_Length = nchar(Peptide),
#     Spans_Junction = sapply(Peptide, is_junction_spanning)
#   ) %>%
#   filter(Peptide_Length >= 8 & Peptide_Length <= 12)  # Only include 8-12mers
# 
# # Extract detected fusion-spanning peptides
# # FIXED: First group by Peptide AND SampleID to avoid double-counting
# sample_level_detection <- combined_data %>%
#   filter(Spans_Junction == TRUE) %>%
#   group_by(Peptide, SampleID) %>%
#   summarize(
#     detected_2cv_this_sample = any(detected_2cv == TRUE, na.rm = TRUE),
#     detected_3cv_this_sample = any(detected_3cv == TRUE, na.rm = TRUE),
#     detected_both_this_sample = any(detected_both == TRUE, na.rm = TRUE),
#     max_final_intensity = max(final_intensity, na.rm = TRUE),
#     max_intensity_2cv = max(Intensity_2cv, na.rm = TRUE),
#     max_intensity_3cv = max(Intensity_3cv, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # Now aggregate by peptide to get correct counts
# detected_junction_peptides <- sample_level_detection %>%
#   group_by(Peptide) %>%
#   summarize(
#     Peptide_Length = nchar(Peptide[1]),
#     Total_Samples = n(),
#     Detected_2CV_Count = sum(detected_2cv_this_sample),
#     Detected_3CV_Count = sum(detected_3cv_this_sample),
#     Detected_Both_Count = sum(detected_both_this_sample),
#     Avg_Final_Intensity = mean(max_final_intensity[is.finite(max_final_intensity)], na.rm = TRUE),
#     Avg_Intensity_2CV = mean(max_intensity_2cv[is.finite(max_intensity_2cv)], na.rm = TRUE),
#     Avg_Intensity_3CV = mean(max_intensity_3cv[is.finite(max_intensity_3cv)], na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(Peptide_Length, Peptide)
# 
# # Get summary statistics about the dataset
# total_samples <- length(unique(combined_data$SampleID))
# samples_with_2cv <- sum(!is.na(combined_data$detected_2cv) & combined_data$detected_2cv == TRUE)
# samples_with_3cv <- sum(!is.na(combined_data$detected_3cv) & combined_data$detected_3cv == TRUE)
# 
# cat("Dataset summary:\n")
# cat("  Total unique samples:", total_samples, "\n")
# cat("  Detections in 2CV:", samples_with_2cv, "\n")
# cat("  Detections in 3CV:", samples_with_3cv, "\n\n")
# 
# if (nrow(detected_junction_peptides) > 0) {
#   # Create detailed analysis of each detected junction peptide
#   detected_details <- data.frame()
#   
#   for (i in 1:nrow(detected_junction_peptides)) {
#     peptide <- detected_junction_peptides$Peptide[i]
#     start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
#     end_pos <- start_pos + nchar(peptide) - 1
#     
#     # Calculate how many residues come from each protein
#     dnajb1_residues <- max(0, min(junction_position - start_pos + 1, nchar(peptide)))
#     prkaca_residues <- max(0, min(end_pos - junction_position, nchar(peptide)))
#     
#     # Visualization string (D for DNAJB1, P for PRKACA)
#     vis_string <- paste0(
#       paste(rep("D", dnajb1_residues), collapse = ""),
#       paste(rep("P", prkaca_residues), collapse = "")
#     )
#     
#     # Add to details dataframe
#     detected_details <- rbind(detected_details, data.frame(
#       Peptide = peptide,
#       Length = nchar(peptide),
#       Start_Position = start_pos,
#       End_Position = end_pos,
#       DNAJB1_Residues = dnajb1_residues,
#       PRKACA_Residues = prkaca_residues,
#       Visualization = vis_string,
#       Total_Samples = detected_junction_peptides$Total_Samples[i],
#       Detected_2CV_Count = detected_junction_peptides$Detected_2CV_Count[i],
#       Detected_3CV_Count = detected_junction_peptides$Detected_3CV_Count[i],
#       Detected_Both_Count = detected_junction_peptides$Detected_Both_Count[i],
#       Avg_Final_Intensity = detected_junction_peptides$Avg_Final_Intensity[i],
#       Avg_Intensity_2CV = detected_junction_peptides$Avg_Intensity_2CV[i],
#       Avg_Intensity_3CV = detected_junction_peptides$Avg_Intensity_3CV[i],
#       stringsAsFactors = FALSE
#     ))
#   }
#   
#   # Create heatmaps showing detection patterns
#   
#   # Prepare data for heatmaps - focusing on detected junction peptides
#   junction_data <- combined_data %>%
#     filter(Spans_Junction == TRUE) %>%
#     select(Peptide, SampleID, detected_2cv, detected_3cv, detected_both, 
#            final_intensity, Intensity_2cv, Intensity_3cv)
#   
#   # Create presence matrices for 2CV and 3CV
#   presence_2cv_data <- junction_data %>%
#     filter(!is.na(detected_2cv)) %>%
#     select(Peptide, SampleID, detected_2cv) %>%
#     pivot_wider(names_from = SampleID, values_from = detected_2cv, values_fill = FALSE) %>%
#     column_to_rownames("Peptide")
#   
#   presence_3cv_data <- junction_data %>%
#     filter(!is.na(detected_3cv)) %>%
#     select(Peptide, SampleID, detected_3cv) %>%
#     pivot_wider(names_from = SampleID, values_from = detected_3cv, values_fill = FALSE) %>%
#     column_to_rownames("Peptide")
#   
#   # Convert to numeric matrices
#   if(nrow(presence_2cv_data) > 0 && ncol(presence_2cv_data) > 0) {
#     presence_2cv_matrix <- as.matrix(presence_2cv_data)
#     presence_2cv_numeric <- matrix(as.numeric(presence_2cv_matrix), 
#                                    nrow = nrow(presence_2cv_matrix),
#                                    dimnames = dimnames(presence_2cv_matrix))
#     
#     # Create 2CV presence heatmap
#     pdf(file.path(results_dir, "junction_peptides_presence_2CV.pdf"), width = 12, height = 8)
#     pheatmap(
#       presence_2cv_numeric,
#       main = "Presence of Junction Peptides in 2CV Samples",
#       color = c("white", "steelblue"),
#       cluster_rows = FALSE,
#       cluster_cols = TRUE,
#       fontsize_row = 10,
#       fontsize_col = 8,
#       display_numbers = TRUE,
#       number_format = "%.0f"
#     )
#     dev.off()
#     
#     png(file.path(results_dir, "junction_peptides_presence_2CV.png"), width = 1000, height = 600, res = 100)
#     pheatmap(
#       presence_2cv_numeric,
#       main = "Presence of Junction Peptides in 2CV Samples",
#       color = c("white", "steelblue"),
#       cluster_rows = FALSE,
#       cluster_cols = TRUE,
#       fontsize_row = 10,
#       fontsize_col = 8,
#       display_numbers = TRUE,
#       number_format = "%.0f"
#     )
#     dev.off()
#   }
#   
#   if(nrow(presence_3cv_data) > 0 && ncol(presence_3cv_data) > 0) {
#     presence_3cv_matrix <- as.matrix(presence_3cv_data)
#     presence_3cv_numeric <- matrix(as.numeric(presence_3cv_matrix), 
#                                    nrow = nrow(presence_3cv_matrix),
#                                    dimnames = dimnames(presence_3cv_matrix))
#     
#     # Create 3CV presence heatmap
#     pdf(file.path(results_dir, "junction_peptides_presence_3CV.pdf"), width = 12, height = 8)
#     pheatmap(
#       presence_3cv_numeric,
#       main = "Presence of Junction Peptides in 3CV Samples",
#       color = c("white", "steelblue"),
#       cluster_rows = FALSE,
#       cluster_cols = TRUE,
#       fontsize_row = 10,
#       fontsize_col = 8,
#       display_numbers = TRUE,
#       number_format = "%.0f"
#     )
#     dev.off()
#     
#     png(file.path(results_dir, "junction_peptides_presence_3CV.png"), width = 1000, height = 600, res = 100)
#     pheatmap(
#       presence_3cv_numeric,
#       main = "Presence of Junction Peptides in 3CV Samples",
#       color = c("white", "steelblue"),
#       cluster_rows = FALSE,
#       cluster_cols = TRUE,
#       fontsize_row = 10,
#       fontsize_col = 8,
#       display_numbers = TRUE,
#       number_format = "%.0f"
#     )
#     dev.off()
#   }
#   
#   # Create intensity heatmaps
#   intensity_2cv_data <- junction_data %>%
#     filter(!is.na(Intensity_2cv) & Intensity_2cv > 0) %>%
#     select(Peptide, SampleID, Intensity_2cv) %>%
#     pivot_wider(names_from = SampleID, values_from = Intensity_2cv, values_fill = 0) %>%
#     column_to_rownames("Peptide")
#   
#   intensity_3cv_data <- junction_data %>%
#     filter(!is.na(Intensity_3cv) & Intensity_3cv > 0) %>%
#     select(Peptide, SampleID, Intensity_3cv) %>%
#     pivot_wider(names_from = SampleID, values_from = Intensity_3cv, values_fill = 0) %>%
#     column_to_rownames("Peptide")
#   
#   if(nrow(intensity_2cv_data) > 0 && ncol(intensity_2cv_data) > 0) {
#     log_intensity_2cv <- log10(as.matrix(intensity_2cv_data) + 1)
#     
#     pdf(file.path(results_dir, "junction_peptides_intensity_2CV.pdf"), width = 12, height = 8)
#     pheatmap(
#       log_intensity_2cv,
#       main = "Intensity of Junction Peptides in 2CV Samples (log10)",
#       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
#       cluster_rows = FALSE,
#       cluster_cols = TRUE,
#       fontsize_row = 10,
#       fontsize_col = 8,
#       display_numbers = TRUE,
#       number_format = "%.1f"
#     )
#     dev.off()
#     
#     png(file.path(results_dir, "junction_peptides_intensity_2CV.png"), width = 1000, height = 600, res = 100)
#     pheatmap(
#       log_intensity_2cv,
#       main = "Intensity of Junction Peptides in 2CV Samples (log10)",
#       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
#       cluster_rows = FALSE,
#       cluster_cols = TRUE,
#       fontsize_row = 10,
#       fontsize_col = 8,
#       display_numbers = TRUE,
#       number_format = "%.1f"
#     )
#     dev.off()
#   }
#   
#   if(nrow(intensity_3cv_data) > 0 && ncol(intensity_3cv_data) > 0) {
#     log_intensity_3cv <- log10(as.matrix(intensity_3cv_data) + 1)
#     
#     pdf(file.path(results_dir, "junction_peptides_intensity_3CV.pdf"), width = 12, height = 8)
#     pheatmap(
#       log_intensity_3cv,
#       main = "Intensity of Junction Peptides in 3CV Samples (log10)",
#       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
#       cluster_rows = FALSE,
#       cluster_cols = TRUE,
#       fontsize_row = 10,
#       fontsize_col = 8,
#       display_numbers = TRUE,
#       number_format = "%.1f"
#     )
#     dev.off()
#     
#     png(file.path(results_dir, "junction_peptides_intensity_3CV.png"), width = 1000, height = 600, res = 100)
#     pheatmap(
#       log_intensity_3cv,
#       main = "Intensity of Junction Peptides in 3CV Samples (log10)",
#       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
#       cluster_rows = FALSE,
#       cluster_cols = TRUE,
#       fontsize_row = 10,
#       fontsize_col = 8,
#       display_numbers = TRUE,
#       number_format = "%.1f"
#     )
#     dev.off()
#   }
#   
#   # Compare 2CV vs 3CV detection efficiency for junction peptides
#   cv_comparison_plot <- ggplot(detected_details, aes(x = Peptide, y = Length)) +
#     geom_point(aes(size = Total_Samples, color = Detected_3CV_Count / (Detected_2CV_Count + 0.001))) +
#     scale_color_gradient2(
#       low = "blue", 
#       mid = "white", 
#       high = "red", 
#       midpoint = 1,
#       name = "3CV/2CV\nDetection Ratio"
#     ) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(
#       title = "Fusion Junction Peptide Detection: 2CV vs 3CV",
#       x = "Peptide Sequence",
#       y = "Peptide Length",
#       size = "Total\nSamples"
#     )
#   
#   ggsave(file.path(results_dir, "junction_peptide_cv_comparison.pdf"), cv_comparison_plot, width = 12, height = 8)
#   ggsave(file.path(results_dir, "junction_peptide_cv_comparison.png"), cv_comparison_plot, width = 12, height = 8)
#   
#   # Create a visualization showing peptide coverage across the fusion protein
#   peptide_coverage_plot <- ggplot(detected_details, 
#                                   aes(x = Start_Position, xend = End_Position, 
#                                       y = reorder(Peptide, Start_Position), yend = reorder(Peptide, Start_Position))) +
#     geom_segment(aes(color = Total_Samples), linewidth = 5) +
#     geom_vline(xintercept = junction_position, linetype = "dashed", color = "red") +
#     annotate("text", x = junction_position - 2, y = nrow(detected_details) + 1, 
#              label = "DNAJB1", color = "darkgreen", fontface = "bold") +
#     annotate("text", x = junction_position + 2, y = nrow(detected_details) + 1, 
#              label = "PRKACA", color = "purple", fontface = "bold") +
#     scale_color_gradient(low = "lightsteelblue", high = "steelblue") +
#     theme_minimal() +
#     labs(
#       title = "Coverage of Fusion Junction by Detected Peptides",
#       x = "Position in Fusion Protein",
#       y = "Peptide",
#       color = "Total\nSamples"
#     )
#   
#   ggsave(file.path(results_dir, "junction_peptide_coverage.pdf"), peptide_coverage_plot, width = 12, height = 8)
#   ggsave(file.path(results_dir, "junction_peptide_coverage.png"), peptide_coverage_plot, width = 12, height = 8)
#   
#   #===============================#
#   # Create OLD SCRIPT STYLE combined heatmap for direct comparison
#   #===============================#
#   
#   # Mimic the old script's approach: create matrices and fill them manually
#   # This matches exactly what the old script does for paired samples
#   
#   # First, simulate the old script's paired sample logic using our data
#   # Get all unique samples that have junction peptides
#   samples_with_junction_peptides <- junction_data %>%
#     filter(!is.na(SampleID)) %>%
#     pull(SampleID) %>%
#     unique()
#   
#   cat("Samples with junction peptides:", length(samples_with_junction_peptides), "\n")
#   cat("Sample IDs:", paste(samples_with_junction_peptides, collapse = ", "), "\n")
#   
#   if(length(samples_with_junction_peptides) > 0 && nrow(detected_details) > 0) {
#     
#     # Create the exact same structure as the old script's paired heatmap
#     # This simulates: paired_data <- combined_data %>% filter(Sample_ID %in% paired_samples)
#     old_style_data_2cv <- junction_data %>%
#       filter(!is.na(detected_2cv) & detected_2cv == TRUE) %>%
#       mutate(Sample_CV = paste0(SampleID, "_2CV")) %>%
#       select(Peptide, Sample_CV, detected = detected_2cv, Intensity = Intensity_2cv)
#     
#     old_style_data_3cv <- junction_data %>%
#       filter(!is.na(detected_3cv) & detected_3cv == TRUE) %>%
#       mutate(Sample_CV = paste0(SampleID, "_3CV")) %>%
#       select(Peptide, Sample_CV, detected = detected_3cv, Intensity = Intensity_3cv)
#     
#     # Combine both 2CV and 3CV data
#     old_style_combined <- bind_rows(old_style_data_2cv, old_style_data_3cv)
#     
#     # Create all possible combinations (like the old script does)
#     all_sample_cv_combinations <- expand.grid(
#       Peptide = detected_details$Peptide,
#       SampleID = samples_with_junction_peptides,
#       CV_Type = c("2CV", "3CV"),
#       stringsAsFactors = FALSE
#     ) %>%
#       mutate(Sample_CV = paste0(SampleID, "_", CV_Type))
#     
#     # Merge with actual detection data (mimicking old script's left_join)
#     old_style_matrix_data <- all_sample_cv_combinations %>%
#       left_join(
#         old_style_combined,
#         by = c("Peptide", "Sample_CV")
#       ) %>%
#       mutate(
#         detected = ifelse(is.na(detected), FALSE, detected),
#         Intensity = ifelse(is.na(Intensity), 0, Intensity)
#       )
#     
#     # Create presence matrix (matching old script exactly)
#     old_style_presence <- old_style_matrix_data %>%
#       select(Peptide, Sample_CV, present = detected) %>%
#       pivot_wider(
#         names_from = Sample_CV,
#         values_from = present,
#         values_fill = FALSE
#       ) %>%
#       column_to_rownames("Peptide")
#     
#     # Create intensity matrix (matching old script exactly)
#     old_style_intensity <- old_style_matrix_data %>%
#       select(Peptide, Sample_CV, Intensity) %>%
#       pivot_wider(
#         names_from = Sample_CV,
#         values_from = Intensity,
#         values_fill = 0
#       ) %>%
#       column_to_rownames("Peptide")
#     
#     if(nrow(old_style_presence) > 0 && ncol(old_style_presence) > 0) {
#       
#       # Convert to numeric (exactly like old script)
#       old_style_presence_mat <- as.matrix(old_style_presence)
#       old_style_presence_numeric <- matrix(as.numeric(old_style_presence_mat), 
#                                            nrow = nrow(old_style_presence_mat),
#                                            dimnames = dimnames(old_style_presence_mat))
#       
#       # Log transform intensities (exactly like old script)
#       log_old_style_intensity <- log10(as.matrix(old_style_intensity) + 1)
#       
#       # Create column annotations (exactly like old script)
#       old_style_column_ann <- data.frame(
#         CV_Type = str_extract(colnames(old_style_presence), "2CV|3CV"),
#         Sample = gsub("_[23]CV$", "", colnames(old_style_presence))
#       )
#       rownames(old_style_column_ann) <- colnames(old_style_presence)
#       
#       # Define colors (exactly like old script)
#       old_style_ann_colors <- list(
#         CV_Type = c("2CV" = "steelblue", "3CV" = "tomato"),
#         Sample = setNames(
#           rainbow(length(samples_with_junction_peptides)),
#           samples_with_junction_peptides
#         )
#       )
#       
#       # Create OLD SCRIPT STYLE heatmaps for direct comparison
#       pdf(file.path(results_dir, "OLD_STYLE_junction_peptides_presence_comparison.pdf"), width = 16, height = 8)
#       pheatmap(
#         old_style_presence_numeric,
#         main = "OLD SCRIPT STYLE: Presence of Junction Peptides (2CV vs 3CV)",
#         color = c("white", "steelblue"),
#         cluster_rows = FALSE,
#         cluster_cols = FALSE,  # No clustering like old script
#         fontsize_row = 10,
#         fontsize_col = 8,
#         display_numbers = TRUE,
#         number_format = "%.0f",
#         annotation_col = old_style_column_ann,
#         annotation_colors = old_style_ann_colors
#       )
#       dev.off()
#       
#       png(file.path(results_dir, "OLD_STYLE_junction_peptides_presence_comparison.png"), width = 1400, height = 600, res = 100)
#       pheatmap(
#         old_style_presence_numeric,
#         main = "OLD SCRIPT STYLE: Presence of Junction Peptides (2CV vs 3CV)",
#         color = c("white", "steelblue"),
#         cluster_rows = FALSE,
#         cluster_cols = FALSE,
#         fontsize_row = 10,
#         fontsize_col = 8,
#         display_numbers = TRUE,
#         number_format = "%.0f",
#         annotation_col = old_style_column_ann,
#         annotation_colors = old_style_ann_colors
#       )
#       dev.off()
#       
#       pdf(file.path(results_dir, "OLD_STYLE_junction_peptides_intensity_comparison.pdf"), width = 16, height = 8)
#       pheatmap(
#         log_old_style_intensity,
#         main = "OLD SCRIPT STYLE: Intensity of Junction Peptides (2CV vs 3CV, log10)",
#         color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
#         cluster_rows = FALSE,
#         cluster_cols = FALSE,
#         fontsize_row = 10,
#         fontsize_col = 8,
#         display_numbers = TRUE,
#         number_format = "%.1f",
#         annotation_col = old_style_column_ann,
#         annotation_colors = old_style_ann_colors
#       )
#       dev.off()
#       
#       png(file.path(results_dir, "OLD_STYLE_junction_peptides_intensity_comparison.png"), width = 1400, height = 600, res = 100)
#       pheatmap(
#         log_old_style_intensity,
#         main = "OLD SCRIPT STYLE: Intensity of Junction Peptides (2CV vs 3CV, log10)",
#         color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
#         cluster_rows = FALSE,
#         cluster_cols = FALSE,
#         fontsize_row = 10,
#         fontsize_col = 8,
#         display_numbers = TRUE,
#         number_format = "%.1f",
#         annotation_col = old_style_column_ann,
#         annotation_colors = old_style_ann_colors
#       )
#       dev.off()
#       
#       # Print comparison information
#       cat("\n=== COMPARISON INFORMATION ===\n")
#       cat("OLD SCRIPT STYLE Results:\n")
#       cat("  Peptides in heatmap:", nrow(old_style_presence_numeric), "\n")
#       cat("  Samples in heatmap:", ncol(old_style_presence_numeric), "\n")  
#       cat("  Sample-CV combinations:", paste(colnames(old_style_presence_numeric), collapse = ", "), "\n")
#       cat("  Peptides detected:", paste(rownames(old_style_presence_numeric), collapse = ", "), "\n\n")
#     }
#   }
#   
#   # Create intensity comparison plot
#   intensity_data <- detected_details %>%
#     filter(!is.na(Avg_Intensity_2CV) & !is.na(Avg_Intensity_3CV) & 
#              Avg_Intensity_2CV > 0 & Avg_Intensity_3CV > 0)
#   
#   if(nrow(intensity_data) > 0) {
#     intensity_comparison_plot <- ggplot(intensity_data, 
#                                         aes(x = Avg_Intensity_2CV + 1, y = Avg_Intensity_3CV + 1)) +
#       geom_point(aes(size = Total_Samples, color = Length)) +
#       geom_text_repel(aes(label = Peptide), size = 3) +
#       scale_x_log10() +
#       scale_y_log10() +
#       geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
#       theme_minimal() +
#       labs(
#         title = "Average Intensity of Junction Peptides: 2CV vs 3CV",
#         x = "Average Intensity in 2CV (log scale)",
#         y = "Average Intensity in 3CV (log scale)",
#         color = "Peptide\nLength",
#         size = "Total\nSamples"
#       )
#     
#     ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.pdf"), intensity_comparison_plot, width = 10, height = 8)
#     ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.png"), intensity_comparison_plot, width = 10, height = 8)
#   }
#   
#   #===============================#
#   # Create Excel Output
#   #===============================#
#   
#   # Sheet with theoretical junction-spanning peptides
#   excel_theoretical <- theoretical_junction_peptides %>%
#     arrange(Length, Start_Position)
#   
#   # Sheet with detected junction peptides details
#   excel_detected <- detected_details %>%
#     arrange(Length, Start_Position)
#   
#   # Sheet with sample-level detection of junction peptides
#   excel_sample_detection <- combined_data %>%
#     filter(Spans_Junction == TRUE) %>%
#     select(SampleID, Peptide, Peptide_Length, detected_2cv, detected_3cv, detected_both, 
#            final_intensity, Intensity_2cv, Intensity_3cv, Spans_Junction) %>%
#     arrange(SampleID, Peptide)
#   
#   # Sheet with 2CV vs 3CV detection statistics
#   excel_cv_comparison <- detected_details %>%
#     select(
#       Peptide,
#       Length,
#       DNAJB1_Residues,
#       PRKACA_Residues,
#       Visualization,
#       Total_Samples,
#       Detected_2CV_Count,
#       Detected_3CV_Count,
#       Detected_Both_Count,
#       Avg_Final_Intensity,
#       Avg_Intensity_2CV,
#       Avg_Intensity_3CV
#     ) %>%
#     mutate(
#       `2CV Detection %` = (Detected_2CV_Count / Total_Samples) * 100,
#       `3CV Detection %` = (Detected_3CV_Count / Total_Samples) * 100,
#       `Both Detection %` = (Detected_Both_Count / Total_Samples) * 100,
#       `3CV/2CV Detection Ratio` = (Detected_3CV_Count + 0.001) / (Detected_2CV_Count + 0.001),
#       `3CV/2CV Intensity Ratio` = (Avg_Intensity_3CV + 1) / (Avg_Intensity_2CV + 1)
#     ) %>%
#     arrange(Length, Peptide)
#   
#   # Create a list of sheets for the Excel file
#   excel_sheets <- list(
#     "Theoretical_Junction_Peptides" = excel_theoretical,
#     "Detected_Junction_Peptides" = excel_detected,
#     "Sample_Level_Detection" = excel_sample_detection,
#     "2CV_vs_3CV_Comparison" = excel_cv_comparison
#   )
#   
#   # Write Excel file with multiple sheets
#   write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
#   
# } else {
#   cat("No junction-spanning peptides were detected in the dataset.\n")
#   
#   # Create Excel with just theoretical peptides
#   excel_sheets <- list(
#     "Theoretical_Junction_Peptides" = theoretical_junction_peptides
#   )
#   
#   write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
# }
# 
# # Print summary information
# cat("\nAnalysis complete! Results saved to:", results_dir, "\n")
# cat("\nThe following files were generated:\n")
# cat("1. Fusion_Junction_Peptides_Analysis.xlsx - Excel file with detailed analysis\n")
# 
# if (nrow(detected_junction_peptides) > 0) {
#   cat("2. junction_peptides_presence_2CV.pdf/png - Heatmap showing presence of junction peptides in 2CV samples\n")
#   cat("3. junction_peptides_presence_3CV.pdf/png - Heatmap showing presence of junction peptides in 3CV samples\n")
#   cat("4. junction_peptides_intensity_2CV.pdf/png - Heatmap showing intensity of junction peptides in 2CV samples\n")
#   cat("5. junction_peptides_intensity_3CV.pdf/png - Heatmap showing intensity of junction peptides in 3CV samples\n")
#   cat("6. OLD_STYLE_junction_peptides_presence_comparison.pdf/png - OLD SCRIPT STYLE heatmap for direct comparison\n")
#   cat("7. OLD_STYLE_junction_peptides_intensity_comparison.pdf/png - OLD SCRIPT STYLE intensity heatmap for comparison\n")
#   cat("8. junction_peptide_cv_comparison.pdf/png - Plot comparing detection in 2CV vs 3CV\n")
#   cat("9. junction_peptide_coverage.pdf/png - Plot showing coverage of fusion junction by detected peptides\n")
#   cat("10. junction_peptide_intensity_comparison.pdf/png - Plot comparing peptide intensities between 2CV and 3CV\n")
# }
# 
# cat("\nSummary of detected junction-spanning peptides:\n")
# if (nrow(detected_junction_peptides) > 0) {
#   for (i in 1:nrow(detected_junction_peptides)) {
#     peptide_info <- detected_junction_peptides[i,]
#     cat(sprintf("  %s (Length: %d) - 2CV: %d, 3CV: %d, Both: %d\n", 
#                 peptide_info$Peptide, peptide_info$Peptide_Length,
#                 peptide_info$Detected_2CV_Count, peptide_info$Detected_3CV_Count,
#                 peptide_info$Detected_Both_Count))
#   }
# } else {
#   cat("  No junction-spanning peptides detected.\n")
# }
# 
# 
# # # Streamlined Script to analyze DNAJB1-PRKACA fusion protein peptides in 2CV vs 3CV samples
# # # Looking for 8-12mers that span the fusion junction
# # # Uses pre-processed combined dataset
# # 
# # # Setting directory
# # setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/")
# # 
# # # Load required packages
# # library(tidyverse)
# # library(ggplot2)
# # library(pheatmap)
# # library(writexl)
# # library(stringr)
# # # Try to load ggrepel - install first if not available
# # if (!requireNamespace("ggrepel", quietly = TRUE)) {
# #   install.packages("ggrepel")
# # }
# # library(ggrepel)
# # 
# # # Define the fusion protein sequence
# # fusion_protein <- "RKREIFDRYGEEVKEFLAKAKEDF"
# # 
# # # Identify the fusion junction - between E (DNAJB1) and V (PRKACA)
# # dnajb1_part <- "RKREIFDRYGEE"  # Ends with YGEE (full DNAJB1 part)
# # prkaca_part <- "VKEFLAKAKEDF"  # Starts with VKEF (full PRKACA part)
# # junction_position <- nchar(dnajb1_part)  # Position after the last E in YGEE
# # 
# # cat("DNAJB1 part:", dnajb1_part, "\n")
# # cat("PRKACA part:", prkaca_part, "\n")
# # cat("Junction position:", junction_position, "\n")
# # cat("Fusion protein:", fusion_protein, "\n")
# # cat("Checking junction: character at junction position is", substr(fusion_protein, junction_position, junction_position), "\n")
# # cat("Checking junction: next character is", substr(fusion_protein, junction_position+1, junction_position+1), "\n\n")
# # 
# # # Generate theoretical junction-spanning peptides (8-12mers)
# # theoretical_peptides <- list()
# # peptide_lengths <- 8:12  # Looking for 8-12mers
# # 
# # for (length in peptide_lengths) {
# #   for (start_pos in 1:(nchar(fusion_protein) - length + 1)) {
# #     peptide <- substr(fusion_protein, start_pos, start_pos + length - 1)
# #     
# #     # Check if this peptide spans the junction
# #     # It spans if it includes at least one residue from both proteins
# #     peptide_end_pos <- start_pos + length - 1
# #     spans_junction <- start_pos <= junction_position && peptide_end_pos > junction_position
# #     
# #     if (spans_junction) {
# #       theoretical_peptides[[length(theoretical_peptides) + 1]] <- list(
# #         Peptide = peptide,
# #         Length = length,
# #         Start_Position = start_pos,
# #         End_Position = peptide_end_pos,
# #         DNAJB1_Part = paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
# #         PRKACA_Part = paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = ""),
# #         Visualization = paste0(
# #           paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
# #           paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = "")
# #         )
# #       )
# #     }
# #   }
# # }
# # 
# # # Convert to dataframe
# # theoretical_junction_peptides <- do.call(rbind, lapply(theoretical_peptides, function(x) {
# #   data.frame(
# #     Peptide = x$Peptide,
# #     Length = x$Length,
# #     Start_Position = x$Start_Position,
# #     End_Position = x$End_Position,
# #     DNAJB1_Part = x$DNAJB1_Part,
# #     PRKACA_Part = x$PRKACA_Part,
# #     Visualization = x$Visualization,
# #     stringsAsFactors = FALSE
# #   )
# # }))
# # 
# # cat("Generated", nrow(theoretical_junction_peptides), "theoretical junction-spanning peptides\n")
# # 
# # # Count by length
# # for (length in peptide_lengths) {
# #   count <- sum(theoretical_junction_peptides$Length == length)
# #   cat("  Length", length, ":", count, "peptides\n")
# # }
# # 
# # # Print all theoretical peptides in a neat table
# # cat("\nAll theoretical junction-spanning peptides:\n")
# # for (i in 1:nrow(theoretical_junction_peptides)) {
# #   peptide <- theoretical_junction_peptides[i,]
# #   cat(sprintf("%2d. %s (Length: %d, Start: %d, End: %d, Vis: %s)\n", 
# #               i, peptide$Peptide, peptide$Length, peptide$Start_Position, 
# #               peptide$End_Position, peptide$Visualization))
# # }
# # 
# # # Read the combined dataset
# # cat("Reading combined dataset...\n")
# # data_file <- "unique_peptides_unmodified.tsv"
# # combined_data <- read.delim(data_file, stringsAsFactors = FALSE)
# # 
# # cat("Dataset loaded with", nrow(combined_data), "peptide records\n")
# # cat("Columns available:", paste(colnames(combined_data), collapse = ", "), "\n\n")
# # 
# # # Create output directory for results
# # results_dir <- "Fusion_Junction_Analysis"
# # dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
# # 
# # #===============================#
# # # Function to check if a peptide spans the fusion junction
# # #===============================#
# # 
# # is_junction_spanning <- function(peptide) {
# #   # First check if it's entirely within the fusion protein
# #   if (!grepl(peptide, fusion_protein, fixed = TRUE)) {
# #     return(FALSE)
# #   }
# #   
# #   # Find position in fusion protein
# #   start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
# #   if (start_pos == -1) {
# #     return(FALSE)
# #   }
# #   
# #   end_pos <- start_pos + nchar(peptide) - 1
# #   
# #   # Check if it spans the junction
# #   spans_junction <- start_pos <= junction_position && end_pos > junction_position
# #   
# #   return(spans_junction)
# # }
# # 
# # # Filter data for 8-12mers and identify junction-spanning peptides
# # combined_data <- combined_data %>%
# #   mutate(
# #     Peptide_Length = nchar(Peptide),
# #     Spans_Junction = sapply(Peptide, is_junction_spanning)
# #   ) %>%
# #   filter(Peptide_Length >= 8 & Peptide_Length <= 12)  # Only include 8-12mers
# # 
# # # Extract detected fusion-spanning peptides
# # detected_junction_peptides <- combined_data %>%
# #   filter(Spans_Junction == TRUE) %>%
# #   group_by(Peptide, Peptide_Length) %>%
# #   summarize(
# #     Total_Samples = n(),
# #     Detected_2CV_Count = sum(!is.na(detected_2cv) & detected_2cv == TRUE, na.rm = TRUE),
# #     Detected_3CV_Count = sum(!is.na(detected_3cv) & detected_3cv == TRUE, na.rm = TRUE),
# #     Detected_Both_Count = sum(!is.na(detected_both) & detected_both == TRUE, na.rm = TRUE),
# #     Avg_Final_Intensity = mean(final_intensity, na.rm = TRUE),
# #     Avg_Intensity_2CV = mean(Intensity_2cv, na.rm = TRUE),
# #     Avg_Intensity_3CV = mean(Intensity_3cv, na.rm = TRUE),
# #     .groups = "drop"
# #   ) %>%
# #   arrange(Peptide_Length, Peptide)
# # 
# # cat("Found", nrow(detected_junction_peptides), "detected peptides that span the fusion junction\n")
# # 
# # # Get summary statistics about the dataset
# # total_samples <- length(unique(combined_data$SampleID))
# # samples_with_2cv <- sum(!is.na(combined_data$detected_2cv) & combined_data$detected_2cv == TRUE)
# # samples_with_3cv <- sum(!is.na(combined_data$detected_3cv) & combined_data$detected_3cv == TRUE)
# # 
# # cat("Dataset summary:\n")
# # cat("  Total unique samples:", total_samples, "\n")
# # cat("  Detections in 2CV:", samples_with_2cv, "\n")
# # cat("  Detections in 3CV:", samples_with_3cv, "\n\n")
# # 
# # if (nrow(detected_junction_peptides) > 0) {
# #   # Create detailed analysis of each detected junction peptide
# #   detected_details <- data.frame()
# #   
# #   for (i in 1:nrow(detected_junction_peptides)) {
# #     peptide <- detected_junction_peptides$Peptide[i]
# #     start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
# #     end_pos <- start_pos + nchar(peptide) - 1
# #     
# #     # Calculate how many residues come from each protein
# #     dnajb1_residues <- max(0, min(junction_position - start_pos + 1, nchar(peptide)))
# #     prkaca_residues <- max(0, min(end_pos - junction_position, nchar(peptide)))
# #     
# #     # Visualization string (D for DNAJB1, P for PRKACA)
# #     vis_string <- paste0(
# #       paste(rep("D", dnajb1_residues), collapse = ""),
# #       paste(rep("P", prkaca_residues), collapse = "")
# #     )
# #     
# #     # Add to details dataframe
# #     detected_details <- rbind(detected_details, data.frame(
# #       Peptide = peptide,
# #       Length = nchar(peptide),
# #       Start_Position = start_pos,
# #       End_Position = end_pos,
# #       DNAJB1_Residues = dnajb1_residues,
# #       PRKACA_Residues = prkaca_residues,
# #       Visualization = vis_string,
# #       Total_Samples = detected_junction_peptides$Total_Samples[i],
# #       Detected_2CV_Count = detected_junction_peptides$Detected_2CV_Count[i],
# #       Detected_3CV_Count = detected_junction_peptides$Detected_3CV_Count[i],
# #       Detected_Both_Count = detected_junction_peptides$Detected_Both_Count[i],
# #       Avg_Final_Intensity = detected_junction_peptides$Avg_Final_Intensity[i],
# #       Avg_Intensity_2CV = detected_junction_peptides$Avg_Intensity_2CV[i],
# #       Avg_Intensity_3CV = detected_junction_peptides$Avg_Intensity_3CV[i],
# #       stringsAsFactors = FALSE
# #     ))
# #   }
# #   
# #   # Create heatmaps showing detection patterns
# #   
# #   # Prepare data for heatmaps - focusing on detected junction peptides
# #   junction_data <- combined_data %>%
# #     filter(Spans_Junction == TRUE) %>%
# #     select(Peptide, SampleID, detected_2cv, detected_3cv, detected_both, 
# #            final_intensity, Intensity_2cv, Intensity_3cv)
# #   
# #   # Create presence matrices for 2CV and 3CV
# #   presence_2cv_data <- junction_data %>%
# #     filter(!is.na(detected_2cv)) %>%
# #     select(Peptide, SampleID, detected_2cv) %>%
# #     pivot_wider(names_from = SampleID, values_from = detected_2cv, values_fill = FALSE) %>%
# #     column_to_rownames("Peptide")
# #   
# #   presence_3cv_data <- junction_data %>%
# #     filter(!is.na(detected_3cv)) %>%
# #     select(Peptide, SampleID, detected_3cv) %>%
# #     pivot_wider(names_from = SampleID, values_from = detected_3cv, values_fill = FALSE) %>%
# #     column_to_rownames("Peptide")
# #   
# #   # Convert to numeric matrices
# #   if(nrow(presence_2cv_data) > 0 && ncol(presence_2cv_data) > 0) {
# #     presence_2cv_matrix <- as.matrix(presence_2cv_data)
# #     presence_2cv_numeric <- matrix(as.numeric(presence_2cv_matrix), 
# #                                    nrow = nrow(presence_2cv_matrix),
# #                                    dimnames = dimnames(presence_2cv_matrix))
# #     
# #     # Create 2CV presence heatmap
# #     pdf(file.path(results_dir, "junction_peptides_presence_2CV.pdf"), width = 12, height = 8)
# #     pheatmap(
# #       presence_2cv_numeric,
# #       main = "Presence of Junction Peptides in 2CV Samples",
# #       color = c("white", "steelblue"),
# #       cluster_rows = FALSE,
# #       cluster_cols = TRUE,
# #       fontsize_row = 10,
# #       fontsize_col = 8,
# #       display_numbers = TRUE,
# #       number_format = "%.0f"
# #     )
# #     dev.off()
# #     
# #     png(file.path(results_dir, "junction_peptides_presence_2CV.png"), width = 1000, height = 600, res = 100)
# #     pheatmap(
# #       presence_2cv_numeric,
# #       main = "Presence of Junction Peptides in 2CV Samples",
# #       color = c("white", "steelblue"),
# #       cluster_rows = FALSE,
# #       cluster_cols = TRUE,
# #       fontsize_row = 10,
# #       fontsize_col = 8,
# #       display_numbers = TRUE,
# #       number_format = "%.0f"
# #     )
# #     dev.off()
# #   }
# #   
# #   if(nrow(presence_3cv_data) > 0 && ncol(presence_3cv_data) > 0) {
# #     presence_3cv_matrix <- as.matrix(presence_3cv_data)
# #     presence_3cv_numeric <- matrix(as.numeric(presence_3cv_matrix), 
# #                                    nrow = nrow(presence_3cv_matrix),
# #                                    dimnames = dimnames(presence_3cv_matrix))
# #     
# #     # Create 3CV presence heatmap
# #     pdf(file.path(results_dir, "junction_peptides_presence_3CV.pdf"), width = 12, height = 8)
# #     pheatmap(
# #       presence_3cv_numeric,
# #       main = "Presence of Junction Peptides in 3CV Samples",
# #       color = c("white", "steelblue"),
# #       cluster_rows = FALSE,
# #       cluster_cols = TRUE,
# #       fontsize_row = 10,
# #       fontsize_col = 8,
# #       display_numbers = TRUE,
# #       number_format = "%.0f"
# #     )
# #     dev.off()
# #     
# #     png(file.path(results_dir, "junction_peptides_presence_3CV.png"), width = 1000, height = 600, res = 100)
# #     pheatmap(
# #       presence_3cv_numeric,
# #       main = "Presence of Junction Peptides in 3CV Samples",
# #       color = c("white", "steelblue"),
# #       cluster_rows = FALSE,
# #       cluster_cols = TRUE,
# #       fontsize_row = 10,
# #       fontsize_col = 8,
# #       display_numbers = TRUE,
# #       number_format = "%.0f"
# #     )
# #     dev.off()
# #   }
# #   
# #   # Create intensity heatmaps
# #   intensity_2cv_data <- junction_data %>%
# #     filter(!is.na(Intensity_2cv) & Intensity_2cv > 0) %>%
# #     select(Peptide, SampleID, Intensity_2cv) %>%
# #     pivot_wider(names_from = SampleID, values_from = Intensity_2cv, values_fill = 0) %>%
# #     column_to_rownames("Peptide")
# #   
# #   intensity_3cv_data <- junction_data %>%
# #     filter(!is.na(Intensity_3cv) & Intensity_3cv > 0) %>%
# #     select(Peptide, SampleID, Intensity_3cv) %>%
# #     pivot_wider(names_from = SampleID, values_from = Intensity_3cv, values_fill = 0) %>%
# #     column_to_rownames("Peptide")
# #   
# #   if(nrow(intensity_2cv_data) > 0 && ncol(intensity_2cv_data) > 0) {
# #     log_intensity_2cv <- log10(as.matrix(intensity_2cv_data) + 1)
# #     
# #     pdf(file.path(results_dir, "junction_peptides_intensity_2CV.pdf"), width = 12, height = 8)
# #     pheatmap(
# #       log_intensity_2cv,
# #       main = "Intensity of Junction Peptides in 2CV Samples (log10)",
# #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# #       cluster_rows = FALSE,
# #       cluster_cols = TRUE,
# #       fontsize_row = 10,
# #       fontsize_col = 8,
# #       display_numbers = TRUE,
# #       number_format = "%.1f"
# #     )
# #     dev.off()
# #     
# #     png(file.path(results_dir, "junction_peptides_intensity_2CV.png"), width = 1000, height = 600, res = 100)
# #     pheatmap(
# #       log_intensity_2cv,
# #       main = "Intensity of Junction Peptides in 2CV Samples (log10)",
# #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# #       cluster_rows = FALSE,
# #       cluster_cols = TRUE,
# #       fontsize_row = 10,
# #       fontsize_col = 8,
# #       display_numbers = TRUE,
# #       number_format = "%.1f"
# #     )
# #     dev.off()
# #   }
# #   
# #   if(nrow(intensity_3cv_data) > 0 && ncol(intensity_3cv_data) > 0) {
# #     log_intensity_3cv <- log10(as.matrix(intensity_3cv_data) + 1)
# #     
# #     pdf(file.path(results_dir, "junction_peptides_intensity_3CV.pdf"), width = 12, height = 8)
# #     pheatmap(
# #       log_intensity_3cv,
# #       main = "Intensity of Junction Peptides in 3CV Samples (log10)",
# #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# #       cluster_rows = FALSE,
# #       cluster_cols = TRUE,
# #       fontsize_row = 10,
# #       fontsize_col = 8,
# #       display_numbers = TRUE,
# #       number_format = "%.1f"
# #     )
# #     dev.off()
# #     
# #     png(file.path(results_dir, "junction_peptides_intensity_3CV.png"), width = 1000, height = 600, res = 100)
# #     pheatmap(
# #       log_intensity_3cv,
# #       main = "Intensity of Junction Peptides in 3CV Samples (log10)",
# #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# #       cluster_rows = FALSE,
# #       cluster_cols = TRUE,
# #       fontsize_row = 10,
# #       fontsize_col = 8,
# #       display_numbers = TRUE,
# #       number_format = "%.1f"
# #     )
# #     dev.off()
# #   }
# #   
# #   # Compare 2CV vs 3CV detection efficiency for junction peptides
# #   cv_comparison_plot <- ggplot(detected_details, aes(x = Peptide, y = Length)) +
# #     geom_point(aes(size = Total_Samples, color = Detected_3CV_Count / (Detected_2CV_Count + 0.001))) +
# #     scale_color_gradient2(
# #       low = "blue", 
# #       mid = "white", 
# #       high = "red", 
# #       midpoint = 1,
# #       name = "3CV/2CV\nDetection Ratio"
# #     ) +
# #     theme_minimal() +
# #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# #     labs(
# #       title = "Fusion Junction Peptide Detection: 2CV vs 3CV",
# #       x = "Peptide Sequence",
# #       y = "Peptide Length",
# #       size = "Total\nSamples"
# #     )
# #   
# #   ggsave(file.path(results_dir, "junction_peptide_cv_comparison.pdf"), cv_comparison_plot, width = 12, height = 8)
# #   ggsave(file.path(results_dir, "junction_peptide_cv_comparison.png"), cv_comparison_plot, width = 12, height = 8)
# #   
# #   # Create a visualization showing peptide coverage across the fusion protein
# #   peptide_coverage_plot <- ggplot(detected_details, 
# #                                   aes(x = Start_Position, xend = End_Position, 
# #                                       y = reorder(Peptide, Start_Position), yend = reorder(Peptide, Start_Position))) +
# #     geom_segment(aes(color = Total_Samples), linewidth = 5) +
# #     geom_vline(xintercept = junction_position, linetype = "dashed", color = "red") +
# #     annotate("text", x = junction_position - 2, y = nrow(detected_details) + 1, 
# #              label = "DNAJB1", color = "darkgreen", fontface = "bold") +
# #     annotate("text", x = junction_position + 2, y = nrow(detected_details) + 1, 
# #              label = "PRKACA", color = "purple", fontface = "bold") +
# #     scale_color_gradient(low = "lightsteelblue", high = "steelblue") +
# #     theme_minimal() +
# #     labs(
# #       title = "Coverage of Fusion Junction by Detected Peptides",
# #       x = "Position in Fusion Protein",
# #       y = "Peptide",
# #       color = "Total\nSamples"
# #     )
# #   
# #   ggsave(file.path(results_dir, "junction_peptide_coverage.pdf"), peptide_coverage_plot, width = 12, height = 8)
# #   ggsave(file.path(results_dir, "junction_peptide_coverage.png"), peptide_coverage_plot, width = 12, height = 8)
# #   
# #   #===============================#
# #   # Create OLD SCRIPT STYLE combined heatmap for direct comparison
# #   #===============================#
# #   
# #   # Mimic the old script's approach: create matrices and fill them manually
# #   # This matches exactly what the old script does for paired samples
# #   
# #   # First, simulate the old script's paired sample logic using our data
# #   # Get all unique samples that have junction peptides
# #   samples_with_junction_peptides <- junction_data %>%
# #     filter(!is.na(SampleID)) %>%
# #     pull(SampleID) %>%
# #     unique()
# #   
# #   cat("Samples with junction peptides:", length(samples_with_junction_peptides), "\n")
# #   cat("Sample IDs:", paste(samples_with_junction_peptides, collapse = ", "), "\n")
# #   
# #   if(length(samples_with_junction_peptides) > 0 && nrow(detected_details) > 0) {
# #     
# #     # Create the exact same structure as the old script's paired heatmap
# #     # This simulates: paired_data <- combined_data %>% filter(Sample_ID %in% paired_samples)
# #     old_style_data_2cv <- junction_data %>%
# #       filter(!is.na(detected_2cv) & detected_2cv == TRUE) %>%
# #       mutate(Sample_CV = paste0(SampleID, "_2CV")) %>%
# #       select(Peptide, Sample_CV, detected = detected_2cv, Intensity = Intensity_2cv)
# #     
# #     old_style_data_3cv <- junction_data %>%
# #       filter(!is.na(detected_3cv) & detected_3cv == TRUE) %>%
# #       mutate(Sample_CV = paste0(SampleID, "_3CV")) %>%
# #       select(Peptide, Sample_CV, detected = detected_3cv, Intensity = Intensity_3cv)
# #     
# #     # Combine both 2CV and 3CV data
# #     old_style_combined <- bind_rows(old_style_data_2cv, old_style_data_3cv)
# #     
# #     # Create all possible combinations (like the old script does)
# #     all_sample_cv_combinations <- expand.grid(
# #       Peptide = detected_details$Peptide,
# #       SampleID = samples_with_junction_peptides,
# #       CV_Type = c("2CV", "3CV"),
# #       stringsAsFactors = FALSE
# #     ) %>%
# #       mutate(Sample_CV = paste0(SampleID, "_", CV_Type))
# #     
# #     # Merge with actual detection data (mimicking old script's left_join)
# #     old_style_matrix_data <- all_sample_cv_combinations %>%
# #       left_join(
# #         old_style_combined,
# #         by = c("Peptide", "Sample_CV")
# #       ) %>%
# #       mutate(
# #         detected = ifelse(is.na(detected), FALSE, detected),
# #         Intensity = ifelse(is.na(Intensity), 0, Intensity)
# #       )
# #     
# #     # Create presence matrix (matching old script exactly)
# #     old_style_presence <- old_style_matrix_data %>%
# #       select(Peptide, Sample_CV, present = detected) %>%
# #       pivot_wider(
# #         names_from = Sample_CV,
# #         values_from = present,
# #         values_fill = FALSE
# #       ) %>%
# #       column_to_rownames("Peptide")
# #     
# #     # Create intensity matrix (matching old script exactly)
# #     old_style_intensity <- old_style_matrix_data %>%
# #       select(Peptide, Sample_CV, Intensity) %>%
# #       pivot_wider(
# #         names_from = Sample_CV,
# #         values_from = Intensity,
# #         values_fill = 0
# #       ) %>%
# #       column_to_rownames("Peptide")
# #     
# #     if(nrow(old_style_presence) > 0 && ncol(old_style_presence) > 0) {
# #       
# #       # Convert to numeric (exactly like old script)
# #       old_style_presence_mat <- as.matrix(old_style_presence)
# #       old_style_presence_numeric <- matrix(as.numeric(old_style_presence_mat), 
# #                                            nrow = nrow(old_style_presence_mat),
# #                                            dimnames = dimnames(old_style_presence_mat))
# #       
# #       # Log transform intensities (exactly like old script)
# #       log_old_style_intensity <- log10(as.matrix(old_style_intensity) + 1)
# #       
# #       # Create column annotations (exactly like old script)
# #       old_style_column_ann <- data.frame(
# #         CV_Type = str_extract(colnames(old_style_presence), "2CV|3CV"),
# #         Sample = gsub("_[23]CV$", "", colnames(old_style_presence))
# #       )
# #       rownames(old_style_column_ann) <- colnames(old_style_presence)
# #       
# #       # Define colors (exactly like old script)
# #       old_style_ann_colors <- list(
# #         CV_Type = c("2CV" = "steelblue", "3CV" = "tomato"),
# #         Sample = setNames(
# #           rainbow(length(samples_with_junction_peptides)),
# #           samples_with_junction_peptides
# #         )
# #       )
# #       
# #       # Create OLD SCRIPT STYLE heatmaps for direct comparison
# #       pdf(file.path(results_dir, "OLD_STYLE_junction_peptides_presence_comparison.pdf"), width = 16, height = 8)
# #       pheatmap(
# #         old_style_presence_numeric,
# #         main = "OLD SCRIPT STYLE: Presence of Junction Peptides (2CV vs 3CV)",
# #         color = c("white", "steelblue"),
# #         cluster_rows = FALSE,
# #         cluster_cols = FALSE,  # No clustering like old script
# #         fontsize_row = 10,
# #         fontsize_col = 8,
# #         display_numbers = TRUE,
# #         number_format = "%.0f",
# #         annotation_col = old_style_column_ann,
# #         annotation_colors = old_style_ann_colors
# #       )
# #       dev.off()
# #       
# #       png(file.path(results_dir, "OLD_STYLE_junction_peptides_presence_comparison.png"), width = 1400, height = 600, res = 100)
# #       pheatmap(
# #         old_style_presence_numeric,
# #         main = "OLD SCRIPT STYLE: Presence of Junction Peptides (2CV vs 3CV)",
# #         color = c("white", "steelblue"),
# #         cluster_rows = FALSE,
# #         cluster_cols = FALSE,
# #         fontsize_row = 10,
# #         fontsize_col = 8,
# #         display_numbers = TRUE,
# #         number_format = "%.0f",
# #         annotation_col = old_style_column_ann,
# #         annotation_colors = old_style_ann_colors
# #       )
# #       dev.off()
# #       
# #       pdf(file.path(results_dir, "OLD_STYLE_junction_peptides_intensity_comparison.pdf"), width = 16, height = 8)
# #       pheatmap(
# #         log_old_style_intensity,
# #         main = "OLD SCRIPT STYLE: Intensity of Junction Peptides (2CV vs 3CV, log10)",
# #         color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# #         cluster_rows = FALSE,
# #         cluster_cols = FALSE,
# #         fontsize_row = 10,
# #         fontsize_col = 8,
# #         display_numbers = TRUE,
# #         number_format = "%.1f",
# #         annotation_col = old_style_column_ann,
# #         annotation_colors = old_style_ann_colors
# #       )
# #       dev.off()
# #       
# #       png(file.path(results_dir, "OLD_STYLE_junction_peptides_intensity_comparison.png"), width = 1400, height = 600, res = 100)
# #       pheatmap(
# #         log_old_style_intensity,
# #         main = "OLD SCRIPT STYLE: Intensity of Junction Peptides (2CV vs 3CV, log10)",
# #         color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# #         cluster_rows = FALSE,
# #         cluster_cols = FALSE,
# #         fontsize_row = 10,
# #         fontsize_col = 8,
# #         display_numbers = TRUE,
# #         number_format = "%.1f",
# #         annotation_col = old_style_column_ann,
# #         annotation_colors = old_style_ann_colors
# #       )
# #       dev.off()
# #       
# #       # Print comparison information
# #       cat("\n=== COMPARISON INFORMATION ===\n")
# #       cat("OLD SCRIPT STYLE Results:\n")
# #       cat("  Peptides in heatmap:", nrow(old_style_presence_numeric), "\n")
# #       cat("  Samples in heatmap:", ncol(old_style_presence_numeric), "\n")  
# #       cat("  Sample-CV combinations:", paste(colnames(old_style_presence_numeric), collapse = ", "), "\n")
# #       cat("  Peptides detected:", paste(rownames(old_style_presence_numeric), collapse = ", "), "\n\n")
# #     }
# #   }
# #   
# #   # Create intensity comparison plot
# #   intensity_data <- detected_details %>%
# #     filter(!is.na(Avg_Intensity_2CV) & !is.na(Avg_Intensity_3CV) & 
# #              Avg_Intensity_2CV > 0 & Avg_Intensity_3CV > 0)
# #   
# #   if(nrow(intensity_data) > 0) {
# #     intensity_comparison_plot <- ggplot(intensity_data, 
# #                                         aes(x = Avg_Intensity_2CV + 1, y = Avg_Intensity_3CV + 1)) +
# #       geom_point(aes(size = Total_Samples, color = Length)) +
# #       geom_text_repel(aes(label = Peptide), size = 3) +
# #       scale_x_log10() +
# #       scale_y_log10() +
# #       geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
# #       theme_minimal() +
# #       labs(
# #         title = "Average Intensity of Junction Peptides: 2CV vs 3CV",
# #         x = "Average Intensity in 2CV (log scale)",
# #         y = "Average Intensity in 3CV (log scale)",
# #         color = "Peptide\nLength",
# #         size = "Total\nSamples"
# #       )
# #     
# #     ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.pdf"), intensity_comparison_plot, width = 10, height = 8)
# #     ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.png"), intensity_comparison_plot, width = 10, height = 8)
# #   }
# #   
# #   #===============================#
# #   # Create Excel Output
# #   #===============================#
# #   
# #   # Sheet with theoretical junction-spanning peptides
# #   excel_theoretical <- theoretical_junction_peptides %>%
# #     arrange(Length, Start_Position)
# #   
# #   # Sheet with detected junction peptides details
# #   excel_detected <- detected_details %>%
# #     arrange(Length, Start_Position)
# #   
# #   # Sheet with sample-level detection of junction peptides
# #   excel_sample_detection <- combined_data %>%
# #     filter(Spans_Junction == TRUE) %>%
# #     select(SampleID, Peptide, Peptide_Length, detected_2cv, detected_3cv, detected_both, 
# #            final_intensity, Intensity_2cv, Intensity_3cv, Spans_Junction) %>%
# #     arrange(SampleID, Peptide)
# #   
# #   # Sheet with 2CV vs 3CV detection statistics
# #   excel_cv_comparison <- detected_details %>%
# #     select(
# #       Peptide,
# #       Length,
# #       DNAJB1_Residues,
# #       PRKACA_Residues,
# #       Visualization,
# #       Total_Samples,
# #       Detected_2CV_Count,
# #       Detected_3CV_Count,
# #       Detected_Both_Count,
# #       Avg_Final_Intensity,
# #       Avg_Intensity_2CV,
# #       Avg_Intensity_3CV
# #     ) %>%
# #     mutate(
# #       `2CV Detection %` = (Detected_2CV_Count / Total_Samples) * 100,
# #       `3CV Detection %` = (Detected_3CV_Count / Total_Samples) * 100,
# #       `Both Detection %` = (Detected_Both_Count / Total_Samples) * 100,
# #       `3CV/2CV Detection Ratio` = (Detected_3CV_Count + 0.001) / (Detected_2CV_Count + 0.001),
# #       `3CV/2CV Intensity Ratio` = (Avg_Intensity_3CV + 1) / (Avg_Intensity_2CV + 1)
# #     ) %>%
# #     arrange(Length, Peptide)
# #   
# #   # Create a list of sheets for the Excel file
# #   excel_sheets <- list(
# #     "Theoretical_Junction_Peptides" = excel_theoretical,
# #     "Detected_Junction_Peptides" = excel_detected,
# #     "Sample_Level_Detection" = excel_sample_detection,
# #     "2CV_vs_3CV_Comparison" = excel_cv_comparison
# #   )
# #   
# #   # Write Excel file with multiple sheets
# #   write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
# #   
# # } else {
# #   cat("No junction-spanning peptides were detected in the dataset.\n")
# #   
# #   # Create Excel with just theoretical peptides
# #   excel_sheets <- list(
# #     "Theoretical_Junction_Peptides" = theoretical_junction_peptides
# #   )
# #   
# #   write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
# # }
# # 
# # # Print summary information
# # cat("\nAnalysis complete! Results saved to:", results_dir, "\n")
# # cat("\nThe following files were generated:\n")
# # cat("1. Fusion_Junction_Peptides_Analysis.xlsx - Excel file with detailed analysis\n")
# # 
# # if (nrow(detected_junction_peptides) > 0) {
# #   cat("2. junction_peptides_presence_2CV.pdf/png - Heatmap showing presence of junction peptides in 2CV samples\n")
# #   cat("3. junction_peptides_presence_3CV.pdf/png - Heatmap showing presence of junction peptides in 3CV samples\n")
# #   cat("4. junction_peptides_intensity_2CV.pdf/png - Heatmap showing intensity of junction peptides in 2CV samples\n")
# #   cat("5. junction_peptides_intensity_3CV.pdf/png - Heatmap showing intensity of junction peptides in 3CV samples\n")
# #   cat("6. OLD_STYLE_junction_peptides_presence_comparison.pdf/png - OLD SCRIPT STYLE heatmap for direct comparison\n")
# #   cat("7. OLD_STYLE_junction_peptides_intensity_comparison.pdf/png - OLD SCRIPT STYLE intensity heatmap for comparison\n")
# #   cat("8. junction_peptide_cv_comparison.pdf/png - Plot comparing detection in 2CV vs 3CV\n")
# #   cat("9. junction_peptide_coverage.pdf/png - Plot showing coverage of fusion junction by detected peptides\n")
# #   cat("10. junction_peptide_intensity_comparison.pdf/png - Plot comparing peptide intensities between 2CV and 3CV\n")
# # }
# # 
# # cat("\nSummary of detected junction-spanning peptides:\n")
# # if (nrow(detected_junction_peptides) > 0) {
# #   for (i in 1:nrow(detected_junction_peptides)) {
# #     peptide_info <- detected_junction_peptides[i,]
# #     cat(sprintf("  %s (Length: %d) - 2CV: %d, 3CV: %d, Both: %d\n", 
# #                 peptide_info$Peptide, peptide_info$Peptide_Length,
# #                 peptide_info$Detected_2CV_Count, peptide_info$Detected_3CV_Count,
# #                 peptide_info$Detected_Both_Count))
# #   }
# # } else {
# #   cat("  No junction-spanning peptides detected.\n")
# # }
# # 
# # # # Streamlined Script to analyze DNAJB1-PRKACA fusion protein peptides in 2CV vs 3CV samples
# # # # Looking for 8-12mers that span the fusion junction
# # # # Uses pre-processed combined dataset
# # # 
# # # # Setting directory
# # # setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/")
# # # 
# # # # Load required packages
# # # library(tidyverse)
# # # library(ggplot2)
# # # library(pheatmap)
# # # library(writexl)
# # # library(stringr)
# # # # Try to load ggrepel - install first if not available
# # # if (!requireNamespace("ggrepel", quietly = TRUE)) {
# # #   install.packages("ggrepel")
# # # }
# # # library(ggrepel)
# # # 
# # # # Define the fusion protein sequence
# # # fusion_protein <- "RKREIFDRYGEEVKEFLAKAKEDF"
# # # 
# # # # Identify the fusion junction - between E (DNAJB1) and V (PRKACA)
# # # dnajb1_part <- "RKREIFDRYGEE"  # Ends with YGEE (full DNAJB1 part)
# # # prkaca_part <- "VKEFLAKAKEDF"  # Starts with VKEF (full PRKACA part)
# # # junction_position <- nchar(dnajb1_part)  # Position after the last E in YGEE
# # # 
# # # cat("DNAJB1 part:", dnajb1_part, "\n")
# # # cat("PRKACA part:", prkaca_part, "\n")
# # # cat("Junction position:", junction_position, "\n")
# # # cat("Fusion protein:", fusion_protein, "\n")
# # # cat("Checking junction: character at junction position is", substr(fusion_protein, junction_position, junction_position), "\n")
# # # cat("Checking junction: next character is", substr(fusion_protein, junction_position+1, junction_position+1), "\n\n")
# # # 
# # # # Generate theoretical junction-spanning peptides (8-12mers)
# # # theoretical_peptides <- list()
# # # peptide_lengths <- 8:12  # Looking for 8-12mers
# # # 
# # # for (length in peptide_lengths) {
# # #   for (start_pos in 1:(nchar(fusion_protein) - length + 1)) {
# # #     peptide <- substr(fusion_protein, start_pos, start_pos + length - 1)
# # #     
# # #     # Check if this peptide spans the junction
# # #     # It spans if it includes at least one residue from both proteins
# # #     peptide_end_pos <- start_pos + length - 1
# # #     spans_junction <- start_pos <= junction_position && peptide_end_pos > junction_position
# # #     
# # #     if (spans_junction) {
# # #       theoretical_peptides[[length(theoretical_peptides) + 1]] <- list(
# # #         Peptide = peptide,
# # #         Length = length,
# # #         Start_Position = start_pos,
# # #         End_Position = peptide_end_pos,
# # #         DNAJB1_Part = paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
# # #         PRKACA_Part = paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = ""),
# # #         Visualization = paste0(
# # #           paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
# # #           paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = "")
# # #         )
# # #       )
# # #     }
# # #   }
# # # }
# # # 
# # # # Convert to dataframe
# # # theoretical_junction_peptides <- do.call(rbind, lapply(theoretical_peptides, function(x) {
# # #   data.frame(
# # #     Peptide = x$Peptide,
# # #     Length = x$Length,
# # #     Start_Position = x$Start_Position,
# # #     End_Position = x$End_Position,
# # #     DNAJB1_Part = x$DNAJB1_Part,
# # #     PRKACA_Part = x$PRKACA_Part,
# # #     Visualization = x$Visualization,
# # #     stringsAsFactors = FALSE
# # #   )
# # # }))
# # # 
# # # cat("Generated", nrow(theoretical_junction_peptides), "theoretical junction-spanning peptides\n")
# # # 
# # # # Count by length
# # # for (length in peptide_lengths) {
# # #   count <- sum(theoretical_junction_peptides$Length == length)
# # #   cat("  Length", length, ":", count, "peptides\n")
# # # }
# # # 
# # # # Print all theoretical peptides in a neat table
# # # cat("\nAll theoretical junction-spanning peptides:\n")
# # # for (i in 1:nrow(theoretical_junction_peptides)) {
# # #   peptide <- theoretical_junction_peptides[i,]
# # #   cat(sprintf("%2d. %s (Length: %d, Start: %d, End: %d, Vis: %s)\n", 
# # #               i, peptide$Peptide, peptide$Length, peptide$Start_Position, 
# # #               peptide$End_Position, peptide$Visualization))
# # # }
# # # 
# # # # Read the combined dataset
# # # cat("Reading combined dataset...\n")
# # # data_file <- "unique_peptides_unmodified.tsv"
# # # combined_data <- read.delim(data_file, stringsAsFactors = FALSE)
# # # 
# # # cat("Dataset loaded with", nrow(combined_data), "peptide records\n")
# # # cat("Columns available:", paste(colnames(combined_data), collapse = ", "), "\n\n")
# # # 
# # # # Create output directory for results
# # # results_dir <- "Fusion_Junction_Analysis"
# # # dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
# # # 
# # # #===============================#
# # # # Function to check if a peptide spans the fusion junction
# # # #===============================#
# # # 
# # # is_junction_spanning <- function(peptide) {
# # #   # First check if it's entirely within the fusion protein
# # #   if (!grepl(peptide, fusion_protein, fixed = TRUE)) {
# # #     return(FALSE)
# # #   }
# # #   
# # #   # Find position in fusion protein
# # #   start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
# # #   if (start_pos == -1) {
# # #     return(FALSE)
# # #   }
# # #   
# # #   end_pos <- start_pos + nchar(peptide) - 1
# # #   
# # #   # Check if it spans the junction
# # #   spans_junction <- start_pos <= junction_position && end_pos > junction_position
# # #   
# # #   return(spans_junction)
# # # }
# # # 
# # # # Filter data for 8-12mers and identify junction-spanning peptides
# # # combined_data <- combined_data %>%
# # #   mutate(
# # #     Peptide_Length = nchar(Peptide),
# # #     Spans_Junction = sapply(Peptide, is_junction_spanning)
# # #   ) %>%
# # #   filter(Peptide_Length >= 8 & Peptide_Length <= 12)  # Only include 8-12mers
# # # 
# # # # Extract detected fusion-spanning peptides
# # # detected_junction_peptides <- combined_data %>%
# # #   filter(Spans_Junction == TRUE) %>%
# # #   group_by(Peptide, Peptide_Length) %>%
# # #   summarize(
# # #     Total_Samples = n(),
# # #     Detected_2CV_Count = sum(!is.na(detected_2cv) & detected_2cv == TRUE, na.rm = TRUE),
# # #     Detected_3CV_Count = sum(!is.na(detected_3cv) & detected_3cv == TRUE, na.rm = TRUE),
# # #     Detected_Both_Count = sum(!is.na(detected_both) & detected_both == TRUE, na.rm = TRUE),
# # #     Avg_Final_Intensity = mean(final_intensity, na.rm = TRUE),
# # #     Avg_Intensity_2CV = mean(Intensity_2cv, na.rm = TRUE),
# # #     Avg_Intensity_3CV = mean(Intensity_3cv, na.rm = TRUE),
# # #     .groups = "drop"
# # #   ) %>%
# # #   arrange(Peptide_Length, Peptide)
# # # 
# # # cat("Found", nrow(detected_junction_peptides), "detected peptides that span the fusion junction\n")
# # # 
# # # # Get summary statistics about the dataset
# # # total_samples <- length(unique(combined_data$SampleID))
# # # samples_with_2cv <- sum(!is.na(combined_data$detected_2cv) & combined_data$detected_2cv == TRUE)
# # # samples_with_3cv <- sum(!is.na(combined_data$detected_3cv) & combined_data$detected_3cv == TRUE)
# # # 
# # # cat("Dataset summary:\n")
# # # cat("  Total unique samples:", total_samples, "\n")
# # # cat("  Detections in 2CV:", samples_with_2cv, "\n")
# # # cat("  Detections in 3CV:", samples_with_3cv, "\n\n")
# # # 
# # # if (nrow(detected_junction_peptides) > 0) {
# # #   # Create detailed analysis of each detected junction peptide
# # #   detected_details <- data.frame()
# # #   
# # #   for (i in 1:nrow(detected_junction_peptides)) {
# # #     peptide <- detected_junction_peptides$Peptide[i]
# # #     start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
# # #     end_pos <- start_pos + nchar(peptide) - 1
# # #     
# # #     # Calculate how many residues come from each protein
# # #     dnajb1_residues <- max(0, min(junction_position - start_pos + 1, nchar(peptide)))
# # #     prkaca_residues <- max(0, min(end_pos - junction_position, nchar(peptide)))
# # #     
# # #     # Visualization string (D for DNAJB1, P for PRKACA)
# # #     vis_string <- paste0(
# # #       paste(rep("D", dnajb1_residues), collapse = ""),
# # #       paste(rep("P", prkaca_residues), collapse = "")
# # #     )
# # #     
# # #     # Add to details dataframe
# # #     detected_details <- rbind(detected_details, data.frame(
# # #       Peptide = peptide,
# # #       Length = nchar(peptide),
# # #       Start_Position = start_pos,
# # #       End_Position = end_pos,
# # #       DNAJB1_Residues = dnajb1_residues,
# # #       PRKACA_Residues = prkaca_residues,
# # #       Visualization = vis_string,
# # #       Total_Samples = detected_junction_peptides$Total_Samples[i],
# # #       Detected_2CV_Count = detected_junction_peptides$Detected_2CV_Count[i],
# # #       Detected_3CV_Count = detected_junction_peptides$Detected_3CV_Count[i],
# # #       Detected_Both_Count = detected_junction_peptides$Detected_Both_Count[i],
# # #       Avg_Final_Intensity = detected_junction_peptides$Avg_Final_Intensity[i],
# # #       Avg_Intensity_2CV = detected_junction_peptides$Avg_Intensity_2CV[i],
# # #       Avg_Intensity_3CV = detected_junction_peptides$Avg_Intensity_3CV[i],
# # #       stringsAsFactors = FALSE
# # #     ))
# # #   }
# # #   
# # #   # Create heatmaps showing detection patterns
# # #   
# # #   # Prepare data for heatmaps - focusing on detected junction peptides
# # #   junction_data <- combined_data %>%
# # #     filter(Spans_Junction == TRUE) %>%
# # #     select(Peptide, SampleID, detected_2cv, detected_3cv, detected_both, 
# # #            final_intensity, Intensity_2cv, Intensity_3cv)
# # #   
# # #   # Create presence matrices for 2CV and 3CV
# # #   presence_2cv_data <- junction_data %>%
# # #     filter(!is.na(detected_2cv)) %>%
# # #     select(Peptide, SampleID, detected_2cv) %>%
# # #     pivot_wider(names_from = SampleID, values_from = detected_2cv, values_fill = FALSE) %>%
# # #     column_to_rownames("Peptide")
# # #   
# # #   presence_3cv_data <- junction_data %>%
# # #     filter(!is.na(detected_3cv)) %>%
# # #     select(Peptide, SampleID, detected_3cv) %>%
# # #     pivot_wider(names_from = SampleID, values_from = detected_3cv, values_fill = FALSE) %>%
# # #     column_to_rownames("Peptide")
# # #   
# # #   # Convert to numeric matrices
# # #   if(nrow(presence_2cv_data) > 0 && ncol(presence_2cv_data) > 0) {
# # #     presence_2cv_matrix <- as.matrix(presence_2cv_data)
# # #     presence_2cv_numeric <- matrix(as.numeric(presence_2cv_matrix), 
# # #                                    nrow = nrow(presence_2cv_matrix),
# # #                                    dimnames = dimnames(presence_2cv_matrix))
# # #     
# # #     # Create 2CV presence heatmap
# # #     pdf(file.path(results_dir, "junction_peptides_presence_2CV.pdf"), width = 12, height = 8)
# # #     pheatmap(
# # #       presence_2cv_numeric,
# # #       main = "Presence of Junction Peptides in 2CV Samples",
# # #       color = c("white", "steelblue"),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.0f"
# # #     )
# # #     dev.off()
# # #     
# # #     png(file.path(results_dir, "junction_peptides_presence_2CV.png"), width = 1000, height = 600, res = 100)
# # #     pheatmap(
# # #       presence_2cv_numeric,
# # #       main = "Presence of Junction Peptides in 2CV Samples",
# # #       color = c("white", "steelblue"),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.0f"
# # #     )
# # #     dev.off()
# # #   }
# # #   
# # #   if(nrow(presence_3cv_data) > 0 && ncol(presence_3cv_data) > 0) {
# # #     presence_3cv_matrix <- as.matrix(presence_3cv_data)
# # #     presence_3cv_numeric <- matrix(as.numeric(presence_3cv_matrix), 
# # #                                    nrow = nrow(presence_3cv_matrix),
# # #                                    dimnames = dimnames(presence_3cv_matrix))
# # #     
# # #     # Create 3CV presence heatmap
# # #     pdf(file.path(results_dir, "junction_peptides_presence_3CV.pdf"), width = 12, height = 8)
# # #     pheatmap(
# # #       presence_3cv_numeric,
# # #       main = "Presence of Junction Peptides in 3CV Samples",
# # #       color = c("white", "steelblue"),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.0f"
# # #     )
# # #     dev.off()
# # #     
# # #     png(file.path(results_dir, "junction_peptides_presence_3CV.png"), width = 1000, height = 600, res = 100)
# # #     pheatmap(
# # #       presence_3cv_numeric,
# # #       main = "Presence of Junction Peptides in 3CV Samples",
# # #       color = c("white", "steelblue"),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.0f"
# # #     )
# # #     dev.off()
# # #   }
# # #   
# # #   # Create intensity heatmaps
# # #   intensity_2cv_data <- junction_data %>%
# # #     filter(!is.na(Intensity_2cv) & Intensity_2cv > 0) %>%
# # #     select(Peptide, SampleID, Intensity_2cv) %>%
# # #     pivot_wider(names_from = SampleID, values_from = Intensity_2cv, values_fill = 0) %>%
# # #     column_to_rownames("Peptide")
# # #   
# # #   intensity_3cv_data <- junction_data %>%
# # #     filter(!is.na(Intensity_3cv) & Intensity_3cv > 0) %>%
# # #     select(Peptide, SampleID, Intensity_3cv) %>%
# # #     pivot_wider(names_from = SampleID, values_from = Intensity_3cv, values_fill = 0) %>%
# # #     column_to_rownames("Peptide")
# # #   
# # #   if(nrow(intensity_2cv_data) > 0 && ncol(intensity_2cv_data) > 0) {
# # #     log_intensity_2cv <- log10(as.matrix(intensity_2cv_data) + 1)
# # #     
# # #     pdf(file.path(results_dir, "junction_peptides_intensity_2CV.pdf"), width = 12, height = 8)
# # #     pheatmap(
# # #       log_intensity_2cv,
# # #       main = "Intensity of Junction Peptides in 2CV Samples (log10)",
# # #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.1f"
# # #     )
# # #     dev.off()
# # #     
# # #     png(file.path(results_dir, "junction_peptides_intensity_2CV.png"), width = 1000, height = 600, res = 100)
# # #     pheatmap(
# # #       log_intensity_2cv,
# # #       main = "Intensity of Junction Peptides in 2CV Samples (log10)",
# # #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.1f"
# # #     )
# # #     dev.off()
# # #   }
# # #   
# # #   if(nrow(intensity_3cv_data) > 0 && ncol(intensity_3cv_data) > 0) {
# # #     log_intensity_3cv <- log10(as.matrix(intensity_3cv_data) + 1)
# # #     
# # #     pdf(file.path(results_dir, "junction_peptides_intensity_3CV.pdf"), width = 12, height = 8)
# # #     pheatmap(
# # #       log_intensity_3cv,
# # #       main = "Intensity of Junction Peptides in 3CV Samples (log10)",
# # #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.1f"
# # #     )
# # #     dev.off()
# # #     
# # #     png(file.path(results_dir, "junction_peptides_intensity_3CV.png"), width = 1000, height = 600, res = 100)
# # #     pheatmap(
# # #       log_intensity_3cv,
# # #       main = "Intensity of Junction Peptides in 3CV Samples (log10)",
# # #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.1f"
# # #     )
# # #     dev.off()
# # #   }
# # #   
# # #   # Compare 2CV vs 3CV detection efficiency for junction peptides
# # #   cv_comparison_plot <- ggplot(detected_details, aes(x = Peptide, y = Length)) +
# # #     geom_point(aes(size = Total_Samples, color = Detected_3CV_Count / (Detected_2CV_Count + 0.001))) +
# # #     scale_color_gradient2(
# # #       low = "blue", 
# # #       mid = "white", 
# # #       high = "red", 
# # #       midpoint = 1,
# # #       name = "3CV/2CV\nDetection Ratio"
# # #     ) +
# # #     theme_minimal() +
# # #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# # #     labs(
# # #       title = "Fusion Junction Peptide Detection: 2CV vs 3CV",
# # #       x = "Peptide Sequence",
# # #       y = "Peptide Length",
# # #       size = "Total\nSamples"
# # #     )
# # #   
# # #   ggsave(file.path(results_dir, "junction_peptide_cv_comparison.pdf"), cv_comparison_plot, width = 12, height = 8)
# # #   ggsave(file.path(results_dir, "junction_peptide_cv_comparison.png"), cv_comparison_plot, width = 12, height = 8)
# # #   
# # #   # Create a visualization showing peptide coverage across the fusion protein
# # #   peptide_coverage_plot <- ggplot(detected_details, 
# # #                                   aes(x = Start_Position, xend = End_Position, 
# # #                                       y = reorder(Peptide, Start_Position), yend = reorder(Peptide, Start_Position))) +
# # #     geom_segment(aes(color = Total_Samples), linewidth = 5) +
# # #     geom_vline(xintercept = junction_position, linetype = "dashed", color = "red") +
# # #     annotate("text", x = junction_position - 2, y = nrow(detected_details) + 1, 
# # #              label = "DNAJB1", color = "darkgreen", fontface = "bold") +
# # #     annotate("text", x = junction_position + 2, y = nrow(detected_details) + 1, 
# # #              label = "PRKACA", color = "purple", fontface = "bold") +
# # #     scale_color_gradient(low = "lightsteelblue", high = "steelblue") +
# # #     theme_minimal() +
# # #     labs(
# # #       title = "Coverage of Fusion Junction by Detected Peptides",
# # #       x = "Position in Fusion Protein",
# # #       y = "Peptide",
# # #       color = "Total\nSamples"
# # #     )
# # #   
# # #   ggsave(file.path(results_dir, "junction_peptide_coverage.pdf"), peptide_coverage_plot, width = 12, height = 8)
# # #   ggsave(file.path(results_dir, "junction_peptide_coverage.png"), peptide_coverage_plot, width = 12, height = 8)
# # #   
# # #   # Create combined 2CV + 3CV heatmap with annotation
# # #   # Prepare data for combined heatmap
# # #   combined_intensity_data <- junction_data %>%
# # #     filter(!is.na(final_intensity) & final_intensity > 0) %>%
# # #     # Create sample-CV combinations
# # #     mutate(Sample_CV = paste0(SampleID, "_", ifelse(!is.na(detected_2cv) & detected_2cv, "2CV", ""))) %>%
# # #     mutate(Sample_CV = ifelse(Sample_CV == paste0(SampleID, "_"), 
# # #                               paste0(SampleID, "_", ifelse(!is.na(detected_3cv) & detected_3cv, "3CV", "")), 
# # #                               Sample_CV)) %>%
# # #     filter(Sample_CV != paste0(SampleID, "_")) %>%  # Remove entries with no CV type
# # #     select(Peptide, Sample_CV, final_intensity)
# # #   
# # #   # Also handle cases where peptides are detected in both CV types
# # #   combined_2cv_data <- junction_data %>%
# # #     filter(!is.na(detected_2cv) & detected_2cv == TRUE & !is.na(Intensity_2cv) & Intensity_2cv > 0) %>%
# # #     mutate(Sample_CV = paste0(SampleID, "_2CV")) %>%
# # #     select(Peptide, Sample_CV, Intensity = Intensity_2cv)
# # #   
# # #   combined_3cv_data <- junction_data %>%
# # #     filter(!is.na(detected_3cv) & detected_3cv == TRUE & !is.na(Intensity_3cv) & Intensity_3cv > 0) %>%
# # #     mutate(Sample_CV = paste0(SampleID, "_3CV")) %>%
# # #     select(Peptide, Sample_CV, Intensity = Intensity_3cv)
# # #   
# # #   # Combine both datasets
# # #   all_combined_data <- bind_rows(combined_2cv_data, combined_3cv_data)
# # #   
# # #   if(nrow(all_combined_data) > 0) {
# # #     # Create combined intensity matrix
# # #     combined_intensity_matrix <- all_combined_data %>%
# # #       pivot_wider(names_from = Sample_CV, values_from = Intensity, values_fill = 0) %>%
# # #       column_to_rownames("Peptide") %>%
# # #       as.matrix()
# # #     
# # #     # Log transform
# # #     log_combined_intensity <- log10(combined_intensity_matrix + 1)
# # #     
# # #     # Create column annotations for CV type
# # #     column_ann <- data.frame(
# # #       CV_Type = str_extract(colnames(log_combined_intensity), "2CV|3CV"),
# # #       Sample = gsub("_[23]CV$", "", colnames(log_combined_intensity))
# # #     )
# # #     rownames(column_ann) <- colnames(log_combined_intensity)
# # #     
# # #     # Define colors for annotation
# # #     ann_colors <- list(
# # #       CV_Type = c("2CV" = "steelblue", "3CV" = "tomato")
# # #     )
# # #     
# # #     # Create combined heatmap
# # #     pdf(file.path(results_dir, "junction_peptides_intensity_combined.pdf"), width = 16, height = 8)
# # #     pheatmap(
# # #       log_combined_intensity,
# # #       main = "Intensity of Fusion Junction Peptides (2CV vs 3CV, log10)",
# # #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.1f",
# # #       annotation_col = column_ann,
# # #       annotation_colors = ann_colors
# # #     )
# # #     dev.off()
# # #     
# # #     png(file.path(results_dir, "junction_peptides_intensity_combined.png"), width = 1200, height = 600, res = 100)
# # #     pheatmap(
# # #       log_combined_intensity,
# # #       main = "Intensity of Fusion Junction Peptides (2CV vs 3CV, log10)",
# # #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.1f",
# # #       annotation_col = column_ann,
# # #       annotation_colors = ann_colors
# # #     )
# # #     dev.off()
# # #     
# # #     # Also create combined presence heatmap
# # #     combined_presence_2cv <- junction_data %>%
# # #       filter(!is.na(detected_2cv) & detected_2cv == TRUE) %>%
# # #       mutate(Sample_CV = paste0(SampleID, "_2CV")) %>%
# # #       select(Peptide, Sample_CV, detected = detected_2cv)
# # #     
# # #     combined_presence_3cv <- junction_data %>%
# # #       filter(!is.na(detected_3cv) & detected_3cv == TRUE) %>%
# # #       mutate(Sample_CV = paste0(SampleID, "_3CV")) %>%
# # #       select(Peptide, Sample_CV, detected = detected_3cv)
# # #     
# # #     all_combined_presence <- bind_rows(combined_presence_2cv, combined_presence_3cv)
# # #     
# # #     if(nrow(all_combined_presence) > 0) {
# # #       combined_presence_matrix <- all_combined_presence %>%
# # #         pivot_wider(names_from = Sample_CV, values_from = detected, values_fill = FALSE) %>%
# # #         column_to_rownames("Peptide") %>%
# # #         as.matrix()
# # #       
# # #       # Convert to numeric
# # #       combined_presence_numeric <- matrix(as.numeric(combined_presence_matrix), 
# # #                                           nrow = nrow(combined_presence_matrix),
# # #                                           dimnames = dimnames(combined_presence_matrix))
# # #       
# # #       # Create column annotations
# # #       presence_column_ann <- data.frame(
# # #         CV_Type = str_extract(colnames(combined_presence_numeric), "2CV|3CV"),
# # #         Sample = gsub("_[23]CV$", "", colnames(combined_presence_numeric))
# # #       )
# # #       rownames(presence_column_ann) <- colnames(combined_presence_numeric)
# # #       
# # #       pdf(file.path(results_dir, "junction_peptides_presence_combined.pdf"), width = 16, height = 8)
# # #       pheatmap(
# # #         combined_presence_numeric,
# # #         main = "Presence of Fusion Junction Peptides (2CV vs 3CV)",
# # #         color = c("white", "steelblue"),
# # #         cluster_rows = FALSE,
# # #         cluster_cols = TRUE,
# # #         fontsize_row = 10,
# # #         fontsize_col = 8,
# # #         display_numbers = TRUE,
# # #         number_format = "%.0f",
# # #         annotation_col = presence_column_ann,
# # #         annotation_colors = ann_colors
# # #       )
# # #       dev.off()
# # #       
# # #       png(file.path(results_dir, "junction_peptides_presence_combined.png"), width = 1200, height = 600, res = 100)
# # #       pheatmap(
# # #         combined_presence_numeric,
# # #         main = "Presence of Fusion Junction Peptides (2CV vs 3CV)",
# # #         color = c("white", "steelblue"),
# # #         cluster_rows = FALSE,
# # #         cluster_cols = TRUE,
# # #         fontsize_row = 10,
# # #         fontsize_col = 8,
# # #         display_numbers = TRUE,
# # #         number_format = "%.0f",
# # #         annotation_col = presence_column_ann,
# # #         annotation_colors = ann_colors
# # #       )
# # #       dev.off()
# # #     }
# # #   }
# # #   
# # #   # Create intensity comparison plot
# # #   intensity_data <- detected_details %>%
# # #     filter(!is.na(Avg_Intensity_2CV) & !is.na(Avg_Intensity_3CV) & 
# # #              Avg_Intensity_2CV > 0 & Avg_Intensity_3CV > 0)
# # #   
# # #   if(nrow(intensity_data) > 0) {
# # #     intensity_comparison_plot <- ggplot(intensity_data, 
# # #                                         aes(x = Avg_Intensity_2CV + 1, y = Avg_Intensity_3CV + 1)) +
# # #       geom_point(aes(size = Total_Samples, color = Length)) +
# # #       geom_text_repel(aes(label = Peptide), size = 3) +
# # #       scale_x_log10() +
# # #       scale_y_log10() +
# # #       geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
# # #       theme_minimal() +
# # #       labs(
# # #         title = "Average Intensity of Junction Peptides: 2CV vs 3CV",
# # #         x = "Average Intensity in 2CV (log scale)",
# # #         y = "Average Intensity in 3CV (log scale)",
# # #         color = "Peptide\nLength",
# # #         size = "Total\nSamples"
# # #       )
# # #     
# # #     ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.pdf"), intensity_comparison_plot, width = 10, height = 8)
# # #     ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.png"), intensity_comparison_plot, width = 10, height = 8)
# # #   }
# # #   
# # #   #===============================#
# # #   # Create Excel Output
# # #   #===============================#
# # #   
# # #   # Sheet with theoretical junction-spanning peptides
# # #   excel_theoretical <- theoretical_junction_peptides %>%
# # #     arrange(Length, Start_Position)
# # #   
# # #   # Sheet with detected junction peptides details
# # #   excel_detected <- detected_details %>%
# # #     arrange(Length, Start_Position)
# # #   
# # #   # Sheet with sample-level detection of junction peptides
# # #   excel_sample_detection <- combined_data %>%
# # #     filter(Spans_Junction == TRUE) %>%
# # #     select(SampleID, Peptide, Peptide_Length, detected_2cv, detected_3cv, detected_both, 
# # #            final_intensity, Intensity_2cv, Intensity_3cv, Spans_Junction) %>%
# # #     arrange(SampleID, Peptide)
# # #   
# # #   # Sheet with 2CV vs 3CV detection statistics
# # #   excel_cv_comparison <- detected_details %>%
# # #     select(
# # #       Peptide,
# # #       Length,
# # #       DNAJB1_Residues,
# # #       PRKACA_Residues,
# # #       Visualization,
# # #       Total_Samples,
# # #       Detected_2CV_Count,
# # #       Detected_3CV_Count,
# # #       Detected_Both_Count,
# # #       Avg_Final_Intensity,
# # #       Avg_Intensity_2CV,
# # #       Avg_Intensity_3CV
# # #     ) %>%
# # #     mutate(
# # #       `2CV Detection %` = (Detected_2CV_Count / Total_Samples) * 100,
# # #       `3CV Detection %` = (Detected_3CV_Count / Total_Samples) * 100,
# # #       `Both Detection %` = (Detected_Both_Count / Total_Samples) * 100,
# # #       `3CV/2CV Detection Ratio` = (Detected_3CV_Count + 0.001) / (Detected_2CV_Count + 0.001),
# # #       `3CV/2CV Intensity Ratio` = (Avg_Intensity_3CV + 1) / (Avg_Intensity_2CV + 1)
# # #     ) %>%
# # #     arrange(Length, Peptide)
# # #   
# # #   # Create a list of sheets for the Excel file
# # #   excel_sheets <- list(
# # #     "Theoretical_Junction_Peptides" = excel_theoretical,
# # #     "Detected_Junction_Peptides" = excel_detected,
# # #     "Sample_Level_Detection" = excel_sample_detection,
# # #     "2CV_vs_3CV_Comparison" = excel_cv_comparison
# # #   )
# # #   
# # #   # Write Excel file with multiple sheets
# # #   write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
# # #   
# # # } else {
# # #   cat("No junction-spanning peptides were detected in the dataset.\n")
# # #   
# # #   # Create Excel with just theoretical peptides
# # #   excel_sheets <- list(
# # #     "Theoretical_Junction_Peptides" = theoretical_junction_peptides
# # #   )
# # #   
# # #   write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
# # # }
# # # 
# # # # Print summary information
# # # cat("\nAnalysis complete! Results saved to:", results_dir, "\n")
# # # cat("\nThe following files were generated:\n")
# # # cat("1. Fusion_Junction_Peptides_Analysis.xlsx - Excel file with detailed analysis\n")
# # # 
# # # if (nrow(detected_junction_peptides) > 0) {
# # #   cat("2. junction_peptides_presence_2CV.pdf/png - Heatmap showing presence of junction peptides in 2CV samples\n")
# # #   cat("3. junction_peptides_presence_3CV.pdf/png - Heatmap showing presence of junction peptides in 3CV samples\n")
# # #   cat("4. junction_peptides_intensity_2CV.pdf/png - Heatmap showing intensity of junction peptides in 2CV samples\n")
# # #   cat("5. junction_peptides_intensity_3CV.pdf/png - Heatmap showing intensity of junction peptides in 3CV samples\n")
# # #   cat("6. junction_peptides_intensity_combined.pdf/png - Combined heatmap showing intensity with 2CV/3CV annotation\n")
# # #   cat("7. junction_peptides_presence_combined.pdf/png - Combined heatmap showing presence with 2CV/3CV annotation\n")
# # #   cat("8. junction_peptide_cv_comparison.pdf/png - Plot comparing detection in 2CV vs 3CV\n")
# # #   cat("9. junction_peptide_coverage.pdf/png - Plot showing coverage of fusion junction by detected peptides\n")
# # #   cat("10. junction_peptide_intensity_comparison.pdf/png - Plot comparing peptide intensities between 2CV and 3CV\n")
# # # }
# # # 
# # # cat("\nSummary of detected junction-spanning peptides:\n")
# # # if (nrow(detected_junction_peptides) > 0) {
# # #   for (i in 1:nrow(detected_junction_peptides)) {
# # #     peptide_info <- detected_junction_peptides[i,]
# # #     cat(sprintf("  %s (Length: %d) - 2CV: %d, 3CV: %d, Both: %d\n", 
# # #                 peptide_info$Peptide, peptide_info$Peptide_Length,
# # #                 peptide_info$Detected_2CV_Count, peptide_info$Detected_3CV_Count,
# # #                 peptide_info$Detected_Both_Count))
# # #   }
# # # } else {
# # #   cat("  No junction-spanning peptides detected.\n")
# # # }
# # # # Streamlined Script to analyze DNAJB1-PRKACA fusion protein peptides in 2CV vs 3CV samples
# # # # Looking for 8-12mers that span the fusion junction
# # # # Uses pre-processed combined dataset
# # # 
# # # # Setting directory
# # # setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/")
# # # 
# # # # Load required packages
# # # library(tidyverse)
# # # library(ggplot2)
# # # library(pheatmap)
# # # library(writexl)
# # # library(stringr)
# # # # Try to load ggrepel - install first if not available
# # # if (!requireNamespace("ggrepel", quietly = TRUE)) {
# # #   install.packages("ggrepel")
# # # }
# # # library(ggrepel)
# # # 
# # # # Define the fusion protein sequence
# # # fusion_protein <- "RKREIFDRYGEEVKEFLAKAKEDF"
# # # 
# # # # Identify the fusion junction - between E (DNAJB1) and V (PRKACA)
# # # dnajb1_part <- "RKREIFDRYGEE"  # Ends with YGEE (full DNAJB1 part)
# # # prkaca_part <- "VKEFLAKAKEDF"  # Starts with VKEF (full PRKACA part)
# # # junction_position <- nchar(dnajb1_part)  # Position after the last E in YGEE
# # # 
# # # cat("DNAJB1 part:", dnajb1_part, "\n")
# # # cat("PRKACA part:", prkaca_part, "\n")
# # # cat("Junction position:", junction_position, "\n")
# # # cat("Fusion protein:", fusion_protein, "\n")
# # # cat("Checking junction: character at junction position is", substr(fusion_protein, junction_position, junction_position), "\n")
# # # cat("Checking junction: next character is", substr(fusion_protein, junction_position+1, junction_position+1), "\n\n")
# # # 
# # # # Generate theoretical junction-spanning peptides (8-12mers)
# # # theoretical_peptides <- list()
# # # peptide_lengths <- 8:12  # Looking for 8-12mers
# # # 
# # # for (length in peptide_lengths) {
# # #   for (start_pos in 1:(nchar(fusion_protein) - length + 1)) {
# # #     peptide <- substr(fusion_protein, start_pos, start_pos + length - 1)
# # #     
# # #     # Check if this peptide spans the junction
# # #     # It spans if it includes at least one residue from both proteins
# # #     peptide_end_pos <- start_pos + length - 1
# # #     spans_junction <- start_pos <= junction_position && peptide_end_pos > junction_position
# # #     
# # #     if (spans_junction) {
# # #       theoretical_peptides[[length(theoretical_peptides) + 1]] <- list(
# # #         Peptide = peptide,
# # #         Length = length,
# # #         Start_Position = start_pos,
# # #         End_Position = peptide_end_pos,
# # #         DNAJB1_Part = paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
# # #         PRKACA_Part = paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = ""),
# # #         Visualization = paste0(
# # #           paste(rep("D", min(junction_position - start_pos + 1, length)), collapse = ""),
# # #           paste(rep("P", max(0, peptide_end_pos - junction_position)), collapse = "")
# # #         )
# # #       )
# # #     }
# # #   }
# # # }
# # # 
# # # # Convert to dataframe
# # # theoretical_junction_peptides <- do.call(rbind, lapply(theoretical_peptides, function(x) {
# # #   data.frame(
# # #     Peptide = x$Peptide,
# # #     Length = x$Length,
# # #     Start_Position = x$Start_Position,
# # #     End_Position = x$End_Position,
# # #     DNAJB1_Part = x$DNAJB1_Part,
# # #     PRKACA_Part = x$PRKACA_Part,
# # #     Visualization = x$Visualization,
# # #     stringsAsFactors = FALSE
# # #   )
# # # }))
# # # 
# # # cat("Generated", nrow(theoretical_junction_peptides), "theoretical junction-spanning peptides\n")
# # # 
# # # # Count by length
# # # for (length in peptide_lengths) {
# # #   count <- sum(theoretical_junction_peptides$Length == length)
# # #   cat("  Length", length, ":", count, "peptides\n")
# # # }
# # # 
# # # # Print all theoretical peptides in a neat table
# # # cat("\nAll theoretical junction-spanning peptides:\n")
# # # for (i in 1:nrow(theoretical_junction_peptides)) {
# # #   peptide <- theoretical_junction_peptides[i,]
# # #   cat(sprintf("%2d. %s (Length: %d, Start: %d, End: %d, Vis: %s)\n", 
# # #               i, peptide$Peptide, peptide$Length, peptide$Start_Position, 
# # #               peptide$End_Position, peptide$Visualization))
# # # }
# # # 
# # # # Read the combined dataset
# # # cat("Reading combined dataset...\n")
# # # data_file <- "unique_peptides_unmodified.tsv"
# # # combined_data <- read.delim(data_file, stringsAsFactors = FALSE)
# # # 
# # # cat("Dataset loaded with", nrow(combined_data), "peptide records\n")
# # # cat("Columns available:", paste(colnames(combined_data), collapse = ", "), "\n\n")
# # # 
# # # # Create output directory for results
# # # results_dir <- "Fusion_Junction_Analysis"
# # # dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
# # # 
# # # #===============================#
# # # # Function to check if a peptide spans the fusion junction
# # # #===============================#
# # # 
# # # is_junction_spanning <- function(peptide) {
# # #   # First check if it's entirely within the fusion protein
# # #   if (!grepl(peptide, fusion_protein, fixed = TRUE)) {
# # #     return(FALSE)
# # #   }
# # #   
# # #   # Find position in fusion protein
# # #   start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
# # #   if (start_pos == -1) {
# # #     return(FALSE)
# # #   }
# # #   
# # #   end_pos <- start_pos + nchar(peptide) - 1
# # #   
# # #   # Check if it spans the junction
# # #   spans_junction <- start_pos <= junction_position && end_pos > junction_position
# # #   
# # #   return(spans_junction)
# # # }
# # # 
# # # # Filter data for 8-12mers and identify junction-spanning peptides
# # # combined_data <- combined_data %>%
# # #   mutate(
# # #     Peptide_Length = nchar(Peptide),
# # #     Spans_Junction = sapply(Peptide, is_junction_spanning)
# # #   ) %>%
# # #   filter(Peptide_Length >= 8 & Peptide_Length <= 12)  # Only include 8-12mers
# # # 
# # # # Extract detected fusion-spanning peptides
# # # detected_junction_peptides <- combined_data %>%
# # #   filter(Spans_Junction == TRUE) %>%
# # #   group_by(Peptide, Peptide_Length) %>%
# # #   summarize(
# # #     Total_Samples = n(),
# # #     Detected_2CV_Count = sum(!is.na(detected_2cv) & detected_2cv == TRUE, na.rm = TRUE),
# # #     Detected_3CV_Count = sum(!is.na(detected_3cv) & detected_3cv == TRUE, na.rm = TRUE),
# # #     Detected_Both_Count = sum(!is.na(detected_both) & detected_both == TRUE, na.rm = TRUE),
# # #     Avg_Final_Intensity = mean(final_intensity, na.rm = TRUE),
# # #     Avg_Intensity_2CV = mean(Intensity_2cv, na.rm = TRUE),
# # #     Avg_Intensity_3CV = mean(Intensity_3cv, na.rm = TRUE),
# # #     .groups = "drop"
# # #   ) %>%
# # #   arrange(Peptide_Length, Peptide)
# # # 
# # # cat("Found", nrow(detected_junction_peptides), "detected peptides that span the fusion junction\n")
# # # 
# # # # Get summary statistics about the dataset
# # # total_samples <- length(unique(combined_data$SampleID))
# # # samples_with_2cv <- sum(!is.na(combined_data$detected_2cv) & combined_data$detected_2cv == TRUE)
# # # samples_with_3cv <- sum(!is.na(combined_data$detected_3cv) & combined_data$detected_3cv == TRUE)
# # # 
# # # cat("Dataset summary:\n")
# # # cat("  Total unique samples:", total_samples, "\n")
# # # cat("  Detections in 2CV:", samples_with_2cv, "\n")
# # # cat("  Detections in 3CV:", samples_with_3cv, "\n\n")
# # # 
# # # if (nrow(detected_junction_peptides) > 0) {
# # #   # Create detailed analysis of each detected junction peptide
# # #   detected_details <- data.frame()
# # #   
# # #   for (i in 1:nrow(detected_junction_peptides)) {
# # #     peptide <- detected_junction_peptides$Peptide[i]
# # #     start_pos <- regexpr(peptide, fusion_protein, fixed = TRUE)[1]
# # #     end_pos <- start_pos + nchar(peptide) - 1
# # #     
# # #     # Calculate how many residues come from each protein
# # #     dnajb1_residues <- max(0, min(junction_position - start_pos + 1, nchar(peptide)))
# # #     prkaca_residues <- max(0, min(end_pos - junction_position, nchar(peptide)))
# # #     
# # #     # Visualization string (D for DNAJB1, P for PRKACA)
# # #     vis_string <- paste0(
# # #       paste(rep("D", dnajb1_residues), collapse = ""),
# # #       paste(rep("P", prkaca_residues), collapse = "")
# # #     )
# # #     
# # #     # Add to details dataframe
# # #     detected_details <- rbind(detected_details, data.frame(
# # #       Peptide = peptide,
# # #       Length = nchar(peptide),
# # #       Start_Position = start_pos,
# # #       End_Position = end_pos,
# # #       DNAJB1_Residues = dnajb1_residues,
# # #       PRKACA_Residues = prkaca_residues,
# # #       Visualization = vis_string,
# # #       Total_Samples = detected_junction_peptides$Total_Samples[i],
# # #       Detected_2CV_Count = detected_junction_peptides$Detected_2CV_Count[i],
# # #       Detected_3CV_Count = detected_junction_peptides$Detected_3CV_Count[i],
# # #       Detected_Both_Count = detected_junction_peptides$Detected_Both_Count[i],
# # #       Avg_Final_Intensity = detected_junction_peptides$Avg_Final_Intensity[i],
# # #       Avg_Intensity_2CV = detected_junction_peptides$Avg_Intensity_2CV[i],
# # #       Avg_Intensity_3CV = detected_junction_peptides$Avg_Intensity_3CV[i],
# # #       stringsAsFactors = FALSE
# # #     ))
# # #   }
# # #   
# # #   # Create heatmaps showing detection patterns
# # #   
# # #   # Prepare data for heatmaps - focusing on detected junction peptides
# # #   junction_data <- combined_data %>%
# # #     filter(Spans_Junction == TRUE) %>%
# # #     select(Peptide, SampleID, detected_2cv, detected_3cv, detected_both, 
# # #            final_intensity, Intensity_2cv, Intensity_3cv)
# # #   
# # #   # Create presence matrices for 2CV and 3CV
# # #   presence_2cv_data <- junction_data %>%
# # #     filter(!is.na(detected_2cv)) %>%
# # #     select(Peptide, SampleID, detected_2cv) %>%
# # #     pivot_wider(names_from = SampleID, values_from = detected_2cv, values_fill = FALSE) %>%
# # #     column_to_rownames("Peptide")
# # #   
# # #   presence_3cv_data <- junction_data %>%
# # #     filter(!is.na(detected_3cv)) %>%
# # #     select(Peptide, SampleID, detected_3cv) %>%
# # #     pivot_wider(names_from = SampleID, values_from = detected_3cv, values_fill = FALSE) %>%
# # #     column_to_rownames("Peptide")
# # #   
# # #   # Convert to numeric matrices
# # #   if(nrow(presence_2cv_data) > 0 && ncol(presence_2cv_data) > 0) {
# # #     presence_2cv_matrix <- as.matrix(presence_2cv_data)
# # #     presence_2cv_numeric <- matrix(as.numeric(presence_2cv_matrix), 
# # #                                    nrow = nrow(presence_2cv_matrix),
# # #                                    dimnames = dimnames(presence_2cv_matrix))
# # #     
# # #     # Create 2CV presence heatmap
# # #     pdf(file.path(results_dir, "junction_peptides_presence_2CV.pdf"), width = 12, height = 8)
# # #     pheatmap(
# # #       presence_2cv_numeric,
# # #       main = "Presence of Junction Peptides in 2CV Samples",
# # #       color = c("white", "steelblue"),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.0f"
# # #     )
# # #     dev.off()
# # #     
# # #     png(file.path(results_dir, "junction_peptides_presence_2CV.png"), width = 1000, height = 600, res = 100)
# # #     pheatmap(
# # #       presence_2cv_numeric,
# # #       main = "Presence of Junction Peptides in 2CV Samples",
# # #       color = c("white", "steelblue"),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.0f"
# # #     )
# # #     dev.off()
# # #   }
# # #   
# # #   if(nrow(presence_3cv_data) > 0 && ncol(presence_3cv_data) > 0) {
# # #     presence_3cv_matrix <- as.matrix(presence_3cv_data)
# # #     presence_3cv_numeric <- matrix(as.numeric(presence_3cv_matrix), 
# # #                                    nrow = nrow(presence_3cv_matrix),
# # #                                    dimnames = dimnames(presence_3cv_matrix))
# # #     
# # #     # Create 3CV presence heatmap
# # #     pdf(file.path(results_dir, "junction_peptides_presence_3CV.pdf"), width = 12, height = 8)
# # #     pheatmap(
# # #       presence_3cv_numeric,
# # #       main = "Presence of Junction Peptides in 3CV Samples",
# # #       color = c("white", "steelblue"),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.0f"
# # #     )
# # #     dev.off()
# # #     
# # #     png(file.path(results_dir, "junction_peptides_presence_3CV.png"), width = 1000, height = 600, res = 100)
# # #     pheatmap(
# # #       presence_3cv_numeric,
# # #       main = "Presence of Junction Peptides in 3CV Samples",
# # #       color = c("white", "steelblue"),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.0f"
# # #     )
# # #     dev.off()
# # #   }
# # #   
# # #   # Create intensity heatmaps
# # #   intensity_2cv_data <- junction_data %>%
# # #     filter(!is.na(Intensity_2cv) & Intensity_2cv > 0) %>%
# # #     select(Peptide, SampleID, Intensity_2cv) %>%
# # #     pivot_wider(names_from = SampleID, values_from = Intensity_2cv, values_fill = 0) %>%
# # #     column_to_rownames("Peptide")
# # #   
# # #   intensity_3cv_data <- junction_data %>%
# # #     filter(!is.na(Intensity_3cv) & Intensity_3cv > 0) %>%
# # #     select(Peptide, SampleID, Intensity_3cv) %>%
# # #     pivot_wider(names_from = SampleID, values_from = Intensity_3cv, values_fill = 0) %>%
# # #     column_to_rownames("Peptide")
# # #   
# # #   if(nrow(intensity_2cv_data) > 0 && ncol(intensity_2cv_data) > 0) {
# # #     log_intensity_2cv <- log10(as.matrix(intensity_2cv_data) + 1)
# # #     
# # #     pdf(file.path(results_dir, "junction_peptides_intensity_2CV.pdf"), width = 12, height = 8)
# # #     pheatmap(
# # #       log_intensity_2cv,
# # #       main = "Intensity of Junction Peptides in 2CV Samples (log10)",
# # #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.1f"
# # #     )
# # #     dev.off()
# # #     
# # #     png(file.path(results_dir, "junction_peptides_intensity_2CV.png"), width = 1000, height = 600, res = 100)
# # #     pheatmap(
# # #       log_intensity_2cv,
# # #       main = "Intensity of Junction Peptides in 2CV Samples (log10)",
# # #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.1f"
# # #     )
# # #     dev.off()
# # #   }
# # #   
# # #   if(nrow(intensity_3cv_data) > 0 && ncol(intensity_3cv_data) > 0) {
# # #     log_intensity_3cv <- log10(as.matrix(intensity_3cv_data) + 1)
# # #     
# # #     pdf(file.path(results_dir, "junction_peptides_intensity_3CV.pdf"), width = 12, height = 8)
# # #     pheatmap(
# # #       log_intensity_3cv,
# # #       main = "Intensity of Junction Peptides in 3CV Samples (log10)",
# # #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.1f"
# # #     )
# # #     dev.off()
# # #     
# # #     png(file.path(results_dir, "junction_peptides_intensity_3CV.png"), width = 1000, height = 600, res = 100)
# # #     pheatmap(
# # #       log_intensity_3cv,
# # #       main = "Intensity of Junction Peptides in 3CV Samples (log10)",
# # #       color = colorRampPalette(c("white", "lightsteelblue", "steelblue"))(100),
# # #       cluster_rows = FALSE,
# # #       cluster_cols = TRUE,
# # #       fontsize_row = 10,
# # #       fontsize_col = 8,
# # #       display_numbers = TRUE,
# # #       number_format = "%.1f"
# # #     )
# # #     dev.off()
# # #   }
# # #   
# # #   # Compare 2CV vs 3CV detection efficiency for junction peptides
# # #   cv_comparison_plot <- ggplot(detected_details, aes(x = Peptide, y = Length)) +
# # #     geom_point(aes(size = Total_Samples, color = Detected_3CV_Count / (Detected_2CV_Count + 0.001))) +
# # #     scale_color_gradient2(
# # #       low = "blue", 
# # #       mid = "white", 
# # #       high = "red", 
# # #       midpoint = 1,
# # #       name = "3CV/2CV\nDetection Ratio"
# # #     ) +
# # #     theme_minimal() +
# # #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# # #     labs(
# # #       title = "Fusion Junction Peptide Detection: 2CV vs 3CV",
# # #       x = "Peptide Sequence",
# # #       y = "Peptide Length",
# # #       size = "Total\nSamples"
# # #     )
# # #   
# # #   ggsave(file.path(results_dir, "junction_peptide_cv_comparison.pdf"), cv_comparison_plot, width = 12, height = 8)
# # #   ggsave(file.path(results_dir, "junction_peptide_cv_comparison.png"), cv_comparison_plot, width = 12, height = 8)
# # #   
# # #   # Create a visualization showing peptide coverage across the fusion protein
# # #   peptide_coverage_plot <- ggplot(detected_details, 
# # #                                   aes(x = Start_Position, xend = End_Position, 
# # #                                       y = reorder(Peptide, Start_Position), yend = reorder(Peptide, Start_Position))) +
# # #     geom_segment(aes(color = Total_Samples), size = 5) +
# # #     geom_vline(xintercept = junction_position, linetype = "dashed", color = "red") +
# # #     annotate("text", x = junction_position - 2, y = nrow(detected_details) + 1, 
# # #              label = "DNAJB1", color = "darkgreen", fontface = "bold") +
# # #     annotate("text", x = junction_position + 2, y = nrow(detected_details) + 1, 
# # #              label = "PRKACA", color = "purple", fontface = "bold") +
# # #     scale_color_gradient(low = "lightsteelblue", high = "steelblue") +
# # #     theme_minimal() +
# # #     labs(
# # #       title = "Coverage of Fusion Junction by Detected Peptides",
# # #       x = "Position in Fusion Protein",
# # #       y = "Peptide",
# # #       color = "Total\nSamples"
# # #     )
# # #   
# # #   ggsave(file.path(results_dir, "junction_peptide_coverage.pdf"), peptide_coverage_plot, width = 12, height = 8)
# # #   ggsave(file.path(results_dir, "junction_peptide_coverage.png"), peptide_coverage_plot, width = 12, height = 8)
# # #   
# # #   # Create intensity comparison plot
# # #   intensity_data <- detected_details %>%
# # #     filter(!is.na(Avg_Intensity_2CV) & !is.na(Avg_Intensity_3CV) & 
# # #              Avg_Intensity_2CV > 0 & Avg_Intensity_3CV > 0)
# # #   
# # #   if(nrow(intensity_data) > 0) {
# # #     intensity_comparison_plot <- ggplot(intensity_data, 
# # #                                         aes(x = Avg_Intensity_2CV + 1, y = Avg_Intensity_3CV + 1)) +
# # #       geom_point(aes(size = Total_Samples, color = Length)) +
# # #       geom_text_repel(aes(label = Peptide), size = 3) +
# # #       scale_x_log10() +
# # #       scale_y_log10() +
# # #       geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
# # #       theme_minimal() +
# # #       labs(
# # #         title = "Average Intensity of Junction Peptides: 2CV vs 3CV",
# # #         x = "Average Intensity in 2CV (log scale)",
# # #         y = "Average Intensity in 3CV (log scale)",
# # #         color = "Peptide\nLength",
# # #         size = "Total\nSamples"
# # #       )
# # #     
# # #     ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.pdf"), intensity_comparison_plot, width = 10, height = 8)
# # #     ggsave(file.path(results_dir, "junction_peptide_intensity_comparison.png"), intensity_comparison_plot, width = 10, height = 8)
# # #   }
# # #   
# # #   #===============================#
# # #   # Create Excel Output
# # #   #===============================#
# # #   
# # #   # Sheet with theoretical junction-spanning peptides
# # #   excel_theoretical <- theoretical_junction_peptides %>%
# # #     arrange(Length, Start_Position)
# # #   
# # #   # Sheet with detected junction peptides details
# # #   excel_detected <- detected_details %>%
# # #     arrange(Length, Start_Position)
# # #   
# # #   # Sheet with sample-level detection of junction peptides
# # #   excel_sample_detection <- combined_data %>%
# # #     filter(Spans_Junction == TRUE) %>%
# # #     select(SampleID, Peptide, Peptide_Length, detected_2cv, detected_3cv, detected_both, 
# # #            final_intensity, Intensity_2cv, Intensity_3cv, Spans_Junction) %>%
# # #     arrange(SampleID, Peptide)
# # #   
# # #   # Sheet with 2CV vs 3CV detection statistics
# # #   excel_cv_comparison <- detected_details %>%
# # #     select(
# # #       Peptide,
# # #       Length,
# # #       DNAJB1_Residues,
# # #       PRKACA_Residues,
# # #       Visualization,
# # #       Total_Samples,
# # #       Detected_2CV_Count,
# # #       Detected_3CV_Count,
# # #       Detected_Both_Count,
# # #       Avg_Final_Intensity,
# # #       Avg_Intensity_2CV,
# # #       Avg_Intensity_3CV
# # #     ) %>%
# # #     mutate(
# # #       `2CV Detection %` = (Detected_2CV_Count / Total_Samples) * 100,
# # #       `3CV Detection %` = (Detected_3CV_Count / Total_Samples) * 100,
# # #       `Both Detection %` = (Detected_Both_Count / Total_Samples) * 100,
# # #       `3CV/2CV Detection Ratio` = (Detected_3CV_Count + 0.001) / (Detected_2CV_Count + 0.001),
# # #       `3CV/2CV Intensity Ratio` = (Avg_Intensity_3CV + 1) / (Avg_Intensity_2CV + 1)
# # #     ) %>%
# # #     arrange(Length, Peptide)
# # #   
# # #   # Create a list of sheets for the Excel file
# # #   excel_sheets <- list(
# # #     "Theoretical_Junction_Peptides" = excel_theoretical,
# # #     "Detected_Junction_Peptides" = excel_detected,
# # #     "Sample_Level_Detection" = excel_sample_detection,
# # #     "2CV_vs_3CV_Comparison" = excel_cv_comparison
# # #   )
# # #   
# # #   # Write Excel file with multiple sheets
# # #   write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
# # #   
# # # } else {
# # #   cat("No junction-spanning peptides were detected in the dataset.\n")
# # #   
# # #   # Create Excel with just theoretical peptides
# # #   excel_sheets <- list(
# # #     "Theoretical_Junction_Peptides" = theoretical_junction_peptides
# # #   )
# # #   
# # #   write_xlsx(excel_sheets, path = file.path(results_dir, "Fusion_Junction_Peptides_Analysis.xlsx"))
# # # }
# # # 
# # # # Print summary information
# # # cat("\nAnalysis complete! Results saved to:", results_dir, "\n")
# # # cat("\nThe following files were generated:\n")
# # # cat("1. Fusion_Junction_Peptides_Analysis.xlsx - Excel file with detailed analysis\n")
# # # 
# # # if (nrow(detected_junction_peptides) > 0) {
# # #   cat("2. junction_peptides_presence_2CV.pdf/png - Heatmap showing presence of junction peptides in 2CV samples\n")
# # #   cat("3. junction_peptides_presence_3CV.pdf/png - Heatmap showing presence of junction peptides in 3CV samples\n")
# # #   cat("4. junction_peptides_intensity_2CV.pdf/png - Heatmap showing intensity of junction peptides in 2CV samples\n")
# # #   cat("5. junction_peptides_intensity_3CV.pdf/png - Heatmap showing intensity of junction peptides in 3CV samples\n")
# # #   cat("6. junction_peptide_cv_comparison.pdf/png - Plot comparing detection in 2CV vs 3CV\n")
# # #   cat("7. junction_peptide_coverage.pdf/png - Plot showing coverage of fusion junction by detected peptides\n")
# # #   cat("8. junction_peptide_intensity_comparison.pdf/png - Plot comparing peptide intensities between 2CV and 3CV\n")
# # # }
# # # 
# # # cat("\nSummary of detected junction-spanning peptides:\n")
# # # if (nrow(detected_junction_peptides) > 0) {
# # #   for (i in 1:nrow(detected_junction_peptides)) {
# # #     peptide_info <- detected_junction_peptides[i,]
# # #     cat(sprintf("  %s (Length: %d) - 2CV: %d, 3CV: %d, Both: %d\n", 
# # #                 peptide_info$Peptide, peptide_info$Peptide_Length,
# # #                 peptide_info$Detected_2CV_Count, peptide_info$Detected_3CV_Count,
# # #                 peptide_info$Detected_Both_Count))
# # #   }
# # # } else {
# # #   cat("  No junction-spanning peptides detected.\n")
# # # }