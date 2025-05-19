# Create Realistic Test Data for Peptide Analysis with Modifications
# This script creates test data files with some peptides having modifications

# Load required packages
library(tidyverse)
library(readr)

# Set up directory for test data
test_dir <- "realistic_test_data"
dir.create(test_dir, showWarnings = FALSE)

# Create sample data with consistent IDs between 2CV and 3CV
# Sample 101
create_sample_101_files <- function() {
  # 2CV data for sample 101
  sample_101_2cv <- tibble(
    Peptide = c("AAAAAVRVC", "AAAALVLKA", "DDDDIAKRM", "MMMMMKVTS"),
    `Prev AA` = c("A", "S", "P", "K"),
    `Next AA` = c("A", "K", "A", "L"),
    Start = c(41, 695, 302, 450),
    End = c(49, 703, 310, 458),
    `Peptide Length` = c(9, 9, 9, 9),
    Charges = c(2, 2, 3, 2),
    Probability = c(0.9638, 0.9966, 0.9892, 0.9788),
    `Spectral Count` = c(1, 2, 3, 1),
    Intensity = c(10000.5, 20000.0, 30000.0, 22000.0),
    `Match Type` = c("MS/MS", "MS/MS", "MS/MS", "MS/MS"),
    `Assigned Modifications` = c("", "", "", "Oxidation@M1; Oxidation@M3"),
    `Observed Modifications` = c("", "", "", "Oxidation@M1; Oxidation@M3"),
    Protein = c("sp|Q9NZV5-2|SELN_HUMAN", "sp|Q9Y490-2|TLN1_HUMAN", 
                "sp|P05556|ITB1_HUMAN", "sp|P02545|LMNA_HUMAN"),
    `Protein ID` = c("Q9NZV5-2", "Q9Y490-2", "P05556", "P02545"),
    `Entry Name` = c("SELN_HUMAN", "TLN1_HUMAN", "ITB1_HUMAN", "LMNA_HUMAN"),
    Gene = c("SELENON", "TLN1", "ITGB1", "LMNA"),
    `Protein Description` = c("Isoform 2 of Selenoprotein N", "Isoform 2 of Talin-1", 
                              "Integrin beta-1", "Prelamin-A/C"),
    `Mapped Genes` = c("", "", "", ""),
    `Mapped Proteins` = c("sp|Q9NZV5|SELN_HUMAN", "sp|Q9Y490|TLN1_HUMAN", "", "")
  )
  
  # 3CV data for sample 101 (contains some of the same peptides with different intensities)
  sample_101_3cv <- tibble(
    Peptide = c("AAAAAVRVC", "AAAALVLKA", "MMMMMKVTS"),
    `Prev AA` = c("A", "S", "K"),
    `Next AA` = c("A", "K", "L"),
    Start = c(41, 695, 450),
    End = c(49, 703, 458),
    `Peptide Length` = c(9, 9, 9),
    Charges = c(2, 2, 2),
    Probability = c(0.9745, 0.9977, 0.9834),
    `Spectral Count` = c(2, 3, 2),
    Intensity = c(15000.0, 25000.0, 28000.0),
    `Match Type` = c("MS/MS", "MS/MS", "MS/MS"),
    `Assigned Modifications` = c("", "", "Oxidation@M1; Oxidation@M3"),
    `Observed Modifications` = c("", "", "Oxidation@M1; Oxidation@M3"),
    Protein = c("sp|Q9NZV5-2|SELN_HUMAN", "sp|Q9Y490-2|TLN1_HUMAN", "sp|P02545|LMNA_HUMAN"),
    `Protein ID` = c("Q9NZV5-2", "Q9Y490-2", "P02545"),
    `Entry Name` = c("SELN_HUMAN", "TLN1_HUMAN", "LMNA_HUMAN"),
    Gene = c("SELENON", "TLN1", "LMNA"),
    `Protein Description` = c("Isoform 2 of Selenoprotein N", "Isoform 2 of Talin-1", "Prelamin-A/C"),
    `Mapped Genes` = c("", "", ""),
    `Mapped Proteins` = c("sp|Q9NZV5|SELN_HUMAN", "sp|Q9Y490|TLN1_HUMAN", "")
  )
  
  # Add source file and sample ID information
  sample_101_2cv$SourceFile <- "DDA_2CV_101_peptides.tsv"
  sample_101_2cv$SampleID <- 101
  sample_101_2cv$CVType <- "2cv"
  
  sample_101_3cv$SourceFile <- "DDA_3CV_101_peptides.tsv"
  sample_101_3cv$SampleID <- 101
  sample_101_3cv$CVType <- "3cv"
  
  # Write files
  write_tsv(sample_101_2cv, file.path(test_dir, "DDA_2CV_101_peptides.tsv"))
  write_tsv(sample_101_3cv, file.path(test_dir, "DDA_3CV_101_peptides.tsv"))
  
  return(list(
    cv2 = sample_101_2cv,
    cv3 = sample_101_3cv
  ))
}

# Sample 102
create_sample_102_files <- function() {
  # 2CV data for sample 102
  sample_102_2cv <- tibble(
    Peptide = c("EEEEFNAAF", "BBBBDASLK", "CCCCCTGEK"),
    `Prev AA` = c("F", "R", "L"),
    `Next AA` = c("T", "L", "V"),
    Start = c(225, 117, 210),
    End = c(233, 125, 218),
    `Peptide Length` = c(9, 9, 9),
    Charges = c(2, 2, 2),
    Probability = c(0.9810, 0.9999, 0.9842),
    `Spectral Count` = c(2, 5, 3),
    Intensity = c(40000.0, 50000.0, 35000.0),
    `Match Type` = c("MS/MS", "MS/MS", "MS/MS"),
    `Assigned Modifications` = c("", "", "Carbamidomethyl@C5"),
    `Observed Modifications` = c("", "", "Carbamidomethyl@C5"),
    Protein = c("sp|P04406|G3P_HUMAN", "sp|P68104|EF1A1_HUMAN", "sp|P04075|ALDOA_HUMAN"),
    `Protein ID` = c("P04406", "P68104", "P04075"),
    `Entry Name` = c("G3P_HUMAN", "EF1A1_HUMAN", "ALDOA_HUMAN"),
    Gene = c("GAPDH", "EEF1A1", "ALDOA"),
    `Protein Description` = c("Glyceraldehyde-3-phosphate dehydrogenase", 
                              "Elongation factor 1-alpha 1",
                              "Fructose-bisphosphate aldolase A"),
    `Mapped Genes` = c("", "", ""),
    `Mapped Proteins` = c("", "sp|Q5VTE0|EF1A3_HUMAN", "")
  )
  
  # 3CV data for sample 102 (contains some of the same peptides plus a new one)
  sample_102_3cv <- tibble(
    Peptide = c("EEEEFNAAF", "BBBBDASLK", "FFFFFQYLA", "CCCCCTGEK"),
    `Prev AA` = c("F", "R", "V", "L"),
    `Next AA` = c("T", "L", "E", "V"),
    Start = c(225, 117, 180, 210),
    End = c(233, 125, 188, 218),
    `Peptide Length` = c(9, 9, 9, 9),
    Charges = c(2, 2, 2, 2),
    Probability = c(0.9856, 0.9990, 0.9912, 0.9778),
    `Spectral Count` = c(3, 4, 2, 1),
    Intensity = c(45000.0, 55000.0, 60000.0, 32000.0),
    `Match Type` = c("MS/MS", "MS/MS", "MS/MS", "MS/MS"),
    `Assigned Modifications` = c("", "", "", "Carbamidomethyl@C5"),
    `Observed Modifications` = c("", "", "", "Carbamidomethyl@C5"),
    Protein = c("sp|P04406|G3P_HUMAN", "sp|P68104|EF1A1_HUMAN", 
                "sp|P27348|1433T_HUMAN", "sp|P04075|ALDOA_HUMAN"),
    `Protein ID` = c("P04406", "P68104", "P27348", "P04075"),
    `Entry Name` = c("G3P_HUMAN", "EF1A1_HUMAN", "1433T_HUMAN", "ALDOA_HUMAN"),
    Gene = c("GAPDH", "EEF1A1", "YWHAQ", "ALDOA"),
    `Protein Description` = c("Glyceraldehyde-3-phosphate dehydrogenase", 
                              "Elongation factor 1-alpha 1", 
                              "14-3-3 protein theta",
                              "Fructose-bisphosphate aldolase A"),
    `Mapped Genes` = c("", "", "", ""),
    `Mapped Proteins` = c("", "sp|Q5VTE0|EF1A3_HUMAN", "", "")
  )
  
  # Add source file and sample ID information
  sample_102_2cv$SourceFile <- "DDA_2CV_102_peptides.tsv"
  sample_102_2cv$SampleID <- 102
  sample_102_2cv$CVType <- "2cv"
  
  sample_102_3cv$SourceFile <- "DDA_3CV_102_peptides.tsv"
  sample_102_3cv$SampleID <- 102
  sample_102_3cv$CVType <- "3cv"
  
  # Write files
  write_tsv(sample_102_2cv, file.path(test_dir, "DDA_2CV_102_peptides.tsv"))
  write_tsv(sample_102_3cv, file.path(test_dir, "DDA_3CV_102_peptides.tsv"))
  
  return(list(
    cv2 = sample_102_2cv,
    cv3 = sample_102_3cv
  ))
}

# Sample 103
create_sample_103_files <- function() {
  # 2CV data for sample 103
  sample_103_2cv <- tibble(
    Peptide = c("VVVVVSPQR", "GGGGGSLLA", "DDDDDASNFK"),
    `Prev AA` = c("K", "R", "A"),
    `Next AA` = c("S", "T", "L"),
    Start = c(125, 311, 278),
    End = c(133, 319, 287),
    `Peptide Length` = c(9, 9, 10),
    Charges = c(2, 2, 2),
    Probability = c(0.9725, 0.9814, 0.9603),
    `Spectral Count` = c(1, 3, 2),
    Intensity = c(22000.0, 33000.0, 18000.0),
    `Match Type` = c("MS/MS", "MS/MS", "MS/MS"),
    `Assigned Modifications` = c("", "", "Phospho@S7"),
    `Observed Modifications` = c("", "", "Phospho@S7"),
    Protein = c("sp|P62258|1433E_HUMAN", "sp|P35579|MYH9_HUMAN", "sp|P06748|NPM_HUMAN"),
    `Protein ID` = c("P62258", "P35579", "P06748"),
    `Entry Name` = c("1433E_HUMAN", "MYH9_HUMAN", "NPM_HUMAN"),
    Gene = c("YWHAE", "MYH9", "NPM1"),
    `Protein Description` = c("14-3-3 protein epsilon", "Myosin-9", "Nucleophosmin"),
    `Mapped Genes` = c("", "", ""),
    `Mapped Proteins` = c("", "", "")
  )
  
  # 3CV data for sample 103
  sample_103_3cv <- tibble(
    Peptide = c("VVVVVSPQR", "GGGGGSLLA", "GGGGGFPSL", "DDDDDASNFK"),
    `Prev AA` = c("K", "R", "L", "A"),
    `Next AA` = c("S", "T", "L", "L"),
    Start = c(125, 311, 401, 278),
    End = c(133, 319, 409, 287),
    `Peptide Length` = c(9, 9, 9, 10),
    Charges = c(2, 2, 3, 2),
    Probability = c(0.9802, 0.9891, 0.9845, 0.9712),
    `Spectral Count` = c(2, 4, 1, 3),
    Intensity = c(28000.0, 38000.0, 70000.0, 26000.0),
    `Match Type` = c("MS/MS", "MS/MS", "MS/MS", "MS/MS"),
    `Assigned Modifications` = c("", "", "", "Phospho@S7"),
    `Observed Modifications` = c("", "", "", "Phospho@S7"),
    Protein = c("sp|P62258|1433E_HUMAN", "sp|P35579|MYH9_HUMAN", 
                "sp|P07900|HS90A_HUMAN", "sp|P06748|NPM_HUMAN"),
    `Protein ID` = c("P62258", "P35579", "P07900", "P06748"),
    `Entry Name` = c("1433E_HUMAN", "MYH9_HUMAN", "HS90A_HUMAN", "NPM_HUMAN"),
    Gene = c("YWHAE", "MYH9", "HSP90AA1", "NPM1"),
    `Protein Description` = c("14-3-3 protein epsilon", "Myosin-9", 
                              "Heat shock protein HSP 90-alpha", "Nucleophosmin"),
    `Mapped Genes` = c("", "", "", ""),
    `Mapped Proteins` = c("", "", "", "")
  )
  
  # Add source file and sample ID information
  sample_103_2cv$SourceFile <- "DDA_2CV_103_peptides.tsv"
  sample_103_2cv$SampleID <- 103
  sample_103_2cv$CVType <- "2cv"
  
  sample_103_3cv$SourceFile <- "DDA_3CV_103_peptides.tsv"
  sample_103_3cv$SampleID <- 103
  sample_103_3cv$CVType <- "3cv"
  
  # Write files
  write_tsv(sample_103_2cv, file.path(test_dir, "DDA_2CV_103_peptides.tsv"))
  write_tsv(sample_103_3cv, file.path(test_dir, "DDA_3CV_103_peptides.tsv"))
  
  return(list(
    cv2 = sample_103_2cv,
    cv3 = sample_103_3cv
  ))
}

# Create all sample files
sample_101 <- create_sample_101_files()
sample_102 <- create_sample_102_files()
sample_103 <- create_sample_103_files()

# Print test files creation confirmation
cat("Realistic test files with modifications created in directory:", test_dir, "\n")
cat("Files:\n")
cat("- DDA_2CV_101_peptides.tsv:", nrow(sample_101$cv2), "peptides\n")
cat("- DDA_3CV_101_peptides.tsv:", nrow(sample_101$cv3), "peptides\n")
cat("- DDA_2CV_102_peptides.tsv:", nrow(sample_102$cv2), "peptides\n")
cat("- DDA_3CV_102_peptides.tsv:", nrow(sample_102$cv3), "peptides\n")
cat("- DDA_2CV_103_peptides.tsv:", nrow(sample_103$cv2), "peptides\n")
cat("- DDA_3CV_103_peptides.tsv:", nrow(sample_103$cv3), "peptides\n")

# Analyze the modified peptides
cat("\nPeptides with modifications:\n")

# Function to count modified peptides in a dataset
count_modified <- function(data) {
  modified <- data %>% 
    filter(`Assigned Modifications` != "" | `Observed Modifications` != "") %>%
    pull(Peptide) %>%
    unique()
  
  return(list(
    count = length(modified),
    peptides = modified
  ))
}

# Count modified peptides in each file
mod_2cv_101 <- count_modified(sample_101$cv2)
mod_3cv_101 <- count_modified(sample_101$cv3)
mod_2cv_102 <- count_modified(sample_102$cv2)
mod_3cv_102 <- count_modified(sample_102$cv3)
mod_2cv_103 <- count_modified(sample_103$cv2)
mod_3cv_103 <- count_modified(sample_103$cv3)

# Print summary
cat("- DDA_2CV_101_peptides.tsv:", mod_2cv_101$count, "modified peptides -", paste(mod_2cv_101$peptides, collapse=", "), "\n")
cat("- DDA_3CV_101_peptides.tsv:", mod_3cv_101$count, "modified peptides -", paste(mod_3cv_101$peptides, collapse=", "), "\n")
cat("- DDA_2CV_102_peptides.tsv:", mod_2cv_102$count, "modified peptides -", paste(mod_2cv_102$peptides, collapse=", "), "\n")
cat("- DDA_3CV_102_peptides.tsv:", mod_3cv_102$count, "modified peptides -", paste(mod_3cv_102$peptides, collapse=", "), "\n")
cat("- DDA_2CV_103_peptides.tsv:", mod_2cv_103$count, "modified peptides -", paste(mod_2cv_103$peptides, collapse=", "), "\n")
cat("- DDA_3CV_103_peptides.tsv:", mod_3cv_103$count, "modified peptides -", paste(mod_3cv_103$peptides, collapse=", "), "\n")

# Get all unique modified peptides
all_modified <- unique(c(
  mod_2cv_101$peptides, mod_3cv_101$peptides,
  mod_2cv_102$peptides, mod_3cv_102$peptides,
  mod_2cv_103$peptides, mod_3cv_103$peptides
))

cat("\nTotal unique modified peptides:", length(all_modified), "-", paste(all_modified, collapse=", "), "\n")

# 
# # Create Realistic Test Data for Peptide Analysis
# 
# # Load required packages
# library(tidyverse)
# library(readr)
# 
# setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516")
# 
# # Set up directory for test data
# test_dir <- "realistic_test_data"
# dir.create(test_dir, showWarnings = FALSE)
# 
# # Create sample data with consistent IDs between 2CV and 3CV
# # Sample 101
# create_sample_101_files <- function() {
#   # 2CV data for sample 101
#   sample_101_2cv <- tibble(
#     Peptide = c("AAAAAVRVC", "AAAALVLKA", "DDDDIAKRM"),
#     `Prev AA` = c("A", "S", "P"),
#     `Next AA` = c("A", "K", "A"),
#     Start = c(41, 695, 302),
#     End = c(49, 703, 310),
#     `Peptide Length` = c(9, 9, 9),
#     Charges = c(2, 2, 3),
#     Probability = c(0.9638, 0.9966, 0.9892),
#     `Spectral Count` = c(1, 2, 3),
#     Intensity = c(10000.5, 20000.0, 30000.0),
#     `Match Type` = c("MS/MS", "MS/MS", "MS/MS"),
#     `Assigned Modifications` = c("", "", ""),
#     `Observed Modifications` = c("", "", ""),
#     Protein = c("sp|Q9NZV5-2|SELN_HUMAN", "sp|Q9Y490-2|TLN1_HUMAN", "sp|P05556|ITB1_HUMAN"),
#     `Protein ID` = c("Q9NZV5-2", "Q9Y490-2", "P05556"),
#     `Entry Name` = c("SELN_HUMAN", "TLN1_HUMAN", "ITB1_HUMAN"),
#     Gene = c("SELENON", "TLN1", "ITGB1"),
#     `Protein Description` = c("Isoform 2 of Selenoprotein N", "Isoform 2 of Talin-1", "Integrin beta-1"),
#     `Mapped Genes` = c("", "", ""),
#     `Mapped Proteins` = c("sp|Q9NZV5|SELN_HUMAN", "sp|Q9Y490|TLN1_HUMAN", "")
#   )
#   
#   # 3CV data for sample 101 (contains some of the same peptides with different intensities)
#   sample_101_3cv <- tibble(
#     Peptide = c("AAAAAVRVC", "AAAALVLKA"),
#     `Prev AA` = c("A", "S"),
#     `Next AA` = c("A", "K"),
#     Start = c(41, 695),
#     End = c(49, 703),
#     `Peptide Length` = c(9, 9),
#     Charges = c(2, 2),
#     Probability = c(0.9745, 0.9977),
#     `Spectral Count` = c(2, 3),
#     Intensity = c(15000.0, 25000.0),
#     `Match Type` = c("MS/MS", "MS/MS"),
#     `Assigned Modifications` = c("", ""),
#     `Observed Modifications` = c("", ""),
#     Protein = c("sp|Q9NZV5-2|SELN_HUMAN", "sp|Q9Y490-2|TLN1_HUMAN"),
#     `Protein ID` = c("Q9NZV5-2", "Q9Y490-2"),
#     `Entry Name` = c("SELN_HUMAN", "TLN1_HUMAN"),
#     Gene = c("SELENON", "TLN1"),
#     `Protein Description` = c("Isoform 2 of Selenoprotein N", "Isoform 2 of Talin-1"),
#     `Mapped Genes` = c("", ""),
#     `Mapped Proteins` = c("sp|Q9NZV5|SELN_HUMAN", "sp|Q9Y490|TLN1_HUMAN")
#   )
#   
#   # Write files
#   write_tsv(sample_101_2cv, file.path(test_dir, "DDA_2CV_101_peptides.tsv"))
#   write_tsv(sample_101_3cv, file.path(test_dir, "DDA_3CV_101_peptides.tsv"))
#   
#   return(list(
#     cv2 = sample_101_2cv,
#     cv3 = sample_101_3cv
#   ))
# }
# 
# # Sample 102
# create_sample_102_files <- function() {
#   # 2CV data for sample 102
#   sample_102_2cv <- tibble(
#     Peptide = c("EEEEFNAAF", "BBBBDASLK"),
#     `Prev AA` = c("F", "R"),
#     `Next AA` = c("T", "L"),
#     Start = c(225, 117),
#     End = c(233, 125),
#     `Peptide Length` = c(9, 9),
#     Charges = c(2, 2),
#     Probability = c(0.9810, 0.9999),
#     `Spectral Count` = c(2, 5),
#     Intensity = c(40000.0, 50000.0),
#     `Match Type` = c("MS/MS", "MS/MS"),
#     `Assigned Modifications` = c("", ""),
#     `Observed Modifications` = c("", ""),
#     Protein = c("sp|P04406|G3P_HUMAN", "sp|P68104|EF1A1_HUMAN"),
#     `Protein ID` = c("P04406", "P68104"),
#     `Entry Name` = c("G3P_HUMAN", "EF1A1_HUMAN"),
#     Gene = c("GAPDH", "EEF1A1"),
#     `Protein Description` = c("Glyceraldehyde-3-phosphate dehydrogenase", "Elongation factor 1-alpha 1"),
#     `Mapped Genes` = c("", ""),
#     `Mapped Proteins` = c("", "sp|Q5VTE0|EF1A3_HUMAN")
#   )
#   
#   # 3CV data for sample 102 (contains some of the same peptides plus a new one)
#   sample_102_3cv <- tibble(
#     Peptide = c("EEEEFNAAF", "BBBBDASLK", "FFFFFQYLA"),
#     `Prev AA` = c("F", "R", "V"),
#     `Next AA` = c("T", "L", "E"),
#     Start = c(225, 117, 180),
#     End = c(233, 125, 188),
#     `Peptide Length` = c(9, 9, 9),
#     Charges = c(2, 2, 2),
#     Probability = c(0.9856, 0.9990, 0.9912),
#     `Spectral Count` = c(3, 4, 2),
#     Intensity = c(45000.0, 55000.0, 60000.0),
#     `Match Type` = c("MS/MS", "MS/MS", "MS/MS"),
#     `Assigned Modifications` = c("", "", ""),
#     `Observed Modifications` = c("", "", ""),
#     Protein = c("sp|P04406|G3P_HUMAN", "sp|P68104|EF1A1_HUMAN", "sp|P27348|1433T_HUMAN"),
#     `Protein ID` = c("P04406", "P68104", "P27348"),
#     `Entry Name` = c("G3P_HUMAN", "EF1A1_HUMAN", "1433T_HUMAN"),
#     Gene = c("GAPDH", "EEF1A1", "YWHAQ"),
#     `Protein Description` = c("Glyceraldehyde-3-phosphate dehydrogenase", "Elongation factor 1-alpha 1", "14-3-3 protein theta"),
#     `Mapped Genes` = c("", "", ""),
#     `Mapped Proteins` = c("", "sp|Q5VTE0|EF1A3_HUMAN", "")
#   )
#   
#   # Write files
#   write_tsv(sample_102_2cv, file.path(test_dir, "DDA_2CV_102_peptides.tsv"))
#   write_tsv(sample_102_3cv, file.path(test_dir, "DDA_3CV_102_peptides.tsv"))
#   
#   return(list(
#     cv2 = sample_102_2cv,
#     cv3 = sample_102_3cv
#   ))
# }
# 
# # Sample 103 (an additional sample)
# create_sample_103_files <- function() {
#   # 2CV data for sample 103
#   sample_103_2cv <- tibble(
#     Peptide = c("VVVVVSPQR", "GGGGGSLLA"),
#     `Prev AA` = c("K", "R"),
#     `Next AA` = c("S", "T"),
#     Start = c(125, 311),
#     End = c(133, 319),
#     `Peptide Length` = c(9, 9),
#     Charges = c(2, 2),
#     Probability = c(0.9725, 0.9814),
#     `Spectral Count` = c(1, 3),
#     Intensity = c(22000.0, 33000.0),
#     `Match Type` = c("MS/MS", "MS/MS"),
#     `Assigned Modifications` = c("", ""),
#     `Observed Modifications` = c("", ""),
#     Protein = c("sp|P62258|1433E_HUMAN", "sp|P35579|MYH9_HUMAN"),
#     `Protein ID` = c("P62258", "P35579"),
#     `Entry Name` = c("1433E_HUMAN", "MYH9_HUMAN"),
#     Gene = c("YWHAE", "MYH9"),
#     `Protein Description` = c("14-3-3 protein epsilon", "Myosin-9"),
#     `Mapped Genes` = c("", ""),
#     `Mapped Proteins` = c("", "")
#   )
#   
#   # 3CV data for sample 103
#   sample_103_3cv <- tibble(
#     Peptide = c("VVVVVSPQR", "GGGGGSLLA", "GGGGGFPSL"),
#     `Prev AA` = c("K", "R", "L"),
#     `Next AA` = c("S", "T", "L"),
#     Start = c(125, 311, 401),
#     End = c(133, 319, 409),
#     `Peptide Length` = c(9, 9, 9),
#     Charges = c(2, 2, 3),
#     Probability = c(0.9802, 0.9891, 0.9845),
#     `Spectral Count` = c(2, 4, 1),
#     Intensity = c(28000.0, 38000.0, 70000.0),
#     `Match Type` = c("MS/MS", "MS/MS", "MS/MS"),
#     `Assigned Modifications` = c("", "", ""),
#     `Observed Modifications` = c("", "", ""),
#     Protein = c("sp|P62258|1433E_HUMAN", "sp|P35579|MYH9_HUMAN", "sp|P07900|HS90A_HUMAN"),
#     `Protein ID` = c("P62258", "P35579", "P07900"),
#     `Entry Name` = c("1433E_HUMAN", "MYH9_HUMAN", "HS90A_HUMAN"),
#     Gene = c("YWHAE", "MYH9", "HSP90AA1"),
#     `Protein Description` = c("14-3-3 protein epsilon", "Myosin-9", "Heat shock protein HSP 90-alpha"),
#     `Mapped Genes` = c("", "", ""),
#     `Mapped Proteins` = c("", "", "")
#   )
#   
#   # Write files
#   write_tsv(sample_103_2cv, file.path(test_dir, "DDA_2CV_103_peptides.tsv"))
#   write_tsv(sample_103_3cv, file.path(test_dir, "DDA_3CV_103_peptides.tsv"))
#   
#   return(list(
#     cv2 = sample_103_2cv,
#     cv3 = sample_103_3cv
#   ))
# }
# 
# # Create all sample files
# sample_101 <- create_sample_101_files()
# sample_102 <- create_sample_102_files()
# sample_103 <- create_sample_103_files()
# 
# # Combine all 2CV data
# all_2cv_data <- bind_rows(
#   sample_101$cv2,
#   sample_102$cv2,
#   sample_103$cv2
# )
# 
# # Combine all 3CV data
# all_3cv_data <- bind_rows(
#   sample_101$cv3,
#   sample_102$cv3,
#   sample_103$cv3
# )
# 
# # Print test files creation confirmation
# cat("Realistic test files created in directory:", test_dir, "\n")
# cat("2CV files:\n")
# cat("- DDA_2CV_101_peptides.tsv:", nrow(sample_101$cv2), "peptides\n")
# cat("- DDA_2CV_102_peptides.tsv:", nrow(sample_102$cv2), "peptides\n")
# cat("- DDA_2CV_103_peptides.tsv:", nrow(sample_103$cv2), "peptides\n")
# cat("3CV files:\n")
# cat("- DDA_3CV_101_peptides.tsv:", nrow(sample_101$cv3), "peptides\n")
# cat("- DDA_3CV_102_peptides.tsv:", nrow(sample_102$cv3), "peptides\n")
# cat("- DDA_3CV_103_peptides.tsv:", nrow(sample_103$cv3), "peptides\n")
# 
# # Analyze the overlap by sample
# analyze_overlap_by_sample <- function(sample_name, cv2_data, cv3_data) {
#   all_peptides <- unique(c(cv2_data$Peptide, cv3_data$Peptide))
#   overlap_peptides <- intersect(cv2_data$Peptide, cv3_data$Peptide)
#   
#   cat("\nSample", sample_name, "Analysis:\n")
#   cat("- Total unique peptides:", length(all_peptides), "\n")
#   cat("- Peptides in 2CV only:", sum(!cv2_data$Peptide %in% cv3_data$Peptide), "\n")
#   cat("- Peptides in 3CV only:", sum(!cv3_data$Peptide %in% cv2_data$Peptide), "\n")
#   cat("- Peptides in both:", length(overlap_peptides), "\n")
#   
#   if (length(overlap_peptides) > 0) {
#     cat("\nOverlapping peptides and their intensity values for Sample", sample_name, ":\n")
#     for (peptide in overlap_peptides) {
#       intensity_2cv <- cv2_data$Intensity[cv2_data$Peptide == peptide]
#       intensity_3cv <- cv3_data$Intensity[cv3_data$Peptide == peptide]
#       
#       cat("- ", peptide, ":\n")
#       cat("  * 2CV: Intensity =", intensity_2cv, "\n")
#       cat("  * 3CV: Intensity =", intensity_3cv, "\n")
#       cat("  * Avg: Intensity =", (intensity_2cv + intensity_3cv)/2, "\n")
#     }
#   }
# }
# 
# # Analyze overall overlap
# all_peptides <- unique(c(all_2cv_data$Peptide, all_3cv_data$Peptide))
# overlap_peptides <- intersect(all_2cv_data$Peptide, all_3cv_data$Peptide)
# 
# cat("\nOverall Analysis:\n")
# cat("- Total unique peptides:", length(all_peptides), "\n")
# cat("- Peptides in 2CV only:", sum(!all_2cv_data$Peptide %in% all_3cv_data$Peptide), "\n")
# cat("- Peptides in 3CV only:", sum(!all_3cv_data$Peptide %in% all_2cv_data$Peptide), "\n")
# cat("- Peptides in both:", length(overlap_peptides), "\n")
# 
# # Analyze overlap for each sample
# analyze_overlap_by_sample("101", sample_101$cv2, sample_101$cv3)
# analyze_overlap_by_sample("102", sample_102$cv2, sample_102$cv3)
# analyze_overlap_by_sample("103", sample_103$cv2, sample_103$cv3)