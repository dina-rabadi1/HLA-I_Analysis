# Cancer Antigen Atlas Integration with HLA-I Immunopeptidome
# caatlas.R
# 
# Cancer Antigen Atlas Integration with HLA-I Immunopeptidome
# Libraries we'll use
library(readxl)
library(dplyr)       # For data manipulation
library(tidyr)       # For data reshaping
library(ggplot2)     # For visualization
library(RColorBrewer) # For color palettes
library(VennDiagram) # For Venn diagrams
library(grid)        # Required for grid.draw

setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis")

# Create output directory structure
if (!dir.exists("results/cancer_atlas")) {
  dir.create("results/cancer_atlas")
}

# Function to create documentation file
create_documentation <- function(file_path, title, description, files_list) {
  sink(file_path)
  cat("# ", title, "\n\n")
  cat(description, "\n\n")
  cat("## Files Generated\n\n")
  for (i in 1:nrow(files_list)) {
    cat("### ", files_list$filename[i], "\n")
    cat(files_list$description[i], "\n\n")
  }
  sink()
}

# 1. Read and preprocess files
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

# Cancer Antigen Atlas files
ca_list <- read.delim("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/caatlas/CancerAssociatedAntigensList.txt", 
                      stringsAsFactors = FALSE)
ptm_list <- read.delim("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/caatlas/PTMAntigenList_SiteLevel.txt", 
                       stringsAsFactors = FALSE)
ct_list <- read.delim("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/caatlas/CTAntigensList.txt", 
                      stringsAsFactors = FALSE)

# 2. Data validation and sanity checks
# Check structure and dimensions
str(ca_list)
str(ptm_list)
str(ct_list)

# Check for potential issues
summary(ca_list)
summary(ptm_list)
summary(ct_list)

# 3. Process peptide data to extract genes
# Create a data frame with peptide-gene pairs
peptide_genes <- data.frame(Peptide = character(), Gene = character(), stringsAsFactors = FALSE)

# Process genes column
for (i in 1:nrow(peptide_data)) {
  peptide <- peptide_data$Peptide[i]
  genes_str <- peptide_data$genes[i]
  
  # Check if genes_str is not NA
  if (!is.na(genes_str) && genes_str != "") {
    # Split genes by semicolon
    genes_list <- strsplit(as.character(genes_str), ";")[[1]]
    genes_list <- trimws(genes_list)  # Remove spaces
    
    # Add each peptide-gene pair
    for (gene in genes_list) {
      if (!is.na(gene) && gene != "") {
        peptide_genes <- rbind(peptide_genes, 
                               data.frame(Peptide = peptide, Gene = gene, stringsAsFactors = FALSE))
      }
    }
  } else {
    # If no gene is available, still record the peptide with NA as gene
    peptide_genes <- rbind(peptide_genes, 
                           data.frame(Peptide = peptide, Gene = NA_character_, stringsAsFactors = FALSE))
  }
}

# Print summary of the processed data
cat("Total peptides processed:", nrow(peptide_data), "\n")
cat("Total peptide-gene pairs created:", nrow(peptide_genes), "\n")
cat("Number of unique peptides in pairs:", length(unique(peptide_genes$Peptide)), "\n")
cat("Number of unique genes in pairs:", length(unique(na.omit(peptide_genes$Gene))), "\n")

# 4. Cross-reference with cancer-associated antigens
ca_matches <- peptide_genes %>%
  inner_join(ca_list, by = c("Gene" = "GeneName"))

# 5. Check for PTM matches - need to compare peptide sequences
ptm_matches <- peptide_data %>%
  inner_join(ptm_list, by = c("Peptide" = "Peptide_Sequence"))

# 6. Cross-reference with cancer-testis antigens
ct_matches <- peptide_genes %>%
  inner_join(ct_list, by = c("Gene" = "GeneName"))

# Add after line 130, after the cross-reference sections

# Create detailed peptide reports for matches
# For cancer-associated antigens (gene-based matches)
if (nrow(ca_matches) > 0) {
  # Create a more readable report with peptides and their matching cancer genes
  ca_peptide_report <- ca_matches %>%
    select(Peptide, Gene, CancerType, SampleNum) %>%
    arrange(CancerType, Gene, Peptide)
  
  # Save detailed report
  write.csv(ca_peptide_report, "results/cancer_atlas/cancer_associated_peptides.csv", row.names = FALSE)
  
  # Create a summary by gene
  ca_gene_summary <- ca_matches %>%
    group_by(Gene, CancerType) %>%
    summarise(
      PeptideCount = n_distinct(Peptide),
      Peptides = paste(unique(Peptide), collapse = ", ")
    ) %>%
    arrange(desc(PeptideCount))
  
  # Save gene summary
  write.csv(ca_gene_summary, "results/cancer_atlas/cancer_associated_genes.csv", row.names = FALSE)
  
  # Print some examples
  cat("\nExample cancer-associated peptides (gene-based matches):\n")
  print(head(ca_peptide_report, 10))
  
  cat("\nNote: Cancer-associated matches are based on gene names, not peptide sequences.\n")
  cat("The peptides listed are from your dataset and are derived from genes found in the Cancer Atlas.\n")
}

# For PTM matches (peptide sequence-based)
if (nrow(ptm_matches) > 0) {
  # Create a readable report with peptides and their modifications
  ptm_peptide_report <- ptm_matches %>%
    select(Peptide, Modification, Site, AA, GeneName) %>%
    arrange(Modification, Peptide)
  
  # Save detailed report
  write.csv(ptm_peptide_report, "results/cancer_atlas/ptm_peptides.csv", row.names = FALSE)
  
  # Print some examples
  cat("\nExample PTM-modified peptides (sequence-based matches):\n")
  print(head(ptm_peptide_report, 10))
  
  cat("\nNote: PTM matches are based on exact peptide sequence matches.\n")
}

# For cancer-testis antigen matches (gene-based)
if (nrow(ct_matches) > 0) {
  # Create a readable report
  ct_peptide_report <- ct_matches %>%
    select(Peptide, Gene, Index) %>%
    arrange(Index, Gene, Peptide)
  
  # Save detailed report
  write.csv(ct_peptide_report, "results/cancer_atlas/cancer_testis_peptides.csv", row.names = FALSE)
  
  # Print some examples
  cat("\nExample cancer-testis antigen peptides (gene-based matches):\n")
  print(head(ct_peptide_report, 10))
  
  cat("\nNote: Cancer-testis antigen matches are based on gene names, not peptide sequences.\n")
}

# 7. Calculate summary statistics
ca_summary <- ca_matches %>%
  group_by(CancerType) %>%
  summarise(count = n_distinct(Peptide)) %>%
  arrange(desc(count))

ct_summary <- ct_matches %>%
  group_by(Index) %>%
  summarise(count = n_distinct(Peptide)) %>%
  arrange(desc(count))

# 8. Visualization
# Bar chart of cancer-associated antigen matches by cancer type
ca_plot <- ggplot(ca_summary, aes(x = reorder(CancerType, count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Cancer-Associated Antigens by Cancer Type",
       x = "Cancer Type", 
       y = "Number of Matched Peptides") +
  scale_y_continuous(breaks = function(x) seq(0, max(x), by = 1))  # Force integer breaks


print(ca_plot)
ggsave("results/cancer_atlas/cancer_associated_barplot.pdf", ca_plot, width = 10, height = 6)

# Bar chart for cancer-testis antigen matches
ct_plot <- ggplot(ct_summary, aes(x = reorder(Index, count), y = count)) +
  geom_bar(stat = "identity", fill = "darkred") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Cancer-Testis Antigens by Index",
       x = "Index", 
       y = "Number of Matched Peptides")

print(ct_plot)
ggsave("results/cancer_atlas/cancer_testis_barplot.pdf", ct_plot, width = 8, height = 5)

# 9. Create a text-based overlap summary
all_genes <- unique(peptide_genes$Gene)
ca_genes <- unique(ca_list$GeneName)
ct_genes <- unique(ct_list$GeneName)

ca_overlap <- length(intersect(all_genes, ca_genes))
ct_overlap <- length(intersect(all_genes, ct_genes))
both_overlap <- length(intersect(intersect(all_genes, ca_genes), ct_genes))

cat("===============================================\n")
cat("OVERLAP SUMMARY\n")
cat("===============================================\n")
cat("Total genes in peptide dataset:", length(all_genes), "\n")
cat("Total cancer-associated genes:", length(ca_genes), "\n")
cat("Total cancer-testis genes:", length(ct_genes), "\n")
cat("\n")
cat("Genes in peptide dataset that are cancer-associated:", ca_overlap, 
    "(", round(ca_overlap/length(all_genes)*100, 2), "%)\n")
cat("Genes in peptide dataset that are cancer-testis:", ct_overlap, 
    "(", round(ct_overlap/length(all_genes)*100, 2), "%)\n")
cat("Genes in peptide dataset that are both cancer-associated and cancer-testis:", both_overlap, 
    "(", round(both_overlap/length(all_genes)*100, 2), "%)\n")

# 10. Create a Venn diagram to visualize gene overlaps
# Filter out NA values
all_genes_clean <- all_genes[!is.na(all_genes)]
ca_genes_clean <- ca_genes[!is.na(ca_genes)]
ct_genes_clean <- ct_genes[!is.na(ct_genes)]

# Check if we have data for the Venn diagram
if (length(all_genes_clean) > 0 && length(ca_genes_clean) > 0 && length(ct_genes_clean) > 0) {
  venn_list <- list(
    "Peptide Genes" = all_genes_clean,
    "Cancer-Associated" = ca_genes_clean,
    "Cancer-Testis" = ct_genes_clean
  )
  
  # Try to create the Venn diagram
  tryCatch({
    venn_plot <- venn.diagram(
      x = venn_list,
      filename = NULL,  # Don't save to file, return the plot
      fill = c("lightblue", "pink", "lightgreen"),
      alpha = 0.5,
      lwd = 1,
      cex = 1,
      fontfamily = "sans",
      cat.fontfamily = "sans",
      cat.cex = 1,
      margins = c(1, 1)
    )
    
    # Try to open a PDF device to save the Venn diagram
    pdf("results/cancer_atlas/gene_overlap_venn.pdf", width = 8, height = 7)
    grid.draw(venn_plot)
    dev.off()
    
    # Display the Venn diagram in the current device
    grid.draw(venn_plot)
    
  }, error = function(e) {
    cat("Error creating Venn diagram:", e$message, "\n")
    cat("Creating alternative visualization instead.\n")
    
    # Create alternative visualization (3-way table)
    overlap_df <- data.frame(
      Category = c(
        "Peptide Genes Only", 
        "Cancer-Associated Only", 
        "Cancer-Testis Only",
        "Peptide & Cancer-Associated",
        "Peptide & Cancer-Testis",
        "Cancer-Associated & Cancer-Testis",
        "All Three Categories"
      ),
      Count = c(
        length(setdiff(all_genes_clean, union(ca_genes_clean, ct_genes_clean))),
        length(setdiff(ca_genes_clean, union(all_genes_clean, ct_genes_clean))),
        length(setdiff(ct_genes_clean, union(all_genes_clean, ca_genes_clean))),
        length(intersect(all_genes_clean, ca_genes_clean)) - length(intersect(intersect(all_genes_clean, ca_genes_clean), ct_genes_clean)),
        length(intersect(all_genes_clean, ct_genes_clean)) - length(intersect(intersect(all_genes_clean, ca_genes_clean), ct_genes_clean)),
        length(intersect(ca_genes_clean, ct_genes_clean)) - length(intersect(intersect(all_genes_clean, ca_genes_clean), ct_genes_clean)),
        length(intersect(intersect(all_genes_clean, ca_genes_clean), ct_genes_clean))
      )
    )
    
    write.csv(overlap_df, "results/cancer_atlas/gene_overlap_table.csv", row.names = FALSE)
    print(overlap_df)
  })
} else {
  cat("Not enough data for Venn diagram (empty gene sets or NA values).\n")
  cat("Peptide genes available:", length(all_genes_clean), "\n")
  cat("Cancer-associated genes available:", length(ca_genes_clean), "\n")
  cat("Cancer-testis genes available:", length(ct_genes_clean), "\n")
}

# 11. Save results
write.csv(ca_matches, "results/cancer_atlas/cancer_associated_matches.csv", row.names = FALSE)
write.csv(ptm_matches, "results/cancer_atlas/ptm_antigen_matches.csv", row.names = FALSE)
write.csv(ct_matches, "results/cancer_atlas/cancer_testis_matches.csv", row.names = FALSE)

# 12. Create a summary table of findings
summary_table <- data.frame(
  Category = c("Total Peptides", 
               "Total Genes", 
               "Cancer-Associated Antigen Matches (Peptides)",
               "PTM Antigen Matches (Peptides)", 
               "Cancer-Testis Antigen Matches (Peptides)"),
  Count = c(nrow(peptide_data),
            length(all_genes),
            n_distinct(ca_matches$Peptide),
            n_distinct(ptm_matches$Peptide),
            n_distinct(ct_matches$Peptide)),
  Percentage = c(100,
                 100,
                 n_distinct(ca_matches$Peptide) / nrow(peptide_data) * 100,
                 n_distinct(ptm_matches$Peptide) / nrow(peptide_data) * 100,
                 n_distinct(ct_matches$Peptide) / nrow(peptide_data) * 100)
)

# Round percentages to 2 decimal places
summary_table$Percentage <- round(summary_table$Percentage, 2)

# Add % symbol
summary_table$Percentage <- paste0(summary_table$Percentage, "%")

write.csv(summary_table, "results/cancer_atlas/summary_statistics.csv", row.names = FALSE)

# 13. Additional analysis - look at unique modifications in PTM matches
if (nrow(ptm_matches) > 0) {
  ptm_counts <- ptm_matches %>%
    group_by(Modification) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  
  # Plot PTM types
  ptm_plot <- ggplot(ptm_counts, aes(x = reorder(Modification, count), y = count)) +
    geom_bar(stat = "identity", fill = "darkgreen") +
    coord_flip() +
    theme_minimal() +
    labs(title = "PTM Types in Matched Peptides",
         x = "Modification Type", 
         y = "Count")
  
  print(ptm_plot)
  ggsave("results/cancer_atlas/ptm_types_barplot.pdf", ptm_plot, width = 8, height = 5)
}

# At the end of the script
# Update files_list to include the new files
files_list <- data.frame(
  filename = c(
    "cancer_associated_matches.csv",
    "cancer_associated_peptides.csv",  # NEW
    "cancer_associated_genes.csv",     # NEW
    "ptm_antigen_matches.csv",
    "ptm_peptides.csv",                # NEW
    "cancer_testis_matches.csv",
    "cancer_testis_peptides.csv",      # NEW
    "summary_statistics.csv",
    "cancer_associated_barplot.pdf",
    "gene_overlap_table.csv"
  ),
  description = c(
    "Contains all peptides that match with cancer-associated antigens. Each row represents a match between a peptide from our dataset and a gene listed in the Cancer Antigen Atlas. Columns include the peptide sequence, gene name, cancer type, and other metadata from the Cancer Atlas.",
    
    "A simplified report listing peptides from our dataset that are derived from cancer-associated genes. This file focuses on the peptide sequences themselves and the cancer types they're associated with.",
    
    "A summary by gene showing how many peptides in our dataset are derived from each cancer-associated gene, and listing the actual peptide sequences.",
    
    "Lists all peptides that match with post-translationally modified (PTM) antigens from the PTM Antigen List. Each row shows a peptide that matches exactly with a PTM antigen sequence. Columns include information about the modification type, site, and cancer types associated with the PTM.",
    
    "A simplified report of peptides with post-translational modifications, focusing on the modification type, site, and amino acid affected.",
    
    "Contains all peptides that match with cancer-testis antigens. Each row represents a match between a peptide from our dataset and a gene listed as a cancer-testis antigen in the Cancer Atlas.",
    
    "A simplified report listing peptides from our dataset that are derived from cancer-testis antigen genes, which are primarily expressed in testis tissue and various cancers.",
    
    "Provides an overall summary of the analysis, including total counts and percentages of peptides matching different antigen types. This gives a high-level view of how many peptides in our dataset correspond to known cancer antigens.",
    
    "Bar chart showing the number of peptide matches per cancer type. Cancer types are ordered by the number of matching peptides, with the most matches at the top. This visualizes which cancer types have the most representation in our peptide dataset.",
    
    "Table showing the overlap between gene sets (our peptide genes, cancer-associated genes, and cancer-testis genes). This helps quantify how many genes are shared between these categories."
  ),
  stringsAsFactors = FALSE
)

# Update the documentation text to clarify the matching process
create_documentation(
  "results/cancer_atlas/README.md",
  "Cancer Antigen Atlas Integration Results",
  paste(
    "This analysis cross-references the HLA-I peptide dataset with the Cancer Antigen Atlas to identify peptides derived from known cancer antigens.",
    "The integration examines three key areas:",
    "1. Cancer-associated antigens: Peptides derived from genes associated with specific cancer types (GENE-BASED MATCHING)",
    "2. Post-translational modifications (PTMs): Peptides with modifications that may be immunogenic (SEQUENCE-BASED MATCHING)",
    "3. Cancer-testis antigens: Peptides derived from genes with expression restricted to testis and various cancer types (GENE-BASED MATCHING)",
    "",
    "Matching Methods:",
    "- For cancer-associated and cancer-testis antigens: We matched the SOURCE GENES of peptides in our dataset with genes in the cancer antigen databases. This identifies peptides derived from genes known to be associated with cancer.",
    "- For PTM antigens: We matched PEPTIDE SEQUENCES directly, identifying peptides in our dataset that are identical to known modified peptides in the Cancer Atlas.",
    sep = "\n"
  ),
  files_list
)

# Additional detailed explanation of the analysis
cat("\n## Analysis Details\n\n",
    "### Gene Extraction Process\n",
    "The analysis began by extracting gene symbols from the 'genes' column of the peptide dataset. ",
    "Multiple genes associated with a single peptide (separated by semicolons) were processed as separate entries. ",
    "In total, ", length(all_genes_clean), " unique genes were identified from ", nrow(peptide_data), " peptides, ",
    "with some peptides being associated with multiple genes.\n\n",
    
    "### Cancer Association Matching\n",
    "We matched peptide genes with three Cancer Antigen Atlas databases:\n",
    "1. Cancer-Associated Antigens: ", length(ca_genes), " genes associated with specific cancer types\n",
    "2. PTM Antigens: Peptides with post-translational modifications\n",
    "3. Cancer-Testis Antigens: ", length(ct_genes), " genes with expression primarily in testis and various cancers\n\n",
    
    "The matching process identified:\n",
    "- ", n_distinct(ca_matches$Peptide), " peptides (", 
    round(n_distinct(ca_matches$Peptide)/nrow(peptide_data)*100, 2), 
    "%) matching with cancer-associated antigens\n",
    "- ", n_distinct(ptm_matches$Peptide), " peptides (", 
    round(n_distinct(ptm_matches$Peptide)/nrow(peptide_data)*100, 2), 
    "%) matching with PTM antigens\n",
    "- ", n_distinct(ct_matches$Peptide), " peptides (", 
    round(n_distinct(ct_matches$Peptide)/nrow(peptide_data)*100, 2), 
    "%) matching with cancer-testis antigens\n\n",
    
    "### Interpretation of Results\n",
    "The low overlap between our peptide dataset and the Cancer Antigen Atlas suggests several possibilities:\n",
    "1. The sample may not be derived from cancer tissue or may represent normal tissue\n",
    "2. The HLA-I peptidome captured in this dataset may contain antigens not yet characterized in cancer databases\n",
    "3. The specific HLA alleles in this sample may present different peptides than those commonly found in cancer studies\n\n",
    
    "Additional experiments or integration with other datasets would be valuable to further characterize these peptides.",
    file = "results/cancer_atlas/README.md", 
    append = TRUE)