# GO Analysis PDF Generator
# This script generates a comprehensive PDF report from your Gene Ontology analysis

# Install required packages if not already installed
if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown")
if (!requireNamespace("knitr", quietly = TRUE)) install.packages("knitr")
if (!requireNamespace("kableExtra", quietly = TRUE)) install.packages("kableExtra")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tinytex", quietly = TRUE)) install.packages("tinytex")

# Load necessary libraries
library(rmarkdown)
library(knitr)
library(kableExtra)
library(ggplot2)
library(dplyr)

# Set working directory to where your files are located
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/results/topGO")

# Create R Markdown document
rmd_content <- '---
title: "Gene Ontology Analysis of Immunopeptidomics Data"
author: "Your Name"
date: "`r Sys.Date()`"
output:
  pdf_document:
    toc: true
    toc_depth: 3
    fig_caption: true
    fig_height: 6
    fig_width: 8
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(ggplot2)
library(dplyr)
library(knitr)
library(kableExtra)
```

# Executive Summary

This report presents a comprehensive Gene Ontology (GO) analysis of immunopeptidomics data. The analysis focused on peptides that are shared across greater than half of samples, with the goal of identifying biological processes, molecular functions, and cellular components that are enriched in the dataset.

Key findings include:

1. Strong enrichment of translation machinery proteins, including cytoplasmic translation and ribosomal components
2. Significant representation of DNA repair proteins
3. Enrichment of proteins associated with extracellular exosomes and cytosolic components
4. RNA binding emerged as the most significant molecular function

These findings provide insights into which cellular processes contribute most significantly to the immunopeptidome in the samples analyzed, with potential implications for understanding immune system visibility and identifying potential neoantigens.

# Data Overview and Processing

The analysis started with a dataset containing 534 peptides from immunopeptidomics experiments. These peptides were shared across greater than half of the samples. Each peptide had an associated gene or genes, with some peptides mapping to multiple genes.

## Data Processing Steps

1. Extracted 616 total genes after splitting entries with multiple genes
2. Identified 492 unique genes for GO analysis
3. Found that 70 peptide entries were associated with multiple genes (separated by semicolons)

These genes were used as the input for Gene Ontology (GO) analysis, compared against the entire human genome as background.

```{r data-summary-table}
data_summary <- data.frame(
  Metric = c("Total peptides", "Total genes after splitting", "Unique genes", "Peptides with multiple genes"),
  Count = c(534, 616, 492, 70)
)

kable(data_summary, caption = "Summary of Dataset Characteristics") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)
```

# Gene Ontology Enrichment Results

Gene Ontology (GO) enrichment analysis was performed using the topGO R package. Three main ontologies were analyzed:

1. Biological Process (BP): Activities at the molecular or cellular level
2. Molecular Function (MF): Activities of gene products at the molecular level
3. Cellular Component (CC): Locations where gene products are active

## Biological Processes (BP)

The most significantly enriched biological processes are related to protein synthesis, nucleotide metabolism, and DNA repair.

```{r bp-table, echo=FALSE}
# Load the BP results file
bp_results <- read.csv("GO_BP_results.csv")

# Display top 10 rows with selected columns
bp_top10 <- bp_results[1:10, c("Term", "Significant", "Expected", "weightFisher", "log10_weightFisher")]
colnames(bp_top10) <- c("GO Term", "Genes Found", "Genes Expected", "P-value", "-log10(p)")

kable(bp_top10, caption = "Top 10 Enriched Biological Processes") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = TRUE)
```

```{r bp-plot, fig.cap="Bar plot showing the top 15 enriched Biological Processes.", fig.height=8}
# Recreate the BP barplot if needed, otherwise use existing file
if (file.exists("GO_BP_barplot.png")) {
  knitr::include_graphics("GO_BP_barplot.png")
} else {
  # Create plot from data
  p <- ggplot(bp_results[1:15,], aes(x = reorder(Term, log10_weightFisher), y = log10_weightFisher)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Top 15 Enriched Biological Processes",
         x = "GO Term",
         y = "-log10(p-value)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  print(p)
}
```

The most significant biological process is **cytoplasmic translation** (p=6.4e-11) with 26 genes found compared to 3.55 expected by chance. This strong enrichment of translation machinery suggests that proteins involved in protein synthesis are well-represented in the immunopeptidome, which aligns with their high abundance and turnover in cells.

Other significantly enriched processes include:

- **5-phosphoribose 1-diphosphate biosynthesis**: Critical for nucleotide synthesis
- **rRNA processing**: Essential for ribosome biogenesis
- **DNA repair**: Suggesting that nuclear proteins involved in genome maintenance are well-represented

## Molecular Functions (MF)

The molecular function analysis reveals which biochemical activities are overrepresented in the dataset.

```{r mf-table}
# Load the MF results file
mf_results <- read.csv("GO_MF_results.csv")

# Display top 10 rows with selected columns
mf_top10 <- mf_results[1:10, c("Term", "Significant", "Expected", "weightFisher", "log10_weightFisher")]
colnames(mf_top10) <- c("GO Term", "Genes Found", "Genes Expected", "P-value", "-log10(p)")

kable(mf_top10, caption = "Top 10 Enriched Molecular Functions") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = TRUE)
```

```{r mf-plot, fig.cap="Bar plot showing the top 15 enriched Molecular Functions.", fig.height=8}
# Recreate the MF barplot if needed, otherwise use existing file
if (file.exists("GO_MF_barplot.png")) {
  knitr::include_graphics("GO_MF_barplot.png")
} else {
  # Create plot from data
  p_mf <- ggplot(mf_results[1:15,], aes(x = reorder(Term, log10_weightFisher), y = log10_weightFisher)) +
    geom_bar(stat = "identity", fill = "forestgreen") +
    coord_flip() +
    labs(title = "Top 15 Enriched Molecular Functions",
         x = "GO Term",
         y = "-log10(p-value)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  print(p_mf)
}
```

**RNA binding** emerged as the most significantly enriched molecular function (p=1.7e-14), with 114 genes found compared to 45.67 expected by chance. This correlates well with the enriched biological processes related to translation and RNA processing.

Other significant functions include:

- **Cadherin binding**: Related to cell-cell adhesion
- **Structural constituent of cytoskeleton**: Involved in maintaining cellular structure
- **Structural constituent of ribosome**: Consistent with enrichment of translation-related processes

## Cellular Components (CC)

The cellular component analysis indicates which cellular locations are overrepresented in the dataset.

```{r cc-table}
# Load the CC results file
cc_results <- read.csv("GO_CC_results.csv")

# Display top 10 rows with selected columns
cc_top10 <- cc_results[1:10, c("Term", "Significant", "Expected", "weightFisher")]
colnames(cc_top10) <- c("GO Term", "Genes Found", "Genes Expected", "P-value")

kable(cc_top10, caption = "Top 10 Enriched Cellular Components") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = TRUE)
```

```{r cc-plot, fig.cap="Bar plot showing the top 15 enriched Cellular Components.", fig.height=8}
# Recreate the CC barplot if needed, otherwise use existing file
if (file.exists("GO_CC_barplot.png")) {
  knitr::include_graphics("GO_CC_barplot.png")
} else {
  # Create plot from data
  p_cc <- ggplot(cc_results[1:15,], aes(x = reorder(Term, log10_weightFisher), y = log10_weightFisher)) +
    geom_bar(stat = "identity", fill = "darkred") +
    coord_flip() +
    labs(title = "Top 15 Enriched Cellular Components",
         x = "GO Term",
         y = "-log10(p-value)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  print(p_cc)
}
```

**Extracellular exosome** shows extremely significant enrichment (p < 1e-30), with 136 genes found compared to 46.92 expected by chance. This is particularly interesting for immunopeptidomics, as exosomes can carry MHC-I molecules and associated peptides, potentially influencing immune responses.

**Cytosol** is also highly enriched (p=2.4e-25), which is expected given that the MHC-I antigen processing machinery primarily samples from the cytosolic protein pool.

# Integrated Analysis

## Combined Visualization of All Ontologies

```{r dot-plot, fig.cap="Dot plot showing the top 10 terms from each ontology category. Dot size represents the number of genes, while color indicates statistical significance.", fig.height=10}
# Recreate the dot plot if needed, otherwise use existing file
if (file.exists("GO_dotplot.png")) {
  knitr::include_graphics("GO_dotplot.png")
} else {
  # Create combined dotplot
  # First prepare the data
  prepare_dotplot_data <- function(go_table, ontology) {
    go_table$Ontology <- ontology
    go_table$GeneRatio <- go_table$Significant / go_table$Annotated
    return(go_table)
  }
  
  # Combine data from all three ontologies
  dotplot_data <- rbind(
    prepare_dotplot_data(bp_results[1:10,], "Biological Process"),
    prepare_dotplot_data(mf_results[1:10,], "Molecular Function"),
    prepare_dotplot_data(cc_results[1:10,], "Cellular Component")
  )
  
  # Create the dot plot
  p_dot <- ggplot(dotplot_data, aes(x = GeneRatio, y = reorder(Term, GeneRatio))) +
    geom_point(aes(size = Significant, color = log10_weightFisher)) +
    facet_grid(Ontology ~ ., scales = "free_y", space = "free") +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = "GO Enrichment Analysis",
         x = "Gene Ratio (Significant/Annotated)",
         y = "",
         size = "Gene Count",
         color = "-log10(p-value)") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8))
  print(p_dot)
}
```

The dot plot provides a comprehensive view of all three ontologies simultaneously, allowing for comparison across different categories. Key observations include:

1. Extracellular exosome and cytosol (CC terms) show the highest significance
2. RNA binding (MF) shows high significance and affects many genes
3. Cytoplasmic translation (BP) is the most significant biological process

## Multi-functional Genes

```{r multi-func-genes}
# Load the multi-functional genes file if it exists
if (file.exists("multi_functional_genes.csv")) {
  multi_genes <- read.csv("multi_functional_genes.csv")
  
  # Display the top genes
  kable(head(multi_genes, 10), caption = "Top Genes Appearing in Multiple GO Terms") %>%
    kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE)
} else {
  # Create placeholder
  cat("Multi-functional genes analysis data not available. Generate multi_functional_genes.csv file by running the GO analysis script.")
}
```

Genes appearing in multiple GO terms may represent particularly interesting targets for further study, as they contribute to multiple biological functions or processes.

# Immunological Implications

The GO enrichment results provide several insights relevant to immunopeptidomics:

## Source of MHC-I Peptides

The data suggests that the MHC-I peptides in the samples come predominantly from:

- Translation machinery (ribosomes, translation factors)
- Cytoskeletal components
- DNA repair proteins
- Membrane and exosome components

## Visibility to Immune System

The enrichment pattern indicates which cellular processes are most "visible" to the immune system in the samples analyzed. This visibility could be influenced by:

- Protein abundance in the cells
- Protein turnover rates
- Efficiency of processing by the proteasome
- Affinity of resulting peptides for MHC-I molecules

## Potential for Neoantigen Discovery

The enriched genes and their associated peptides provide a foundation for identifying potential neoantigens. Particularly interesting candidates would be peptides that:

1. Are derived from genes with cancer-related functions (e.g., DNA repair)
2. Are shared across multiple patient samples
3. Have sequences distinct from self-peptides

# Conclusions and Next Steps

This comprehensive GO analysis provides valuable insights into the biological context of the immunopeptidomics data. Key conclusions include:

1. The immunopeptidome is enriched for proteins involved in fundamental cellular processes, particularly protein synthesis
2. Extracellular components and membrane proteins are well-represented
3. Multiple cellular compartments contribute to the MHC-I peptide repertoire

## Recommended Next Steps

1. Investigate the most enriched genes for potential neoantigen candidates
2. Analyze patient-specific patterns in peptide presentation
3. Correlate GO enrichment with clinical or experimental variables
4. Explore the specific peptides derived from genes of interest for their immunogenic potential

# Appendix: Methodology

## Gene Ontology Analysis

Gene Ontology analysis was performed using the topGO R package with the following parameters:

- Algorithm: weight01 (accounts for GO hierarchy)
- Statistical test: Fisher\'s exact test
- Background: All human genes (from org.Hs.eg.db)
- Ontologies analyzed: Biological Process (BP), Molecular Function (MF), Cellular Component (CC)

The weight01 algorithm was chosen because it accounts for the hierarchical structure of GO terms, reducing redundancy in the results.

```{r session-info, echo=TRUE}
sessionInfo()
```'

# Write the R Markdown content to a file
writeLines(rmd_content, "GO_Analysis_Report.Rmd")

# Check if TinyTex is installed (needed for PDF rendering)
if (!tinytex::is_tinytex()) {
  message("TinyTeX is not installed. Would you like to install it now? (y/n)")
  choice <- readline()
  if (tolower(choice) == "y") {
    tinytex::install_tinytex()
  } else {
    message("PDF rendering may fail without LaTeX. Please install TinyTeX or another LaTeX distribution.")
  }
}

# Render the document to PDF
message("Rendering PDF... This may take a few moments.")
render("GO_Analysis_Report.Rmd", output_format = "pdf_document")

message("Done! PDF report has been generated: GO_Analysis_Report.pdf")

# Create a script for neoantigen analysis
neoantigen_script <- '# Neoantigen Analysis Script
# This script analyzes peptides, genes, and patient distribution for potential neoantigens

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(readxl)

# Load the peptide data
peptide_data <- read_excel("shared50_filtered_out_148N_manual.xlsx", sheet = "Sheet3")

# 1. Create a summary table: Gene → Peptides → Patient distribution
peptide_gene_summary <- peptide_data %>%
  mutate(
    # Split genes column into multiple rows if needed
    gene_list = strsplit(as.character(genes), "; |;"),
    # Count number of patients/samples
    patient_count = sample_count,
    # Extract patient IDs
    patient_ids = sample_list
  ) %>%
  unnest(gene_list) %>%
  group_by(gene_list) %>%
  summarize(
    total_peptides = n(),
    unique_peptides = n_distinct(Peptide),
    max_patients = max(patient_count),
    min_patients = min(patient_count),
    avg_patients = mean(patient_count),
    peptide_list = paste(unique(Peptide), collapse = "; "),
    all_patients = paste(unique(unlist(strsplit(patient_ids, ", |,"))), collapse = ", ")
  ) %>%
  arrange(desc(max_patients), desc(unique_peptides))

# Save the summary
write.csv(peptide_gene_summary, "gene_peptide_patient_summary.csv", row.names = FALSE)

# 2. Create a detailed peptide-level analysis
peptide_detail <- peptide_data %>%
  mutate(
    gene_list = strsplit(as.character(genes), "; |;")
  ) %>%
  unnest(gene_list) %>%
  select(Peptide, `Peptide Length`, gene_list, sample_list, sample_count) %>%
  rename(
    Gene = gene_list,
    Patient_IDs = sample_list,
    Patient_Count = sample_count
  )

# Save detailed peptide data
write.csv(peptide_detail, "peptide_gene_patient_detail.csv", row.names = FALSE)

# 3. Create a patient-gene matrix (for heatmap visualization)
# First, create a table of patient occurrences per gene
patient_gene_matrix <- peptide_detail %>%
  # Split patient IDs into separate rows
  mutate(patients = strsplit(as.character(Patient_IDs), ", |,")) %>%
  unnest(patients) %>%
  # Count unique peptides for each gene-patient combination
  group_by(Gene, patients) %>%
  summarize(peptide_count = n_distinct(Peptide)) %>%
  # Reshape to wide format for heatmap
  pivot_wider(
    names_from = patients,
    values_from = peptide_count,
    values_fill = 0
  )

# Save the matrix
write.csv(patient_gene_matrix, "patient_gene_matrix.csv", row.names = FALSE)

# 4. Create visualizations

# 4.1 Top genes by patient coverage
top_genes_plot <- ggplot(
  peptide_gene_summary %>% 
    arrange(desc(max_patients), desc(unique_peptides)) %>% 
    head(20),
  aes(x = reorder(gene_list, max_patients), y = max_patients)
) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = unique_peptides), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(
    title = "Top 20 Genes by Patient Coverage",
    subtitle = "Numbers indicate unique peptides per gene",
    x = "Gene",
    y = "Maximum Patient Count"
  ) +
  theme_minimal()

ggsave("top_genes_by_patient.png", top_genes_plot, width = 10, height = 8)

# 4.2 Gene-Patient Heatmap
# Select top genes for better visualization
top_genes <- peptide_gene_summary %>% 
  arrange(desc(max_patients), desc(unique_peptides)) %>% 
  head(30) %>% 
  pull(gene_list)

# Filter matrix for top genes
heatmap_data <- patient_gene_matrix %>%
  filter(Gene %in% top_genes)

# Convert to matrix format for pheatmap
row_names <- heatmap_data$Gene
heatmap_data <- as.matrix(heatmap_data[,-1])
rownames(heatmap_data) <- row_names

# Save heatmap to PNG
png("gene_patient_heatmap.png", width = 12, height = 10, units = "in", res = 300)
pheatmap(
  heatmap_data,
  main = "Gene-Patient Peptide Distribution",
  color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  display_numbers = FALSE
)
dev.off()

# 4.3 Peptide length distribution by top genes
peptide_length_plot <- ggplot(
  peptide_detail %>% 
    filter(Gene %in% top_genes) %>%
    group_by(Gene, `Peptide Length`) %>%
    summarize(count = n()),
  aes(x = `Peptide Length`, y = count, fill = Gene)
) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Peptide Length Distribution by Top Genes",
    x = "Peptide Length",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("peptide_length_by_gene.png", peptide_length_plot, width = 12, height = 8)

# 4.4 Patient overlap network visualization (requires igraph)
# This creates a network showing which patients share the same peptides
if (requireNamespace("igraph", quietly = TRUE)) {
  library(igraph)
  
  # Create patient-patient connections based on shared peptides
  patient_connections <- peptide_data %>%
    # Split patient IDs into separate rows
    mutate(patients = strsplit(as.character(sample_list), ", |,")) %>%
    select(Peptide, patients) %>%
    unnest(patients) %>%
    group_by(Peptide) %>%
    do({
      pts <- unique(.$patients)
      if (length(pts) > 1) {
        expand.grid(from = pts, to = pts, stringsAsFactors = FALSE) %>%
          filter(from < to)  # Avoid self-loops and duplicates
      } else {
        data.frame(from = character(0), to = character(0))
      }
    }) %>%
    ungroup() %>%
    count(from, to) %>%
    rename(weight = n)
  
  # Create graph
  g <- graph_from_data_frame(
    patient_connections,
    directed = FALSE
  )
  
  # Save network visualization
  png("patient_overlap_network.png", width = 10, height = 10, units = "in", res = 300)
  plot(
    g,
    vertex.size = 10,
    vertex.label.cex = 0.8,
    vertex.color = "lightblue",
    edge.width = E(g)$weight / max(E(g)$weight) * 5,
    main = "Patient Overlap Network (Based on Shared Peptides)"
  )
  dev.off()
}

# Print summary information
cat("\nAnalysis complete. Files generated:\n")
cat("1. gene_peptide_patient_summary.csv - Overview of genes, peptides, and patient distribution\n")
cat("2. peptide_gene_patient_detail.csv - Detailed peptide-level information\n")
cat("3. patient_gene_matrix.csv - Matrix showing peptide counts per gene per patient\n")
cat("4. top_genes_by_patient.png - Bar chart of top genes by patient coverage\n")
cat("5. gene_patient_heatmap.png - Heatmap showing peptide distribution across genes and patients\n")
cat("6. peptide_length_by_gene.png - Peptide length distribution for top genes\n")
cat("7. patient_overlap_network.png - Network visualization of patient overlap (if igraph is installed)\n")
'

# Save the neoantigen analysis script
writeLines(neoantigen_script, "Neoantigen_Analysis_Script.R")

message("A script for neoantigen analysis has also been created: Neoantigen_Analysis_Script.R")