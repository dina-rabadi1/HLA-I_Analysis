# # What I used to install
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   + install.packages("BiocManager")
# BiocManager::install()
# 
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   + install.packages("BiocManager")
# BiocManager::install("topGO")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GO.db")
# # 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")
# 
# This script uses topGO to do gene ontology analysis 


library(topGO)
library(GOplot)
library(ggplot2)
library(dplyr)
library(readxl)
library(org.Hs.eg.db)
library(knitr)
library(grid)
library(gridExtra)
library(gtable)

setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis")

# Create output directory structure
if (!dir.exists("results/topGO")) {
  dir.create("results/topGO")
}

# 1. Read and preprocess the Excel file
peptide_data <- read_excel("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/shared50_filtered_out_148N_manual.xlsx", sheet = "Sheet3")

# Check the total number of rows in your dataset
print(paste("Total number of rows in dataset:", nrow(peptide_data)))

# 2. Better extraction of gene list
# Print a few examples of the gene column to see the format
print("Sample gene entries:")
print(head(peptide_data$genes, 20))

# Count genes with ";" separator to see how many entries have multiple genes
multiple_genes_count <- sum(grepl(";", peptide_data$genes))
print(paste("Number of entries with multiple genes:", multiple_genes_count))

# Split genes when there are multiple genes separated by semicolons
all_genes <- unlist(strsplit(as.character(peptide_data$genes), "; |;"))
# Remove any potential whitespace
all_genes <- trimws(all_genes)
# Get unique genes
unique_genes <- unique(all_genes)
print(paste("Number of all genes after splitting:", length(all_genes)))
print(paste("Number of unique genes after splitting:", length(unique_genes)))

# 3. Prepare gene list for topGO - CORRECTED APPROACH
# For topGO, we need to create a named vector with gene names and their "importance"
# First, get all possible gene symbols from the org.Hs.eg.db database
gene_universe <- rep(0, length(unique_genes))
names(gene_universe) <- unique_genes

# Mark all genes in our list as "of interest" (value = 1)
gene_universe[] <- 1

# 4. Now get a list of all genes in the database to use as background
all_human_genes <- keys(org.Hs.eg.db, keytype="SYMBOL")
# Create a full universe vector with all human genes
full_gene_universe <- rep(0, length(all_human_genes))
names(full_gene_universe) <- all_human_genes

# Mark genes that are in our list as "of interest"
full_gene_universe[names(gene_universe)] <- 1

# Convert to factor for topGO
gene_vector <- factor(full_gene_universe)

# Check if we now have a factor with 2 levels
print(paste("Number of levels in gene vector:", nlevels(gene_vector)))
print(table(gene_vector))

# 5. Create topGO data object for Biological Process (BP) ontology
GO_data_BP <- new("topGOdata",
                  ontology = "BP",
                  allGenes = gene_vector,
                  geneSel = function(x) x == 1,
                  annot = annFUN.org,
                  mapping = "org.Hs.eg.db",  # Use human database
                  ID = "symbol")  # Using gene symbols

# Print summary of GO data to check if it was created successfully
print(GO_data_BP)

# 6. Run GO enrichment analysis using different algorithms
# Fisher's exact test
result_fisher_BP <- runTest(GO_data_BP, algorithm = "classic", statistic = "fisher")
# Weight algorithm (accounts for GO hierarchy)
result_weight_BP <- runTest(GO_data_BP, algorithm = "weight01", statistic = "fisher")

# 7. Get results table
results_table_BP <- GenTable(GO_data_BP, 
                             classicFisher = result_fisher_BP,
                             weightFisher = result_weight_BP,
                             orderBy = "weightFisher", 
                             ranksOf = "classicFisher",
                             topNodes = 30)  # Top 30 GO terms

# 8. Add more columns to the results table for better interpretation
results_table_BP$log10_weightFisher <- -log10(as.numeric(results_table_BP$weightFisher))
results_table_BP$log10_classicFisher <- -log10(as.numeric(results_table_BP$classicFisher))

# Print the top results
print(results_table_BP[1:10, ])

# 9. Save results to file
write.csv(results_table_BP, "results/topGO/GO_BP_results.csv", row.names = FALSE)

# 10. Create visualization
# Bar plot for Biological Process
p <- ggplot(results_table_BP[1:15,], aes(x = reorder(Term, log10_weightFisher), y = log10_weightFisher)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 15 Enriched Biological Processes",
       x = "GO Term",
       y = "-log10(p-value)") +
  theme_minimal()
print(p)
ggsave("results/topGO/GO_BP_barplot.png", width = 10, height = 8)

# 11. Get the genes associated with top GO terms (for further analysis)
# Example for top 5 BP terms
top_go_terms <- results_table_BP$GO.ID[1:5]
genes_in_top_terms <- list()

for (term in top_go_terms) {
  genes <- genesInTerm(GO_data_BP, term)
  significant_genes <- intersect(names(which(gene_vector == 1)), genes[[1]])
  genes_in_top_terms[[term]] <- significant_genes
}

# Print genes for each top term
for (i in 1:length(top_go_terms)) {
  cat("\nGenes in", results_table_BP$Term[i], "(", top_go_terms[i], "):\n")
  print(genes_in_top_terms[[top_go_terms[i]]])
}

# Let's continue with the analysis and fix the gene retrieval issue

# First, let's try a different approach to get genes in each GO term
# We'll use the sigGenes function from topGO
sig_genes_BP <- sigGenes(GO_data_BP)
print(paste("Number of significant genes in topGO object:", length(sig_genes_BP)))

# Let's examine the top 5 GO terms and their genes more carefully
for (i in 1:5) {
  term <- results_table_BP$GO.ID[i]
  term_genes <- genesInTerm(GO_data_BP, term)[[1]]
  # Filter to just keep the significant genes
  sig_term_genes <- intersect(term_genes, sig_genes_BP)
  
  cat("\nGenes in", results_table_BP$Term[i], "(", term, "):\n")
  print(sig_term_genes)
  cat("Number of genes:", length(sig_term_genes), "\n")
}

# Let's continue with MF and CC ontologies
# Create topGO data object for Molecular Function (MF) ontology
GO_data_MF <- new("topGOdata",
                  ontology = "MF",
                  allGenes = gene_vector,
                  geneSel = function(x) x == 1,
                  annot = annFUN.org,
                  mapping = "org.Hs.eg.db",
                  ID = "symbol")

# Run Fisher's exact test for MF
result_fisher_MF <- runTest(GO_data_MF, algorithm = "classic", statistic = "fisher")
result_weight_MF <- runTest(GO_data_MF, algorithm = "weight01", statistic = "fisher")

# Get results table for MF
results_table_MF <- GenTable(GO_data_MF, 
                             classicFisher = result_fisher_MF,
                             weightFisher = result_weight_MF,
                             orderBy = "weightFisher", 
                             ranksOf = "classicFisher",
                             topNodes = 30)

# Add log10 columns
results_table_MF$log10_weightFisher <- -log10(as.numeric(results_table_MF$weightFisher))
results_table_MF$log10_classicFisher <- -log10(as.numeric(results_table_MF$classicFisher))

# Create topGO data object for Cellular Component (CC) ontology
GO_data_CC <- new("topGOdata",
                  ontology = "CC",
                  allGenes = gene_vector,
                  geneSel = function(x) x == 1,
                  annot = annFUN.org,
                  mapping = "org.Hs.eg.db",
                  ID = "symbol")

# Run Fisher's exact test for CC
result_fisher_CC <- runTest(GO_data_CC, algorithm = "classic", statistic = "fisher")
result_weight_CC <- runTest(GO_data_CC, algorithm = "weight01", statistic = "fisher")

# Get results table for CC
results_table_CC <- GenTable(GO_data_CC, 
                             classicFisher = result_fisher_CC,
                             weightFisher = result_weight_CC,
                             orderBy = "weightFisher", 
                             ranksOf = "classicFisher",
                             topNodes = 30)

# Add log10 columns
results_table_CC$log10_weightFisher <- -log10(as.numeric(results_table_CC$weightFisher))
results_table_CC$log10_classicFisher <- -log10(as.numeric(results_table_CC$classicFisher))

# Save results to files
write.csv(results_table_MF, "results/topGO/GO_MF_results.csv", row.names = FALSE)
write.csv(results_table_CC, "results/topGO/GO_CC_results.csv", row.names = FALSE)

# Print top results for MF and CC
cat("\nTop 10 Molecular Function GO terms:\n")
print(results_table_MF[1:10, c("GO.ID", "Term", "Significant", "Expected", "weightFisher")])

cat("\nTop 10 Cellular Component GO terms:\n")
print(results_table_CC[1:10, c("GO.ID", "Term", "Significant", "Expected", "weightFisher")])

# Create visualizations for MF and CC
# Bar plot for Molecular Function
p_mf <- ggplot(results_table_MF[1:15,], aes(x = reorder(Term, log10_weightFisher), y = log10_weightFisher)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  coord_flip() +
  labs(title = "Top 15 Enriched Molecular Functions",
       x = "GO Term",
       y = "-log10(p-value)") +
  theme_minimal()
print(p_mf)
ggsave("results/topGO/GO_MF_barplot.png", width = 10, height = 8)

# Bar plot for Cellular Component
p_cc <- ggplot(results_table_CC[1:15,], aes(x = reorder(Term, log10_weightFisher), y = log10_weightFisher)) +
  geom_bar(stat = "identity", fill = "darkred") +
  coord_flip() +
  labs(title = "Top 15 Enriched Cellular Components",
       x = "GO Term",
       y = "-log10(p-value)") +
  theme_minimal()
print(p_cc)
ggsave("results/topGO/GO_CC_barplot.png", width = 10, height = 8)

# Create a combined dot plot showing all three ontologies
prepare_dotplot_data <- function(go_table, ontology) {
  go_table$Ontology <- ontology
  go_table$GeneRatio <- go_table$Significant / go_table$Annotated
  return(go_table)
}

# Select top 10 terms from each ontology
dotplot_data <- rbind(
  prepare_dotplot_data(results_table_BP[1:10,], "Biological Process"),
  prepare_dotplot_data(results_table_MF[1:10,], "Molecular Function"),
  prepare_dotplot_data(results_table_CC[1:10,], "Cellular Component")
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
ggsave("results/topGO/GO_dotplot.png", width = 12, height = 10)

# Define helper function for extracting genes in GO terms
genes_in_GO_terms <- function(go_data, term_index) {
  if (term_index <= 0 || term_index > length(go_data@graph@nodes)) {
    return(character(0))
  }
  
  term <- go_data@graph@nodes[term_index]
  if (is.null(term)) {
    return(character(0))
  }
  
  term_genes <- genesInTerm(go_data, term)[[1]]
  sig_term_genes <- intersect(term_genes, sigGenes(go_data))
  return(sig_term_genes)
}

# Find enriched GO terms common across different ontologies
# This can help identify biological themes from different perspectives
common_terms <- list()
for (bp_term in results_table_BP$Term[1:10]) {
  for (mf_term in results_table_MF$Term[1:10]) {
    common_genes_bp_mf <- intersect(
      genes_in_GO_terms(GO_data_BP, grep(bp_term, results_table_BP$Term, fixed=TRUE)[1]), 
      genes_in_GO_terms(GO_data_MF, grep(mf_term, results_table_MF$Term, fixed=TRUE)[1])
    )
    if (length(common_genes_bp_mf) > 3) {
      common_terms[[paste(bp_term, " & ", mf_term)]] <- common_genes_bp_mf
    }
  }
}

# Create a network visualization of relationships between top GO terms
# First, prepare data for GOplot
id2gene <- data.frame(
  ID = names(gene_vector)[gene_vector == 1],
  genes = names(gene_vector)[gene_vector == 1]
)

# Combine top terms from all three ontologies
top_terms <- c(
  results_table_BP$GO.ID[1:10],
  results_table_MF$GO.ID[1:10],
  results_table_CC$GO.ID[1:10]
)

# Create a data frame for GOplot's enrichment function
enrich_df <- data.frame(
  ID = c(results_table_BP$GO.ID[1:10], results_table_MF$GO.ID[1:10], results_table_CC$GO.ID[1:10]),
  term = c(results_table_BP$Term[1:10], results_table_MF$Term[1:10], results_table_CC$Term[1:10]),
  ontology = c(rep("BP", 10), rep("MF", 10), rep("CC", 10)),
  p.adj = c(results_table_BP$weightFisher[1:10], results_table_MF$weightFisher[1:10], results_table_CC$weightFisher[1:10])
)

# Add count and gene ratio
enrich_df$count <- c(results_table_BP$Significant[1:10], results_table_MF$Significant[1:10], results_table_CC$Significant[1:10])
enrich_df$gene.ratio <- c(
  results_table_BP$Significant[1:10] / results_table_BP$Annotated[1:10],
  results_table_MF$Significant[1:10] / results_table_MF$Annotated[1:10],
  results_table_CC$Significant[1:10] / results_table_CC$Annotated[1:10]
)

# Save this enrichment data for future use
write.csv(enrich_df, "results/topGO/GO_enrichment_combined.csv", row.names = FALSE)

# Create an interactive chord diagram using GOplot if possible
# (This may require additional steps depending on your environment)
# Try to create the chord diagram
tryCatch({
  # Get gene lists for each term
  genes_by_term <- list()
  for (i in 1:nrow(enrich_df)) {
    go_term <- enrich_df$ID[i]
    if (i <= 10) {
      go_data <- GO_data_BP
    } else if (i <= 20) {
      go_data <- GO_data_MF
    } else {
      go_data <- GO_data_CC
    }
    
    term_genes <- genesInTerm(go_data, go_term)[[1]]
    sig_genes <- intersect(term_genes, sigGenes(go_data))
    genes_by_term[[go_term]] <- sig_genes
  }
  
  # Create a matrix for the chord diagram
  chord_data <- data.frame(
    term = rep(enrich_df$ID, sapply(genes_by_term, length)),
    genes = unlist(genes_by_term)
  )
  
  # Save this for potential manual visualization later
  write.csv(chord_data, "results/topGO/GO_chord_data.csv", row.names = FALSE)
  
  # If GOplot's chord diagram function is available, create one
  if ("chord_data" %in% ls() && requireNamespace("GOplot", quietly = TRUE)) {
    GOplot::chord_data(chord_data, enrich_df)
    GOplot::GOChord(chord_data, enrich_df)
    ggsave("results/topGO/GO_chord_diagram.png", width = 12, height = 12)
  }
}, error = function(e) {
  cat("Chord diagram creation failed:", e$message, "\n")
  cat("Saving the chord data for manual visualization later.\n")
})

# Create a summary of the top enriched GO terms and their significance
summary_df <- data.frame(
  Rank = 1:10,
  BP_Term = results_table_BP$Term[1:10],
  BP_PValue = results_table_BP$weightFisher[1:10],
  MF_Term = results_table_MF$Term[1:10],
  MF_PValue = results_table_MF$weightFisher[1:10],
  CC_Term = results_table_CC$Term[1:10],
  CC_PValue = results_table_CC$weightFisher[1:10]
)

# Save this summary table
write.csv(summary_df, "results/topGO/GO_summary_table.csv", row.names = FALSE)

# Find enriched GO terms common across different ontologies - simplified approach
# We'll get genes for the top GO terms in each category and look for overlaps
cat("\nGetting genes for top GO terms in each category:\n")

# For BP
top_bp_genes <- list()
for (i in 1:5) {
  term <- results_table_BP$GO.ID[i]
  term_genes <- genesInTerm(GO_data_BP, term)[[1]]
  sig_term_genes <- intersect(term_genes, sig_genes_BP)
  top_bp_genes[[results_table_BP$Term[i]]] <- sig_term_genes
}

# For MF
sig_genes_MF <- sigGenes(GO_data_MF)
top_mf_genes <- list()
for (i in 1:5) {
  term <- results_table_MF$GO.ID[i]
  term_genes <- genesInTerm(GO_data_MF, term)[[1]]
  sig_term_genes <- intersect(term_genes, sig_genes_MF)
  top_mf_genes[[results_table_MF$Term[i]]] <- sig_term_genes
}

# For CC
sig_genes_CC <- sigGenes(GO_data_CC)
top_cc_genes <- list()
for (i in 1:5) {
  term <- results_table_CC$GO.ID[i]
  term_genes <- genesInTerm(GO_data_CC, term)[[1]]
  sig_term_genes <- intersect(term_genes, sig_genes_CC)
  top_cc_genes[[results_table_CC$Term[i]]] <- sig_term_genes
}

# Find common genes between top BP and MF terms
cat("\nLooking for common genes between top Biological Process and Molecular Function terms:\n")
common_genes_bp_mf <- list()
for (bp_name in names(top_bp_genes)) {
  for (mf_name in names(top_mf_genes)) {
    common <- intersect(top_bp_genes[[bp_name]], top_mf_genes[[mf_name]])
    if (length(common) >= 3) {
      key <- paste(bp_name, "AND", mf_name)
      common_genes_bp_mf[[key]] <- common
      cat("\n", key, ":", length(common), "genes\n")
      print(common)
    }
  }
}

# Find genes that appear in multiple top GO terms (across all categories)
# This can help identify key genes that contribute to multiple biological functions
all_top_genes <- c(
  unlist(top_bp_genes, use.names = FALSE),
  unlist(top_mf_genes, use.names = FALSE),
  unlist(top_cc_genes, use.names = FALSE)
)

gene_counts <- table(all_top_genes)
multi_term_genes <- names(gene_counts[gene_counts >= 3])

cat("\nGenes appearing in 3 or more top GO terms:\n")
print(multi_term_genes)

# For each multi-functional gene, find which terms it belongs to
gene_term_mapping <- list()
for (gene in multi_term_genes) {
  terms_for_gene <- character(0)
  
  # Check BP terms
  for (term_name in names(top_bp_genes)) {
    if (gene %in% top_bp_genes[[term_name]]) {
      terms_for_gene <- c(terms_for_gene, paste("BP:", term_name))
    }
  }
  
  # Check MF terms
  for (term_name in names(top_mf_genes)) {
    if (gene %in% top_mf_genes[[term_name]]) {
      terms_for_gene <- c(terms_for_gene, paste("MF:", term_name))
    }
  }
  
  # Check CC terms
  for (term_name in names(top_cc_genes)) {
    if (gene %in% top_cc_genes[[term_name]]) {
      terms_for_gene <- c(terms_for_gene, paste("CC:", term_name))
    }
  }
  
  gene_term_mapping[[gene]] <- terms_for_gene
}

cat("\nTerms associated with multi-functional genes:\n")
for (gene in names(gene_term_mapping)) {
  cat("\n", gene, "appears in these terms:\n")
  print(gene_term_mapping[[gene]])
}

# Save the multi-functional genes information
multi_gene_df <- data.frame(
  Gene = names(gene_term_mapping),
  NumTerms = sapply(gene_term_mapping, length)
)
multi_gene_df <- multi_gene_df[order(multi_gene_df$NumTerms, decreasing = TRUE), ]
write.csv(multi_gene_df, "results/topGO/multi_functional_genes.csv", row.names = FALSE)

# Create a summary table with the top terms from each category
summary_df <- data.frame(
  Rank = 1:10,
  BP_Term = results_table_BP$Term[1:10],
  BP_PValue = results_table_BP$weightFisher[1:10],
  BP_Genes = sapply(results_table_BP$GO.ID[1:10], function(term) {
    genes <- genesInTerm(GO_data_BP, term)[[1]]
    return(length(intersect(genes, sig_genes_BP)))
  }),
  MF_Term = results_table_MF$Term[1:10],
  MF_PValue = results_table_MF$weightFisher[1:10],
  MF_Genes = sapply(results_table_MF$GO.ID[1:10], function(term) {
    genes <- genesInTerm(GO_data_MF, term)[[1]]
    return(length(intersect(genes, sig_genes_MF)))
  }),
  CC_Term = results_table_CC$Term[1:10],
  CC_PValue = results_table_CC$weightFisher[1:10],
  CC_Genes = sapply(results_table_CC$GO.ID[1:10], function(term) {
    genes <- genesInTerm(GO_data_CC, term)[[1]]
    return(length(intersect(genes, sig_genes_CC)))
  })
)

# Save this summary table
write.csv(summary_df, "results/topGO/GO_summary_table.csv", row.names = FALSE)

# Finally, create a list of all significant genes with their GO annotations
all_sig_genes <- unique(c(sig_genes_BP, sig_genes_MF, sig_genes_CC))
gene_annotations <- list()

# Function to get all GO terms for a gene
get_gene_go_terms <- function(gene, go_data, results_table) {
  gene_terms <- character(0)
  
  for (i in 1:nrow(results_table)) {
    term <- results_table$GO.ID[i]
    term_genes <- genesInTerm(go_data, term)[[1]]
    
    if (gene %in% term_genes && gene %in% sigGenes(go_data)) {
      gene_terms <- c(gene_terms, paste(results_table$Term[i], " (", term, ")", sep = ""))
    }
  }
  
  return(gene_terms)
}

# Get annotations for the top 50 most significant genes
cat("\nGetting GO annotations for top genes...\n")
gene_annotations_df <- data.frame(
  Gene = character(0),
  BP_Terms = character(0),
  MF_Terms = character(0),
  CC_Terms = character(0)
)

for (gene in all_sig_genes[1:min(50, length(all_sig_genes))]) {
  bp_terms <- get_gene_go_terms(gene, GO_data_BP, results_table_BP[1:20,])
  mf_terms <- get_gene_go_terms(gene, GO_data_MF, results_table_MF[1:20,])
  cc_terms <- get_gene_go_terms(gene, GO_data_CC, results_table_CC[1:20,])
  
  if (length(bp_terms) > 0 || length(mf_terms) > 0 || length(cc_terms) > 0) {
    gene_annotations_df <- rbind(gene_annotations_df, data.frame(
      Gene = gene,
      BP_Terms = paste(bp_terms, collapse = "; "),
      MF_Terms = paste(mf_terms, collapse = "; "),
      CC_Terms = paste(cc_terms, collapse = "; ")
    ))
  }
}

# Save annotations
write.csv(gene_annotations_df, "results/topGO/gene_GO_annotations.csv", row.names = FALSE)

cat("\nAnalysis complete. Files saved:\n")
cat("1. GO_BP_results.csv - Biological Process GO terms\n")
cat("2. GO_MF_results.csv - Molecular Function GO terms\n")
cat("3. GO_CC_results.csv - Cellular Component GO terms\n")
cat("4. GO_summary_table.csv - Summary of top terms from all categories\n")
cat("5. multi_functional_genes.csv - Genes appearing in multiple GO terms\n")
cat("6. gene_GO_annotations.csv - GO annotations for top genes\n")
cat("7. GO_BP_barplot.png, GO_MF_barplot.png, GO_CC_barplot.png - Bar plots of top terms\n")
cat("8. GO_dotplot.png - Combined dot plot of all categories\n")

# ===== SAVE ALL VISUALIZATIONS TO PDF =====
pdf("topGO_analysis_report.pdf", width=11, height=8.5)

# Page 1: GO dot plot
title_grob <- textGrob("GO Dotplot", 
                       gp=gpar(fontsize=16, fontface="bold"))

grid.arrange(
  p_dot,
  ncol=1, nrow=1,
  top=title_grob
)

# Page 2: 
title_grob <- textGrob("GO Barplots", 
                       gp=gpar(fontsize=16, fontface="bold"))

grid.arrange(
  p, p_mf, p_cc,
  ncol=1, nrow=3,
  top=title_grob
)

# Page 3: Summary statistics (fully working)

# Create title
title_grob <- textGrob("GO Summary Statistics", 
                       gp = gpar(fontsize = 16, fontface = "bold"))

# Make table grob with formatting
summary_table <- tableGrob(summary_df, rows = NULL, theme = ttheme_minimal(
  core = list(fg_params = list(hjust = 0, x = 0.1, fontsize = 8)),
  colhead = list(fg_params = list(hjust = 0.5, fontface = "bold", fontsize = 9))
))

# Scale to fit page width
table_width <- sum(summary_table$widths)
if (convertWidth(table_width, "inches", valueOnly = TRUE) > 10) {
  summary_table$widths <- summary_table$widths * (10 / convertWidth(table_width, "inches", valueOnly = TRUE))
}

# Arrange title and table
grid.arrange(title_grob, summary_table, 
             ncol = 1,
             heights = unit.c(unit(1, "lines"), unit(1, "npc") - unit(1, "lines")))


dev.off()