# Load required libraries
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)  # For creating multi-sheet Excel files

# Install additional required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")
if (!requireNamespace("factoextra", quietly = TRUE)) install.packages("factoextra")
if (!requireNamespace("FactoMineR", quietly = TRUE)) install.packages("FactoMineR")

# Try to install Bioconductor packages if needed 
# Uncomment these lines if you want to install bioconductor packages
# if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
# if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")

# Load additional libraries
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(grid)
library(FactoMineR)
library(factoextra)

# Create the output directory if it doesn't exist
dir.create("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/RU148_transcriptome", showWarnings = FALSE)

# Set working directory to the raw data folder to read the input file
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata")

# Read the Excel file
gene_data <- read_excel("Normalized_Gene_counts_FLCdb_Panel_1.xlsx")

# Select columns of interest
selected_data <- gene_data %>%
  select(RU148_N, RU148_T8, RU148_T11, Mean.Normal, Mean.Tumor, 
         geneID, symbol, biotype, chromosome, gene_start, gene_end, gene_length, description)

# Calculate log2 fold change for each tumor sample compared to normal
selected_data <- selected_data %>%
  mutate(
    log2FC_T8 = log2((RU148_T8 + 0.1) / (RU148_N + 0.1)),
    log2FC_T11 = log2((RU148_T11 + 0.1) / (RU148_N + 0.1)),
    mean_log2FC = (log2((RU148_T8 + 0.1) / (RU148_N + 0.1)) + 
                     log2((RU148_T11 + 0.1) / (RU148_N + 0.1))) / 2
  )

# Option 1: Using individual sample criteria
upregulated_individual <- selected_data %>%
  filter(
    # Require log2FC > 1 in at least one tumor sample
    (log2FC_T8 > 1 | log2FC_T11 > 1),
    # Ensure expression level > 10 in the corresponding sample
    ((log2FC_T8 > 1 & RU148_T8 > 10) | (log2FC_T11 > 1 & RU148_T11 > 10))
  ) %>%
  arrange(desc(mean_log2FC))

# Option 2: Focus on genes up in both tumor samples compared to normal
upregulated_both_samples <- selected_data %>%
  filter(
    # Up in both samples (less stringent threshold)
    log2FC_T8 > 0.5 & log2FC_T11 > 0.5,
    # Reasonable expression in both
    RU148_T8 > 5 & RU148_T11 > 5
  ) %>%
  arrange(desc(mean_log2FC))

# Option 3: Broader approach focusing on the sample with stronger signal
upregulated_any_sample <- selected_data %>%
  filter(
    # Either sample shows upregulation with higher expression
    (log2FC_T8 > 1 & RU148_T8 > 10) | (log2FC_T11 > 1 & RU148_T11 > 10)
  ) %>%
  # Add metadata to identify which sample showed the stronger effect
  mutate(
    stronger_in_T8 = ifelse(log2FC_T8 > log2FC_T11, "Yes", "No"),
    significantly_upregulated = case_when(
      log2FC_T8 > 1 & log2FC_T11 > 1 ~ "Both samples",
      log2FC_T8 > 1 ~ "T8 only",
      log2FC_T11 > 1 ~ "T11 only",
      TRUE ~ "Neither"
    )
  ) %>%
  arrange(desc(mean_log2FC))

# Additional analysis: overlap between the two tumor samples
sample_comparison <- selected_data %>%
  filter(
    (log2FC_T8 > 1 & RU148_T8 > 10) | (log2FC_T11 > 1 & RU148_T11 > 10)
  ) %>%
  mutate(
    up_in_T8 = log2FC_T8 > 1 & RU148_T8 > 10,
    up_in_T11 = log2FC_T11 > 1 & RU148_T11 > 10,
    up_in_both = up_in_T8 & up_in_T11,
    pattern = case_when(
      up_in_both ~ "Both tumors",
      up_in_T8 ~ "T8 only",
      up_in_T11 ~ "T11 only",
      TRUE ~ "Neither"
    )
  )

# Create summary counts for each pattern
pattern_summary <- sample_comparison %>%
  group_by(pattern) %>%
  summarize(count = n()) %>%
  arrange(desc(count))

# Print counts of genes found
cat("Number of genes upregulated in either sample:", nrow(upregulated_individual), "\n")
cat("Number of genes upregulated in both samples:", nrow(upregulated_both_samples), "\n")
cat("Number of genes using broader approach:", nrow(upregulated_any_sample), "\n")

# Create a threshold sensitivity analysis
thresholds <- expand.grid(
  fc = c(0.5, 1, 1.5, 2),
  expr = c(1, 5, 10, 20)
)

results <- data.frame()
for (i in 1:nrow(thresholds)) {
  fc_thresh <- thresholds$fc[i]
  expr_thresh <- thresholds$expr[i]
  
  # Convert to log2 scale for fold change
  log2fc_thresh <- fc_thresh
  
  count_T8 <- sum(selected_data$log2FC_T8 > log2fc_thresh & 
                    selected_data$RU148_T8 > expr_thresh, na.rm=TRUE)
  count_T11 <- sum(selected_data$log2FC_T11 > log2fc_thresh & 
                     selected_data$RU148_T11 > expr_thresh, na.rm=TRUE)
  count_either <- sum((selected_data$log2FC_T8 > log2fc_thresh & selected_data$RU148_T8 > expr_thresh) | 
                        (selected_data$log2FC_T11 > log2fc_thresh & selected_data$RU148_T11 > expr_thresh), 
                      na.rm=TRUE)
  count_both <- sum((selected_data$log2FC_T8 > log2fc_thresh & selected_data$RU148_T8 > expr_thresh) & 
                      (selected_data$log2FC_T11 > log2fc_thresh & selected_data$RU148_T11 > expr_thresh), 
                    na.rm=TRUE)
  
  results <- rbind(results, data.frame(
    log2FC_threshold = log2fc_thresh,
    expression_threshold = expr_thresh,
    genes_T8 = count_T8,
    genes_T11 = count_T11,
    genes_either_tumor = count_either,
    genes_both_tumors = count_both
  ))
}

# Move to the output directory for saving files
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/RU148_transcriptome")

# Create a workbook for the multi-sheet Excel file
wb <- createWorkbook()

# Create detailed overview sheet with threshold information
overview <- data.frame(
  Sheet_Name = c(
    "Upregulated_either_tumor", 
    "Upregulated_both_tumors", 
    "Comprehensive_analysis", 
    "Sample_comparison", 
    "Pattern_summary", 
    "Threshold_sensitivity",
    "Top_100_genes",
    "Biotype_distribution",
    "Chromosome_distribution"
  ),
  Description = c(
    "Genes with log2FC > 1 and expression > 10 in at least one tumor",
    "Genes with log2FC > 0.5 and expression > 5 in both tumors",
    "Combined analysis with annotations about tumor patterns (log2FC > 1, expr > 10)",
    "All genes with pattern classification (log2FC > 1, expr > 10)",
    "Summary counts for each pattern (log2FC > 1, expr > 10)",
    "Analysis of gene counts with various thresholds",
    "Top 100 upregulated genes with detailed annotation",
    "Distribution of gene biotypes in upregulated genes",
    "Chromosomal distribution of upregulated genes"
  ),
  Log2FC_Threshold = c(
    "> 1 (2-fold) in at least one tumor", 
    "> 0.5 (1.4-fold) in both tumors", 
    "> 1 (2-fold) in at least one tumor", 
    "> 1 (2-fold) in at least one tumor", 
    "> 1 (2-fold) in at least one tumor", 
    "various",
    "> 1 (2-fold) in at least one tumor",
    "> 1 (2-fold) in at least one tumor",
    "> 1 (2-fold) in at least one tumor"
  ),
  Expression_Threshold = c(
    "> 10 in corresponding tumor", 
    "> 5 in both tumors", 
    "> 10 in corresponding tumor", 
    "> 10 in corresponding tumor", 
    "> 10 in corresponding tumor", 
    "various",
    "> 10 in corresponding tumor",
    "> 10 in corresponding tumor",
    "> 10 in corresponding tumor"
  ),
  Gene_Count = c(
    nrow(upregulated_individual),
    nrow(upregulated_both_samples),
    nrow(upregulated_any_sample),
    nrow(sample_comparison),
    nrow(pattern_summary),
    nrow(results),
    100,
    NA,
    NA
  )
)

# Add metadata to each data frame by creating a new column with threshold info
upregulated_individual$threshold_info <- "log2FC > 1, expr > 10"
upregulated_both_samples$threshold_info <- "log2FC > 0.5, expr > 5"
upregulated_any_sample$threshold_info <- "log2FC > 1, expr > 10"
sample_comparison$threshold_info <- "log2FC > 1, expr > 10"
pattern_summary$threshold_info <- "log2FC > 1, expr > 10"

# Add sheets to the workbook
addWorksheet(wb, "Overview", gridLines = TRUE)
writeData(wb, "Overview", overview)
setColWidths(wb, "Overview", cols = 1:5, widths = c(25, 60, 30, 30, 15))

# Now add the data sheets
addWorksheet(wb, "Upregulated_either_tumor")
writeData(wb, "Upregulated_either_tumor", upregulated_individual)

addWorksheet(wb, "Upregulated_both_tumors")
writeData(wb, "Upregulated_both_tumors", upregulated_both_samples)

addWorksheet(wb, "Comprehensive_analysis")
writeData(wb, "Comprehensive_analysis", upregulated_any_sample)

addWorksheet(wb, "Sample_comparison") 
writeData(wb, "Sample_comparison", sample_comparison)

addWorksheet(wb, "Pattern_summary")
writeData(wb, "Pattern_summary", pattern_summary)

addWorksheet(wb, "Threshold_sensitivity")
writeData(wb, "Threshold_sensitivity", results)

# Generate visualizations in both PDF and PNG formats
# ====================================================

# 1. Fold change distribution histograms
# PDF version
pdf("Fold_change_distribution.pdf", width=10, height=8)
par(mfrow = c(2, 1))
hist(selected_data$log2FC_T8, breaks=50, 
     main="Log2 Fold Change Distribution (T8 vs Normal)",
     xlab="Log2 Fold Change",
     col="lightblue")
hist(selected_data$log2FC_T11, breaks=50,
     main="Log2 Fold Change Distribution (T11 vs Normal)",
     xlab="Log2 Fold Change",
     col="lightgreen") 
dev.off()

# PNG version
png("Fold_change_distribution.png", width=1000, height=800, res=120)
par(mfrow = c(2, 1))
hist(selected_data$log2FC_T8, breaks=50, 
     main="Log2 Fold Change Distribution (T8 vs Normal)",
     xlab="Log2 Fold Change",
     col="lightblue")
hist(selected_data$log2FC_T11, breaks=50,
     main="Log2 Fold Change Distribution (T11 vs Normal)",
     xlab="Log2 Fold Change",
     col="lightgreen") 
dev.off()

# 2. Expression scatter plot
# Improved expression scatter plot that adds 0.1 to all values to avoid infinite values from log10(0)
expr_scatter <- ggplot(selected_data, aes(x=RU148_T8 + 0.1, y=RU148_T11 + 0.1, color=RU148_N + 0.1)) +
  geom_point(alpha=0.5) +
  scale_color_gradient(low="blue", high="red", trans="log10") +
  labs(title="Expression in Tumors with Normal Expression Colored",
       subtitle="Blue = low expression in normal tissue, Red = high expression in normal tissue",
       x="RU148_T8 Expression (log10 scale)", 
       y="RU148_T11 Expression (log10 scale)",
       color="Normal Expression (log10)") +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10()

# PDF version
ggsave("Expression_comparison_fixed.pdf", expr_scatter, width=10, height=8)

# PNG version
ggsave("Expression_comparison_fixed.png", expr_scatter, width=10, height=8, dpi=120)

# 3. Pattern summary bar plot
# PDF version
pdf("Pattern_summary.pdf", width=10, height=8)
barplot(pattern_summary$count, 
        names.arg = pattern_summary$pattern,
        col = c("purple", "lightblue", "pink", "gray"),
        main = "Patterns of Upregulation in Tumor Samples (log2FC > 1, expr > 10)",
        ylab = "Number of genes")
dev.off()

# PNG version
png("Pattern_summary.png", width=1000, height=800, res=120)
barplot(pattern_summary$count, 
        names.arg = pattern_summary$pattern,
        col = c("purple", "lightblue", "pink", "gray"),
        main = "Patterns of Upregulation in Tumor Samples (log2FC > 1, expr > 10)",
        ylab = "Number of genes")
dev.off()

# 4. Improved volcano plots
# Create improved volcano plots showing expression level vs fold change with better annotations
# T8 volcano plot
volcano_T8 <- ggplot(selected_data, aes(x=log2FC_T8, y=log10(RU148_T8 + 0.1))) +
  geom_point(aes(color=log2FC_T8 > 1 & RU148_T8 > 10), alpha=0.5) +
  scale_color_manual(values=c("gray", "red"), 
                     labels=c("Not significant", "log2FC > 1 & expr > 10")) +
  geom_vline(xintercept=1, linetype="dashed", color="blue") +
  geom_hline(yintercept=log10(10), linetype="dashed", color="blue") +
  labs(title="Modified Volcano Plot - T8 vs Normal",
       subtitle="Genes in red are upregulated in tumor compared to normal tissue",
       x="Log2 Fold Change (T8/Normal)",
       y="Log10 Expression Level in T8",
       color="Significance") +
  theme_minimal()

# T11 volcano plot
volcano_T11 <- ggplot(selected_data, aes(x=log2FC_T11, y=log10(RU148_T11 + 0.1))) +
  geom_point(aes(color=log2FC_T11 > 1 & RU148_T11 > 10), alpha=0.5) +
  scale_color_manual(values=c("gray", "red"), 
                     labels=c("Not significant", "log2FC > 1 & expr > 10")) +
  geom_vline(xintercept=1, linetype="dashed", color="blue") +
  geom_hline(yintercept=log10(10), linetype="dashed", color="blue") +
  labs(title="Modified Volcano Plot - T11 vs Normal",
       subtitle="Genes in red are upregulated in tumor compared to normal tissue",
       x="Log2 Fold Change (T11/Normal)",
       y="Log10 Expression Level in T11",
       color="Significance") +
  theme_minimal()

# PDF versions
ggsave("Volcano_T8_improved.pdf", volcano_T8, width=10, height=8)
ggsave("Volcano_T11_improved.pdf", volcano_T11, width=10, height=8)

# PNG versions
ggsave("Volcano_T8_improved.png", volcano_T8, width=10, height=8, dpi=120)
ggsave("Volcano_T11_improved.png", volcano_T11, width=10, height=8, dpi=120)

# 5. Combined tumor comparison plot
# Create a comparison of log2FC between the two tumor samples
comparison_plot <- ggplot(selected_data, aes(x=log2FC_T8, y=log2FC_T11)) +
  geom_point(aes(color=RU148_N > 10), alpha=0.5) +
  scale_color_manual(values=c("blue", "red"), 
                     labels=c("Normal expr <= 10", "Normal expr > 10")) +
  geom_hline(yintercept=1, linetype="dashed", color="darkgreen") +
  geom_vline(xintercept=1, linetype="dashed", color="darkgreen") +
  labs(title="Comparison of Log2FC Between Tumor Samples",
       subtitle="Points in different quadrants show tumor-specific expression patterns",
       x="Log2FC in T8",
       y="Log2FC in T11",
       color="Normal Expression") +
  theme_minimal() +
  # Add quadrant labels
  annotate("text", x=3, y=3, label="Both tumors", color="darkgreen", size=5, fontface="bold") +
  annotate("text", x=3, y=-1, label="T8 only", color="darkblue", size=5, fontface="bold") +
  annotate("text", x=-1, y=3, label="T11 only", color="darkred", size=5, fontface="bold") +
  annotate("text", x=-1, y=-1, label="Downregulated", color="gray40", size=5, fontface="bold")

# PDF version
ggsave("Tumor_comparison.pdf", comparison_plot, width=10, height=8)

# PNG version
ggsave("Tumor_comparison.png", comparison_plot, width=10, height=8, dpi=120)

#==============================================================================
# NEW ANALYSES
#==============================================================================

# 1. Heatmap of Top Differentially Expressed Genes
#------------------------------------------------
# Select top 50 upregulated genes based on mean log2FC
top_genes <- upregulated_any_sample %>%
  arrange(desc(mean_log2FC)) %>%
  head(50)

# Create a matrix for the heatmap
heatmap_data <- as.matrix(top_genes[, c("RU148_N", "RU148_T8", "RU148_T11")])
rownames(heatmap_data) <- top_genes$symbol

# Log transform for better visualization
heatmap_data <- log2(heatmap_data + 1)

# Create annotation for samples
sample_annotation <- data.frame(
  Sample_Type = c("Normal", "Tumor", "Tumor"),
  row.names = c("RU148_N", "RU148_T8", "RU148_T11")
)

# Set colors
ann_colors <- list(
  Sample_Type = c(Normal = "blue", Tumor = "red")
)

# Create heatmap
pdf("Top_genes_heatmap.pdf", width=10, height=12)
pheatmap(heatmap_data,
         scale = "row",  # Scale by row to see relative changes
         cluster_cols = FALSE,  # Don't cluster columns (samples)
         annotation_col = sample_annotation,
         annotation_colors = ann_colors,
         main = "Top 50 Upregulated Genes",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         fontsize_row = 8,
         cellwidth = 20,
         cellheight = 10)
dev.off()

# PNG version
png("Top_genes_heatmap.png", width=1000, height=1200, res=120)
pheatmap(heatmap_data,
         scale = "row",
         cluster_cols = FALSE,
         annotation_col = sample_annotation,
         annotation_colors = ann_colors,
         main = "Top 50 Upregulated Genes",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         fontsize_row = 8,
         cellwidth = 20,
         cellheight = 10)
dev.off()

# 2. Tumor-Specific vs Shared Expression Patterns
#------------------------------------------------
# Extract genes specific to each tumor and shared between them
genes_T8_only <- sample_comparison %>% filter(pattern == "T8 only") %>% pull(symbol)
genes_T11_only <- sample_comparison %>% filter(pattern == "T11 only") %>% pull(symbol)
genes_both <- sample_comparison %>% filter(pattern == "Both tumors") %>% pull(symbol)

# Create lists for Venn diagram
gene_lists <- list(
  T8 = genes_T8_only,
  T11 = genes_T11_only,
  Both = genes_both
)

# Generate Venn diagram
venn_colors <- c("skyblue", "pink", "purple")
venn.plot <- venn.diagram(
  x = list(
    "T8 specific" = genes_T8_only,
    "T11 specific" = genes_T11_only,
    "Shared" = genes_both
  ),
  filename = NULL,
  fill = venn_colors,
  alpha = 0.5,
  main = "Distribution of Upregulated Genes",
  sub = "log2FC > 1, expression > 10",
  main.cex = 2,
  cat.cex = 1.5
)

# Save as PDF and PNG
pdf("Upregulated_genes_venn.pdf", width=10, height=8)
grid.draw(venn.plot)
dev.off()

png("Upregulated_genes_venn.png", width=1000, height=800, res=120)
grid.draw(venn.plot)
dev.off()

# 3. PCA Analysis to See Overall Sample Relationships
#---------------------------------------------------
# PCA Analysis
# Prepare data for PCA - use all genes
pca_data <- t(as.matrix(selected_data[, c("RU148_N", "RU148_T8", "RU148_T11")]))

# Run PCA 
pca_result <- PCA(pca_data, scale.unit = TRUE, graph = FALSE)

# Visualize PCA
pdf("PCA_analysis.pdf", width=10, height=8)
fviz_pca_ind(pca_result, 
             col.ind = c("blue", "red", "red"),
             pointsize = 5,
             labelsize = 5,
             title = "PCA - Sample Relationships")
dev.off()

png("PCA_analysis.png", width=1000, height=800, res=120)
fviz_pca_ind(pca_result, 
             col.ind = c("blue", "red", "red"),
             pointsize = 5,
             labelsize = 5,
             title = "PCA - Sample Relationships")
dev.off()

# 4. Biotype Distribution Analysis
#--------------------------------
# Analyze distribution of gene biotypes in upregulated genes
biotype_distribution <- upregulated_any_sample %>%
  group_by(biotype) %>%
  summarize(count = n()) %>%
  arrange(desc(count))

# Plot biotype distribution
pdf("Biotype_distribution.pdf", width=12, height=8)
par(mar=c(10, 4, 4, 2))  # Increase bottom margin for labels
barplot(biotype_distribution$count, 
        names.arg = biotype_distribution$biotype,
        col = rainbow(nrow(biotype_distribution)),
        main = "Biotype Distribution in Upregulated Genes",
        ylab = "Number of genes",
        las = 2)  # Rotate labels
dev.off()

png("Biotype_distribution.png", width=1200, height=800, res=120)
par(mar=c(10, 4, 4, 2))  # Increase bottom margin for labels
barplot(biotype_distribution$count, 
        names.arg = biotype_distribution$biotype,
        col = rainbow(nrow(biotype_distribution)),
        main = "Biotype Distribution in Upregulated Genes",
        ylab = "Number of genes",
        las = 2)  # Rotate labels
dev.off()

# Save biotype distribution to Excel sheet
addWorksheet(wb, "Biotype_distribution")
writeData(wb, "Biotype_distribution", biotype_distribution)

# 5. Chromosomal Distribution of Upregulated Genes
#-----------------------------------------------
# Analyze chromosomal distribution
chrom_distribution <- upregulated_any_sample %>%
  group_by(chromosome) %>%
  summarize(count = n()) %>%
  arrange(desc(count))

# Plot chromosome distribution
pdf("Chromosome_distribution.pdf", width=12, height=8)
barplot(chrom_distribution$count, 
        names.arg = chrom_distribution$chromosome,
        col = topo.colors(nrow(chrom_distribution)),
        main = "Chromosomal Distribution of Upregulated Genes",
        ylab = "Number of genes",
        las = 2)
dev.off()

png("Chromosome_distribution.png", width=1200, height=800, res=120)
barplot(chrom_distribution$count, 
        names.arg = chrom_distribution$chromosome,
        col = topo.colors(nrow(chrom_distribution)),
        main = "Chromosomal Distribution of Upregulated Genes",
        ylab = "Number of genes",
        las = 2)
dev.off()

# Save chromosome distribution to Excel sheet
addWorksheet(wb, "Chromosome_distribution")
writeData(wb, "Chromosome_distribution", chrom_distribution)

# 6. Create Expression Ratio Plot
#------------------------------
# Calculate expression ratios
ratio_plot_data <- selected_data %>%
  mutate(
    T8_to_N_ratio = RU148_T8 / (RU148_N + 0.1),
    T11_to_N_ratio = RU148_T11 / (RU148_N + 0.1),
    significant = (log2FC_T8 > 1 & RU148_T8 > 10) | (log2FC_T11 > 1 & RU148_T11 > 10)
  )

# Plot expression ratios
ratio_plot <- ggplot(ratio_plot_data, aes(x=T8_to_N_ratio, y=T11_to_N_ratio)) +
  geom_point(aes(color=significant), alpha=0.5) +
  scale_color_manual(values=c("gray", "red"), 
                     labels=c("Not significant", "Significant")) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(title="Expression Ratios in Tumor vs Normal",
       subtitle="Diagonal line represents equal fold change in both tumor samples",
       x="T8/Normal Ratio (log10 scale)",
       y="T11/Normal Ratio (log10 scale)",
       color="Significance") +
  theme_minimal()

# Save plots
ggsave("Expression_ratio.pdf", ratio_plot, width=10, height=8)
ggsave("Expression_ratio.png", ratio_plot, width=10, height=8, dpi=120)

# 7. Table of Top Upregulated Genes with Annotation
#-----------------------------------------------
# Create a table of top 100 upregulated genes with detailed annotation
top_100_genes <- upregulated_any_sample %>%
  arrange(desc(mean_log2FC)) %>%
  head(100) %>%
  select(symbol, geneID, log2FC_T8, log2FC_T11, mean_log2FC, 
         RU148_N, RU148_T8, RU148_T11, 
         significantly_upregulated, biotype, chromosome, description)

# Create sheet in Excel file with this data
addWorksheet(wb, "Top_100_genes")
writeData(wb, "Top_100_genes", top_100_genes)

# 8. Expression levels across all three samples
#------------------------------------------
# Create a normalized expression matrix for visualization 
top_50_means <- top_genes %>%
  select(symbol, RU148_N, RU148_T8, RU148_T11) %>%
  mutate(mean_expr = (RU148_N + RU148_T8 + RU148_T11)/3) %>%
  arrange(desc(mean_expr)) %>%
  head(50)

# Reshape for ggplot
top_50_long <- top_50_means %>%
  select(symbol, RU148_N, RU148_T8, RU148_T11) %>%
  tidyr::pivot_longer(cols = c(RU148_N, RU148_T8, RU148_T11), 
                      names_to = "sample", 
                      values_to = "expression")

# Create a grouped bar plot of expression levels
expr_barplot <- ggplot(top_50_long, aes(x = reorder(symbol, -expression), y = expression, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("RU148_N" = "blue", "RU148_T8" = "red", "RU148_T11" = "orange"),
                    labels = c("RU148_N" = "Normal", "RU148_T8" = "Tumor T8", "RU148_T11" = "Tumor T11")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  labs(title = "Expression Levels of Top 50 Genes",
       subtitle = "Comparing normal and tumor samples",
       x = "Gene Symbol",
       y = "Expression Level",
       fill = "Sample")

# Save plots
ggsave("Top_genes_expression.pdf", expr_barplot, width=15, height=8)
ggsave("Top_genes_expression.png", expr_barplot, width=15, height=8, dpi=120)

# 9. Analysis of T8-specific and T11-specific genes
#------------------------------------------------
# Create additional sheets for tumor-specific gene lists
# T8-specific genes
t8_specific <- sample_comparison %>%
  filter(pattern == "T8 only") %>%
  select(symbol, geneID, log2FC_T8, log2FC_T11, mean_log2FC, 
         RU148_N, RU148_T8, RU148_T11, biotype, description) %>%
  arrange(desc(log2FC_T8))

addWorksheet(wb, "T8_specific_genes")
writeData(wb, "T8_specific_genes", t8_specific)

# T11-specific genes
t11_specific <- sample_comparison %>%
  filter(pattern == "T11 only") %>%
  select(symbol, geneID, log2FC_T8, log2FC_T11, mean_log2FC, 
         RU148_N, RU148_T8, RU148_T11, biotype, description) %>%
  arrange(desc(log2FC_T11))

addWorksheet(wb, "T11_specific_genes")
writeData(wb, "T11_specific_genes", t11_specific)

# Shared genes
shared_genes <- sample_comparison %>%
  filter(pattern == "Both tumors") %>%
  select(symbol, geneID, log2FC_T8, log2FC_T11, mean_log2FC, 
         RU148_N, RU148_T8, RU148_T11, biotype, description) %>%
  arrange(desc(mean_log2FC))

addWorksheet(wb, "Shared_upregulated_genes")
writeData(wb, "Shared_upregulated_genes", shared_genes)

# 10. Create a summary statistics sheet
#-----------------------------------
# Generate summary statistics for each sample
summary_stats <- data.frame(
  Statistic = c("Total genes analyzed", 
                "Genes upregulated in T8 only", 
                "Genes upregulated in T11 only",
                "Genes upregulated in both tumors",
                "Mean log2FC in T8",
                "Mean log2FC in T11",
                "Max log2FC in T8",
                "Max log2FC in T11"),
  Value = c(nrow(selected_data),
            nrow(t8_specific),
            nrow(t11_specific),
            nrow(shared_genes),
            mean(selected_data$log2FC_T8, na.rm=TRUE),
            mean(selected_data$log2FC_T11, na.rm=TRUE),
            max(selected_data$log2FC_T8, na.rm=TRUE),
            max(selected_data$log2FC_T11, na.rm=TRUE))
)

addWorksheet(wb, "Summary_Statistics")
writeData(wb, "Summary_Statistics", summary_stats)

# Save the final workbook
saveWorkbook(wb, "RU148_transcriptome_analysis.xlsx", overwrite = TRUE)

# Print completion message
cat("Complete analysis finished. Excel workbook with multiple sheets has been created in the RU148_transcriptome folder,")
cat("along with visualization PDFs and PNGs.\n")
cat("The analyses include:\n")
cat("1. Filter for upregulated genes with various thresholds\n")
cat("2. Multiple visualizations to explore expression patterns\n")
cat("3. Comparisons between the two tumor samples\n")
cat("4. Detailed gene annotations and distributions\n")
cat("5. Top gene rankings and pathway exploration\n")