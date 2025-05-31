

# # Load required libraries
# library(readxl)
# library(writexl)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(openxlsx)  # For creating multi-sheet Excel files
# 
# # Create the output directory if it doesn't exist
# dir.create("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/RU148_transcriptome", showWarnings = FALSE)
# 
# # Set working directory to the raw data folder to read the input file
# setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata")
# 
# # Read the Excel file
# gene_data <- read_excel("Normalized_Gene_counts_FLCdb_Panel_1.xlsx")
# 
# # Select columns of interest
# selected_data <- gene_data %>%
#   select(RU148_N, RU148_T8, RU148_T11, Mean.Normal, Mean.Tumor, 
#          geneID, symbol, biotype, chromosome, gene_start, gene_end, gene_length, description)
# 
# # Calculate log2 fold change for each tumor sample compared to normal
# selected_data <- selected_data %>%
#   mutate(
#     log2FC_T8 = log2((RU148_T8 + 0.1) / (RU148_N + 0.1)),
#     log2FC_T11 = log2((RU148_T11 + 0.1) / (RU148_N + 0.1)),
#     mean_log2FC = (log2((RU148_T8 + 0.1) / (RU148_N + 0.1)) + 
#                      log2((RU148_T11 + 0.1) / (RU148_N + 0.1))) / 2
#   )
# 
# # Option 1: Using individual sample criteria
# upregulated_individual <- selected_data %>%
#   filter(
#     # Require log2FC > 1 in at least one tumor sample
#     (log2FC_T8 > 1 | log2FC_T11 > 1),
#     # Ensure expression level > 10 in the corresponding sample
#     ((log2FC_T8 > 1 & RU148_T8 > 10) | (log2FC_T11 > 1 & RU148_T11 > 10))
#   ) %>%
#   arrange(desc(mean_log2FC))
# 
# # Option 2: Focus on genes up in both tumor samples compared to normal
# upregulated_both_samples <- selected_data %>%
#   filter(
#     # Up in both samples (less stringent threshold)
#     log2FC_T8 > 0.5 & log2FC_T11 > 0.5,
#     # Reasonable expression in both
#     RU148_T8 > 5 & RU148_T11 > 5
#   ) %>%
#   arrange(desc(mean_log2FC))
# 
# # Option 3: Broader approach focusing on the sample with stronger signal
# upregulated_any_sample <- selected_data %>%
#   filter(
#     # Either sample shows upregulation with higher expression
#     (log2FC_T8 > 1 & RU148_T8 > 10) | (log2FC_T11 > 1 & RU148_T11 > 10)
#   ) %>%
#   # Add metadata to identify which sample showed the stronger effect
#   mutate(
#     stronger_in_T8 = ifelse(log2FC_T8 > log2FC_T11, "Yes", "No"),
#     significantly_upregulated = case_when(
#       log2FC_T8 > 1 & log2FC_T11 > 1 ~ "Both samples",
#       log2FC_T8 > 1 ~ "T8 only",
#       log2FC_T11 > 1 ~ "T11 only",
#       TRUE ~ "Neither"
#     )
#   ) %>%
#   arrange(desc(mean_log2FC))
# 
# # Additional analysis: overlap between the two tumor samples
# sample_comparison <- selected_data %>%
#   filter(
#     (log2FC_T8 > 1 & RU148_T8 > 10) | (log2FC_T11 > 1 & RU148_T11 > 10)
#   ) %>%
#   mutate(
#     up_in_T8 = log2FC_T8 > 1 & RU148_T8 > 10,
#     up_in_T11 = log2FC_T11 > 1 & RU148_T11 > 10,
#     up_in_both = up_in_T8 & up_in_T11,
#     pattern = case_when(
#       up_in_both ~ "Both tumors",
#       up_in_T8 ~ "T8 only",
#       up_in_T11 ~ "T11 only",
#       TRUE ~ "Neither"
#     )
#   )
# 
# # Create summary counts for each pattern
# pattern_summary <- sample_comparison %>%
#   group_by(pattern) %>%
#   summarize(count = n()) %>%
#   arrange(desc(count))
# 
# # Print counts of genes found
# cat("Number of genes upregulated in either sample:", nrow(upregulated_individual), "\n")
# cat("Number of genes upregulated in both samples:", nrow(upregulated_both_samples), "\n")
# cat("Number of genes using broader approach:", nrow(upregulated_any_sample), "\n")
# 
# # Create a threshold sensitivity analysis
# thresholds <- expand.grid(
#   fc = c(0.5, 1, 1.5, 2),
#   expr = c(1, 5, 10, 20)
# )
# 
# results <- data.frame()
# for (i in 1:nrow(thresholds)) {
#   fc_thresh <- thresholds$fc[i]
#   expr_thresh <- thresholds$expr[i]
#   
#   # Convert to log2 scale for fold change
#   log2fc_thresh <- fc_thresh
#   
#   count_T8 <- sum(selected_data$log2FC_T8 > log2fc_thresh & 
#                     selected_data$RU148_T8 > expr_thresh, na.rm=TRUE)
#   count_T11 <- sum(selected_data$log2FC_T11 > log2fc_thresh & 
#                      selected_data$RU148_T11 > expr_thresh, na.rm=TRUE)
#   count_either <- sum((selected_data$log2FC_T8 > log2fc_thresh & selected_data$RU148_T8 > expr_thresh) | 
#                         (selected_data$log2FC_T11 > log2fc_thresh & selected_data$RU148_T11 > expr_thresh), 
#                       na.rm=TRUE)
#   count_both <- sum((selected_data$log2FC_T8 > log2fc_thresh & selected_data$RU148_T8 > expr_thresh) & 
#                       (selected_data$log2FC_T11 > log2fc_thresh & selected_data$RU148_T11 > expr_thresh), 
#                     na.rm=TRUE)
#   
#   results <- rbind(results, data.frame(
#     log2FC_threshold = log2fc_thresh,
#     expression_threshold = expr_thresh,
#     genes_T8 = count_T8,
#     genes_T11 = count_T11,
#     genes_either_tumor = count_either,
#     genes_both_tumors = count_both
#   ))
# }
# 
# # Move to the output directory for saving files
# setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/RU148_transcriptome")
# 
# # Create a workbook for the multi-sheet Excel file
# wb <- createWorkbook()
# 
# # Create detailed overview sheet with threshold information
# overview <- data.frame(
#   Sheet_Name = c(
#     "Upregulated_either_tumor", 
#     "Upregulated_both_tumors", 
#     "Comprehensive_analysis", 
#     "Sample_comparison", 
#     "Pattern_summary", 
#     "Threshold_sensitivity"
#   ),
#   Description = c(
#     "Genes with log2FC > 1 and expression > 10 in at least one tumor",
#     "Genes with log2FC > 0.5 and expression > 5 in both tumors",
#     "Combined analysis with annotations about tumor patterns (log2FC > 1, expr > 10)",
#     "All genes with pattern classification (log2FC > 1, expr > 10)",
#     "Summary counts for each pattern (log2FC > 1, expr > 10)",
#     "Analysis of gene counts with various thresholds"
#   ),
#   Log2FC_Threshold = c(
#     "> 1 (2-fold) in at least one tumor", 
#     "> 0.5 (1.4-fold) in both tumors", 
#     "> 1 (2-fold) in at least one tumor", 
#     "> 1 (2-fold) in at least one tumor", 
#     "> 1 (2-fold) in at least one tumor", 
#     "various"
#   ),
#   Expression_Threshold = c(
#     "> 10 in corresponding tumor", 
#     "> 5 in both tumors", 
#     "> 10 in corresponding tumor", 
#     "> 10 in corresponding tumor", 
#     "> 10 in corresponding tumor", 
#     "various"
#   ),
#   Gene_Count = c(
#     nrow(upregulated_individual),
#     nrow(upregulated_both_samples),
#     nrow(upregulated_any_sample),
#     nrow(sample_comparison),
#     nrow(pattern_summary),
#     nrow(results)
#   )
# )
# 
# # Add metadata to each data frame by creating a new column with threshold info
# upregulated_individual$threshold_info <- "log2FC > 1, expr > 10"
# upregulated_both_samples$threshold_info <- "log2FC > 0.5, expr > 5"
# upregulated_any_sample$threshold_info <- "log2FC > 1, expr > 10"
# sample_comparison$threshold_info <- "log2FC > 1, expr > 10"
# pattern_summary$threshold_info <- "log2FC > 1, expr > 10"
# 
# # Add sheets to the workbook
# addWorksheet(wb, "Overview", gridLines = TRUE)
# writeData(wb, "Overview", overview)
# setColWidths(wb, "Overview", cols = 1:5, widths = c(25, 60, 30, 30, 15))
# 
# # Now add the data sheets
# addWorksheet(wb, "Upregulated_either_tumor")
# writeData(wb, "Upregulated_either_tumor", upregulated_individual)
# 
# addWorksheet(wb, "Upregulated_both_tumors")
# writeData(wb, "Upregulated_both_tumors", upregulated_both_samples)
# 
# addWorksheet(wb, "Comprehensive_analysis")
# writeData(wb, "Comprehensive_analysis", upregulated_any_sample)
# 
# addWorksheet(wb, "Sample_comparison") 
# writeData(wb, "Sample_comparison", sample_comparison)
# 
# addWorksheet(wb, "Pattern_summary")
# writeData(wb, "Pattern_summary", pattern_summary)
# 
# addWorksheet(wb, "Threshold_sensitivity")
# writeData(wb, "Threshold_sensitivity", results)
# 
# # Save the workbook
# saveWorkbook(wb, "RU148_transcriptome_analysis.xlsx", overwrite = TRUE)
# 
# # Generate visualizations in both PDF and PNG formats
# # ====================================================
# 
# # 1. Fold change distribution histograms
# # PDF version
# pdf("Fold_change_distribution.pdf", width=10, height=8)
# par(mfrow = c(2, 1))
# hist(selected_data$log2FC_T8, breaks=50, 
#      main="Log2 Fold Change Distribution (T8 vs Normal)",
#      xlab="Log2 Fold Change",
#      col="lightblue")
# hist(selected_data$log2FC_T11, breaks=50,
#      main="Log2 Fold Change Distribution (T11 vs Normal)",
#      xlab="Log2 Fold Change",
#      col="lightgreen") 
# dev.off()
# 
# # PNG version
# png("Fold_change_distribution.png", width=1000, height=800, res=120)
# par(mfrow = c(2, 1))
# hist(selected_data$log2FC_T8, breaks=50, 
#      main="Log2 Fold Change Distribution (T8 vs Normal)",
#      xlab="Log2 Fold Change",
#      col="lightblue")
# hist(selected_data$log2FC_T11, breaks=50,
#      main="Log2 Fold Change Distribution (T11 vs Normal)",
#      xlab="Log2 Fold Change",
#      col="lightgreen") 
# dev.off()
# 
# # 2. Expression scatter plot
# # Improved expression scatter plot that add 0.1 to all values to avoid infinite values from log10(0)
# expr_scatter <- ggplot(selected_data, aes(x=RU148_T8 + 0.1, y=RU148_T11 + 0.1, color=RU148_N + 0.1)) +
#   geom_point(alpha=0.5) +
#   scale_color_gradient(low="blue", high="red", trans="log10") +
#   labs(title="Expression in Tumors with Normal Expression Colored",
#        x="RU148_T8 Expression (log10 scale)", 
#        y="RU148_T11 Expression (log10 scale)",
#        color="Normal Expression (log10)") +
#   theme_minimal() +
#   scale_x_log10() +
#   scale_y_log10()
# 
# # PDF version
# ggsave("Expression_comparison_fixed.pdf", expr_scatter, width=10, height=8)
# 
# # PNG version
# ggsave("Expression_comparison_fixed.png", expr_scatter, width=10, height=8, dpi=120)
# 
# # 3. Pattern summary bar plot
# # PDF version
# pdf("Pattern_summary.pdf", width=10, height=8)
# barplot(pattern_summary$count, 
#         names.arg = pattern_summary$pattern,
#         col = c("purple", "lightblue", "pink", "gray"),
#         main = "Patterns of Upregulation in Tumor Samples (log2FC > 1, expr > 10)",
#         ylab = "Number of genes")
# dev.off()
# 
# # PNG version
# png("Pattern_summary.png", width=1000, height=800, res=120)
# barplot(pattern_summary$count, 
#         names.arg = pattern_summary$pattern,
#         col = c("purple", "lightblue", "pink", "gray"),
#         main = "Patterns of Upregulation in Tumor Samples (log2FC > 1, expr > 10)",
#         ylab = "Number of genes")
# dev.off()
# 
# # 4. Additional visualization: Volcano plot for both tumor samples
# # Create volcano plots showing significance vs fold change
# # T8 volcano plot
# volcano_T8 <- ggplot(selected_data, aes(x=log2FC_T8, y=-log10(0.05))) +
#   geom_point(aes(color=log2FC_T8 > 1 & RU148_T8 > 10), alpha=0.5) +
#   scale_color_manual(values=c("gray", "red"), 
#                      labels=c("Not significant", "log2FC > 1 & expr > 10")) +
#   geom_vline(xintercept=1, linetype="dashed", color="blue") +
#   labs(title="Volcano Plot - T8 vs Normal",
#        x="Log2 Fold Change",
#        y="-Log10(p-value) [Placeholder]",
#        color="Significance") +
#   theme_minimal()
# 
# # T11 volcano plot
# volcano_T11 <- ggplot(selected_data, aes(x=log2FC_T11, y=-log10(0.05))) +
#   geom_point(aes(color=log2FC_T11 > 1 & RU148_T11 > 10), alpha=0.5) +
#   scale_color_manual(values=c("gray", "red"), 
#                      labels=c("Not significant", "log2FC > 1 & expr > 10")) +
#   geom_vline(xintercept=1, linetype="dashed", color="blue") +
#   labs(title="Volcano Plot - T11 vs Normal",
#        x="Log2 Fold Change",
#        y="-Log10(p-value) [Placeholder]",
#        color="Significance") +
#   theme_minimal()
# 
# # PDF versions
# ggsave("Volcano_T8.pdf", volcano_T8, width=10, height=8)
# ggsave("Volcano_T11.pdf", volcano_T11, width=10, height=8)
# 
# # PNG versions
# ggsave("Volcano_T8.png", volcano_T8, width=10, height=8, dpi=120)
# ggsave("Volcano_T11.png", volcano_T11, width=10, height=8, dpi=120)
# 
# # 5. Combined tumor comparison plot
# # Create a comparison of log2FC between the two tumor samples
# comparison_plot <- ggplot(selected_data, aes(x=log2FC_T8, y=log2FC_T11)) +
#   geom_point(aes(color=RU148_N > 10), alpha=0.5) +
#   scale_color_manual(values=c("blue", "red"), 
#                      labels=c("Normal expr <= 10", "Normal expr > 10")) +
#   geom_hline(yintercept=1, linetype="dashed", color="darkgreen") +
#   geom_vline(xintercept=1, linetype="dashed", color="darkgreen") +
#   labs(title="Comparison of Log2FC Between Tumor Samples",
#        x="Log2FC in T8",
#        y="Log2FC in T11",
#        color="Normal Expression") +
#   theme_minimal()
# 
# # PDF version
# ggsave("Tumor_comparison.pdf", comparison_plot, width=10, height=8)
# 
# # PNG version
# ggsave("Tumor_comparison.png", comparison_plot, width=10, height=8, dpi=120)
# 
# cat("Analysis complete. A single Excel file with multiple sheets has been created in the RU148_transcriptome folder,")
# cat("along with visualization PDFs and PNGs.")