library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(reshape2)

# Set working directory - update if needed
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation")

# Define input files
merged_file <- "merged_immunopeptidomics_data.tsv"
chop_file <- "final_CHOP_combined_peptide.tsv"
mskcc_file <- "MSKCC_combined_peptides.tsv"

# Read files
merged_data <- fread(merged_file, sep="\t", header=TRUE, fill=TRUE)
chop_data <- fread(chop_file, sep="\t", header=TRUE, fill=TRUE)
mskcc_data <- fread(mskcc_file, sep="\t", header=TRUE, fill=TRUE)

# Extract spectral count columns
mskcc_sc_cols <- grep("MSKCC_.*_Spectral_Count", colnames(merged_data), value=TRUE)
chop_sc_cols <- grep("FL.*Spectral Count", colnames(merged_data), value=TRUE)

# Extract sample IDs
extract_mskcc_id <- function(col_name) {
  parts <- strsplit(col_name, "_")[[1]]
  if(length(parts) >= 3) {
    return(parts[2])
  } else {
    return(NA)
  }
}

extract_chop_id <- function(col_name) {
  match <- regexpr("FL[0-9]+", col_name)
  if(match > 0) {
    fl_part <- regmatches(col_name, match)
    return(gsub("FL", "", fl_part))
  } else {
    return(NA)
  }
}

mskcc_sample_ids <- sapply(mskcc_sc_cols, extract_mskcc_id)
chop_sample_ids <- sapply(chop_sc_cols, extract_chop_id)

# Find matching samples
matching_samples <- intersect(mskcc_sample_ids, chop_sample_ids)
print(paste("Matching samples:", paste(matching_samples, collapse=", ")))

# Create correlation data
corr_data <- data.frame()
sample_overlap_stats <- list()

# Collect sample-by-sample data
for(sample in matching_samples) {
  mskcc_col <- grep(paste0("MSKCC_", sample, "_Spectral_Count"), colnames(merged_data), value=TRUE)[1]
  chop_col <- grep(paste0("FL", sample), chop_sc_cols, value=TRUE)[1]
  
  if(is.na(mskcc_col) || is.na(chop_col)) next
  
  # Get peptides detected in MSKCC and CHOP
  mskcc_peptides <- merged_data$Peptide[merged_data[[mskcc_col]] > 0]
  chop_peptides <- merged_data$Peptide[merged_data[[chop_col]] > 0]
  
  # Find common peptides
  common_peptides <- intersect(mskcc_peptides, chop_peptides)
  
  # Calculate overlap statistics
  total_unique <- length(union(mskcc_peptides, chop_peptides))
  overlap_pct <- 100 * length(common_peptides) / total_unique
  
  # Store overlap statistics
  sample_overlap_stats[[sample]] <- list(
    sample = sample,
    mskcc_count = length(mskcc_peptides),
    chop_count = length(chop_peptides),
    common_count = length(common_peptides),
    total_unique = total_unique,
    overlap_pct = overlap_pct
  )
  
  # Get peptides detected in both datasets for correlation
  detected_idx <- which(merged_data[[mskcc_col]] > 0 & merged_data[[chop_col]] > 0)
  
  if(length(detected_idx) >= 5) {
    sample_data <- data.frame(
      Sample = rep(sample, length(detected_idx)),
      MSKCC_Count = merged_data[[mskcc_col]][detected_idx],
      CHOP_Count = merged_data[[chop_col]][detected_idx],
      Peptide = merged_data$Peptide[detected_idx],
      stringsAsFactors = FALSE
    )
    corr_data <- rbind(corr_data, sample_data)
  }
}

# Create overlap summary dataframe
overlap_df <- do.call(rbind, lapply(sample_overlap_stats, function(x) {
  data.frame(
    Sample = x$sample,
    MSKCC_Peptides = x$mskcc_count,
    CHOP_Peptides = x$chop_count,
    Common_Peptides = x$common_count,
    Total_Unique = x$total_unique,
    Overlap_Percent = x$overlap_pct
  )
}))

# Calculate correlation statistics by sample
sample_stats <- data.frame()

for(sample in unique(corr_data$Sample)) {
  sample_data <- corr_data[corr_data$Sample == sample, ]
  
  if(nrow(sample_data) >= 5) {
    # Calculate Pearson and Spearman correlations
    pearson_cor <- cor(sample_data$MSKCC_Count, sample_data$CHOP_Count, 
                       method="pearson", use="pairwise.complete.obs")
    spearman_cor <- cor(sample_data$MSKCC_Count, sample_data$CHOP_Count, 
                        method="spearman", use="pairwise.complete.obs")
    
    # Calculate R-squared from linear model
    model <- lm(CHOP_Count ~ MSKCC_Count, data=sample_data)
    r_squared <- summary(model)$r.squared
    
    # Add to stats dataframe
    sample_stats <- rbind(sample_stats, data.frame(
      Sample = sample,
      Pearson_Correlation = pearson_cor,
      Spearman_Correlation = spearman_cor,
      R_Squared = r_squared,
      Peptide_Count = nrow(sample_data)
    ))
  }
}

# Sort samples for consistent visualization
sample_stats$Sample <- factor(sample_stats$Sample, 
                              levels=sample_stats$Sample[order(as.numeric(sample_stats$Sample))])
overlap_df$Sample <- factor(overlap_df$Sample, 
                            levels=overlap_df$Sample[order(as.numeric(overlap_df$Sample))])

# ===== VISUALIZATION 1: SAMPLE OVERLAP COMPARISON =====
# Create overlapping bar chart
p_overlap <- ggplot(overlap_df, aes(x=Sample)) +
  geom_bar(aes(y=MSKCC_Peptides, fill="MSKCC Only"), stat="identity", alpha=0.7) +
  geom_bar(aes(y=CHOP_Peptides, fill="CHOP Only"), stat="identity", alpha=0.7) +
  geom_bar(aes(y=Common_Peptides, fill="Common"), stat="identity") +
  scale_fill_brewer(palette="Set1", name="Dataset") +
  geom_text(aes(y=MSKCC_Peptides, label=MSKCC_Peptides), 
            position=position_stack(vjust=0.5), color="white", size=3.5) +
  geom_text(aes(y=CHOP_Peptides, label=CHOP_Peptides), 
            position=position_stack(vjust=0.5), color="white", size=3.5) +
  geom_text(aes(y=Common_Peptides, label=Common_Peptides), 
            position=position_stack(vjust=0.5), color="white", size=3.5) +
  labs(title="Peptide Detection by Sample", 
       subtitle="Number of peptides detected in each dataset", 
       x="Sample", y="Number of Peptides") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5))

# Create percentage overlap plot
p_pct_overlap <- ggplot(overlap_df, aes(x=Sample, y=Overlap_Percent, fill=Sample)) +
  geom_bar(stat="identity", alpha=0.8) +
  geom_text(aes(label=sprintf("%.1f%%", Overlap_Percent)), vjust=-0.5) +
  ylim(0, max(30, max(overlap_df$Overlap_Percent) * 1.1)) +
  scale_fill_brewer(palette="Set2") +
  labs(title="Peptide Overlap Percentage by Sample",
       x="Sample", y="Overlap Percentage (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        legend.position="none")

# ===== VISUALIZATION 2: CORRELATION STATISTICS =====
# Create correlation statistics plot
p_correlation <- ggplot(sample_stats) +
  geom_bar(aes(x=Sample, y=Pearson_Correlation, fill="Pearson"), 
           stat="identity", position=position_dodge(), alpha=0.7) +
  geom_bar(aes(x=Sample, y=Spearman_Correlation, fill="Spearman"), 
           stat="identity", position=position_dodge(), alpha=0.7) +
  geom_bar(aes(x=Sample, y=R_Squared, fill="R²"), 
           stat="identity", position=position_dodge(), alpha=0.7) +
  geom_text(aes(x=Sample, y=Pearson_Correlation, 
                label=sprintf("%.2f", Pearson_Correlation)), 
            position=position_dodge(width=0.9), vjust=-0.5, size=3) +
  geom_text(aes(x=Sample, y=Spearman_Correlation, 
                label=sprintf("%.2f", Spearman_Correlation)), 
            position=position_dodge(width=0.9), vjust=-0.5, size=3) +
  geom_text(aes(x=Sample, y=R_Squared, 
                label=sprintf("%.2f", R_Squared)), 
            position=position_dodge(width=0.9), vjust=-0.5, size=3) +
  ylim(0, max(1, max(c(sample_stats$Pearson_Correlation, 
                       sample_stats$Spearman_Correlation, 
                       sample_stats$R_Squared)) * 1.2)) +
  scale_fill_brewer(palette="Set1", name="Metric") +
  labs(title="Correlation Metrics by Sample",
       subtitle="Pearson, Spearman, and R² values for overlapping peptides",
       x="Sample", y="Correlation Value") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5))

# Add number of peptides to plot title
p_correlation_n <- ggplot(sample_stats, aes(x=Sample, y=Peptide_Count, fill=Sample)) +
  geom_bar(stat="identity", alpha=0.8) +
  geom_text(aes(label=Peptide_Count), vjust=-0.5) +
  scale_fill_brewer(palette="Set2") +
  labs(title="Number of Common Peptides Used for Correlation",
       x="Sample", y="Peptide Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        legend.position="none")

# ===== VISUALIZATION 3: SAMPLE CORRELATION SCATTER PLOTS =====
# Create scatter plots for each sample
sample_plots <- list()

for(sample in unique(corr_data$Sample)) {
  sample_data <- corr_data[corr_data$Sample == sample, ]
  
  if(nrow(sample_data) >= 5) {
    # Get correlation values
    pearson_cor <- cor(sample_data$MSKCC_Count, sample_data$CHOP_Count, 
                       method="pearson", use="pairwise.complete.obs")
    r_squared <- sample_stats$R_Squared[sample_stats$Sample == sample]
    
    # Create scatter plot
    p <- ggplot(sample_data, aes(x=MSKCC_Count, y=CHOP_Count)) +
      geom_point(alpha=0.5, color="blue") +
      geom_smooth(method="lm", color="red", se=TRUE) +
      labs(title=paste0("Sample ", sample, 
                        " (R² = ", sprintf("%.3f", r_squared), 
                        ", n = ", nrow(sample_data), ")"),
           x="MSKCC Spectral Count", 
           y="CHOP Spectral Count") +
      theme_minimal() +
      theme(plot.title = element_text(hjust=0.5))
    
    sample_plots[[sample]] <- p
  }
}

# ===== VISUALIZATION 4: SUMMARY STATISTICS =====
# Create summary statistics
summary_stats <- data.frame(
  Metric = c(
    "Total MSKCC unique peptides",
    "Total CHOP unique peptides",
    "Common peptides between datasets",
    "Overall overlap percentage",
    "Average Pearson correlation",
    "Average Spearman correlation",
    "Average R-squared value",
    "Number of matching samples"
  ),
  Value = c(
    length(unique(mskcc_data$Peptide)),
    ifelse("Peptide" %in% colnames(chop_data), 
           length(unique(chop_data$Peptide)), 
           length(unique(chop_data$`Peptide Sequence`))),
    sum(merged_data$Peptide %in% intersect(
      unique(mskcc_data$Peptide), 
      unique(ifelse("Peptide" %in% colnames(chop_data), 
                    chop_data$Peptide, 
                    chop_data$`Peptide Sequence`)))),
    round(100 * length(intersect(unique(mskcc_data$Peptide), 
                                 unique(ifelse("Peptide" %in% colnames(chop_data), 
                                               chop_data$Peptide, 
                                               chop_data$`Peptide Sequence`)))) / 
            length(union(unique(mskcc_data$Peptide), 
                         unique(ifelse("Peptide" %in% colnames(chop_data), 
                                       chop_data$Peptide, 
                                       chop_data$`Peptide Sequence`)))), 1),
    round(mean(sample_stats$Pearson_Correlation), 3),
    round(mean(sample_stats$Spearman_Correlation), 3),
    round(mean(sample_stats$R_Squared), 3),
    length(matching_samples)
  )
)

# Format summary statistics table
summary_table <- tableGrob(summary_stats, rows=NULL, theme=ttheme_minimal(
  core=list(fg_params=list(hjust=0, x=0.1), bg_params=list(fill=NA)),
  colhead=list(fg_params=list(hjust=0.5, fontface="bold"))
))

# ===== SAVE ALL VISUALIZATIONS TO PDF =====
pdf("immunopeptidomics_visualization_report.pdf", width=11, height=8.5)

# Page 1: Sample overlap and correlation statistics
grid.arrange(
  p_overlap, p_pct_overlap,
  p_correlation, p_correlation_n,
  ncol=2, nrow=2,
  top=textGrob("Immunopeptidomics Datasets Comparison", 
               gp=gpar(fontsize=16, fontface="bold"))
)

# Page 2: Correlation scatter plots
sample_plot_list <- lapply(names(sample_plots), function(s) sample_plots[[s]])

do.call(grid.arrange, 
        c(sample_plot_list, 
          list(ncol=min(3, length(sample_plot_list)), 
               top=textGrob("Sample-by-Sample Correlation (MSKCC vs CHOP)",
                            gp=gpar(fontsize=16, fontface="bold"))))
)

# Page 3: Summary statistics
grid.arrange(
  summary_table,
  top=textGrob("Summary Statistics", gp=gpar(fontsize=16, fontface="bold"))
)

dev.off()

# Save the main figures as individual PNGs for easy viewing
ggsave("sample_overlap.png", p_overlap, width=8, height=6)
ggsave("overlap_percentage.png", p_pct_overlap, width=8, height=6)
ggsave("correlation_metrics.png", p_correlation, width=8, height=6)
ggsave("common_peptides_count.png", p_correlation_n, width=8, height=6)

# Print completion message
cat("Visualizations created and saved to:\n")
cat("1. immunopeptidomics_visualization_report.pdf (comprehensive report)\n")
cat("2. sample_overlap.png\n")
cat("3. overlap_percentage.png\n")
cat("4. correlation_metrics.png\n")
cat("5. common_peptides_count.png\n")

# Print summary statistics
cat("\nSummary Statistics:\n")
print(summary_stats)

# Print correlation statistics by sample
cat("\nCorrelation Statistics by Sample:\n")
print(sample_stats[order(as.numeric(as.character(sample_stats$Sample))), ])

# Print overlap statistics by sample
cat("\nOverlap Statistics by Sample:\n")
print(overlap_df[order(as.numeric(as.character(overlap_df$Sample))), ])
