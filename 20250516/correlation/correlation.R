# This script correlates CHOP data with our immunopeptidomics data that we ran at MSKCC
# First run the script titled "generate_final_CHOP.R"
# Take the combined_peptides.tsv unmodified from the combine2c3v script
# and rename is to MSKCC_combined_peptides.tsv
# REMEMBER TO UPDATE THIS FILE PATH!!!!!

library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(reshape2)
library(dplyr)

# Set working directory - update if needed
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation")

# Define input files
merged_file <- "merged_immunopeptidomics_data.tsv"
chop_file <- "final_CHOP_combined_peptide.tsv"
mskcc_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation/unique_peptides_unmodified.tsv"

# Read files
cat("Reading merged data file...\n")
merged_data <- fread(merged_file, sep="\t", header=TRUE, fill=TRUE)
cat("Reading CHOP data file...\n")
chop_data <- fread(chop_file, sep="\t", header=TRUE, fill=TRUE)
cat("Reading MSKCC data file...\n")
mskcc_data <- fread(mskcc_file, sep="\t", header=TRUE, fill=TRUE)

# Ensure consistent peptide column names
if("Peptide Sequence" %in% colnames(mskcc_data) && !"Peptide" %in% colnames(mskcc_data)) {
  colnames(mskcc_data)[colnames(mskcc_data) == "Peptide Sequence"] <- "Peptide"
}

if("Peptide Sequence" %in% colnames(chop_data) && !"Peptide" %in% colnames(chop_data)) {
  colnames(chop_data)[colnames(chop_data) == "Peptide Sequence"] <- "Peptide"
}

# Extract spectral count columns
mskcc_sc_cols <- grep("MSKCC_.*_Spectral_Count", colnames(merged_data), value=TRUE)
chop_sc_cols <- grep("FL.*Spectral Count", colnames(merged_data), value=TRUE)

# Extract intensity columns
mskcc_int_cols <- grep("MSKCC_.*_Intensity", colnames(merged_data), value=TRUE)
chop_int_cols <- grep("FL.*Intensity", colnames(merged_data), value=TRUE)

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
  # Check for specific patterns
  if(grepl("FL57Liver|H5_FL57", col_name)) {
    return("57Liver")
  } else if(grepl("FL57Lung|H6_FL57", col_name)) {
    return("57Lung")
  } else if(grepl("FL([0-9]+)", col_name)) {
    match <- regexpr("FL([0-9]+)", col_name)
    fl_part <- regmatches(col_name, match)
    return(gsub("FL", "", fl_part))
  } else {
    match2 <- regexpr("H[0-9]+_FL([0-9]+)", col_name)
    if(match2 > 0) {
      full_match <- regmatches(col_name, match2)
      fl_num <- gsub(".*FL", "", full_match)
      return(fl_num)
    }
    return(NA)
  }
}

mskcc_sample_ids <- sapply(mskcc_sc_cols, extract_mskcc_id)
chop_sample_ids <- sapply(chop_sc_cols, extract_chop_id)

# Handle special case for sample 57 - map it to both 57Liver and 57Lung
if("57" %in% mskcc_sample_ids) {
  # Create additional mappings for 57
  mskcc_col_57 <- grep("MSKCC_57_Spectral_Count", colnames(merged_data), value = TRUE)[1]
  
  # Manually add entries for 57Liver and 57Lung if they exist in CHOP
  if("57Liver" %in% chop_sample_ids) {
    chop_col_57liver <- grep("57Liver", chop_sc_cols, value = TRUE)[1]
    
    # Get peptides detected in MSKCC and CHOP
    mskcc_peptides_57 <- merged_data$Peptide[merged_data[[mskcc_col_57]] > 0]
    chop_peptides_57liver <- merged_data$Peptide[merged_data[[chop_col_57liver]] > 0]
    
    # Find common peptides
    common_peptides_57liver <- intersect(mskcc_peptides_57, chop_peptides_57liver)
    
    # Calculate overlap statistics
    total_unique_57liver <- length(union(mskcc_peptides_57, chop_peptides_57liver))
    overlap_pct_57liver <- 100 * length(common_peptides_57liver) / total_unique_57liver
    
    # Store overlap statistics
    sample_overlap_stats[["57Liver"]] <- list(
      sample = "57Liver",
      mskcc_count = length(mskcc_peptides_57),
      chop_count = length(chop_peptides_57liver),
      common_count = length(common_peptides_57liver),
      total_unique = total_unique_57liver,
      overlap_pct = overlap_pct_57liver
    )
    
    # Get peptides detected in both datasets for correlation
    detected_idx_57liver <- which(merged_data[[mskcc_col_57]] > 0 & merged_data[[chop_col_57liver]] > 0)
    
    if(length(detected_idx_57liver) >= 5) {
      sample_data_57liver <- data.frame(
        Sample = rep("57Liver", length(detected_idx_57liver)),
        MSKCC_Count = merged_data[[mskcc_col_57]][detected_idx_57liver],
        CHOP_Count = merged_data[[chop_col_57liver]][detected_idx_57liver],
        Peptide = merged_data$Peptide[detected_idx_57liver],
        stringsAsFactors = FALSE
      )
      corr_data <- rbind(corr_data, sample_data_57liver)
    }
  }
  
  if("57Lung" %in% chop_sample_ids) {
    chop_col_57lung <- grep("57Lung", chop_sc_cols, value = TRUE)[1]
    
    # Get peptides detected in MSKCC and CHOP
    mskcc_peptides_57 <- merged_data$Peptide[merged_data[[mskcc_col_57]] > 0]
    chop_peptides_57lung <- merged_data$Peptide[merged_data[[chop_col_57lung]] > 0]
    
    # Find common peptides
    common_peptides_57lung <- intersect(mskcc_peptides_57, chop_peptides_57lung)
    
    # Calculate overlap statistics
    total_unique_57lung <- length(union(mskcc_peptides_57, chop_peptides_57lung))
    overlap_pct_57lung <- 100 * length(common_peptides_57lung) / total_unique_57lung
    
    # Store overlap statistics
    sample_overlap_stats[["57Lung"]] <- list(
      sample = "57Lung",
      mskcc_count = length(mskcc_peptides_57),
      chop_count = length(chop_peptides_57lung),
      common_count = length(common_peptides_57lung),
      total_unique = total_unique_57lung,
      overlap_pct = overlap_pct_57lung
    )
    
    # Get peptides detected in both datasets for correlation
    detected_idx_57lung <- which(merged_data[[mskcc_col_57]] > 0 & merged_data[[chop_col_57lung]] > 0)
    
    if(length(detected_idx_57lung) >= 5) {
      sample_data_57lung <- data.frame(
        Sample = rep("57Lung", length(detected_idx_57lung)),
        MSKCC_Count = merged_data[[mskcc_col_57]][detected_idx_57lung],
        CHOP_Count = merged_data[[chop_col_57lung]][detected_idx_57lung],
        Peptide = merged_data$Peptide[detected_idx_57lung],
        stringsAsFactors = FALSE
      )
      corr_data <- rbind(corr_data, sample_data_57lung)
    }
  }
}

# Find matching samples
all_samples <- union(mskcc_sample_ids, chop_sample_ids)
matching_samples <- intersect(mskcc_sample_ids, chop_sample_ids)

# Manually add 57Liver and 57Lung if they exist in CHOP and 57 exists in MSKCC
if("57" %in% mskcc_sample_ids) {
  if("57Liver" %in% chop_sample_ids) {
    matching_samples <- c(matching_samples, "57Liver")
  }
  if("57Lung" %in% chop_sample_ids) {
    matching_samples <- c(matching_samples, "57Lung")
  }
}

print(paste("All samples:", paste(all_samples, collapse=", ")))
print(paste("Matching samples:", paste(matching_samples, collapse=", ")))

# Create correlation data
corr_data <- data.frame()
sample_overlap_stats <- list()

# Create correlation data for intensity
corr_data_intensity <- data.frame()

# Collect sample-by-sample data
for(sample in matching_samples) {
  # Special handling for 57Liver and 57Lung
  if(sample == "57Liver" || sample == "57Lung") {
    mskcc_col <- grep("MSKCC_57_Spectral_Count", colnames(merged_data), value=TRUE)[1]
    # Make sure to find the right CHOP column for Liver or Lung specifically
    if(sample == "57Liver") {
      chop_col <- grep("FL57Liver|H5_FL57", chop_sc_cols, value=TRUE)[1]
    } else { # 57Lung
      chop_col <- grep("FL57Lung|H6_FL57", chop_sc_cols, value=TRUE)[1]
    }
  } else {
    # Normal handling for other samples
    mskcc_col <- grep(paste0("MSKCC_", sample, "_Spectral_Count"), colnames(merged_data), value=TRUE)[1]
    chop_col <- grep(paste0(sample), chop_sc_cols, value=TRUE)[1]
  }
  
  if(is.na(mskcc_col) || is.na(chop_col)) {
    cat("Warning: Could not find columns for sample", sample, "\n")
    next
  }
  
  # Get peptides with intensity values in both datasets for correlation
  if(length(mskcc_int_cols) > 0 && length(chop_int_cols) > 0) {
    if(sample == "57Liver" || sample == "57Lung") {
      mskcc_int_col <- grep("MSKCC_57_Intensity", colnames(merged_data), value=TRUE)[1]
      # Find the right CHOP column for Liver or Lung
      if(sample == "57Liver") {
        chop_int_col <- grep("FL57Liver|H5_FL57.*Intensity", chop_int_cols, value=TRUE)[1]
      } else { # 57Lung
        chop_int_col <- grep("FL57Lung|H6_FL57.*Intensity", chop_int_cols, value=TRUE)[1]
      }
    } else {
      # Normal handling for other samples
      mskcc_int_col <- grep(paste0("MSKCC_", sample, "_Intensity"), colnames(merged_data), value=TRUE)[1]
      chop_int_col <- grep(paste0(sample, ".*Intensity"), chop_int_cols, value=TRUE)[1]
    }
    
    # Skip if columns not found
    if(!is.na(mskcc_int_col) && !is.na(chop_int_col)) {
      detected_idx_int <- which(merged_data[[mskcc_int_col]] > 0 & merged_data[[chop_int_col]] > 0)
      
      if(length(detected_idx_int) >= 5) {
        sample_data_int <- data.frame(
          Sample = rep(sample, length(detected_idx_int)),
          MSKCC_Intensity = merged_data[[mskcc_int_col]][detected_idx_int],
          CHOP_Intensity = merged_data[[chop_int_col]][detected_idx_int],
          Peptide = merged_data$Peptide[detected_idx_int],
          stringsAsFactors = FALSE
        )
        corr_data_intensity <- rbind(corr_data_intensity, sample_data_int)
      }
    }
  }
  
  # Modify this part in your script to focus on unique peptide sequences
  # Get peptides detected in MSKCC and CHOP
  mskcc_peptides <- unique(merged_data$Peptide[merged_data[[mskcc_col]] > 0])
  chop_peptides <- unique(merged_data$Peptide[merged_data[[chop_col]] > 0])
  
  # Find common peptides
  common_peptides <- intersect(mskcc_peptides, chop_peptides)
  
  # Calculate overlap statistics based on unique peptide counts
  mskcc_only <- length(mskcc_peptides) - length(common_peptides)
  chop_only <- length(chop_peptides) - length(common_peptides)
  common_count <- length(common_peptides)
  total_unique <- length(union(mskcc_peptides, chop_peptides))
  overlap_pct <- 100 * common_count / total_unique
  
  # Store overlap statistics
  sample_overlap_stats[[sample]] <- list(
    sample = sample,
    mskcc_count = length(mskcc_peptides),
    chop_count = length(chop_peptides),
    common_count = common_count,
    mskcc_only = mskcc_only,
    chop_only = chop_only,
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
    MSKCC_Only = x$mskcc_only,
    CHOP_Only = x$chop_only,
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

# Calculate correlation statistics by sample for intensity data
sample_stats_intensity <- data.frame()

for(sample in unique(corr_data_intensity$Sample)) {
  sample_data <- corr_data_intensity[corr_data_intensity$Sample == sample, ]
  
  if(nrow(sample_data) >= 5) {
    # Calculate Pearson and Spearman correlations
    pearson_cor <- cor(sample_data$MSKCC_Intensity, sample_data$CHOP_Intensity, 
                       method="pearson", use="pairwise.complete.obs")
    spearman_cor <- cor(sample_data$MSKCC_Intensity, sample_data$CHOP_Intensity, 
                        method="spearman", use="pairwise.complete.obs")
    
    # Calculate R-squared from linear model
    model <- lm(CHOP_Intensity ~ MSKCC_Intensity, data=sample_data)
    r_squared <- summary(model)$r.squared
    
    # Add to stats dataframe
    sample_stats_intensity <- rbind(sample_stats_intensity, data.frame(
      Sample = sample,
      Pearson_Correlation = pearson_cor,
      Spearman_Correlation = spearman_cor,
      R_Squared = r_squared,
      Peptide_Count = nrow(sample_data)
    ))
  }
}

# Order samples consistently
sample_stats_intensity$Sample <- factor(sample_stats_intensity$Sample, 
                                        levels=intersect(sample_order, sample_stats_intensity$Sample))

# Sort samples for consistent visualization
# Custom sort order for better display
sample_order <- c("51", "57", "57Liver", "57Lung", "62", "63", "88", "117", "123")
sample_stats$Sample <- factor(sample_stats$Sample, 
                              levels=intersect(sample_order, sample_stats$Sample))
overlap_df$Sample <- factor(overlap_df$Sample, 
                            levels=intersect(sample_order, overlap_df$Sample))

# Check if there are samples to visualize
if(nrow(sample_stats) > 0 && nrow(overlap_df) > 0) {
  # ===== VISUALIZATION 1: SAMPLE OVERLAP COMPARISON =====
  # Create corrected overlapping bar chart 
  # Prepare data for correct stacking
  overlap_df_stacked <- overlap_df %>%
    mutate(
      MSKCC_Only = MSKCC_Peptides - Common_Peptides,
      CHOP_Only = CHOP_Peptides - Common_Peptides
    ) %>%
    tidyr::pivot_longer(
      cols = c("MSKCC_Only", "CHOP_Only", "Common_Peptides"),
      names_to = "Category",
      values_to = "Count"
    ) %>%
    mutate(
      Category = factor(Category, 
                        levels = c("MSKCC_Only", "CHOP_Only", "Common_Peptides"),
                        labels = c("MSKCC Only", "CHOP Only", "Common"))
    )
  
  # Before creating the plot, make sure to add these calculations to overlap_df
  overlap_df$MSKCC_Only <- overlap_df$MSKCC_Peptides - overlap_df$Common_Peptides
  overlap_df$CHOP_Only <- overlap_df$CHOP_Peptides - overlap_df$Common_Peptides
  
  # Print data to verify
  print(overlap_df[, c("Sample", "MSKCC_Peptides", "CHOP_Peptides", "Common_Peptides", "MSKCC_Only", "CHOP_Only")])
  
  # Create a grouped (not stacked) bar chart that clearly shows each category
  p_overlap <- ggplot(overlap_df) +
    # Position bars side by side, not stacked
    geom_bar(aes(x=as.numeric(Sample)-0.25, y=MSKCC_Only, fill="MSKCC Only"), 
             stat="identity", width=0.2) +
    geom_bar(aes(x=as.numeric(Sample), y=Common_Peptides, fill="Common"), 
             stat="identity", width=0.2) +
    geom_bar(aes(x=as.numeric(Sample)+0.25, y=CHOP_Only, fill="CHOP Only"), 
             stat="identity", width=0.2) +
    
    # Add text labels above each bar
    geom_text(aes(x=as.numeric(Sample)-0.25, y=MSKCC_Only, label=MSKCC_Only),
              vjust=-0.5, size=3) +
    geom_text(aes(x=as.numeric(Sample), y=Common_Peptides, label=Common_Peptides),
              vjust=-0.5, size=3) +
    geom_text(aes(x=as.numeric(Sample)+0.25, y=CHOP_Only, label=CHOP_Only),
              vjust=-0.5, size=3) +
    
    # Use proper color scale
    scale_fill_manual(values=c("MSKCC Only"="forestgreen", 
                               "Common"="#377EB8",
                               "CHOP Only"="#E41A1C"), 
                      name="Dataset") +
    
    # Set x-axis breaks and labels
    scale_x_continuous(breaks=1:length(levels(overlap_df$Sample)), 
                       labels=levels(overlap_df$Sample)) +
    
    # Labels and theme
    labs(title="Peptide Detection by Sample", 
         subtitle="Number of unique peptides detected in each dataset", 
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
  # Create separate correlation statistics plots
  # 1. Pearson Correlation
  p_pearson <- ggplot(sample_stats, aes(x=Sample, y=Pearson_Correlation, fill=Sample)) +
    geom_bar(stat="identity", alpha=0.8) +
    geom_text(aes(label=sprintf("%.2f", Pearson_Correlation)), vjust=-0.5) +
    ylim(0, max(1, max(sample_stats$Pearson_Correlation) * 1.2)) +
    scale_fill_brewer(palette="Set2") +
    labs(title="Pearson Correlation by Sample",
         x="Sample", y="Pearson Correlation") +
    theme_minimal() +
    theme(plot.title = element_text(hjust=0.5, face="bold"),
          legend.position="none")
  
  # 2. Spearman Correlation
  p_spearman <- ggplot(sample_stats, aes(x=Sample, y=Spearman_Correlation, fill=Sample)) +
    geom_bar(stat="identity", alpha=0.8) +
    geom_text(aes(label=sprintf("%.2f", Spearman_Correlation)), vjust=-0.5) +
    ylim(0, max(1, max(sample_stats$Spearman_Correlation) * 1.2)) +
    scale_fill_brewer(palette="Set2") +
    labs(title="Spearman Correlation by Sample",
         x="Sample", y="Spearman Correlation") +
    theme_minimal() +
    theme(plot.title = element_text(hjust=0.5, face="bold"),
          legend.position="none")
  
  # 3. R-squared
  p_rsquared <- ggplot(sample_stats, aes(x=Sample, y=R_Squared, fill=Sample)) +
    geom_bar(stat="identity", alpha=0.8) +
    geom_text(aes(label=sprintf("%.2f", R_Squared)), vjust=-0.5) +
    ylim(0, max(1, max(sample_stats$R_Squared) * 1.2)) +
    scale_fill_brewer(palette="Set2") +
    labs(title="R² Values by Sample",
         x="Sample", y="R-squared") +
    theme_minimal() +
    theme(plot.title = element_text(hjust=0.5, face="bold"),
          legend.position="none")
  
  # ===== VISUALIZATION 3: SAMPLE CORRELATION SCATTER PLOTS =====
  # Create square scatter plots for each sample
  sample_plots <- list()
  
  for(sample in unique(corr_data$Sample)) {
    sample_data <- corr_data[corr_data$Sample == sample, ]
    
    if(nrow(sample_data) >= 5) {
      # Get correlation values
      pearson_cor <- cor(sample_data$MSKCC_Count, sample_data$CHOP_Count, 
                         method="pearson", use="pairwise.complete.obs")
      r_squared <- sample_stats$R_Squared[sample_stats$Sample == sample]
      
      # Create square scatter plot
      p <- ggplot(sample_data, aes(x=MSKCC_Count, y=CHOP_Count)) +
        geom_point(alpha=0.5, color="blue") +
        geom_smooth(method="lm", color="red", se=TRUE) +
        labs(title=paste0("Sample ", sample, 
                          " (R² = ", sprintf("%.3f", r_squared), 
                          ", n = ", nrow(sample_data), ")"),
             x="MSKCC Spectral Count", 
             y="CHOP Spectral Count") +
        theme_minimal() +
        theme(plot.title = element_text(hjust=0.5),
              aspect.ratio = 1)  # Force square aspect ratio
      
      sample_plots[[sample]] <- p
    }
  }
  
  # Create intensity correlation statistics plots
  # 1. Pearson Correlation for intensity
  p_pearson_int <- ggplot(sample_stats_intensity, aes(x=Sample, y=Pearson_Correlation, fill=Sample)) +
    geom_bar(stat="identity", alpha=0.8) +
    geom_text(aes(label=sprintf("%.2f", Pearson_Correlation)), vjust=-0.5) +
    ylim(0, max(1, max(sample_stats_intensity$Pearson_Correlation) * 1.2)) +
    scale_fill_brewer(palette="Set3") +
    labs(title="Pearson Correlation by Sample (Intensity)",
         x="Sample", y="Pearson Correlation") +
    theme_minimal() +
    theme(plot.title = element_text(hjust=0.5, face="bold"),
          legend.position="none")
  
  # 2. Spearman Correlation for intensity
  p_spearman_int <- ggplot(sample_stats_intensity, aes(x=Sample, y=Spearman_Correlation, fill=Sample)) +
    geom_bar(stat="identity", alpha=0.8) +
    geom_text(aes(label=sprintf("%.2f", Spearman_Correlation)), vjust=-0.5) +
    ylim(0, max(1, max(sample_stats_intensity$Spearman_Correlation) * 1.2)) +
    scale_fill_brewer(palette="Set3") +
    labs(title="Spearman Correlation by Sample (Intensity)",
         x="Sample", y="Spearman Correlation") +
    theme_minimal() +
    theme(plot.title = element_text(hjust=0.5, face="bold"),
          legend.position="none")
  
  # 3. R-squared for intensity
  p_rsquared_int <- ggplot(sample_stats_intensity, aes(x=Sample, y=R_Squared, fill=Sample)) +
    geom_bar(stat="identity", alpha=0.8) +
    geom_text(aes(label=sprintf("%.2f", R_Squared)), vjust=-0.5) +
    ylim(0, max(1, max(sample_stats_intensity$R_Squared) * 1.2)) +
    scale_fill_brewer(palette="Set3") +
    labs(title="R² Values by Sample (Intensity)",
         x="Sample", y="R-squared") +
    theme_minimal() +
    theme(plot.title = element_text(hjust=0.5, face="bold"),
          legend.position="none")
  
  # Create scatter plots for each sample's intensity data
  sample_plots_intensity <- list()
  
  for(sample in unique(corr_data_intensity$Sample)) {
    sample_data <- corr_data_intensity[corr_data_intensity$Sample == sample, ]
    
    if(nrow(sample_data) >= 5) {
      # Get correlation values
      pearson_cor <- cor(sample_data$MSKCC_Intensity, sample_data$CHOP_Intensity, 
                         method="pearson", use="pairwise.complete.obs")
      r_squared <- sample_stats_intensity$R_Squared[sample_stats_intensity$Sample == sample]
      
      # Create square scatter plot
      p <- ggplot(sample_data, aes(x=MSKCC_Intensity, y=CHOP_Intensity)) +
        geom_point(alpha=0.5, color="darkgreen") +
        geom_smooth(method="lm", color="red", se=TRUE) +
        labs(title=paste0("Sample ", sample, 
                          " (R² = ", sprintf("%.3f", r_squared), 
                          ", n = ", nrow(sample_data), ")"),
             x="MSKCC Intensity", 
             y="CHOP Intensity") +
        theme_minimal() +
        theme(plot.title = element_text(hjust=0.5),
              aspect.ratio = 1)  # Force square aspect ratio
      
      sample_plots_intensity[[sample]] <- p
    }
  }
  
  # ===== VISUALIZATION 4: SUMMARY STATISTICS =====
  # Create summary statistics
  chop_peptide_col <- ifelse("Peptide" %in% colnames(chop_data), "Peptide", "Peptide Sequence")
  
  # Modify summary statistics to include intensity correlations
  summary_stats <- data.frame(
    Metric = c(
      "Total MSKCC unique peptides",
      "Total CHOP unique peptides",
      "Common peptides between datasets",
      "Overall overlap percentage",
      "Average Spectral Count Pearson correlation",  # Renamed for clarity
      "Average Spectral Count Spearman correlation", # Renamed for clarity
      "Average Spectral Count R-squared value",      # Renamed for clarity
      "Average Intensity Pearson correlation",       # Added
      "Average Intensity Spearman correlation",      # Added
      "Average Intensity R-squared value",           # Added
      "Number of matching samples"
    ),
    Value = c(
      length(unique(mskcc_data$Peptide)),
      length(unique(chop_data[[chop_peptide_col]])),
      length(intersect(unique(mskcc_data$Peptide), unique(chop_data[[chop_peptide_col]]))),
      round(100 * length(intersect(unique(mskcc_data$Peptide), unique(chop_data[[chop_peptide_col]]))) / 
              length(union(unique(mskcc_data$Peptide), unique(chop_data[[chop_peptide_col]]))), 1),
      round(mean(sample_stats$Pearson_Correlation), 3),
      round(mean(sample_stats$Spearman_Correlation), 3),
      round(mean(sample_stats$R_Squared), 3),
      round(mean(sample_stats_intensity$Pearson_Correlation), 3),     # Added
      round(mean(sample_stats_intensity$Spearman_Correlation), 3),    # Added 
      round(mean(sample_stats_intensity$R_Squared), 3),               # Added
      length(unique(matching_samples))
    )
  )
  
  # Format summary statistics table
  summary_table <- tableGrob(summary_stats, rows=NULL, theme=ttheme_minimal(
    core=list(fg_params=list(hjust=0, x=0.1), bg_params=list(fill=NA)),
    colhead=list(fg_params=list(hjust=0.5, fontface="bold"))
  ))
  
  # ===== SAVE ALL VISUALIZATIONS TO PDF =====
  pdf("immunopeptidomics_visualization_report.pdf", width=11, height=8.5)
  
  # Page 1: Sample overlap statistics
  title_grob <- textGrob("Immunopeptidomics Datasets Comparison: Sample Overlap", 
                         gp=gpar(fontsize=16, fontface="bold"))
  
  grid.arrange(
    p_overlap, p_pct_overlap,
    ncol=2, nrow=1,
    top=title_grob
  )
  
  # Page 2: Spectral Count correlation metrics
  title_grob2 <- textGrob("Immunopeptidomics Datasets Comparison: Spectral Count Correlation Metrics", 
                          gp=gpar(fontsize=16, fontface="bold"))
  
  grid.arrange(
    p_pearson, p_spearman, p_rsquared,
    ncol=1, nrow=3,  # Stacked vertically for better readability
    top=title_grob2
  )
  
  # Page 3: Intensity correlation metrics (NEW)
  title_grob_int <- textGrob("Immunopeptidomics Datasets Comparison: Intensity Correlation Metrics", 
                             gp=gpar(fontsize=16, fontface="bold"))
  
  grid.arrange(
    p_pearson_int, p_spearman_int, p_rsquared_int,
    ncol=1, nrow=3,
    top=title_grob_int
  )
  
  # Page 4: Spectral Count scatter plots
  if(length(sample_plots) > 0) {
    sample_plot_list <- lapply(names(sample_plots), function(s) sample_plots[[s]])
    title_grob3 <- textGrob("Sample-by-Sample Correlation (MSKCC vs CHOP)\nn = number of peptides detected in both datasets with spectral count > 0",
                            gp=gpar(fontsize=16, fontface="bold"))
    
    # Arrange plots in a grid with at most 3 plots per row
    plot_grid_ncol <- min(3, length(sample_plot_list))
    
    do.call(grid.arrange, 
            c(sample_plot_list, 
              list(ncol=plot_grid_ncol, 
                   top=title_grob3)
            )
    )
  }
  
  # Page 5: Intensity scatter plots (NEW)
  if(length(sample_plots_intensity) > 0) {
    sample_plot_list_int <- lapply(names(sample_plots_intensity), function(s) sample_plots_intensity[[s]])
    title_grob_int2 <- textGrob("Sample-by-Sample Intensity Correlation (MSKCC vs CHOP)\nn = number of peptides with intensity values > 0 in both datasets",
                                gp=gpar(fontsize=16, fontface="bold"))
    
    # Arrange plots in a grid with at most 3 plots per row
    plot_grid_ncol_int <- min(3, length(sample_plot_list_int))
    
    do.call(grid.arrange, 
            c(sample_plot_list_int, 
              list(ncol=plot_grid_ncol_int, 
                   top=title_grob_int2)
            )
    )
  }
  
  # Page 6: Summary statistics
  title_grob4 <- textGrob("Summary Statistics", 
                          gp=gpar(fontsize=16, fontface="bold"))
  
  # Clean up duplicate intensity statistics (currently appearing twice in your output)
  summary_stats <- summary_stats[!duplicated(summary_stats$Metric),]
  
  grid.arrange(
    tableGrob(summary_stats, rows=NULL, theme=ttheme_minimal(
      core=list(fg_params=list(hjust=0, x=0.1), bg_params=list(fill=NA)),
      colhead=list(fg_params=list(hjust=0.5, fontface="bold"))
    )),
    top=title_grob4
  )
  
  # Close the PDF file only after adding all content
  dev.off()
  
  # Save the main figures as individual PNGs for easy viewing
  ggsave("sample_overlap.png", p_overlap, width=10, height=6)
  ggsave("overlap_percentage.png", p_pct_overlap, width=10, height=6)
  ggsave("pearson_correlation.png", p_pearson, width=10, height=6)
  ggsave("spearman_correlation.png", p_spearman, width=10, height=6)
  ggsave("r_squared.png", p_rsquared, width=10, height=6)
  ggsave("pearson_correlation_intensity.png", p_pearson_int, width=10, height=6)
  ggsave("spearman_correlation_intensity.png", p_spearman_int, width=10, height=6)
  ggsave("r_squared_intensity.png", p_rsquared_int, width=10, height=6)
  
  # Print completion message
  cat("Visualizations created and saved to:\n")
  cat("1. immunopeptidomics_visualization_report.pdf (comprehensive report)\n")
  cat("2. sample_overlap.png\n")
  cat("3. overlap_percentage.png\n")
  cat("4. pearson_correlation.png\n")
  cat("5. spearman_correlation.png\n")
  cat("6. r_squared.png\n")
  cat("7. pearson_correlation_intensity.png\n")
  cat("8. spearman_correlation_intensity.png\n")
  cat("9. r_squared_intensity.png\n")
  
  # Print summary statistics
  cat("\nSummary Statistics:\n")
  print(summary_stats)
  
  # Print correlation statistics by sample
  cat("\nCorrelation Statistics by Sample:\n")
  print(sample_stats)
  
  # Print overlap statistics by sample
  print("Overlap data frame:")
  print(overlap_df[, c("Sample", "MSKCC_Peptides", "CHOP_Peptides", "Common_Peptides", "MSKCC_Only", "CHOP_Only")])
  
  # Print intensity correlation stats
  print("Intensity correlation stats:")
  print(sample_stats_intensity)
  
} else {
  cat("\nWARNING: No matching samples found with sufficient data for correlation analysis.\n")
}

# Optional: Add annotation analysis
if(file.exists("peptide_annotations_all.tsv")) {
  cat("\n===== ANNOTATION ANALYSIS =====\n")
  # Load the annotations file
  annotations_file <- "peptide_annotations_all.tsv"
  annotation_data <- fread(annotations_file, sep="\t", header=TRUE, fill=TRUE)
  
  # Count peptides with multiple annotations
  multi_gene_count <- sum(grepl(";", annotation_data$All_Genes, fixed=TRUE), na.rm=TRUE)
  multi_protein_count <- sum(grepl(";", annotation_data$All_Proteins, fixed=TRUE), na.rm=TRUE)
  
  cat("Peptides with multiple gene annotations:", multi_gene_count, "\n")
  cat("Peptides with multiple protein annotations:", multi_protein_count, "\n")
  
  # Count peptides found in both sources
  both_sources <- sum(grepl("MSKCC.*CHOP|CHOP.*MSKCC", annotation_data$Sources, ignore.case=TRUE), na.rm=TRUE)
  cat("Peptides detected in both MSKCC and CHOP:", both_sources, "\n")
}

# # This script correlates CHOP data with our immunopeptidomics data that we ran at MSKCC
# # First run the script titled "generate_final_CHOP.R"
# # Take the combined_peptides.tsv unmodified from the combine2c3v script
# # and rename is to MSKCC_combined_peptides.tsv
# 
# library(data.table)
# library(ggplot2)
# library(gridExtra)
# library(grid)  # Add this for textGrob
# library(RColorBrewer)
# library(reshape2)
# 
# # Set working directory - update if needed
# setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation")
# 
# # Define input files
# merged_file <- "merged_immunopeptidomics_data.tsv"
# chop_file <- "final_CHOP_combined_peptide.tsv"
# mskcc_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation/unique_peptides_unmodified.tsv"
# 
# # Read files
# cat("Reading merged data file...\n")
# merged_data <- fread(merged_file, sep="\t", header=TRUE, fill=TRUE)
# cat("Reading CHOP data file...\n")
# chop_data <- fread(chop_file, sep="\t", header=TRUE, fill=TRUE)
# cat("Reading MSKCC data file...\n")
# mskcc_data <- fread(mskcc_file, sep="\t", header=TRUE, fill=TRUE)
# 
# # Ensure consistent peptide column names
# if("Peptide Sequence" %in% colnames(mskcc_data) && !"Peptide" %in% colnames(mskcc_data)) {
#   colnames(mskcc_data)[colnames(mskcc_data) == "Peptide Sequence"] <- "Peptide"
# }
# 
# if("Peptide Sequence" %in% colnames(chop_data) && !"Peptide" %in% colnames(chop_data)) {
#   colnames(chop_data)[colnames(chop_data) == "Peptide Sequence"] <- "Peptide"
# }
# 
# # Extract spectral count columns
# mskcc_sc_cols <- grep("MSKCC_.*_Spectral_Count", colnames(merged_data), value=TRUE)
# chop_sc_cols <- grep("FL.*Spectral Count", colnames(merged_data), value=TRUE)
# 
# # Extract sample IDs
# extract_mskcc_id <- function(col_name) {
#   parts <- strsplit(col_name, "_")[[1]]
#   if(length(parts) >= 3) {
#     return(parts[2])
#   } else {
#     return(NA)
#   }
# }
# 
# # Modify the extraction function to handle the lung/liver variants
# extract_chop_id <- function(col_name) {
#   # First check for the specific liver/lung cases
#   if(grepl("FL57Liver", col_name) || grepl("H5_FL57", col_name)) {
#     return("57Liver")
#   } else if(grepl("FL57Lung", col_name) || grepl("H6_FL57", col_name)) {
#     return("57Lung")
#   }
# 
#   # Regular extraction for other samples
#   # More robust regex to handle different formats like "FL51 Spectral Count", "H1_FL62 Spectral Count"
#   match <- regexpr("FL([0-9]+)", col_name)
#   if(match > 0) {
#     fl_part <- regmatches(col_name, match)
#     return(gsub("FL", "", fl_part))
#   } else {
#     # Try to extract sample IDs in the format H1_FL62
#     match2 <- regexpr("H[0-9]+_FL([0-9]+)", col_name)
#     if(match2 > 0) {
#       full_match <- regmatches(col_name, match2)
#       # Extract just the number after FL
#       fl_num <- gsub(".*FL", "", full_match)
#       return(fl_num)
#     }
#     return(NA)
#   }
# }
# 
# # And then add this special mapping for MSKCC sample 57
# # This will create duplicate entries for sample 57 to match with both liver and lung
# mskcc_sample_ids <- sapply(mskcc_sc_cols, extract_mskcc_id)
# # Duplicate the MSKCC 57 sample to match both variants
# if("57" %in% mskcc_sample_ids) {
#   mskcc_sample_ids_expanded <- c(mskcc_sample_ids, "57")
#   names(mskcc_sample_ids_expanded)[length(mskcc_sample_ids_expanded)] <- 
#     names(mskcc_sample_ids)[which(mskcc_sample_ids == "57")]
#   
#   # Replace 57 with 57Lung and 57Liver
#   mskcc_sample_ids_expanded[mskcc_sample_ids_expanded == "57"] <- 
#     c("57Lung", "57Liver")
#   
#   # Use the expanded list for matching
#   mskcc_sample_ids <- mskcc_sample_ids_expanded
# }
# 
# mskcc_sample_ids <- sapply(mskcc_sc_cols, extract_mskcc_id)
# chop_sample_ids <- sapply(chop_sc_cols, extract_chop_id)
# 
# # Find matching samples
# matching_samples <- intersect(mskcc_sample_ids, chop_sample_ids)
# print(paste("Matching samples:", paste(matching_samples, collapse=", ")))
# 
# # Create correlation data
# corr_data <- data.frame()
# sample_overlap_stats <- list()
# 
# # Collect sample-by-sample data
# for(sample in matching_samples) {
#   mskcc_col <- grep(paste0("MSKCC_", sample, "_Spectral_Count"), colnames(merged_data), value=TRUE)[1]
#   chop_col <- grep(paste0("FL", sample), chop_sc_cols, value=TRUE)[1]
#   
#   if(is.na(mskcc_col) || is.na(chop_col)) next
#   
#   # Get peptides detected in MSKCC and CHOP
#   mskcc_peptides <- merged_data$Peptide[merged_data[[mskcc_col]] > 0]
#   chop_peptides <- merged_data$Peptide[merged_data[[chop_col]] > 0]
#   
#   # Find common peptides
#   common_peptides <- intersect(mskcc_peptides, chop_peptides)
#   
#   # Calculate overlap statistics
#   total_unique <- length(union(mskcc_peptides, chop_peptides))
#   overlap_pct <- 100 * length(common_peptides) / total_unique
#   
#   # Store overlap statistics
#   sample_overlap_stats[[sample]] <- list(
#     sample = sample,
#     mskcc_count = length(mskcc_peptides),
#     chop_count = length(chop_peptides),
#     common_count = length(common_peptides),
#     total_unique = total_unique,
#     overlap_pct = overlap_pct
#   )
#   
#   # Get peptides detected in both datasets for correlation
#   detected_idx <- which(merged_data[[mskcc_col]] > 0 & merged_data[[chop_col]] > 0)
#   
#   if(length(detected_idx) >= 5) {
#     sample_data <- data.frame(
#       Sample = rep(sample, length(detected_idx)),
#       MSKCC_Count = merged_data[[mskcc_col]][detected_idx],
#       CHOP_Count = merged_data[[chop_col]][detected_idx],
#       Peptide = merged_data$Peptide[detected_idx],
#       stringsAsFactors = FALSE
#     )
#     corr_data <- rbind(corr_data, sample_data)
#   }
# }
# 
# # Create overlap summary dataframe
# overlap_df <- do.call(rbind, lapply(sample_overlap_stats, function(x) {
#   data.frame(
#     Sample = x$sample,
#     MSKCC_Peptides = x$mskcc_count,
#     CHOP_Peptides = x$chop_count,
#     Common_Peptides = x$common_count,
#     Total_Unique = x$total_unique,
#     Overlap_Percent = x$overlap_pct
#   )
# }))
# 
# # Calculate correlation statistics by sample
# sample_stats <- data.frame()
# 
# for(sample in unique(corr_data$Sample)) {
#   sample_data <- corr_data[corr_data$Sample == sample, ]
#   
#   if(nrow(sample_data) >= 5) {
#     # Calculate Pearson and Spearman correlations
#     pearson_cor <- cor(sample_data$MSKCC_Count, sample_data$CHOP_Count, 
#                        method="pearson", use="pairwise.complete.obs")
#     spearman_cor <- cor(sample_data$MSKCC_Count, sample_data$CHOP_Count, 
#                         method="spearman", use="pairwise.complete.obs")
#     
#     # Calculate R-squared from linear model
#     model <- lm(CHOP_Count ~ MSKCC_Count, data=sample_data)
#     r_squared <- summary(model)$r.squared
#     
#     # Add to stats dataframe
#     sample_stats <- rbind(sample_stats, data.frame(
#       Sample = sample,
#       Pearson_Correlation = pearson_cor,
#       Spearman_Correlation = spearman_cor,
#       R_Squared = r_squared,
#       Peptide_Count = nrow(sample_data)
#     ))
#   }
# }
# 
# # Sort samples for consistent visualization
# sample_stats$Sample <- factor(sample_stats$Sample, 
#                               levels=sample_stats$Sample[order(as.numeric(sample_stats$Sample))])
# overlap_df$Sample <- factor(overlap_df$Sample, 
#                             levels=overlap_df$Sample[order(as.numeric(overlap_df$Sample))])
# 
# # Check if there are samples to visualize
# if(nrow(sample_stats) > 0 && nrow(overlap_df) > 0) {
#   # ===== VISUALIZATION 1: SAMPLE OVERLAP COMPARISON =====
# # Create overlapping bar chart with corrected color mapping
#   # Create overlapping bar chart with corrected color mapping
#   p_overlap <- ggplot(overlap_df) +
#     # First plot MSKCC peptides
#     geom_bar(aes(x=Sample, y=MSKCC_Peptides - Common_Peptides, fill="MSKCC Only"), 
#              stat="identity", alpha=0.7) +
#     # Then plot CHOP peptides
#     geom_bar(aes(x=Sample, y=CHOP_Peptides - Common_Peptides, fill="CHOP Only"), 
#              stat="identity", alpha=0.7) +
#     # Finally plot common peptides
#     geom_bar(aes(x=Sample, y=Common_Peptides, fill="Common"), 
#              stat="identity") +
#     # Set colors explicitly 
#     scale_fill_manual(values=c("MSKCC Only"="forestgreen", "CHOP Only"="#E41A1C", "Common"="#377EB8"), 
#                       name="Dataset") +
#     # Add text labels
#     geom_text(aes(x=Sample, y=MSKCC_Peptides, label=MSKCC_Peptides), 
#               vjust=-0.5, color="black", size=3) +
#     geom_text(aes(x=Sample, y=CHOP_Peptides, label=CHOP_Peptides), 
#               vjust=-0.5, color="black", size=3) +
#     geom_text(aes(x=Sample, y=Common_Peptides/2, label=Common_Peptides), 
#               color="white", size=3) +
#     labs(title="Peptide Detection by Sample", 
#          subtitle="Number of peptides detected in each dataset", 
#          x="Sample", y="Number of Peptides") +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust=0.5, face="bold"),
#           plot.subtitle = element_text(hjust=0.5))
#   
#   # Create percentage overlap plot
#   p_pct_overlap <- ggplot(overlap_df, aes(x=Sample, y=Overlap_Percent, fill=Sample)) +
#     geom_bar(stat="identity", alpha=0.8) +
#     geom_text(aes(label=sprintf("%.1f%%", Overlap_Percent)), vjust=-0.5) +
#     ylim(0, max(30, max(overlap_df$Overlap_Percent) * 1.1)) +
#     scale_fill_brewer(palette="Set2") +
#     labs(title="Peptide Overlap Percentage by Sample",
#          x="Sample", y="Overlap Percentage (%)") +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust=0.5, face="bold"),
#           legend.position="none")
#   
#   # ===== VISUALIZATION 2: CORRELATION STATISTICS =====
#   # Create separate correlation statistics plots
#   # 1. Pearson Correlation
#   p_pearson <- ggplot(sample_stats, aes(x=Sample, y=Pearson_Correlation, fill=Sample)) +
#     geom_bar(stat="identity", alpha=0.8) +
#     geom_text(aes(label=sprintf("%.2f", Pearson_Correlation)), vjust=-0.5) +
#     ylim(0, max(1, max(sample_stats$Pearson_Correlation) * 1.2)) +
#     scale_fill_brewer(palette="Set2") +
#     labs(title="Pearson Correlation by Sample",
#          x="Sample", y="Pearson Correlation") +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust=0.5, face="bold"),
#           legend.position="none")
#   
#   # 2. Spearman Correlation
#   p_spearman <- ggplot(sample_stats, aes(x=Sample, y=Spearman_Correlation, fill=Sample)) +
#     geom_bar(stat="identity", alpha=0.8) +
#     geom_text(aes(label=sprintf("%.2f", Spearman_Correlation)), vjust=-0.5) +
#     ylim(0, max(1, max(sample_stats$Spearman_Correlation) * 1.2)) +
#     scale_fill_brewer(palette="Set2") +
#     labs(title="Spearman Correlation by Sample",
#          x="Sample", y="Spearman Correlation") +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust=0.5, face="bold"),
#           legend.position="none")
#   
#   # 3. R-squared
#   p_rsquared <- ggplot(sample_stats, aes(x=Sample, y=R_Squared, fill=Sample)) +
#     geom_bar(stat="identity", alpha=0.8) +
#     geom_text(aes(label=sprintf("%.2f", R_Squared)), vjust=-0.5) +
#     ylim(0, max(1, max(sample_stats$R_Squared) * 1.2)) +
#     scale_fill_brewer(palette="Set2") +
#     labs(title="R² Values by Sample",
#          x="Sample", y="R-squared") +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust=0.5, face="bold"),
#           legend.position="none")
#   
#   # ===== VISUALIZATION 3: SAMPLE CORRELATION SCATTER PLOTS =====
#   # Create scatter plots for each sample
#   sample_plots <- list()
#   
#   for(sample in unique(corr_data$Sample)) {
#     sample_data <- corr_data[corr_data$Sample == sample, ]
#     
#     if(nrow(sample_data) >= 5) {
#       # Get correlation values
#       pearson_cor <- cor(sample_data$MSKCC_Count, sample_data$CHOP_Count, 
#                          method="pearson", use="pairwise.complete.obs")
#       r_squared <- sample_stats$R_Squared[sample_stats$Sample == sample]
#       
#       # Create scatter plot
#       p <- ggplot(sample_data, aes(x=MSKCC_Count, y=CHOP_Count)) +
#         geom_point(alpha=0.5, color="blue") +
#         geom_smooth(method="lm", color="red", se=TRUE) +
#         labs(title=paste0("Sample ", sample, 
#                           " (R² = ", sprintf("%.3f", r_squared), 
#                           ", n = ", nrow(sample_data), ")"),
#              x="MSKCC Spectral Count", 
#              y="CHOP Spectral Count") +
#         theme_minimal() +
#         theme(plot.title = element_text(hjust=0.5),
#       aspect.ratio = 1)  # Force square aspect ratio
#       
#       sample_plots[[sample]] <- p
#     }
#   }
#   
#   # ===== VISUALIZATION 4: SUMMARY STATISTICS =====
#   # Create summary statistics
#   chop_peptide_col <- ifelse("Peptide" %in% colnames(chop_data), "Peptide", "Peptide Sequence")
#   
#   summary_stats <- data.frame(
#     Metric = c(
#       "Total MSKCC unique peptides",
#       "Total CHOP unique peptides",
#       "Common peptides between datasets",
#       "Overall overlap percentage",
#       "Average Pearson correlation",
#       "Average Spearman correlation",
#       "Average R-squared value",
#       "Number of matching samples"
#     ),
#     Value = c(
#       length(unique(mskcc_data$Peptide)),
#       length(unique(chop_data[[chop_peptide_col]])),
#       length(intersect(unique(mskcc_data$Peptide), unique(chop_data[[chop_peptide_col]]))),
#       round(100 * length(intersect(unique(mskcc_data$Peptide), unique(chop_data[[chop_peptide_col]]))) / 
#               length(union(unique(mskcc_data$Peptide), unique(chop_data[[chop_peptide_col]]))), 1),
#       round(mean(sample_stats$Pearson_Correlation), 3),
#       round(mean(sample_stats$Spearman_Correlation), 3),
#       round(mean(sample_stats$R_Squared), 3),
#       length(matching_samples)
#     )
#   )
#   
#   # Format summary statistics table
#   summary_table <- tableGrob(summary_stats, rows=NULL, theme=ttheme_minimal(
#     core=list(fg_params=list(hjust=0, x=0.1), bg_params=list(fill=NA)),
#     colhead=list(fg_params=list(hjust=0.5, fontface="bold"))
#   ))
#   
#   # ===== SAVE ALL VISUALIZATIONS TO PDF =====
#   pdf("immunopeptidomics_visualization_report.pdf", width=11, height=8.5)
#   
#   # Page 1: Sample overlap statistics
#   title_grob <- textGrob("Immunopeptidomics Datasets Comparison: Sample Overlap", 
#                          gp=gpar(fontsize=16, fontface="bold"))
#   
#   grid.arrange(
#     p_overlap, p_pct_overlap,
#     ncol=2, nrow=1,
#     top=title_grob
#   )
#   
#   # New Page: Correlation metrics as three separate bar graphs
#   title_grob2 <- textGrob("Immunopeptidomics Datasets Comparison: Correlation Metrics", 
#                           gp=gpar(fontsize=16, fontface="bold"))
#   
#   grid.arrange(
#     p_pearson, p_spearman, p_rsquared,
#     ncol=1, nrow=3,  # Stacked vertically for better readability
#     top=title_grob2
#   )
#   
#   # Page 2: Correlation scatter plots
#   if(length(sample_plots) > 0) {
#     sample_plot_list <- lapply(names(sample_plots), function(s) sample_plots[[s]])
#     title_grob3 <- textGrob("Sample-by-Sample Correlation (MSKCC vs CHOP)\nn = number of peptides detected in both datasets with spectral count > 0",
#                             gp=gpar(fontsize=16, fontface="bold"))
#     
#     do.call(grid.arrange, 
#             c(sample_plot_list, 
#               list(ncol=min(2, length(sample_plot_list)), 
#                    top=title_grob3)
#             )
#     )
#   }
#   
#   # Page 3: Summary statistics
#   title_grob4 <- textGrob("Summary Statistics", 
#                           gp=gpar(fontsize=16, fontface="bold"))
#   
#   grid.arrange(
#     summary_table,
#     top=title_grob4
#   )
#   
#   dev.off()
#   
#   # Save the main figures as individual PNGs for easy viewing
#   ggsave("sample_overlap.png", p_overlap, width=8, height=6)
#   ggsave("overlap_percentage.png", p_pct_overlap, width=8, height=6)
#   ggsave("pearson_correlation.png", p_pearson, width=8, height=6)
#   ggsave("spearman_correlation.png", p_spearman, width=8, height=6)
#   ggsave("r_squared.png", p_rsquared, width=8, height=6)
#   
#   # Print completion message
#   cat("Visualizations created and saved to:\n")
#   cat("1. immunopeptidomics_visualization_report.pdf (comprehensive report)\n")
#   cat("2. sample_overlap.png\n")
#   cat("3. overlap_percentage.png\n")
#   cat("4. pearson_correlation.png\n")
#   cat("5. spearman_correlation.png\n")
#   cat("6. r_squared.png\n")
#   
#   # Print summary statistics
#   cat("\nSummary Statistics:\n")
#   print(summary_stats)
#   
#   # Print correlation statistics by sample
#   cat("\nCorrelation Statistics by Sample:\n")
#   print(sample_stats[order(as.numeric(as.character(sample_stats$Sample))), ])
#   
#   # Print overlap statistics by sample
#   cat("\nOverlap Statistics by Sample:\n")
#   print(overlap_df[order(as.numeric(as.character(overlap_df$Sample))), ])
#   
# } else {
#   cat("\nWARNING: No matching samples found with sufficient data for correlation analysis.\n")
# }
# 
# # Optional: Add annotation analysis
# if(file.exists("peptide_annotations_all.tsv")) {
#   cat("\n===== ANNOTATION ANALYSIS =====\n")
#   # Load the annotations file
#   annotations_file <- "peptide_annotations_all.tsv"
#   annotation_data <- fread(annotations_file, sep="\t", header=TRUE, fill=TRUE)
#   
#   # Count peptides with multiple annotations
#   multi_gene_count <- sum(grepl(";", annotation_data$All_Genes, fixed=TRUE), na.rm=TRUE)
#   multi_protein_count <- sum(grepl(";", annotation_data$All_Proteins, fixed=TRUE), na.rm=TRUE)
#   
#   cat("Peptides with multiple gene annotations:", multi_gene_count, "\n")
#   cat("Peptides with multiple protein annotations:", multi_protein_count, "\n")
#   
#   # Count peptides found in both sources
#   both_sources <- sum(grepl("MSKCC.*CHOP|CHOP.*MSKCC", annotation_data$Sources, ignore.case=TRUE), na.rm=TRUE)
#   cat("Peptides detected in both MSKCC and CHOP:", both_sources, "\n")
# }
# # # This script correlates CHOP data with our immunopeptidomics data that we ran at MSKCC
# # # First run the script titled "generate_final_CHOP.R"
# # # Uses the unique_peptides_unmodified.tsv from the combine2c3v script
# # 
# # library(data.table)
# # library(ggplot2)
# # library(grid)
# # library(gridExtra)
# # library(RColorBrewer)
# # library(reshape2)
# # 
# # # Set working directory - update if needed
# # setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation")
# # 
# # # Define input files
# # merged_file <- "merged_immunopeptidomics_data.tsv"
# # chop_file <- "final_CHOP_combined_peptide.tsv"
# # mskcc_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation/unique_peptides_unmodified.tsv"
# # 
# # # Read files
# # merged_data <- fread(merged_file, sep="\t", header=TRUE, fill=TRUE)
# # chop_data <- fread(chop_file, sep="\t", header=TRUE, fill=TRUE)
# # mskcc_data <- fread(mskcc_file, sep="\t", header=TRUE, fill=TRUE)
# # 
# # # After reading files, ensure consistent peptide column names
# # if("Peptide Sequence" %in% colnames(mskcc_data) && !"Peptide" %in% colnames(mskcc_data)) {
# #   colnames(mskcc_data)[colnames(mskcc_data) == "Peptide Sequence"] <- "Peptide"
# # }
# # 
# # if("Peptide Sequence" %in% colnames(chop_data) && !"Peptide" %in% colnames(chop_data)) {
# #   colnames(chop_data)[colnames(chop_data) == "Peptide Sequence"] <- "Peptide"
# # }
# # 
# # # Extract spectral count columns
# # mskcc_sc_cols <- grep("MSKCC_.*_Spectral_Count", colnames(merged_data), value=TRUE)
# # chop_sc_cols <- grep("FL.*Spectral Count", colnames(merged_data), value=TRUE)
# # 
# # # Extract sample IDs
# # extract_mskcc_id <- function(col_name) {
# #   parts <- strsplit(col_name, "_")[[1]]
# #   if(length(parts) >= 3) {
# #     return(parts[2])
# #   } else {
# #     return(NA)
# #   }
# # }
# # 
# # extract_chop_id <- function(col_name) {
# #   # More robust regex to handle different formats like "FL51 Spectral Count", "H1_FL62 Spectral Count"
# #   match <- regexpr("FL([0-9]+)", col_name)
# #   if(match > 0) {
# #     fl_part <- regmatches(col_name, match)
# #     return(gsub("FL", "", fl_part))
# #   } else {
# #     # Try to extract sample IDs in the format H1_FL62
# #     match2 <- regexpr("H[0-9]+_FL([0-9]+)", col_name)
# #     if(match2 > 0) {
# #       full_match <- regmatches(col_name, match2)
# #       # Extract just the number after FL
# #       fl_num <- gsub(".*FL", "", full_match)
# #       return(fl_num)
# #     }
# #     return(NA)
# #   }
# # }
# # 
# # mskcc_sample_ids <- sapply(mskcc_sc_cols, extract_mskcc_id)
# # chop_sample_ids <- sapply(chop_sc_cols, extract_chop_id)
# # 
# # # Find matching samples
# # matching_samples <- intersect(mskcc_sample_ids, chop_sample_ids)
# # print(paste("Matching samples:", paste(matching_samples, collapse=", ")))
# # 
# # # Create correlation data
# # corr_data <- data.frame()
# # sample_overlap_stats <- list()
# # 
# # # Collect sample-by-sample data
# # for(sample in matching_samples) {
# #   mskcc_col <- grep(paste0("MSKCC_", sample, "_Spectral_Count"), colnames(merged_data), value=TRUE)[1]
# #   chop_col <- grep(paste0("FL", sample), chop_sc_cols, value=TRUE)[1]
# #   
# #   if(is.na(mskcc_col) || is.na(chop_col)) next
# #   
# #   # Get peptides detected in MSKCC and CHOP
# #   mskcc_peptides <- merged_data$Peptide[merged_data[[mskcc_col]] > 0]
# #   chop_peptides <- merged_data$Peptide[merged_data[[chop_col]] > 0]
# #   
# #   # Find common peptides
# #   common_peptides <- intersect(mskcc_peptides, chop_peptides)
# #   
# #   # Calculate overlap statistics
# #   total_unique <- length(union(mskcc_peptides, chop_peptides))
# #   overlap_pct <- 100 * length(common_peptides) / total_unique
# #   
# #   # Store overlap statistics
# #   sample_overlap_stats[[sample]] <- list(
# #     sample = sample,
# #     mskcc_count = length(mskcc_peptides),
# #     chop_count = length(chop_peptides),
# #     common_count = length(common_peptides),
# #     total_unique = total_unique,
# #     overlap_pct = overlap_pct
# #   )
# #   
# #   # Get peptides detected in both datasets for correlation
# #   detected_idx <- which(merged_data[[mskcc_col]] > 0 & merged_data[[chop_col]] > 0)
# #   
# #   if(length(detected_idx) >= 5) {
# #     sample_data <- data.frame(
# #       Sample = rep(sample, length(detected_idx)),
# #       MSKCC_Count = merged_data[[mskcc_col]][detected_idx],
# #       CHOP_Count = merged_data[[chop_col]][detected_idx],
# #       Peptide = merged_data$Peptide[detected_idx],
# #       stringsAsFactors = FALSE
# #     )
# #     corr_data <- rbind(corr_data, sample_data)
# #   }
# # }
# # 
# # # Create overlap summary dataframe
# # overlap_df <- do.call(rbind, lapply(sample_overlap_stats, function(x) {
# #   data.frame(
# #     Sample = x$sample,
# #     MSKCC_Peptides = x$mskcc_count,
# #     CHOP_Peptides = x$chop_count,
# #     Common_Peptides = x$common_count,
# #     Total_Unique = x$total_unique,
# #     Overlap_Percent = x$overlap_pct
# #   )
# # }))
# # 
# # # Calculate correlation statistics by sample
# # sample_stats <- data.frame()
# # 
# # for(sample in unique(corr_data$Sample)) {
# #   sample_data <- corr_data[corr_data$Sample == sample, ]
# #   
# #   if(nrow(sample_data) >= 5) {
# #     # Calculate Pearson and Spearman correlations
# #     pearson_cor <- cor(sample_data$MSKCC_Count, sample_data$CHOP_Count, 
# #                        method="pearson", use="pairwise.complete.obs")
# #     spearman_cor <- cor(sample_data$MSKCC_Count, sample_data$CHOP_Count, 
# #                         method="spearman", use="pairwise.complete.obs")
# #     
# #     # Calculate R-squared from linear model
# #     model <- lm(CHOP_Count ~ MSKCC_Count, data=sample_data)
# #     r_squared <- summary(model)$r.squared
# #     
# #     # Add to stats dataframe
# #     sample_stats <- rbind(sample_stats, data.frame(
# #       Sample = sample,
# #       Pearson_Correlation = pearson_cor,
# #       Spearman_Correlation = spearman_cor,
# #       R_Squared = r_squared,
# #       Peptide_Count = nrow(sample_data)
# #     ))
# #   }
# # }
# # 
# # # Sort samples for consistent visualization
# # sample_stats$Sample <- factor(sample_stats$Sample, 
# #                               levels=sample_stats$Sample[order(as.numeric(sample_stats$Sample))])
# # overlap_df$Sample <- factor(overlap_df$Sample, 
# #                             levels=overlap_df$Sample[order(as.numeric(overlap_df$Sample))])
# # 
# # # ===== VISUALIZATION 1: SAMPLE OVERLAP COMPARISON =====
# # # Create overlapping bar chart
# # p_overlap <- ggplot(overlap_df, aes(x=Sample)) +
# #   geom_bar(aes(y=MSKCC_Peptides, fill="MSKCC Only"), stat="identity", alpha=0.7) +
# #   geom_bar(aes(y=CHOP_Peptides, fill="CHOP Only"), stat="identity", alpha=0.7) +
# #   geom_bar(aes(y=Common_Peptides, fill="Common"), stat="identity") +
# #   scale_fill_manual(values=c("MSKCC Only"="forestgreen", "CHOP Only"="#E41A1C", "Common"="#377EB8"), name="Dataset") +
# #   geom_text(aes(y=MSKCC_Peptides, label=MSKCC_Peptides), 
# #             position=position_stack(vjust=0.5), color="white", size=3.5) +
# #   geom_text(aes(y=CHOP_Peptides, label=CHOP_Peptides), 
# #             position=position_stack(vjust=0.5), color="white", size=3.5) +
# #   geom_text(aes(y=Common_Peptides, label=Common_Peptides), 
# #             position=position_stack(vjust=0.5), color="white", size=3.5) +
# #   labs(title="Peptide Detection by Sample", 
# #        subtitle="Number of peptides detected in each dataset", 
# #        x="Sample", y="Number of Peptides") +
# #   theme_minimal() +
# #   theme(plot.title = element_text(hjust=0.5, face="bold"),
# #         plot.subtitle = element_text(hjust=0.5))
# # 
# # # Create percentage overlap plot
# # p_pct_overlap <- ggplot(overlap_df, aes(x=Sample, y=Overlap_Percent, fill=Sample)) +
# #   geom_bar(stat="identity", alpha=0.8) +
# #   geom_text(aes(label=sprintf("%.1f%%", Overlap_Percent)), vjust=-0.5) +
# #   ylim(0, max(30, max(overlap_df$Overlap_Percent) * 1.1)) +
# #   scale_fill_brewer(palette="Set2") +
# #   labs(title="Peptide Overlap Percentage by Sample",
# #        x="Sample", y="Overlap Percentage (%)") +
# #   theme_minimal() +
# #   theme(plot.title = element_text(hjust=0.5, face="bold"),
# #         legend.position="none")
# # 
# # # ===== VISUALIZATION 2: CORRELATION STATISTICS =====
# # # Create correlation statistics plot
# # # Create separate correlation statistics plots
# # # 1. Pearson Correlation
# # p_pearson <- ggplot(sample_stats, aes(x=Sample, y=Pearson_Correlation, fill=Sample)) +
# #   geom_bar(stat="identity", alpha=0.8) +
# #   geom_text(aes(label=sprintf("%.2f", Pearson_Correlation)), vjust=-0.5) +
# #   ylim(0, max(1, max(sample_stats$Pearson_Correlation) * 1.2)) +
# #   scale_fill_brewer(palette="Set2") +
# #   labs(title="Pearson Correlation by Sample",
# #        x="Sample", y="Pearson Correlation") +
# #   theme_minimal() +
# #   theme(plot.title = element_text(hjust=0.5, face="bold"),
# #         legend.position="none")
# # 
# # # 2. Spearman Correlation
# # p_spearman <- ggplot(sample_stats, aes(x=Sample, y=Spearman_Correlation, fill=Sample)) +
# #   geom_bar(stat="identity", alpha=0.8) +
# #   geom_text(aes(label=sprintf("%.2f", Spearman_Correlation)), vjust=-0.5) +
# #   ylim(0, max(1, max(sample_stats$Spearman_Correlation) * 1.2)) +
# #   scale_fill_brewer(palette="Set2") +
# #   labs(title="Spearman Correlation by Sample",
# #        x="Sample", y="Spearman Correlation") +
# #   theme_minimal() +
# #   theme(plot.title = element_text(hjust=0.5, face="bold"),
# #         legend.position="none")
# # 
# # # 3. R-squared
# # p_rsquared <- ggplot(sample_stats, aes(x=Sample, y=R_Squared, fill=Sample)) +
# #   geom_bar(stat="identity", alpha=0.8) +
# #   geom_text(aes(label=sprintf("%.2f", R_Squared)), vjust=-0.5) +
# #   ylim(0, max(1, max(sample_stats$R_Squared) * 1.2)) +
# #   scale_fill_brewer(palette="Set2") +
# #   labs(title="R² Values by Sample",
# #        x="Sample", y="R-squared") +
# #   theme_minimal() +
# #   theme(plot.title = element_text(hjust=0.5, face="bold"),
# #         legend.position="none")
# # 
# # # p_correlation <- ggplot(sample_stats) +
# # #   geom_bar(aes(x=Sample, y=Pearson_Correlation, fill="Pearson"), 
# # #            stat="identity", position=position_dodge(), alpha=0.7) +
# # #   geom_bar(aes(x=Sample, y=Spearman_Correlation, fill="Spearman"), 
# # #            stat="identity", position=position_dodge(), alpha=0.7) +
# # #   geom_bar(aes(x=Sample, y=R_Squared, fill="R²"), 
# # #            stat="identity", position=position_dodge(), alpha=0.7) +
# # #   geom_text(aes(x=Sample, y=Pearson_Correlation, 
# # #                 label=sprintf("%.2f", Pearson_Correlation)), 
# # #             position=position_dodge(width=0.9), vjust=-0.5, size=3) +
# # #   geom_text(aes(x=Sample, y=Spearman_Correlation, 
# # #                 label=sprintf("%.2f", Spearman_Correlation)), 
# # #             position=position_dodge(width=0.9), vjust=-0.5, size=3) +
# # #   geom_text(aes(x=Sample, y=R_Squared, 
# # #                 label=sprintf("%.2f", R_Squared)), 
# # #             position=position_dodge(width=0.9), vjust=-0.5, size=3) +
# # #   ylim(0, max(1, max(c(sample_stats$Pearson_Correlation, 
# # #                        sample_stats$Spearman_Correlation, 
# # #                        sample_stats$R_Squared)) * 1.2)) +
# # #   scale_fill_brewer(palette="Set1", name="Metric") +
# # #   labs(title="Correlation Metrics by Sample",
# # #        subtitle="Pearson, Spearman, and R² values for overlapping peptides",
# # #        x="Sample", y="Correlation Value") +
# # #   theme_minimal() +
# # #   theme(plot.title = element_text(hjust=0.5, face="bold"),
# # #         plot.subtitle = element_text(hjust=0.5))
# # # 
# # # # Add number of peptides to plot title
# # # p_correlation_n <- ggplot(sample_stats, aes(x=Sample, y=Peptide_Count, fill=Sample)) +
# # #   geom_bar(stat="identity", alpha=0.8) +
# # #   geom_text(aes(label=Peptide_Count), vjust=-0.5) +
# # #   scale_fill_brewer(palette="Set2") +
# # #   labs(title="Number of Common Peptides Used for Correlation",
# # #        x="Sample", y="Peptide Count") +
# # #   theme_minimal() +
# # #   theme(plot.title = element_text(hjust=0.5, face="bold"),
# # #         legend.position="none")
# # 
# # # ===== VISUALIZATION 3: SAMPLE CORRELATION SCATTER PLOTS =====
# # # Create scatter plots for each sample
# # sample_plots <- list()
# # 
# # for(sample in unique(corr_data$Sample)) {
# #   sample_data <- corr_data[corr_data$Sample == sample, ]
# #   
# #   if(nrow(sample_data) >= 5) {
# #     # Get correlation values
# #     pearson_cor <- cor(sample_data$MSKCC_Count, sample_data$CHOP_Count, 
# #                        method="pearson", use="pairwise.complete.obs")
# #     r_squared <- sample_stats$R_Squared[sample_stats$Sample == sample]
# #     
# #     # Create scatter plot
# #     p <- ggplot(sample_data, aes(x=MSKCC_Count, y=CHOP_Count)) +
# #       geom_point(alpha=0.5, color="blue") +
# #       geom_smooth(method="lm", color="red", se=TRUE) +
# #       labs(title=paste0("Sample ", sample, 
# #                         " (R² = ", sprintf("%.3f", r_squared), 
# #                         ", n = ", nrow(sample_data), ")"),
# #            x="MSKCC Spectral Count", 
# #            y="CHOP Spectral Count") +
# #       theme_minimal() +
# #       theme(plot.title = element_text(hjust=0.5))
# #     
# #     sample_plots[[sample]] <- p
# #   }
# # }
# # 
# # # ===== VISUALIZATION 4: SUMMARY STATISTICS =====
# # # Create summary statistics
# # summary_stats <- data.frame(
# #   Metric = c(
# #     "Total MSKCC unique peptides",
# #     "Total CHOP unique peptides",
# #     "Common peptides between datasets",
# #     "Overall overlap percentage",
# #     "Average Pearson correlation",
# #     "Average Spearman correlation",
# #     "Average R-squared value",
# #     "Number of matching samples"
# #   ),
# #   Value = c(
# #     length(unique(mskcc_data$Peptide)),
# #     length(unique(if("Peptide" %in% colnames(chop_data)) chop_data$Peptide else chop_data$`Peptide Sequence`)),
# #     length(intersect(
# #       unique(mskcc_data$Peptide), 
# #       unique(if("Peptide" %in% colnames(chop_data)) chop_data$Peptide else chop_data$`Peptide Sequence`)
# #     )),
# #     round(100 * length(intersect(
# #       unique(mskcc_data$Peptide), 
# #       unique(if("Peptide" %in% colnames(chop_data)) chop_data$Peptide else chop_data$`Peptide Sequence`)
# #     )) / length(union(
# #       unique(mskcc_data$Peptide), 
# #       unique(if("Peptide" %in% colnames(chop_data)) chop_data$Peptide else chop_data$`Peptide Sequence`)
# #     )), 1),
# #     round(mean(sample_stats$Pearson_Correlation), 3),
# #     round(mean(sample_stats$Spearman_Correlation), 3),
# #     round(mean(sample_stats$R_Squared), 3),
# #     length(matching_samples)
# #   )
# # )
# # 
# # # Format summary statistics table
# # summary_table <- tableGrob(summary_stats, rows=NULL, theme=ttheme_minimal(
# #   core=list(fg_params=list(hjust=0, x=0.1), bg_params=list(fill=NA)),
# #   colhead=list(fg_params=list(hjust=0.5, fontface="bold"))
# # ))
# # 
# # # ===== SAVE ALL VISUALIZATIONS TO PDF =====
# # pdf("immunopeptidomics_visualization_report.pdf", width=11, height=8.5)
# # 
# # # Page 1: Sample overlap statistics
# # grid.arrange(
# #   p_overlap, p_pct_overlap,
# #   ncol=2, nrow=1,
# #   top=textGrob("Immunopeptidomics Datasets Comparison: Sample Overlap", 
# #                gp=gpar(fontsize=16, fontface="bold"))
# # )
# # 
# # # New Page: Correlation metrics as three separate bar graphs
# # grid.arrange(
# #   p_pearson, p_spearman, p_rsquared,
# #   ncol=1, nrow=3,  # Stacked vertically for better readability
# #   top=textGrob("Immunopeptidomics Datasets Comparison: Correlation Metrics", 
# #                gp=gpar(fontsize=16, fontface="bold"))
# # )
# # 
# # # Page 2: Correlation scatter plots
# # sample_plot_list <- lapply(names(sample_plots), function(s) sample_plots[[s]])
# # 
# # do.call(grid.arrange, 
# #         c(sample_plot_list, 
# #           list(ncol=min(3, length(sample_plot_list)), 
# #                top=textGrob("Sample-by-Sample Correlation (MSKCC vs CHOP)",
# #                             gp=gpar(fontsize=16, fontface="bold"))))
# # )
# # 
# # # Page 3: Summary statistics
# # grid.arrange(
# #   summary_table,
# #   top=textGrob("Summary Statistics", gp=gpar(fontsize=16, fontface="bold"))
# # )
# # 
# # dev.off()
# # 
# # # Save the main figures as individual PNGs for easy viewing
# # ggsave("sample_overlap.png", p_overlap, width=8, height=6)
# # ggsave("overlap_percentage.png", p_pct_overlap, width=8, height=6)
# # ggsave("pearson_correlation.png", p_pearson, width=8, height=6)
# # ggsave("spearman_correlation.png", p_spearman, width=8, height=6)
# # ggsave("r_squared.png", p_rsquared, width=8, height=6)
# # 
# # # Print completion message
# # cat("Visualizations created and saved to:\n")
# # cat("1. immunopeptidomics_visualization_report.pdf (comprehensive report)\n")
# # cat("2. sample_overlap.png\n")
# # cat("3. overlap_percentage.png\n")
# # cat("4. pearson_correlation.png\n")
# # cat("5. spearman_correlation.png\n")
# # cat("6. r_squared.png\n")
# # 
# # # Print summary statistics
# # cat("\nSummary Statistics:\n")
# # print(summary_stats)
# # 
# # # Print correlation statistics by sample
# # cat("\nCorrelation Statistics by Sample:\n")
# # print(sample_stats[order(as.numeric(as.character(sample_stats$Sample))), ])
# # 
# # # Print overlap statistics by sample
# # cat("\nOverlap Statistics by Sample:\n")
# # print(overlap_df[order(as.numeric(as.character(overlap_df$Sample))), ])
# # 
# # # Add this section to include annotation analysis
# # cat("\n===== ANNOTATION ANALYSIS =====\n")
# # # Load the annotations file
# # annotations_file <- "peptide_annotations_all.tsv"
# # annotation_data <- fread(annotations_file, sep="\t", header=TRUE, fill=TRUE)
# # 
# # # Count peptides with multiple annotations
# # multi_gene_count <- sum(grepl(";", annotation_data$All_Genes, fixed=TRUE), na.rm=TRUE)
# # multi_protein_count <- sum(grepl(";", annotation_data$All_Proteins, fixed=TRUE), na.rm=TRUE)
# # 
# # cat("Peptides with multiple gene annotations:", multi_gene_count, "\n")
# # cat("Peptides with multiple protein annotations:", multi_protein_count, "\n")
# # 
# # # Count peptides found in both sources
# # both_sources <- sum(grepl("MSKCC.*CHOP|CHOP.*MSKCC", annotation_data$Sources, ignore.case=TRUE), na.rm=TRUE)
# # cat("Peptides detected in both MSKCC and CHOP:", both_sources, "\n")
# # 
# # # For overlapping samples, check if peptides detected in same sample from both sources
# # cat("\nPeptides detected in same sample from both sources:\n")
# # for(sample in matching_samples) {
# #   # Count peptides detected in this sample in both sources
# #   sample_pattern <- paste0(";.*", sample, ";|;", sample, ";|^", sample, ";|;", sample, "$|^", sample, "$")
# #   sample_in_mskcc <- sum(grepl(sample_pattern, annotation_data$MSKCC_SampleIDs), na.rm=TRUE)
# #   
# #   # Peptides detected in this sample in MSKCC and also found in CHOP
# #   sample_mskcc_and_chop <- sum(grepl(sample_pattern, annotation_data$MSKCC_SampleIDs) & 
# #                                  grepl("CHOP", annotation_data$Sources), na.rm=TRUE)
# #   
# #   cat("Sample", sample, ":", sample_mskcc_and_chop, "peptides detected in MSKCC and also found in CHOP\n")
# # }