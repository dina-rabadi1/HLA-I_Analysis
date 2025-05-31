#!/usr/bin/env Rscript
#' HLA-I Analysis - Multi-Sample Comparison
#' analyze_multi_sample.R
#' Master script for analyzing immunopeptidome data across multiple samples
#' @author Your Name
#' @version 1.0

# Load helper functions
source("peptide_core_utils.R")
source("peptide_data_processing.R")
source("peptide_visualizations.R")
source("peptide_integration.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Configure run settings
config <- list(
  data_path = if(length(args) >= 1) args[1] else "/path/to/data",
  output_name = if(length(args) >= 2) args[2] else "multi_sample_analysis",
  sample_pattern = ".*_([0-9]+[A-Za-z]?)_.*", # Pattern to extract sample IDs
  fusion_parts = list(
    part1 = "RKREIFDRYGEE",  
    part2 = "VKEFLAKAKEDF" 
  ),
  fusion_sequence = "RKREIFDRYGEEVKEFLAKAKEDF",
  junction_position = 12,
  peptide_length_filter = c(8, 12),
  spiked_peptides = NULL,  # Set this to a vector of spiked peptide sequences if applicable
  spiked_sample = NULL,    # Set this to the sample ID that was spiked
  generate_interactive = TRUE
)

# Load required packages
load_required_packages(c("tidyverse", "openxlsx", "ggplot2", "pheatmap", 
                         "plotly", "VennDiagram", "RColorBrewer", "UpSetR"))

# Create output directories
dirs <- create_output_directories(config$data_path, config$output_name)

# 1. Load and process immunopeptidome data
cat("\n## 1. Loading and processing immunopeptidome data...\n")

# Get list of immunopeptidome TSV files
immuno_files <- list.files(path = config$data_path, 
                           pattern = ".*_peptides\\.tsv$|.*_untargeted_peptide\\.tsv$", 
                           full.names = TRUE)

if (length(immuno_files) == 0) {
  stop("No peptide files found at the specified path")
}

# Process the immunopeptidome data
immuno_data <- process_immunopeptidome_data(
  immuno_files, 
  experiment_type = "multi_sample",
  sample_pattern = config$sample_pattern,
  filter_peptide_length = config$peptide_length_filter
)

# Get unique sample IDs
sample_ids <- unique(immuno_data$Sample_ID)
cat("Found", length(sample_ids), "unique samples:", paste(sample_ids, collapse = ", "), "\n")

# Create peptide summary and matrix
peptide_matrix_data <- create_peptide_sample_matrix(
  immuno_data, 
  value_col = "Intensity", 
  id_col = "Sample_ID", 
  peptide_col = "Peptide"
)

# Perform multi-sample analysis
multi_sample_analysis <- process_multi_sample_data(
  peptide_matrix_data$matrix, 
  samples = sample_ids
)

# 2. Handle spiked peptides if specified
if (!is.null(config$spiked_peptides) && !is.null(config$spiked_sample)) {
  cat("\n## 2. Analyzing spiked peptides in sample", config$spiked_sample, "...\n")
  
  # Filter data for just the spiked peptides
  spiked_peptide_data <- immuno_data %>%
    filter(Peptide %in% config$spiked_peptides) %>%
    mutate(
      is_spiked = ifelse(Sample_ID == config$spiked_sample, "Spiked", "Natural")
    )
  
  # Create summary of spiked peptide detection
  spiked_summary <- spiked_peptide_data %>%
    group_by(Sample_ID, Peptide, is_spiked) %>%
    summarize(
      detected = TRUE,
      spectral_count = sum(Spectral.Count),
      total_intensity = sum(Intensity),
      .groups = "drop"
    ) %>%
    arrange(Peptide, Sample_ID)
  
  # Create matrix of all peptide-sample combinations
  all_combinations <- expand.grid(
    Peptide = config$spiked_peptides,
    Sample_ID = sample_ids,
    stringsAsFactors = FALSE
  )
  
  # Merge with data
  spiked_matrix_data <- all_combinations %>%
    left_join(spiked_summary, by = c("Peptide", "Sample_ID")) %>%
    mutate(
      detected = ifelse(is.na(detected), FALSE, detected),
      spectral_count = ifelse(is.na(spectral_count), 0, spectral_count),
      total_intensity = ifelse(is.na(total_intensity), 0, total_intensity),
      is_spiked = ifelse(is.na(is_spiked), "Not Detected", is_spiked)
    )
  
  cat("Analysis of", length(config$spiked_peptides), "spiked peptides complete\n")
}

# 3. Identify fusion peptides if parameters are provided
cat("\n## 3. Identifying fusion peptides...\n")

if (!is.null(config$fusion_sequence) && !is.null(config$fusion_parts) && !is.null(config$junction_position)) {
  fusion_results <- identify_fusion_peptides(
    multi_sample_analysis,
    fusion_sequence = config$fusion_sequence,
    fusion_parts = config$fusion_parts,
    junction_position = config$junction_position
  )
  
  # Update analysis data with fusion information
  multi_sample_analysis <- fusion_results$all_data_with_fusion
  
  # Extract fusion peptides for reporting
  fusion_peptides <- fusion_results$fusion_peptides
  
  cat("Found", nrow(fusion_peptides), "fusion-derived peptides\n")
  if ("spans_junction" %in% colnames(fusion_peptides)) {
    cat("of which", sum(fusion_peptides$spans_junction), "span the fusion junction\n")
  }
}

# 4. Generate visualizations
cat("\n## 4. Generating visualizations...\n")

# Prepare visualization directory
viz_dir <- dirs$viz_dir

# a) Distribution of peptides by sample count
peptide_counts <- multi_sample_analysis %>%
  count(detected_samples) %>%
  mutate(percentage = n / sum(n) * 100)

count_plot <- ggplot(peptide_counts, aes(x = factor(detected_samples), y = n,
                                         text = paste0("Sample count: ", detected_samples,
                                                       "\nPeptides: ", n,
                                                       "\nPercentage: ", round(percentage, 1), "%"))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Peptide Distribution by Sample Count",
    x = "Number of Samples",
    y = "Number of Peptides"
  ) +
  theme_minimal()

save_plot(count_plot, file.path(viz_dir, "peptide_sample_count_distribution"))

# b) Private vs Shared peptides pie chart
private_shared <- data.frame(
  category = c("Private", "Shared"),
  count = c(
    sum(peptide_counts$n[peptide_counts$detected_samples == 1]),
    sum(peptide_counts$n[peptide_counts$detected_samples > 1])
  )
) %>%
  mutate(percentage = count / sum(count) * 100)

pie_plot <- ggplot(private_shared, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(
    title = "Private vs Shared Peptides",
    fill = "Category"
  ) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

save_plot(pie_plot, file.path(viz_dir, "private_vs_shared_pie"))

# c) Peptide Length Distribution
length_dist <- multi_sample_analysis %>%
  count(peptide_length) %>%
  mutate(percentage = n / sum(n) * 100)

length_plot <- ggplot(length_dist, aes(x = factor(peptide_length), y = n,
                                       text = paste0("Length: ", peptide_length, " amino acids",
                                                     "\nCount: ", n,
                                                     "\nPercentage: ", round(percentage, 1), "%"))) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  labs(
    title = "Peptide Length Distribution",
    x = "Peptide Length (amino acids)",
    y = "Count"
  ) +
  theme_minimal()

save_plot(length_plot, file.path(viz_dir, "peptide_length_distribution"))

# d) Heatmap of peptide presence/absence across samples
# Get the top shared peptides
if (length(sample_ids) > 1) {
  # Get intensity columns for all samples
  intensity_cols <- paste0("total_intensity_", sample_ids)
  
  # Determine how many peptides to show in the heatmap (max 100, but fewer if fewer shared peptides)
  shared_peptides <- multi_sample_analysis %>%
    filter(is_shared == TRUE) %>%
    arrange(desc(detected_samples)) %>%
    head(100)  # Limit to top 100 to keep heatmap readable
  
  if (nrow(shared_peptides) > 0) {
    heatmap_presence <- create_peptide_heatmap(
      shared_peptides,
      value_cols = intensity_cols,
      is_intensity = FALSE,  # Binary presence/absence
      peptide_col = "Peptide",
      annotation_cols = c("detected_samples", "sharing_category"),
      title = paste0("Peptide Presence Across ", length(sample_ids), " Samples"),
      cluster_rows = TRUE,
      cluster_cols = TRUE
    )
    
    pdf(file.path(viz_dir, "peptide_presence_heatmap.pdf"), width = 10, height = min(25, max(8, nrow(shared_peptides)/4)))
    print(heatmap_presence)
    dev.off()
    
    png(file.path(viz_dir, "peptide_presence_heatmap.png"), width = 800, height = min(2000, max(600, nrow(shared_peptides)*20)), res = 100)
    print(heatmap_presence)
    dev.off()
  }
  
  # e) Heatmap of peptide intensities across samples
  # Get top peptides with highest total intensity
  top_intensity_peptides <- multi_sample_analysis %>%
    rowwise() %>%
    mutate(
      total_intensity = sum(across(all_of(intensity_cols), ~ifelse(is.na(.), 0, .)))
    ) %>%
    arrange(desc(total_intensity)) %>%
    head(50)  # Top 50 by intensity
  
  if (nrow(top_intensity_peptides) > 0) {
    heatmap_intensity <- create_peptide_heatmap(
      top_intensity_peptides,
      value_cols = intensity_cols,
      is_intensity = TRUE,  # Intensity values
      log_transform = TRUE,
      peptide_col = "Peptide",
      annotation_cols = c("detected_samples", "sharing_category"),
      title = paste0("Top 50 Peptides by Intensity Across ", length(sample_ids), " Samples"),
      cluster_rows = TRUE,
      cluster_cols = TRUE
    )
    
    pdf(file.path(viz_dir, "peptide_intensity_heatmap.pdf"), width = 10, height = min(20, max(8, nrow(top_intensity_peptides)/3)))
    print(heatmap_intensity)
    dev.off()
    
    png(file.path(viz_dir, "peptide_intensity_heatmap.png"), width = 800, height = min(1600, max(600, nrow(top_intensity_peptides)*25)), res = 100)
    print(heatmap_intensity)
    dev.off()
  }
  
  # f) Create UpSet plot for sample intersections if more than 2 samples
  if (length(sample_ids) >= 2 && length(sample_ids) <= 10) {
    # Create binary matrix for upset plot
    presence_matrix <- matrix(0, nrow = nrow(multi_sample_analysis), ncol = length(sample_ids))
    colnames(presence_matrix) <- sample_ids
    
    # Populate the matrix
    for (i in 1:length(sample_ids)) {
      col_name <- paste0("total_intensity_", sample_ids[i])
      presence_matrix[, i] <- ifelse(multi_sample_analysis[[col_name]] > 0, 1, 0)
    }
    
    # Convert to data frame for UpSetR
    presence_df <- as.data.frame(presence_matrix)
    
    # Create UpSet plot
    upset_file <- file.path(viz_dir, "sample_intersections_upset.pdf")
    pdf(upset_file, width = 10, height = 8)
    upset_plot <- upset(
      presence_df,
      nsets = length(sample_ids),
      nintersects = min(30, 2^length(sample_ids) - 1),  # Limit to 30 intersections or fewer
      mb.ratio = c(0.5, 0.5),
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "darkred",
      keep.order = FALSE,
      set_size.show = TRUE,
      text.scale = 1.2,
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size"
    )
    dev.off()
    
    # Also save as PNG
    png(file.path(viz_dir, "sample_intersections_upset.png"), width = 1000, height = 800, res = 100)
    upset(
      presence_df,
      nsets = length(sample_ids),
      nintersects = min(30, 2^length(sample_ids) - 1),
      mb.ratio = c(0.5, 0.5),
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "darkred",
      keep.order = FALSE,
      set_size.show = TRUE,
      text.scale = 1.2,
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size"
    )
    dev.off()
    
    cat("Created UpSet plot for sample intersections\n")
  }
  
  # g) Create Venn diagram if 2-5 samples
  if (length(sample_ids) >= 2 && length(sample_ids) <= 5) {
    venn_plot <- create_venn_diagram(
      multi_sample_analysis,
      sample_cols = paste0("total_intensity_", sample_ids),
      labels = sample_ids,
      title = "Peptide Overlap Between Samples",
      output_file = file.path(viz_dir, "sample_venn_diagram.pdf")
    )
    
    # Also save as PNG
    png(file.path(viz_dir, "sample_venn_diagram.png"), width = 800, height = 800, res = 100)
    grid.draw(venn_plot)
    dev.off()
    
    cat("Created Venn diagram for peptide sample overlaps\n")
  }
}

# h) Fusion peptide visualizations (if available)
if (exists("fusion_peptides") && nrow(fusion_peptides) > 0) {
  # Get intensity columns for all samples
  intensity_cols <- paste0("total_intensity_", sample_ids)
  
  # Create heatmap of fusion peptides
  fusion_heatmap <- create_peptide_heatmap(
    fusion_peptides,
    value_cols = intensity_cols,
    is_intensity = TRUE,
    log_transform = TRUE,
    peptide_col = "Peptide",
    annotation_cols = c("fusion_peptide_type", "spans_junction"),
    title = "Fusion Peptides Across Samples",
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  
  pdf(file.path(viz_dir, "fusion_peptides_heatmap.pdf"), width = 10, height = max(8, nrow(fusion_peptides)/3))
  print(fusion_heatmap)
  dev.off()
  
  png(file.path(viz_dir, "fusion_peptides_heatmap.png"), width = 800, height = max(600, nrow(fusion_peptides)*40), res = 100)
  print(fusion_heatmap)
  dev.off()
  
  # Create bar plot of fusion peptide types
  if ("fusion_peptide_type" %in% colnames(fusion_peptides)) {
    fusion_type_counts <- fusion_peptides %>%
      count(fusion_peptide_type) %>%
      arrange(desc(n))
    
    fusion_type_plot <- ggplot(fusion_type_counts, aes(x = reorder(fusion_peptide_type, -n), y = n, fill = fusion_peptide_type)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = n), vjust = -0.5) +
      labs(
        title = "Fusion Peptide Types",
        x = "Peptide Type",
        y = "Count"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    
    save_plot(fusion_type_plot, file.path(viz_dir, "fusion_peptide_types"))
  }
}

# i) Spiked peptide visualizations (if available)
if (exists("spiked_matrix_data") && nrow(spiked_matrix_data) > 0) {
  # Create bar plot of spiked peptide detection
  detection_summary <- spiked_matrix_data %>%
    group_by(Peptide) %>%
    summarize(
      detected_samples = sum(detected),
      spiked_sample_detected = sum(detected & (Sample_ID == config$spiked_sample)),
      .groups = "drop"
    ) %>%
    mutate(
      detection_status = case_when(
        spiked_sample_detected > 0 ~ "Detected in spiked sample",
        detected_samples > 0 ~ "Detected in other samples only",
        TRUE ~ "Not detected"
      )
    )
  
  detection_plot <- ggplot(detection_summary, aes(x = reorder(Peptide, -detected_samples), y = detected_samples, fill = detection_status)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = detected_samples), vjust = -0.5) +
    scale_fill_manual(values = c("Detected in spiked sample" = "darkgreen", 
                                 "Detected in other samples only" = "steelblue",
                                 "Not detected" = "gray80")) +
    labs(
      title = paste0("Detection of Spiked Peptides Across ", length(sample_ids), " Samples"),
      x = "Peptide",
      y = "Number of Samples Detected In",
      fill = "Detection Status"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  save_plot(detection_plot, file.path(viz_dir, "spiked_peptide_detection"))
  
  # Create intensity barplot for detected peptides
  intensity_data <- spiked_matrix_data %>%
    filter(detected) %>%
    mutate(
      label_text = paste0(Sample_ID, ifelse(is_spiked == "Spiked", " (Spiked)", ""))
    )
  
  if (nrow(intensity_data) > 0) {
    intensity_plot <- ggplot(intensity_data, aes(x = interaction(Peptide, Sample_ID), y = total_intensity, fill = is_spiked)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Spiked" = "darkred", "Natural" = "steelblue", "Not Detected" = "gray80")) +
      labs(
        title = "Intensity of Detected Spiked Peptides",
        x = "Peptide-Sample",
        y = "Total Intensity",
        fill = "Type"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    save_plot(intensity_plot, file.path(viz_dir, "spiked_peptide_intensity"))
  }
}

# 5. Generate interactive visualizations if requested
if (config$generate_interactive) {
  cat("\n## 5. Creating interactive visualizations...\n")
  
  # Convert appropriate plots to interactive versions
  interactive_viz_list <- list()
  
  # Sample count distribution
  if (exists("count_plot")) {
    count_interactive <- ggplotly(count_plot)
    htmlwidgets::saveWidget(count_interactive, 
                            file.path(viz_dir, "interactive_sample_count.html"), 
                            selfcontained = TRUE)
    
    interactive_viz_list[[length(interactive_viz_list) + 1]] <- list(
      path = file.path(viz_dir, "interactive_sample_count.html"),
      title = "Peptide Distribution by Sample Count"
    )
  }
  
  # Peptide length distribution
  if (exists("length_plot")) {
    length_interactive <- ggplotly(length_plot)
    htmlwidgets::saveWidget(length_interactive, 
                            file.path(viz_dir, "interactive_length_dist.html"), 
                            selfcontained = TRUE)
    
    interactive_viz_list[[length(interactive_viz_list) + 1]] <- list(
      path = file.path(viz_dir, "interactive_length_dist.html"),
      title = "Peptide Length Distribution"
    )
  }
  
  # Fusion peptide type plot
  if (exists("fusion_type_plot")) {
    fusion_type_interactive <- ggplotly(fusion_type_plot)
    htmlwidgets::saveWidget(fusion_type_interactive, 
                            file.path(viz_dir, "interactive_fusion_types.html"), 
                            selfcontained = TRUE)
    
    interactive_viz_list[[length(interactive_viz_list) + 1]] <- list(
      path = file.path(viz_dir, "interactive_fusion_types.html"),
      title = "Fusion Peptide Types"
    )
  }
  
  # Spiked peptide detection
  if (exists("detection_plot")) {
    detection_interactive <- ggplotly(detection_plot)
    htmlwidgets::saveWidget(detection_interactive, 
                            file.path(viz_dir, "interactive_spiked_detection.html"), 
                            selfcontained = TRUE)
    
    interactive_viz_list[[length(interactive_viz_list) + 1]] <- list(
      path = file.path(viz_dir, "interactive_spiked_detection.html"),
      title = "Spiked Peptide Detection"
    )
  }
  
  # Create dashboard if we have interactive visualizations
  if (length(interactive_viz_list) > 0) {
    cat("Creating interactive dashboard...\n")
    
    # Generate appropriate summary text
    summary_text <- paste0(
      "Analysis of peptides across ", length(sample_ids), " samples. ",
      "A total of ", nrow(multi_sample_analysis), " peptides were analyzed, ",
      "of which ", sum(multi_sample_analysis$is_private), " are private and ",
      sum(multi_sample_analysis$is_shared), " are shared between samples. ",
      if(exists("fusion_peptides")) paste0("Found ", nrow(fusion_peptides), " fusion-derived peptides",
                                           if("spans_junction" %in% colnames(fusion_peptides)) 
                                             paste0(", of which ", sum(fusion_peptides$spans_junction), " span the fusion junction. ") 
                                           else ". ") else ""
    )
    
    # Create the dashboard
    create_interactive_dashboard(
      interactive_viz_list,
      title = paste0("Multi-Sample Analysis: ", length(sample_ids), " Samples"),
      output_file = file.path(viz_dir, "interactive_dashboard.html"),
      summary_text = summary_text
    )
  }
}

# 6. Generate Excel reports
cat("\n## 6. Generating Excel reports...\n")

# Prepare Excel data sheets
excel_sheets <- prepare_excel_data(multi_sample_analysis, analysis_type = "multi_sample")

# Add fusion peptides sheet if available
if (exists("fusion_peptides") && !("Fusion_Peptides" %in% names(excel_sheets))) {
  excel_sheets[["Fusion_Peptides"]] <- fusion_peptides
}

# Add spiked peptides sheet if available
if (exists("spiked_matrix_data")) {
  excel_sheets[["Spiked_Peptides"]] <- spiked_matrix_data
}

# Generate Excel report
excel_output_file <- file.path(dirs$excel_dir, 
                               paste0(config$output_name, "_results.xlsx"))
generate_excel_report(excel_sheets, excel_output_file)

# 7. Print summary information
cat("\n## 7. Analysis summary:\n")
cat("\nAnalysis complete! Results saved to:", dirs$main_dir, "\n")

# Summary counts
cat("\nTotal peptides analyzed:", nrow(multi_sample_analysis), "\n")
cat("Private peptides (found in 1 sample):", sum(multi_sample_analysis$is_private), "\n")
cat("Shared peptides (found in >1 sample):", sum(multi_sample_analysis$is_shared), "\n")

# Sample sharing distribution
sample_dist <- multi_sample_analysis %>%
  count(detected_samples) %>%
  mutate(percentage = round(n / sum(n) * 100, 1))

cat("\nPeptide distribution by sample count:\n")
print(sample_dist)

# Fusion peptide summary
if (exists("fusion_peptides") && nrow(fusion_peptides) > 0) {
  cat("\nFusion peptide summary:\n")
  cat("Total fusion-derived peptides:", nrow(fusion_peptides), "\n")
  
  if ("spans_junction" %in% colnames(fusion_peptides)) {
    cat("Junction-spanning peptides:", sum(fusion_peptides$spans_junction), "\n")
  }
  
  if ("fusion_peptide_type" %in% colnames(fusion_peptides)) {
    type_counts <- fusion_peptides %>%
      count(fusion_peptide_type) %>%
      mutate(percentage = round(n / sum(n) * 100, 1))
    
    print(type_counts)
  }
}

# Spiked peptide summary
if (exists("spiked_matrix_data")) {
  spiked_summary <- spiked_matrix_data %>%
    group_by(Peptide) %>%
    summarize(
      detected_count = sum(detected),
      spiked_detected = sum(detected & (Sample_ID == config$spiked_sample)),
      .groups = "drop"
    )
  
  cat("\nSpiked peptide summary:\n")
  cat("Total spiked peptides:", length(config$spiked_peptides), "\n")
  cat("Spiked peptides detected in any sample:", sum(spiked_summary$detected_count > 0), "\n")
  cat("Spiked peptides detected in spiked sample (", config$spiked_sample, "):", 
      sum(spiked_summary$spiked_detected > 0), "\n")
}

cat("\nOutput files generated in the following locations:\n")
cat("- Excel reports:", dirs$excel_dir, "\n")
cat("- Visualizations:", dirs$viz_dir, "\n")
cat("- Processed data:", dirs$data_dir, "\n")

# Save the processed data objects for future use
save(multi_sample_analysis, file = file.path(dirs$data_dir, "multi_sample_analysis.RData"))

if (exists("fusion_peptides")) {
  save(fusion_peptides, file = file.path(dirs$data_dir, "fusion_peptides.RData"))
}

if (exists("spiked_matrix_data")) {
  save(spiked_matrix_data, file = file.path(dirs$data_dir, "spiked_peptides.RData"))
}

cat("\nAnalysis complete!\n")