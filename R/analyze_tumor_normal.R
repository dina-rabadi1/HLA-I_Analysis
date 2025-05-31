#!/usr/bin/env Rscript
#' HLA-I Analysis - Tumor-Normal Comparison
#' analyze_tumor_normal.R
#' Modified to work with the configuration system
#' @author Your Name
#' @version 1.0

#' Main function for tumor-normal analysis
#' @param config Configuration object
#' @param dirs Directory structure created by create_output_directories
#' @return Invisibly returns a list of results
analyze_tumor_normal <- function(config, dirs = NULL) {
  # Create output directories if not provided
  if (is.null(dirs)) {
    dirs <- create_output_directories(config$data_path, config$output_name)
  }
  
  # 1. Load and process immunopeptidome data
  cat("\n## 1. Loading and processing immunopeptidome data...\n")
  
  # Get list of immunopeptidome TSV files
  immuno_files <- list.files(path = config$data_path, 
                             pattern = ".*_peptides\\.tsv$|.*_untargeted_peptide\\.tsv$", 
                             full.names = TRUE)
  
  # Filter for tumor and normal files
  tumor_normal_files <- immuno_files[grepl(paste0(config$tumor_id, "|", config$normal_id), 
                                           immuno_files)]
  
  if (length(tumor_normal_files) == 0) {
    stop(paste("No", config$tumor_id, "or", config$normal_id, "peptide files found"))
  }
  
  # Process the immunopeptidome data
  immuno_data <- process_immunopeptidome_data(
    tumor_normal_files, 
    experiment_type = "tumor_normal",
    sample_pattern = paste0(".*_(", config$tumor_id, "|", config$normal_id, ").*"),
    filter_peptide_length = config$peptide_length_filter
  )
  
  # Create peptide summary and matrix
  peptide_matrix_data <- create_peptide_sample_matrix(
    immuno_data, 
    value_col = "Intensity", 
    id_col = "Sample_ID", 
    peptide_col = "Peptide"
  )
  
  # Perform differential expression analysis
  diff_expr_analysis <- process_tumor_normal_data(
    peptide_matrix_data$matrix, 
    tumor_id = config$tumor_id, 
    normal_id = config$normal_id
  )
  
  # 2. Identify fusion peptides if parameters are provided
  cat("\n## 2. Identifying fusion peptides...\n")
  
  if (!is.null(config$fusion_sequence) && !is.null(config$fusion_parts) && !is.null(config$junction_position)) {
    fusion_results <- identify_fusion_peptides(
      diff_expr_analysis,
      fusion_sequence = config$fusion_sequence,
      fusion_parts = config$fusion_parts,
      junction_position = config$junction_position
    )
    
    # Update analysis data with fusion information
    diff_expr_analysis <- fusion_results$all_data_with_fusion
    
    # Extract fusion peptides for reporting
    fusion_peptides <- fusion_results$fusion_peptides
    
    cat("Found", nrow(fusion_peptides), "fusion-derived peptides\n")
    if ("spans_junction" %in% colnames(fusion_peptides)) {
      cat("of which", sum(fusion_peptides$spans_junction), "span the fusion junction\n")
    }
  }
  
  # 3. Process transcriptome data if available
  if (!is.null(config$transcriptome_path) && file.exists(config$transcriptome_path)) {
    cat("\n## 3. Processing transcriptome data...\n")
    
    transcriptome_data <- process_transcriptome_data(
      config$transcriptome_path,
      tumor_col = NULL,  # Auto-detect
      normal_col = NULL,  # Auto-detect
      gene_col = "symbol"
    )
  } else {
    transcriptome_data <- NULL
  }
  
  # 4. Process proteomics data if available
  lfq_data <- NULL
  tmt_data <- NULL
  
  if (!is.null(config$lfq_path) && file.exists(config$lfq_path)) {
    cat("\n## 4a. Processing LFQ proteomics data...\n")
    
    lfq_data <- process_proteomics_data(
      config$lfq_path,
      data_type = "LFQ",
      sheet_name = "Significant and 1.5x_2"
    )
  }
  
  if (!is.null(config$tmt_path) && file.exists(config$tmt_path)) {
    cat("\n## 4b. Processing TMT proteomics data...\n")
    
    tmt_data <- process_proteomics_data(
      config$tmt_path,
      data_type = "TMT",
      sheet_name = "Significant and 1.5x_2"
    )
  }
  
  # 5. Perform multi-omics integration
  cat("\n## 5. Integrating multi-omics data...\n")
  
  multi_omics_data <- integrate_multi_omics(
    diff_expr_analysis,
    transcriptome_data = transcriptome_data,
    lfq_data = lfq_data,
    tmt_data = tmt_data
  )
  
  # 6. Identify potential public neoantigens if multiple omics datasets are available
  if ((!is.null(transcriptome_data) || !is.null(lfq_data) || !is.null(tmt_data))) {
    cat("\n## 6. Identifying potential public neoantigens...\n")
    
    neoantigen_results <- identify_public_neoantigens(multi_omics_data)
    
    # Update with neoantigen information
    multi_omics_data <- neoantigen_results$all_data
    public_neoantigens <- neoantigen_results$public_neoantigens
    
    cat("Found", nrow(public_neoantigens), "potential public neoantigens\n")
  }
  
  # 7. Generate visualizations
  cat("\n## 7. Generating visualizations...\n")
  
  # Prepare visualization directory
  viz_dir <- dirs$viz_dir
  
  # Create basic immunopeptidome visualizations
  
  # a) Volcano plot of tumor vs normal
  volcano_plot <- create_volcano_plot(
    multi_omics_data,
    fc_col = "log2_fold_change",
    label_col = "primary_gene",
    highlight_col = if ("from_fusion" %in% colnames(multi_omics_data)) "from_fusion" else NULL,
    title = paste0("Peptide Differential Expression: ", config$tumor_id, " vs ", config$normal_id)
  )
  
  save_plot(volcano_plot, file.path(viz_dir, "tumor_normal_volcano"))
  
  # b) Heatmap of peptide intensities
  # Extract intensity columns
  intensity_cols <- c(
    paste0("total_intensity_", config$tumor_id),
    paste0("total_intensity_", config$normal_id)
  )
  
  # Create heatmap for top differentially expressed peptides
  top_peptides <- multi_omics_data %>%
    arrange(desc(abs(log2_fold_change))) %>%
    head(50)
  
  heatmap_expr <- create_peptide_heatmap(
    top_peptides,
    value_cols = intensity_cols,
    is_intensity = TRUE,
    log_transform = TRUE,
    peptide_col = "Peptide",
    annotation_cols = c("detection_status", "expression_category"),
    title = paste0("Top 50 Differential Peptides: ", config$tumor_id, " vs ", config$normal_id),
    cluster_rows = TRUE
  )
  
  pdf(file.path(viz_dir, "top_peptides_heatmap.pdf"), width = 10, height = 12)
  print(heatmap_expr)
  dev.off()
  
  png(file.path(viz_dir, "top_peptides_heatmap.png"), width = 800, height = 1000, res = 100)
  print(heatmap_expr)
  dev.off()
  
  # c) Scatter plots for multi-omics comparisons (if available)
  if (!is.null(transcriptome_data)) {
    immuno_trans_scatter <- create_comparison_scatter(
      multi_omics_data, 
      "log2_fold_change_transcriptome", 
      "log2_fold_change",
      color_col = "expression_category",
      highlight_col = if ("from_fusion" %in% colnames(multi_omics_data)) "from_fusion" else NULL,
      x_label = "Log2 Fold Change Transcriptome (Tumor/Normal)",
      y_label = "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
      title = paste0("Comparison of ", config$tumor_id, "/", config$normal_id, " Fold Changes"),
      subtitle = "Immunopeptidome vs Transcriptome"
    )
    
    save_plot(immuno_trans_scatter, file.path(viz_dir, "immuno_vs_transcriptome_scatter"))
  }
  
  if (!is.null(lfq_data)) {
    immuno_lfq_scatter <- create_comparison_scatter(
      multi_omics_data, 
      "log2_fold_change_lfq", 
      "log2_fold_change",
      color_col = "expression_category_lfq",
      highlight_col = if ("from_fusion" %in% colnames(multi_omics_data)) "from_fusion" else NULL,
      x_label = "Log2 Fold Change LFQ Proteome (Tumor/Normal)",
      y_label = "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
      title = paste0("Comparison of ", config$tumor_id, "/", config$normal_id, " Fold Changes"),
      subtitle = "Immunopeptidome vs LFQ Proteome"
    )
    
    save_plot(immuno_lfq_scatter, file.path(viz_dir, "immuno_vs_lfq_scatter"))
  }
  
  if (!is.null(tmt_data)) {
    immuno_tmt_scatter <- create_comparison_scatter(
      multi_omics_data, 
      "log2_fold_change_tmt", 
      "log2_fold_change",
      color_col = "expression_category_tmt",
      highlight_col = if ("from_fusion" %in% colnames(multi_omics_data)) "from_fusion" else NULL,
      x_label = "Log2 Fold Change TMT Proteome (Tumor/Normal)",
      y_label = "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
      title = paste0("Comparison of ", config$tumor_id, "/", config$normal_id, " Fold Changes"),
      subtitle = "Immunopeptidome vs TMT Proteome"
    )
    
    save_plot(immuno_tmt_scatter, file.path(viz_dir, "immuno_vs_tmt_scatter"))
  }
  
  # d) Fusion peptide visualizations (if available)
  if (exists("fusion_peptides") && nrow(fusion_peptides) > 0) {
    # Create heatmap of fusion peptides
    fusion_heatmap <- create_peptide_heatmap(
      fusion_peptides,
      value_cols = intensity_cols,
      is_intensity = TRUE,
      log_transform = TRUE,
      peptide_col = "Peptide",
      annotation_cols = c("fusion_peptide_type", "spans_junction"),
      title = paste0("Fusion Peptides in ", config$tumor_id, " vs ", config$normal_id),
      cluster_rows = FALSE
    )
    
    pdf(file.path(viz_dir, "fusion_peptides_heatmap.pdf"), width = 10, height = max(8, nrow(fusion_peptides)/3))
    print(fusion_heatmap)
    dev.off()
    
    png(file.path(viz_dir, "fusion_peptides_heatmap.png"), width = 800, height = max(600, nrow(fusion_peptides)*40), res = 100)
    print(fusion_heatmap)
    dev.off()
    
    # Create multi-omics fusion peptide heatmap if data is available
    if (exists("public_neoantigens") && any(public_neoantigens$from_fusion)) {
      fusion_public_neoantigens <- public_neoantigens %>%
        filter(from_fusion)
      
      if (nrow(fusion_public_neoantigens) > 0) {
        # Create matrix for the heatmap
        fusion_omics_cols <- c("log2_fold_change")
        if ("log2_fold_change_transcriptome" %in% colnames(fusion_public_neoantigens)) {
          fusion_omics_cols <- c(fusion_omics_cols, "log2_fold_change_transcriptome")
        }
        if ("log2_fold_change_lfq" %in% colnames(fusion_public_neoantigens)) {
          fusion_omics_cols <- c(fusion_omics_cols, "log2_fold_change_lfq")
        }
        if ("log2_fold_change_tmt" %in% colnames(fusion_public_neoantigens)) {
          fusion_omics_cols <- c(fusion_omics_cols, "log2_fold_change_tmt")
        }
        
        fusion_omics_heatmap <- create_peptide_heatmap(
          fusion_public_neoantigens,
          value_cols = fusion_omics_cols,
          is_intensity = FALSE,
          log_transform = FALSE,
          peptide_col = "Peptide",
          annotation_cols = c("fusion_peptide_type", "spans_junction", "public_neoantigen_classification"),
          title = "Fusion-Derived Public Neoantigens: Fold Changes Across Omics Datasets",
          cluster_rows = FALSE
        )
        
        pdf(file.path(viz_dir, "fusion_neoantigens_heatmap.pdf"), width = 12, height = max(8, nrow(fusion_public_neoantigens)/3))
        print(fusion_omics_heatmap)
        dev.off()
        
        png(file.path(viz_dir, "fusion_neoantigens_heatmap.png"), width = 1000, height = max(600, nrow(fusion_public_neoantigens)*40), res = 100)
        print(fusion_omics_heatmap)
        dev.off()
      }
    }
  }
  
  # e) Public neoantigen visualizations (if available)
  if (exists("public_neoantigens") && nrow(public_neoantigens) > 0) {
    # Create barplot of public neoantigen tiers
    tier_summary <- public_neoantigens %>%
      group_by(public_neoantigen_classification) %>%
      summarise(count = n(), .groups = "drop")
    
    tier_plot <- ggplot(tier_summary, 
                        aes(x = reorder(public_neoantigen_classification, -count), 
                            y = count, 
                            fill = public_neoantigen_classification)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = count), vjust = -0.5) +
      theme_minimal() +
      scale_fill_brewer(palette = "Set1") +
      labs(
        title = "Distribution of Potential Public Neoantigens",
        x = "Classification Tier",
        y = "Count",
        fill = "Classification"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    save_plot(tier_plot, file.path(viz_dir, "public_neoantigen_tiers"))
    
    # Create heatmap of public neoantigens if multiple datasets
    if (length(grep("log2_fold_change_", colnames(public_neoantigens))) > 1) {
      # Create matrix for the heatmap
      neoantigen_omics_cols <- c("log2_fold_change")
      col_labels <- c("Immunopeptidome")
      
      if ("log2_fold_change_transcriptome" %in% colnames(public_neoantigens)) {
        neoantigen_omics_cols <- c(neoantigen_omics_cols, "log2_fold_change_transcriptome")
        col_labels <- c(col_labels, "Transcriptome")
      }
      if ("log2_fold_change_lfq" %in% colnames(public_neoantigens)) {
        neoantigen_omics_cols <- c(neoantigen_omics_cols, "log2_fold_change_lfq")
        col_labels <- c(col_labels, "LFQ Proteome")
      }
      if ("log2_fold_change_tmt" %in% colnames(public_neoantigens)) {
        neoantigen_omics_cols <- c(neoantigen_omics_cols, "log2_fold_change_tmt")
        col_labels <- c(col_labels, "TMT Proteome")
      }
      
      # Adjust height calculation based on number of rows
      max_height <- min(40, max(8, nrow(public_neoantigens)/2))
      
      # Create diverging color palette
      div_colors <- colorRampPalette(c("blue", "white", "red"))(100)
      
      # Create unique row names with gene and peptide info
      public_neoantigens$row_label <- paste0(public_neoantigens$primary_gene, " (", 
                                             public_neoantigens$Peptide, ")")
      
      # Convert to matrix for heatmap
      neoantigen_matrix <- public_neoantigens %>%
        select(all_of(neoantigen_omics_cols), row_label) %>%
        mutate(across(all_of(neoantigen_omics_cols), ~ ifelse(is.na(.), 0, .))) %>%
        column_to_rownames("row_label") %>%
        as.matrix()
      
      # Create row annotations
      row_annotation <- data.frame(
        Classification = public_neoantigens$public_neoantigen_classification,
        row.names = rownames(neoantigen_matrix)
      )
      
      # Create the heatmap
      pdf(file.path(viz_dir, "public_neoantigens_heatmap.pdf"), width = 12, height = max_height)
      pheatmap(
        neoantigen_matrix,
        main = "Potential Public Neoantigens: Log2 Fold Changes Across Omics Datasets",
        color = div_colors,
        breaks = seq(-3, 3, length.out = 101),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        annotation_row = row_annotation,
        display_numbers = nrow(neoantigen_matrix) <= 50,
        number_format = "%.1f",
        fontsize_row = max(4, min(10, 300/nrow(neoantigen_matrix))),
        fontsize_col = 10,
        labels_col = col_labels
      )
      dev.off()
      
      png(file.path(viz_dir, "public_neoantigens_heatmap.png"), width = 1200, height = min(2000, max(800, nrow(neoantigen_matrix)*50)), res = 120)
      pheatmap(
        neoantigen_matrix,
        main = "Potential Public Neoantigens: Log2 Fold Changes Across Omics Datasets",
        color = div_colors,
        breaks = seq(-3, 3, length.out = 101),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        annotation_row = row_annotation,
        display_numbers = nrow(neoantigen_matrix) <= 50,
        number_format = "%.1f",
        fontsize_row = max(4, min(10, 300/nrow(neoantigen_matrix))),
        fontsize_col = 10,
        labels_col = col_labels
      )
      dev.off()
    }
  }
  # 8. Generate interactive visualizations if requested
  if (config$generate_interactive) {
    cat("\n## 8. Creating interactive visualizations...\n")
    
    # Convert appropriate plots to interactive versions
    interactive_viz_list <- list()
    
    # Volcano plot
    if (exists("volcano_plot")) {
      volcano_interactive <- ggplotly(volcano_plot)
      htmlwidgets::saveWidget(volcano_interactive, 
                              file.path(viz_dir, "interactive_volcano.html"), 
                              selfcontained = TRUE)
      
      interactive_viz_list[[length(interactive_viz_list) + 1]] <- list(
        path = file.path(viz_dir, "interactive_volcano.html"),
        title = "Tumor vs Normal Volcano Plot"
      )
    }
    
    # Multi-omics comparison plots
    if (exists("immuno_trans_scatter")) {
      immuno_trans_interactive <- ggplotly(immuno_trans_scatter)
      htmlwidgets::saveWidget(immuno_trans_interactive, 
                              file.path(viz_dir, "interactive_immuno_vs_trans.html"), 
                              selfcontained = TRUE)
      
      interactive_viz_list[[length(interactive_viz_list) + 1]] <- list(
        path = file.path(viz_dir, "interactive_immuno_vs_trans.html"),
        title = "Immunopeptidome vs Transcriptome"
      )
    }
    
    if (exists("immuno_lfq_scatter")) {
      immuno_lfq_interactive <- ggplotly(immuno_lfq_scatter)
      htmlwidgets::saveWidget(immuno_lfq_interactive, 
                              file.path(viz_dir, "interactive_immuno_vs_lfq.html"), 
                              selfcontained = TRUE)
      
      interactive_viz_list[[length(interactive_viz_list) + 1]] <- list(
        path = file.path(viz_dir, "interactive_immuno_vs_lfq.html"),
        title = "Immunopeptidome vs LFQ Proteome"
      )
    }
    
    if (exists("immuno_tmt_scatter")) {
      immuno_tmt_interactive <- ggplotly(immuno_tmt_scatter)
      htmlwidgets::saveWidget(immuno_tmt_interactive, 
                              file.path(viz_dir, "interactive_immuno_vs_tmt.html"), 
                              selfcontained = TRUE)
      
      interactive_viz_list[[length(interactive_viz_list) + 1]] <- list(
        path = file.path(viz_dir, "interactive_immuno_vs_tmt.html"),
        title = "Immunopeptidome vs TMT Proteome"
      )
    }
    
    # Create public neoantigen tier plot if available
    if (exists("tier_plot")) {
      tier_interactive <- ggplotly(tier_plot)
      htmlwidgets::saveWidget(tier_interactive, 
                              file.path(viz_dir, "interactive_neoantigen_tiers.html"), 
                              selfcontained = TRUE)
      
      interactive_viz_list[[length(interactive_viz_list) + 1]] <- list(
        path = file.path(viz_dir, "interactive_neoantigen_tiers.html"),
        title = "Public Neoantigen Tiers"
      )
    }
    
    # Create dashboard if we have interactive visualizations
    if (length(interactive_viz_list) > 0) {
      cat("Creating interactive dashboard...\n")
      
      # Generate appropriate summary text
      summary_text <- paste0(
        "Analysis of tumor vs normal peptides for samples ", config$tumor_id, " and ", config$normal_id, 
        ". A total of ", nrow(multi_omics_data), " peptides were analyzed. ",
        if(exists("public_neoantigens")) paste0("Found ", nrow(public_neoantigens), " potential public neoantigens. ") else "",
        if(exists("fusion_peptides")) paste0("Found ", nrow(fusion_peptides), " fusion-derived peptides",
                                             if("spans_junction" %in% colnames(fusion_peptides)) 
                                               paste0(", of which ", sum(fusion_peptides$spans_junction), " span the fusion junction. ") 
                                             else ". ") else ""
      )
      
      # Create the dashboard
      create_interactive_dashboard(
        interactive_viz_list,
        title = paste0("Tumor-Normal Analysis: ", config$tumor_id, " vs ", config$normal_id),
        output_file = file.path(viz_dir, "interactive_dashboard.html"),
        summary_text = summary_text
      )
    }
  }
  
  # 9. Generate Excel reports
  cat("\n## 9. Generating Excel reports...\n")
  
  # Prepare Excel data sheets
  excel_sheets <- prepare_excel_data(multi_omics_data, analysis_type = "tumor_normal")
  
  # Generate Excel report
  excel_output_file <- file.path(dirs$excel_dir, 
                                 paste0(config$output_name, "_results.xlsx"))
  generate_excel_report(excel_sheets, excel_output_file)
  
  # 10. Print summary information
  cat("\n## 10. Analysis summary:\n")
  cat("\nAnalysis complete! Results saved to:", dirs$main_dir, "\n")
  
  # Summary counts
  cat("\nTotal peptides analyzed:", nrow(multi_omics_data), "\n")
  
  if (!is.null(transcriptome_data)) {
    matching_trans <- sum(!is.na(multi_omics_data$log2_fold_change_transcriptome))
    cat("Peptides with matching transcriptome data:", matching_trans, 
        "(", round(matching_trans/nrow(multi_omics_data)*100, 1), "%)\n")
  }
  
  if (!is.null(lfq_data)) {
    matching_lfq <- sum(!is.na(multi_omics_data$log2_fold_change_lfq))
    cat("Peptides with matching LFQ proteome data:", matching_lfq, 
        "(", round(matching_lfq/nrow(multi_omics_data)*100, 1), "%)\n")
  }
  
  if (!is.null(tmt_data)) {
    matching_tmt <- sum(!is.na(multi_omics_data$log2_fold_change_tmt))
    cat("Peptides with matching TMT proteome data:", matching_tmt, 
        "(", round(matching_tmt/nrow(multi_omics_data)*100, 1), "%)\n")
  }
  
  # Detection status summary
  if ("detection_status" %in% colnames(multi_omics_data)) {
    status_counts <- multi_omics_data %>%
      count(detection_status) %>%
      mutate(percentage = round(n / sum(n) * 100, 1))
    
    cat("\nDetection status summary:\n")
    print(status_counts)
  }
  
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
  
  # Public neoantigen summary
  if (exists("public_neoantigens") && nrow(public_neoantigens) > 0) {
    cat("\nPublic neoantigen summary:\n")
    cat("Total potential public neoantigens:", nrow(public_neoantigens), "\n")
    
    if ("public_neoantigen_classification" %in% colnames(public_neoantigens)) {
      tier_counts <- public_neoantigens %>%
        count(public_neoantigen_classification) %>%
        mutate(percentage = round(n / sum(n) * 100, 1))
      
      print(tier_counts)
    }
  }
  
  cat("\nOutput files generated in the following locations:\n")
  cat("- Excel reports:", dirs$excel_dir, "\n")
  cat("- Visualizations:", dirs$viz_dir, "\n")
  cat("- Processed data:", dirs$data_dir, "\n")
  
  # Save the processed data objects for future use
  save(multi_omics_data, file = file.path(dirs$data_dir, "multi_omics_data.RData"))
  
  if (exists("fusion_peptides")) {
    save(fusion_peptides, file = file.path(dirs$data_dir, "fusion_peptides.RData"))
  }
  
  if (exists("public_neoantigens")) {
    save(public_neoantigens, file = file.path(dirs$data_dir, "public_neoantigens.RData"))
  }
  
  # Return results invisibly
  invisible(list(
    multi_omics_data = multi_omics_data,
    fusion_peptides = if(exists("fusion_peptides")) fusion_peptides else NULL,
    public_neoantigens = if(exists("public_neoantigens")) public_neoantigens else NULL,
    directories = dirs
  ))
}
