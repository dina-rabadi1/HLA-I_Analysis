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
      tumor_col = config$tumor_id,  
      normal_col = config$normal_id,
      gene_col = "symbol",
      sample_pattern = NULL  # Use direct tumor/normal IDs instead of pattern
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
    
    # Create expression category barplot if available
    if ("expression_category" %in% colnames(multi_omics_data)) {
      expression_summary <- multi_omics_data %>%
        group_by(expression_category) %>%
        summarise(count = n(), .groups = "drop") %>%
        arrange(desc(count))
      
      expr_category_plot <- ggplot(expression_summary, 
                                   aes(x = reorder(expression_category, -count), 
                                       y = count, 
                                       fill = expression_category)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = count), vjust = -0.5) +
        theme_minimal() +
        labs(
          title = "Distribution of Expression Categories",
          x = "Category",
          y = "Count",
          fill = "Expression Category"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      save_plot(expr_category_plot, file.path(viz_dir, "expression_category_counts"))
    }
  
  # Create ranked bar plot for tumor-specific peptides
  tumor_specific_plot <- create_tumor_specific_plot(
    data = multi_omics_data,
    tumor_specific_col = "detection_status",
    intensity_col = paste0("total_intensity_", config$tumor_id),
    tumor_value = paste0(config$tumor_id, "-specific"),  # Using tumor_id variable rather than hardcoded
    label_col = "Peptide",
    title = paste0("Peptides Exclusive to ", config$tumor_id)
  )
  save_plot(tumor_specific_plot, file.path(viz_dir, "tumor_specific_peptides"))
  
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
  
  # Extract tumor and normal intensity columns
  intensity_cols <- c(
    paste0("total_intensity_", config$tumor_id),
    paste0("total_intensity_", config$normal_id)
  )
  
  # Ensure peptide_category is added to the data if it doesn't exist
  if (!("peptide_category" %in% colnames(multi_omics_data))) {
    multi_omics_data <- multi_omics_data %>%
      mutate(
        peptide_category = case_when(
          log2_fold_change > 1 ~ "Up in Tumor (FC > 2)",
          log2_fold_change < -1 ~ "Down in Tumor (FC < 0.5)",
          TRUE ~ "Similar (-1 < log2FC < 1)"
        )
      )
    cat("Added peptide_category column to multi_omics_data\n")
  }
  
  # a) Volcano plot
  volcano_plot <- create_volcano_plot(
    multi_omics_data,
    fc_col = "log2_fold_change",
    label_col = "primary_gene",
    highlight_col = if ("from_fusion" %in% colnames(multi_omics_data)) "from_fusion" else NULL,
    title = paste0("Peptide Differential Expression: ", config$tumor_id, " vs ", config$normal_id)
  )
  save_plot(volcano_plot, file.path(viz_dir, "tumor_normal_volcano"))
  
  # b) Fold change histogram 
  hist_plot <- create_fold_change_histogram(
    multi_omics_data,
    fc_col = "log2_fold_change",
    category_col = "peptide_category",
    title = paste0("Distribution of Peptide Fold Changes in ", config$tumor_id, " vs ", config$normal_id)
  )
  save_plot(hist_plot, file.path(viz_dir, "tumor_normal_fold_change_histogram"))
  
  # c) Detection status barplot
  detection_plot <- create_detection_status_plot(
    multi_omics_data,
    status_col = "detection_status",
    title = paste0("Peptide Detection Status in ", config$tumor_id, " vs ", config$normal_id)
  )
  save_plot(detection_plot, file.path(viz_dir, "tumor_normal_detection_status"))
  
  # d) Exclusive peptide analysis
  exclusive_peptides <- create_exclusive_peptide_analysis(
    multi_omics_data,
    tumor_id = config$tumor_id,
    normal_id = config$normal_id
  )
  
  # Add exclusive peptides to Excel sheet list - initialize excel_sheets if not existing
  excel_sheets <- list()
  excel_sheets[["Tumor_Exclusive_Peptides"]] <- exclusive_peptides$tumor_exclusive
  excel_sheets[["Normal_Exclusive_Peptides"]] <- exclusive_peptides$normal_exclusive
  
  # If there are fusion-specific exclusive peptides, add them too
  if (!is.null(exclusive_peptides$tumor_exclusive_fusion)) {
    excel_sheets[["Tumor_Exclusive_Fusion"]] <- exclusive_peptides$tumor_exclusive_fusion
  } else if ("from_fusion" %in% colnames(multi_omics_data) && exists("fusion_peptides")) {
    # Try to extract fusion peptides from the exclusive peptides if they exist
    tumor_exclusive_fusion <- exclusive_peptides$tumor_exclusive %>%
      filter(from_fusion) %>%
      arrange(desc(if("spans_junction" %in% colnames(.)) spans_junction else TRUE), 
              desc(!!sym(paste0("total_intensity_", config$tumor_id))))
    
    normal_exclusive_fusion <- exclusive_peptides$normal_exclusive %>%
      filter(from_fusion) %>%
      arrange(desc(if("spans_junction" %in% colnames(.)) spans_junction else TRUE), 
              desc(!!sym(paste0("total_intensity_", config$normal_id))))
    
    if(nrow(tumor_exclusive_fusion) > 0) {
      excel_sheets[["Tumor_Exclusive_Fusion"]] <- tumor_exclusive_fusion
    }
    
    if(nrow(normal_exclusive_fusion) > 0) {
      excel_sheets[["Normal_Exclusive_Fusion"]] <- normal_exclusive_fusion
    }
  }
  
  # Create ranked bar plot for tumor-specific peptides
  tumor_specific_plot <- create_tumor_specific_plot(
    data = multi_omics_data,
    tumor_specific_col = "detection_status",
    intensity_col = paste0("total_intensity_", config$tumor_id),
    tumor_value = paste0(config$tumor_id, "-specific"),  # Using dynamic value based on tumor_id
    label_col = "Peptide",
    title = paste0("Peptides Exclusive to ", config$tumor_id)
  )
  save_plot(tumor_specific_plot, file.path(viz_dir, "tumor_specific_peptides"))
  
  # e) Heatmap of peptide intensities
  # Create heatmap for top differentially expressed peptides
  top_peptides <- multi_omics_data %>%
    arrange(desc(abs(log2_fold_change))) %>%
    head(50)
  
  annotation_cols <- "detection_status"
  if ("expression_category" %in% colnames(top_peptides)) {
    annotation_cols <- c(annotation_cols, "expression_category")
  } else if ("peptide_category" %in% colnames(top_peptides)) {
    annotation_cols <- c(annotation_cols, "peptide_category")
  }
  
  heatmap_expr <- create_peptide_heatmap(
    top_peptides,
    value_cols = intensity_cols,
    is_intensity = TRUE,
    log_transform = TRUE,
    peptide_col = "Peptide",
    annotation_cols = annotation_cols,
    title = paste0("Top 50 Differential Peptides: ", config$tumor_id, " vs ", config$normal_id),
    cluster_rows = TRUE
  )
  
  pdf(file.path(viz_dir, "top_peptides_heatmap.pdf"), width = 10, height = 12)
  print(heatmap_expr)
  dev.off()
  
  png(file.path(viz_dir, "top_peptides_heatmap.png"), width = 800, height = 1000, res = 100)
  print(heatmap_expr)
  dev.off()
  
  # f) Scatter plots for multi-omics comparisons (if available)
  if (!is.null(transcriptome_data)) {
    immuno_trans_scatter <- create_comparison_scatter(
      multi_omics_data, 
      "log2_fold_change_transcriptome", 
      "log2_fold_change",
      color_col = if ("expression_category" %in% colnames(multi_omics_data)) "expression_category" else NULL,
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
      color_col = if ("expression_category_lfq" %in% colnames(multi_omics_data)) "expression_category_lfq" else NULL,
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
      color_col = if ("expression_category_tmt" %in% colnames(multi_omics_data)) "expression_category_tmt" else NULL,
      highlight_col = if ("from_fusion" %in% colnames(multi_omics_data)) "from_fusion" else NULL,
      x_label = "Log2 Fold Change TMT Proteome (Tumor/Normal)",
      y_label = "Log2 Fold Change Immunopeptidome (Tumor/Normal)",
      title = paste0("Comparison of ", config$tumor_id, "/", config$normal_id, " Fold Changes"),
      subtitle = "Immunopeptidome vs TMT Proteome"
    )
    
    save_plot(immuno_tmt_scatter, file.path(viz_dir, "immuno_vs_tmt_scatter"))
  }
  
  # g) Fusion peptide visualizations (if available)
  if (exists("fusion_peptides") && nrow(fusion_peptides) > 0) {
    # Create heatmap of fusion peptides
    fusion_annotation_cols <- NULL
    if ("fusion_peptide_type" %in% colnames(fusion_peptides)) {
      fusion_annotation_cols <- c("fusion_peptide_type")
      if ("spans_junction" %in% colnames(fusion_peptides)) {
        fusion_annotation_cols <- c(fusion_annotation_cols, "spans_junction")
      }
    }
    
    fusion_heatmap <- create_peptide_heatmap(
      fusion_peptides,
      value_cols = intensity_cols,
      is_intensity = TRUE,
      log_transform = TRUE,
      peptide_col = "Peptide",
      annotation_cols = fusion_annotation_cols,
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
    if (exists("public_neoantigens") && "from_fusion" %in% colnames(public_neoantigens) && any(public_neoantigens$from_fusion)) {
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
        
        # Build annotation columns
        fusion_neo_annotation <- c()
        if ("fusion_peptide_type" %in% colnames(fusion_public_neoantigens)) {
          fusion_neo_annotation <- c(fusion_neo_annotation, "fusion_peptide_type")
        }
        if ("spans_junction" %in% colnames(fusion_public_neoantigens)) {
          fusion_neo_annotation <- c(fusion_neo_annotation, "spans_junction")
        }
        if ("public_neoantigen_classification" %in% colnames(fusion_public_neoantigens)) {
          fusion_neo_annotation <- c(fusion_neo_annotation, "public_neoantigen_classification")
        }
        
        fusion_omics_heatmap <- create_peptide_heatmap(
          fusion_public_neoantigens,
          value_cols = fusion_omics_cols,
          is_intensity = FALSE,
          log_transform = FALSE,
          peptide_col = "Peptide",
          annotation_cols = fusion_neo_annotation,
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
  
  # h) Public neoantigen visualizations (if available)
  if (exists("public_neoantigens") && nrow(public_neoantigens) > 0) {
    # Create barplot of public neoantigen tiers
    if ("public_neoantigen_classification" %in% colnames(public_neoantigens)) {
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
    }
    
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
      
      # Create row annotations if the classification column exists
      row_annotation <- NULL
      if ("public_neoantigen_classification" %in% colnames(public_neoantigens)) {
        row_annotation <- data.frame(
          Classification = public_neoantigens$public_neoantigen_classification,
          row.names = rownames(neoantigen_matrix)
        )
      }
      
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
  
  # i) Run sanity checks
  sanity_results <- create_sanity_check_report(
    multi_omics_data,
    tumor_id = config$tumor_id,
    normal_id = config$normal_id
  )
  
  # Save sanity check results to a text file
  sanity_report_file <- file.path(viz_dir, "sanity_check_report.txt")
  sink(sanity_report_file)
  cat("----- SANITY CHECK REPORT -----\n\n")
  
  cat("PASSED CHECKS:\n")
  for(passed in sanity_results$passed) {
    cat("✓ ", passed, "\n")
  }
  
  cat("\nFAILED CHECKS:\n")
  if(length(sanity_results$failed) == 0) {
    cat("No failed checks - all tests passed!\n")
  } else {
    for(check_name in names(sanity_results$failed)) {
      cat("✗ Failed check: ", check_name, "\n")
      print(sanity_results$failed[[check_name]])
      cat("\n")
    }
  }
  
  cat("\nBASIC STATISTICS:\n")
  print(sanity_results$stats)
  sink()
  
  cat("Sanity check report saved to:", sanity_report_file, "\n")
  
  # 8. Generate interactive visualizations if requested
  if (config$generate_interactive) {
    cat("\n## 8. Creating interactive visualizations...\n")
    
    # Try to create interactive plots without requiring Pandoc
    tryCatch({
      # Volcano plot
      if (exists("volcano_plot")) {
        # Simplify plot by removing geom_text_repel
        volcano_plot_simple <- volcano_plot + 
          theme(legend.position = "right")
        
        # For complicated plots like volcano with text repel, convert manually
        if ("ggplot" %in% class(volcano_plot_simple)) {
          # Remove geom_text_repel layers that cause problems with plotly
          layers_to_keep <- sapply(volcano_plot_simple$layers, function(l) {
            !("GeomTextRepel" %in% class(l$geom))
          })
          
          if (any(!layers_to_keep)) {
            volcano_plot_simple$layers <- volcano_plot_simple$layers[layers_to_keep]
          }
        }
        
        # Convert to plotly and export
        volcano_interactive <- ggplotly(volcano_plot_simple)
        
        # Use htmlwidgets with selfcontained = FALSE (doesn't require Pandoc)
        htmlwidgets::saveWidget(
          volcano_interactive, 
          file.path(viz_dir, "interactive_volcano.html"),
          selfcontained = FALSE
        )
        
        cat("Created interactive volcano plot\n")
      }
      
      # Fold change histogram
      if (exists("hist_plot")) {
        hist_interactive <- ggplotly(hist_plot)
        
        htmlwidgets::saveWidget(
          hist_interactive, 
          file.path(viz_dir, "interactive_fold_change_histogram.html"),
          selfcontained = FALSE
        )
        
        cat("Created interactive fold change histogram\n")
      }
      
      # Detection status plot
      if (exists("detection_plot")) {
        detection_interactive <- ggplotly(detection_plot)
        
        htmlwidgets::saveWidget(
          detection_interactive, 
          file.path(viz_dir, "interactive_detection_status.html"),
          selfcontained = FALSE
        )
        
        cat("Created interactive detection status plot\n")
      }
      
      # Create a simple HTML index file that links to all the interactive plots
      index_html <- paste0(
        "<!DOCTYPE html>
      <html>
      <head>
        <title>Interactive Visualizations</title>
        <style>
          body { font-family: Arial, sans-serif; margin: 20px; }
          h1 { text-align: center; }
          .links { margin: 20px; }
          .links a { display: block; margin: 10px 0; }
        </style>
      </head>
      <body>
        <h1>Interactive Visualizations</h1>
        <div class='links'>")
      
      # Add links to each interactive plot
      if (file.exists(file.path(viz_dir, "interactive_volcano.html"))) {
        index_html <- paste0(index_html, 
                             "<a href='interactive_volcano.html'>Tumor vs Normal Volcano Plot</a>")
      }
      
      if (file.exists(file.path(viz_dir, "interactive_fold_change_histogram.html"))) {
        index_html <- paste0(index_html, 
                             "<a href='interactive_fold_change_histogram.html'>Peptide Fold Change Distribution</a>")
      }
      
      if (file.exists(file.path(viz_dir, "interactive_detection_status.html"))) {
        index_html <- paste0(index_html, 
                             "<a href='interactive_detection_status.html'>Peptide Detection Status</a>")
      }
      
      # Close the HTML
      index_html <- paste0(index_html, 
                           "</div>
      <footer>
        <p>Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</p>
      </footer>
      </body>
      </html>")
      
      # Write to file
      writeLines(index_html, file.path(viz_dir, "interactive_index.html"))
      cat("Created interactive visualization index at:", file.path(viz_dir, "interactive_index.html"), "\n")
      
    }, error = function(e) {
      cat("Warning: Error creating interactive visualizations:", conditionMessage(e), "\n")
      cat("Continuing with the analysis...\n")
    })
  }
  
  # 9. Generate Excel reports
  cat("\n## 9. Generating Excel reports...\n")
  
  # Prepare Excel data sheets
  main_excel_sheets <- prepare_excel_data(multi_omics_data, analysis_type = "tumor_normal")
  
  # Merge with any existing excel_sheets
  if (!exists("excel_sheets")) {
    excel_sheets <- list()
  }
  
  for (sheet_name in names(main_excel_sheets)) {
    excel_sheets[[sheet_name]] <- main_excel_sheets[[sheet_name]]
  }
  
  # Add specific sheets for public neoantigens
  if (exists("public_neoantigens") && nrow(public_neoantigens) > 0) {
    excel_sheets[["Public_Neoantigens"]] <- public_neoantigens
    
    # Add filtered sheets if classification exists
    if ("public_neoantigen_classification" %in% colnames(public_neoantigens)) {
      # Get unique classifications
      classifications <- unique(public_neoantigens$public_neoantigen_classification)
      
      # Create separate sheets for each classification
      for (classification in classifications) {
        sheet_name <- paste0("Neoantigen_", gsub("[^A-Za-z0-9]", "_", classification))
        excel_sheets[[sheet_name]] <- public_neoantigens %>%
          filter(public_neoantigen_classification == classification)
      }
    }
  }
  
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