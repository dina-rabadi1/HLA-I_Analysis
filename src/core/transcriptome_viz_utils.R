#' transcriptome_viz_utils.R

#' Create a scatter plot of peptide intensity vs transcriptome expression
#' 
#' @param integrated_data Data frame with peptide intensity and transcriptome data
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @param sample_id Sample ID for the plot
#' @return Path to saved visualization file
create_intensity_expression_scatter <- function(integrated_data, viz_dir, config, sample_id) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Check if we have intensity data
  if (!"Intensity" %in% colnames(integrated_data)) {
    # Try to find intensity column - common alternatives
    intensity_col <- NULL
    for (col in c("Intensity", "total_intensity", "Total Intensity", "intensity")) {
      if (col %in% colnames(integrated_data)) {
        intensity_col <- col
        break
      }
    }
    
    if (is.null(intensity_col)) {
      warning("No intensity column found for scatter plot")
      return(NULL)
    }
    
    # Rename to standardized column name
    integrated_data$Intensity <- integrated_data[[intensity_col]]
  }
  
  # Check if we have expression data
  if (!"sample_expr" %in% colnames(integrated_data)) {
    warning("No expression data found for scatter plot")
    return(NULL)
  }
  
  # Filter out rows with NA or zero in either column
  plot_data <- integrated_data %>%
    dplyr::filter(!is.na(sample_expr), !is.na(Intensity), Intensity > 0)
  
  if (nrow(plot_data) == 0) {
    warning("No valid data available for intensity-expression scatter plot")
    return(NULL)
  }
  
  # Create the plot
  p <- ggplot2::ggplot(plot_data, 
                       ggplot2::aes(x = sample_expr, y = Intensity)) +
    ggplot2::geom_point(alpha = 0.7, color = "steelblue") +
    # Add smooth trend line
    ggplot2::geom_smooth(method = "loess", se = TRUE, color = "darkred", alpha = 0.2) +
    # Use log scale for intensity
    ggplot2::scale_y_log10() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Peptide Intensity vs Gene Expression in", sample_id),
      x = "Gene Expression Level",
      y = "Peptide Intensity (log scale)"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold")
    )
  
  # Create file paths for each format
  file_paths <- list()
  
  # Base file name
  base_name <- paste0(gsub("[^a-zA-Z0-9]", "_", sample_id), "_intensity_vs_expression")
  
  # Save in each requested format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0(base_name, ".", format))
    
    if (format == "pdf") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height)
    } else if (format == "png") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, dpi = png_dpi)
    }
    
    file_paths[[format]] <- file_path
  }
  
  # Return the file paths
  return(file_paths)
}

#' Create a boxplot comparing expression levels of shared vs unique peptides
#' 
#' @param peptide_data Data frame with peptide and gene expression data
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization file
create_shared_vs_unique_expression_boxplot <- function(peptide_data, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Make sure we have the required columns
  if (!all(c("sample_count", "sample_expr") %in% colnames(peptide_data))) {
    warning("Missing required columns for shared vs unique expression boxplot")
    return(NULL)
  }
  
  # Create sharing category
  plot_data <- peptide_data %>%
    dplyr::mutate(
      sharing_category = dplyr::case_when(
        sample_count == 1 ~ "Unique (1 sample)",
        TRUE ~ "Shared (2+ samples)"
      )
    )
  
  # Create the plot
  p <- ggplot2::ggplot(plot_data, 
                       ggplot2::aes(x = sharing_category, y = sample_expr, fill = sharing_category)) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Gene Expression Levels: Shared vs Unique Peptides",
      x = "",
      y = "Gene Expression Level",
      fill = "Peptide Category"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold")
    )
  
  # Add significance test if we have enough data
  if (sum(plot_data$sharing_category == "Unique (1 sample)") > 3 && 
      sum(plot_data$sharing_category == "Shared (2+ samples)") > 3) {
    # Test for significance
    t_test_result <- t.test(
      sample_expr ~ sharing_category, 
      data = plot_data
    )
    
    p_value <- t_test_result$p.value
    
    # Add p-value annotation
    p <- p + 
      ggplot2::labs(
        subtitle = paste0("p-value: ", format.pval(p_value, digits = 3))
      )
  }
  
  # Create file paths for each format
  file_paths <- list()
  
  # Save in each requested format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0("shared_vs_unique_expression.", format))
    
    if (format == "pdf") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height)
    } else if (format == "png") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, dpi = png_dpi)
    }
    
    file_paths[[format]] <- file_path
  }
  
  # Return the file paths
  return(file_paths)
}

#' Create a correlation heatmap between samples using transcriptome data
#' 
#' @param transcriptome_data Data frame with transcriptome data across samples
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization file
create_transcriptome_correlation_heatmap <- function(transcriptome_data, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Identify sample columns
  sample_cols <- grep("^RU", colnames(transcriptome_data), value = TRUE)
  
  if (length(sample_cols) < 2) {
    warning("Not enough sample columns for correlation heatmap")
    return(NULL)
  }
  
  # Extract expression matrix
  expr_matrix <- transcriptome_data %>%
    dplyr::select(dplyr::all_of(sample_cols)) %>%
    as.matrix()
  
  # Remove rows with NA or 0 in all samples
  valid_rows <- rowSums(!is.na(expr_matrix) & expr_matrix > 0) > 0
  expr_matrix <- expr_matrix[valid_rows, ]
  
  # Calculate correlation matrix
  cor_matrix <- cor(expr_matrix, method = "spearman", use = "pairwise.complete.obs")
  
  # File paths for each format
  file_paths <- list()
  
  # Create heatmap for each format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0("transcriptome_correlation_heatmap.", format))
    
    if (format == "pdf") {
      pdf(file_path, width = plot_width, height = plot_height)
    } else if (format == "png") {
      png(file_path, width = plot_width * png_dpi, height = plot_height * png_dpi, 
          res = png_dpi)
    }
    
    # Create the heatmap
    pheatmap::pheatmap(
      cor_matrix,
      main = "Sample Correlation Based on Transcriptome Data",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-1, 1, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      display_numbers = TRUE,
      number_format = "%.2f",
      fontsize_number = 8,
      fontsize_row = 10,
      fontsize_col = 10
    )
    
    dev.off()
    
    file_paths[[format]] <- file_path
  }
  
  return(file_paths)
}

#' Create a dendrogram of samples based on transcriptome data
#' 
#' @param transcriptome_data Data frame with transcriptome data across samples
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization file
create_transcriptome_dendrogram <- function(transcriptome_data, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Identify sample columns
  sample_cols <- grep("^RU", colnames(transcriptome_data), value = TRUE)
  
  if (length(sample_cols) < 2) {
    warning("Not enough sample columns for dendrogram")
    return(NULL)
  }
  
  # Extract expression matrix
  expr_matrix <- transcriptome_data %>%
    dplyr::select(dplyr::all_of(sample_cols)) %>%
    as.matrix()
  
  # Remove rows with NA or 0 in all samples
  valid_rows <- rowSums(!is.na(expr_matrix) & expr_matrix > 0) > 0
  expr_matrix <- expr_matrix[valid_rows, ]
  
  # Calculate distance matrix
  dist_matrix <- dist(t(expr_matrix), method = "euclidean")
  
  # Calculate hierarchical clustering
  hc <- hclust(dist_matrix, method = "ward.D2")
  
  # File paths for each format
  file_paths <- list()
  
  # Create dendrogram for each format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0("transcriptome_dendrogram.", format))
    
    if (format == "pdf") {
      pdf(file_path, width = plot_width, height = plot_height)
    } else if (format == "png") {
      png(file_path, width = plot_width * png_dpi, height = plot_height * png_dpi, 
          res = png_dpi)
    }
    
    # Plot dendrogram
    plot(hc, main = "Sample Clustering Based on Transcriptome Data",
         xlab = "", sub = "", cex = 0.9)
    
    # Add colored rectangles for tumor/normal
    cols <- rep("black", length(sample_cols))
    cols[grep("T", sample_cols)] <- "red"   # Tumor samples
    cols[grep("N", sample_cols)] <- "blue"  # Normal samples
    
    # Draw colored labels
    text_offset <- par("usr")[3] - 0.1 * diff(par("usr")[3:4])
    text(1:length(sample_cols), rep(text_offset, length(sample_cols)), 
         labels = sample_cols, srt = 45, col = cols, xpd = TRUE)
    
    # Add legend
    legend("topright", 
           legend = c("Tumor", "Normal"),
           fill = c("red", "blue"),
           bty = "n")
    
    dev.off()
    
    file_paths[[format]] <- file_path
  }
  
  return(file_paths)
}


#' Generate all visualizations for transcriptome analysis
#' 
#' @param transcriptome_results List of transcriptome analysis results
#' @param viz_dir Directory to save visualizations
#' @param config Configuration settings
#' @return List of paths to all generated visualizations
generate_transcriptome_visualizations <- function(transcriptome_results, viz_dir, config) {
  cat("\nGenerating transcriptome visualizations...\n")
  
  # Make sure output directory exists
  if (!dir.exists(viz_dir)) {
    dir.create(viz_dir, recursive = TRUE)
  }
  
  # Initialize list to store visualization paths
  viz_paths <- list()
  
  # 1. Generate paired tumor-normal visualizations if available
  if (!is.null(transcriptome_results$paired_analysis)) {
    paired_viz <- list()
    for (pair_name in names(transcriptome_results$paired_analysis)) {
      pair_data <- transcriptome_results$paired_analysis[[pair_name]]
      
      # Create volcano plot
      tryCatch({
        paired_viz[[paste0(pair_name, "_volcano")]] <- create_transcriptome_volcano(
          pair_data, viz_dir, config, pair_name
        )
        cat("Created volcano plot for", pair_name, "\n")
      }, error = function(e) {
        cat("Error creating volcano plot for", pair_name, ":", e$message, "\n")
      })
    }
    viz_paths$paired <- paired_viz
  }
  
  # 2. Generate sample-specific visualizations if available
  if (!is.null(transcriptome_results$sample_analysis) && 
      !is.null(transcriptome_results$sample_analysis$sample_specific)) {
    
    sample_viz <- list()
    for (sample_id in names(transcriptome_results$sample_analysis$sample_specific)) {
      sample_data <- transcriptome_results$sample_analysis$sample_specific[[sample_id]]
      
      # Expression level bar chart
      if (!is.null(sample_data$summary) && !is.null(sample_data$summary$expression_levels)) {
        tryCatch({
          sample_viz[[paste0(sample_id, "_expression_levels")]] <- create_expression_level_bar(
            sample_data$summary$expression_levels,
            viz_dir,
            config,
            paste("Expression Levels in", sample_id)
          )
          cat("Created expression level bar chart for", sample_id, "\n")
        }, error = function(e) {
          cat("Error creating expression level bar chart for", sample_id, ":", e$message, "\n")
        })
      }
    }
    viz_paths$sample_specific <- sample_viz
  }
  
  # 3. Generate shared peptide visualizations if available
  if (!is.null(transcriptome_results$sample_analysis) && 
      !is.null(transcriptome_results$sample_analysis$shared_analysis) &&
      !is.null(transcriptome_results$sample_analysis$shared_analysis$threshold_results)) {
    
    shared_viz <- list()
    
    # Get threshold results
    threshold_results <- transcriptome_results$sample_analysis$shared_analysis$threshold_results
    
    for (threshold_name in names(threshold_results)) {
      threshold_data <- threshold_results[[threshold_name]]
      threshold_value <- as.numeric(gsub("shared_", "", threshold_name))
      
      # Expression levels bar chart for shared peptides
      if (!is.null(threshold_data$summary) && !is.null(threshold_data$summary$expression_levels)) {
        tryCatch({
          shared_viz[[paste0(threshold_name, "_levels")]] <- create_expression_level_bar(
            threshold_data$summary$expression_levels,
            viz_dir,
            config,
            paste("Expression Levels of Peptides Shared in ≥", threshold_value, "Samples")
          )
          cat("Created expression level bar chart for shared peptides (threshold:", threshold_value, ")\n")
        }, error = function(e) {
          cat("Error creating expression level bar chart for shared peptides:", e$message, "\n")
        })
      }
    }
    viz_paths$shared_peptides <- shared_viz
  }
  
  return(viz_paths)
}

#' Create a volcano plot of transcriptome data
#' 
#' @param diff_expr Differential expression data frame
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @param pair_name Name of tumor-normal pair (for file naming)
#' @return Path to saved visualization file
create_transcriptome_volcano <- function(diff_expr, viz_dir, config, pair_name = NULL) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Create volcano plot
  p <- ggplot2::ggplot(diff_expr, 
                       ggplot2::aes(x = log2_fold_change, 
                                    y = 1, 
                                    color = expression_status)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +
    ggplot2::scale_color_manual(values = c("Up in Tumor (FC > 2)" = "red", 
                                           "Down in Tumor (FC < 0.5)" = "blue", 
                                           "Similar (-1 < log2FC < 1)" = "gray")) +
    ggplot2::scale_y_continuous(breaks = NULL) +  # Remove y-axis ticks since we don't have p-values
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = ifelse(is.null(pair_name), 
                     "Transcriptome Volcano Plot", 
                     paste(pair_name, "Transcriptome Volcano Plot")),
      subtitle = "Red: Upregulated in Tumor, Blue: Downregulated in Tumor",
      x = "Log2 Fold Change (Tumor/Normal)",
      y = "Gene density",
      color = "Expression Status"
    ) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12)
    )
  
  # Create file paths for each format
  file_paths <- list()
  
  # Base file name
  base_name <- ifelse(is.null(pair_name), 
                      "transcriptome_volcano", 
                      paste0(gsub("[^a-zA-Z0-9]", "_", pair_name), "_volcano"))
  
  # Save in each requested format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0(base_name, ".", format))
    
    if (format == "pdf") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height)
    } else if (format == "png") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, dpi = png_dpi)
    }
    
    file_paths[[format]] <- file_path
  }
  
  # Return the file paths
  return(file_paths)
}

#' Create a bar chart of transcriptome expression levels for peptide-associated genes
#' 
#' @param expression_data Data frame with expression level categories
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @param title Plot title
#' @return Path to saved visualization file
create_expression_level_bar <- function(expression_data, viz_dir, config, title = NULL) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Default title if not provided
  if (is.null(title)) {
    title <- "Expression Levels of Peptide-Associated Genes"
  }
  
  # Create the plot
  p <- ggplot2::ggplot(expression_data, 
                       ggplot2::aes(x = expression_level, y = n, fill = expression_level)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = paste0(n, " (", round(percentage, 1), "%)")), 
                       vjust = -0.5) +
    ggplot2::scale_fill_brewer(palette = "Set3") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = "Expression Level",
      y = "Number of Peptides",
      fill = "Expression Level"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Create file paths for each format
  file_paths <- list()
  
  # Base file name (create from title)
  base_name <- gsub("[^a-zA-Z0-9]", "_", tolower(substr(title, 1, 20)))
  base_name <- paste0("expression_levels_", base_name)
  
  # Save in each requested format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0(base_name, ".", format))
    
    if (format == "pdf") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height)
    } else if (format == "png") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, dpi = png_dpi)
    }
    
    file_paths[[format]] <- file_path
  }
  
  # Return the file paths
  return(file_paths)
}

#' Create a heatmap of public neoantigens
#' 
#' @param neoantigens Data frame with public neoantigen data
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization file
create_public_neoantigens_heatmap <- function(neoantigens, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Check if we have data to plot
  if (nrow(neoantigens) == 0) {
    warning("No public neoantigens to plot")
    return(NULL)
  }
  
  # Determine which omics data sources we have
  has_transcriptome <- "log2_fold_change_transcriptome" %in% colnames(neoantigens)
  has_lfq <- "log2_fold_change_lfq" %in% colnames(neoantigens)
  has_tmt <- "log2_fold_change_tmt" %in% colnames(neoantigens)
  
  # Collect available fold change columns
  fc_cols <- c()
  if (has_transcriptome) fc_cols <- c(fc_cols, "log2_fold_change_transcriptome")
  if (has_lfq) fc_cols <- c(fc_cols, "log2_fold_change_lfq")
  if (has_tmt) fc_cols <- c(fc_cols, "log2_fold_change_tmt")
  
  # Need at least one fold change column
  if (length(fc_cols) == 0) {
    warning("No fold change columns found for public neoantigens heatmap")
    return(NULL)
  }
  
  # Create a matrix for the heatmap
  heatmap_data <- neoantigens %>%
    dplyr::select(primary_gene, Peptide, dplyr::all_of(fc_cols))
  
  # Create row labels
  row_labels <- paste0(heatmap_data$primary_gene, " (", heatmap_data$Peptide, ")")
  
  # Create the matrix for visualization
  heatmap_matrix <- heatmap_data %>%
    dplyr::select(dplyr::all_of(fc_cols)) %>%
    as.matrix()
  
  # Replace NA with 0 for visualization (this is just for display)
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  
  # Make sure row names are unique by adding a sequence number for duplicates
  make_unique_rownames <- function(labels) {
    result <- character(length(labels))
    counts <- table(labels)
    for (name in names(counts)) {
      if (counts[name] == 1) {
        # If there's only one occurrence, use the original name
        result[labels == name] <- name
      } else {
        # If there are multiple occurrences, add a counter
        counter <- 1
        for (i in which(labels == name)) {
          result[i] <- paste0(name, "_", counter)
          counter <- counter + 1
        }
      }
    }
    return(result)
  }
  
  # Create unique row names
  unique_rownames <- make_unique_rownames(row_labels)
  rownames(heatmap_matrix) <- unique_rownames
  
  # Create row annotations
  row_annotation <- data.frame(
    Classification = neoantigens$public_neoantigen_classification,
    row.names = unique_rownames
  )
  
  # Column labels
  col_labels <- c()
  if (has_transcriptome) col_labels <- c(col_labels, "Transcriptome")
  if (has_lfq) col_labels <- c(col_labels, "LFQ Proteome")
  if (has_tmt) col_labels <- c(col_labels, "TMT Proteome")
  
  # File paths for each format
  file_paths <- list()
  
  # Create heatmap for each format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0("public_neoantigens_heatmap.", format))
    
    if (format == "pdf") {
      pdf(file_path, width = plot_width, height = max(8, nrow(heatmap_matrix)/3))
    } else if (format == "png") {
      png(file_path, width = plot_width * png_dpi, height = max(800, nrow(heatmap_matrix)*40), 
          res = png_dpi)
    }
    
    # Create the heatmap
    pheatmap::pheatmap(
      heatmap_matrix,
      main = "Potential Public Neoantigens: Log2 Fold Changes Across Omics Datasets",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-3, 3, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_row = row_annotation,
      display_numbers = TRUE,
      number_format = "%.1f",
      fontsize_row = 8,
      fontsize_col = 10,
      labels_col = col_labels
    )
    
    dev.off()
    
    file_paths[[format]] <- file_path
  }
  
  return(file_paths)
}

#' Create a bar chart showing expression correlation patterns across datasets
#' 
#' @param integrated_data Data frame with fold changes from multiple omics datasets
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization file
create_expression_correlation_plot <- function(integrated_data, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Initialize file paths list
  file_paths <- list()
  
  # Check which omics data we have
  has_transcriptome <- "log2_fold_change_transcriptome" %in% colnames(integrated_data)
  has_lfq <- "log2_fold_change_lfq" %in% colnames(integrated_data)
  has_tmt <- "log2_fold_change_tmt" %in% colnames(integrated_data)
  
  # Need to have at least immunopeptidome + one other dataset
  if (!has_transcriptome && !has_lfq && !has_tmt) {
    return(NULL)
  }
  
  # Debug: Print the column names and first few rows
  cat("Columns in integrated_data:", paste(colnames(integrated_data), collapse=", "), "\n")
  
  # Count rows that have data for both datasets
  if (has_transcriptome) {
    cat("Immunopeptidome vs. Transcriptome: ", 
        sum(!is.na(integrated_data$log2_fold_change_immuno) & 
              !is.na(integrated_data$log2_fold_change_transcriptome)), " rows\n")
  }
  if (has_lfq) {
    cat("Immunopeptidome vs. LFQ: ", 
        sum(!is.na(integrated_data$log2_fold_change_immuno) & 
              !is.na(integrated_data$log2_fold_change_lfq)), " rows\n")
  }
  if (has_tmt) {
    cat("Immunopeptidome vs. TMT: ", 
        sum(!is.na(integrated_data$log2_fold_change_immuno) & 
              !is.na(integrated_data$log2_fold_change_tmt)), " rows\n")
  }
  
  # Create correlation categories
  if (has_transcriptome) {
    # Check if we have immunopeptidome fold change data
    if ("log2_fold_change_immuno" %in% colnames(integrated_data)) {
      # Generate transcriptome correlation
      integrated_data <- integrated_data %>%
        dplyr::mutate(
          immuno_transcriptome_pattern = dplyr::case_when(
            is.na(log2_fold_change_transcriptome) ~ "No transcriptome data",
            log2_fold_change_immuno > 1 & log2_fold_change_transcriptome > 1 ~ "Up in both",
            log2_fold_change_immuno < -1 & log2_fold_change_transcriptome < -1 ~ "Down in both",
            log2_fold_change_immuno > 1 & log2_fold_change_transcriptome < -1 ~ "Up in immunopeptidome, down in transcriptome",
            log2_fold_change_immuno < -1 & log2_fold_change_transcriptome > 1 ~ "Down in immunopeptidome, up in transcriptome",
            TRUE ~ "No significant change"
          )
        )
    } else {
      # Alternative approach if we don't have explicit immunopeptidome fold change
      cat("No log2_fold_change_immuno column found. Using simplified correlation analysis.\n")
      integrated_data <- integrated_data %>%
        dplyr::mutate(
          transcriptome_status = dplyr::case_when(
            is.na(log2_fold_change_transcriptome) ~ "No data",
            log2_fold_change_transcriptome > 1 ~ "Upregulated",
            log2_fold_change_transcriptome < -1 ~ "Downregulated",
            TRUE ~ "Unchanged"
          )
        )
      
      # Count by status
      trans_data <- integrated_data %>%
        dplyr::filter(!is.na(log2_fold_change_transcriptome)) %>%
        dplyr::count(transcriptome_status) %>%
        dplyr::mutate(
          percentage = n / sum(n) * 100,
          dataset = "Transcriptome"
        )
    }
    
    # Create transcriptome correlation plot
    trans_data <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_transcriptome)) %>%
      dplyr::count(immuno_transcriptome_pattern) %>%
      dplyr::mutate(
        percentage = n / sum(n) * 100,
        dataset = "Transcriptome"
      )
    
    # Debug: Print the prepared data
    cat("Transcriptome correlation data:\n")
    print(trans_data)
    
    # Create bar plot
    p_trans <- ggplot2::ggplot(trans_data, 
                               ggplot2::aes(x = reorder(immuno_transcriptome_pattern, -n), 
                                            y = n, fill = immuno_transcriptome_pattern)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_text(ggplot2::aes(label = paste0(n, " (", round(percentage, 1), "%)")), 
                         vjust = -0.5) +
      ggplot2::scale_fill_brewer(palette = "Set3") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Correlation Patterns: Immunopeptidome vs Transcriptome",
        x = "",
        y = "Number of Peptides",
        fill = "Pattern"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    # Make sure the directory exists
    if (!dir.exists(viz_dir)) {
      dir.create(viz_dir, recursive = TRUE)
    }
    
    # Save the plot
    for (format in output_formats) {
      file_path <- file.path(viz_dir, paste0("immuno_transcriptome_correlation.", format))
      cat("Saving correlation plot to:", file_path, "\n")
      
      if (format == "pdf") {
        ggplot2::ggsave(file_path, p_trans, width = plot_width, height = plot_height)
      } else if (format == "png") {
        ggplot2::ggsave(file_path, p_trans, width = plot_width, height = plot_height, dpi = png_dpi)
      }
      
      file_paths[[paste0("transcriptome_", format)]] <- file_path
    }
  }
  
  # Generate LFQ correlation plot
  if (has_lfq) {
    integrated_data <- integrated_data %>%
      dplyr::mutate(
        immuno_lfq_pattern = dplyr::case_when(
          is.na(log2_fold_change_lfq) ~ "No LFQ data",
          log2_fold_change_immuno > 1 & log2_fold_change_lfq > 1 ~ "Up in both",
          log2_fold_change_immuno < -1 & log2_fold_change_lfq < -1 ~ "Down in both",
          log2_fold_change_immuno > 1 & log2_fold_change_lfq < -1 ~ "Up in immunopeptidome, down in LFQ",
          log2_fold_change_immuno < -1 & log2_fold_change_lfq > 1 ~ "Down in immunopeptidome, up in LFQ",
          TRUE ~ "No significant change"
        )
      )
    
    # Create LFQ correlation plot
    lfq_data <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_lfq)) %>%
      dplyr::count(immuno_lfq_pattern) %>%
      dplyr::mutate(
        percentage = n / sum(n) * 100,
        dataset = "LFQ Proteome"
      )
    
    # Debug: Print the prepared data
    cat("LFQ correlation data:\n")
    print(lfq_data)
    
    # Create bar plot
    p_lfq <- ggplot2::ggplot(lfq_data, 
                             ggplot2::aes(x = reorder(immuno_lfq_pattern, -n), 
                                          y = n, fill = immuno_lfq_pattern)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_text(ggplot2::aes(label = paste0(n, " (", round(percentage, 1), "%)")), 
                         vjust = -0.5) +
      ggplot2::scale_fill_brewer(palette = "Set3") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Correlation Patterns: Immunopeptidome vs LFQ Proteome",
        x = "",
        y = "Number of Peptides",
        fill = "Pattern"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    # Save the plot
    for (format in output_formats) {
      file_path <- file.path(viz_dir, paste0("immuno_lfq_correlation.", format))
      cat("Saving LFQ correlation plot to:", file_path, "\n")
      
      if (format == "pdf") {
        ggplot2::ggsave(file_path, p_lfq, width = plot_width, height = plot_height)
      } else if (format == "png") {
        ggplot2::ggsave(file_path, p_lfq, width = plot_width, height = plot_height, dpi = png_dpi)
      }
      
      file_paths[[paste0("lfq_", format)]] <- file_path
    }
  }
  
  # Generate TMT correlation plot
  if (has_tmt) {
    integrated_data <- integrated_data %>%
      dplyr::mutate(
        immuno_tmt_pattern = dplyr::case_when(
          is.na(log2_fold_change_tmt) ~ "No TMT data",
          log2_fold_change_immuno > 1 & log2_fold_change_tmt > 1 ~ "Up in both",
          log2_fold_change_immuno < -1 & log2_fold_change_tmt < -1 ~ "Down in both",
          log2_fold_change_immuno > 1 & log2_fold_change_tmt < -1 ~ "Up in immunopeptidome, down in TMT",
          log2_fold_change_immuno < -1 & log2_fold_change_tmt > 1 ~ "Down in immunopeptidome, up in TMT",
          TRUE ~ "No significant change"
        )
      )
    
    # Create TMT correlation plot
    tmt_data <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_tmt)) %>%
      dplyr::count(immuno_tmt_pattern) %>%
      dplyr::mutate(
        percentage = n / sum(n) * 100,
        dataset = "TMT Proteome"
      )
    
    # Debug: Print the prepared data
    cat("TMT correlation data:\n")
    print(tmt_data)
    
    # Create bar plot
    p_tmt <- ggplot2::ggplot(tmt_data, 
                             ggplot2::aes(x = reorder(immuno_tmt_pattern, -n), 
                                          y = n, fill = immuno_tmt_pattern)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_text(ggplot2::aes(label = paste0(n, " (", round(percentage, 1), "%)")), 
                         vjust = -0.5) +
      ggplot2::scale_fill_brewer(palette = "Set3") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Correlation Patterns: Immunopeptidome vs TMT Proteome",
        x = "",
        y = "Number of Peptides",
        fill = "Pattern"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    # Save the plot
    for (format in output_formats) {
      file_path <- file.path(viz_dir, paste0("immuno_tmt_correlation.", format))
      cat("Saving TMT correlation plot to:", file_path, "\n")
      
      if (format == "pdf") {
        ggplot2::ggsave(file_path, p_tmt, width = plot_width, height = plot_height)
      } else if (format == "png") {
        ggplot2::ggsave(file_path, p_tmt, width = plot_width, height = plot_height, dpi = png_dpi)
      }
      
      file_paths[[paste0("tmt_", format)]] <- file_path
    }
  }
  
  # Create a combined visualization if multiple datasets are available
  combined_datasets <- c()
  combined_data <- data.frame()
  
  if (has_transcriptome && exists("trans_data") && nrow(trans_data) > 0) {
    combined_datasets <- c(combined_datasets, "Transcriptome")
    combined_data <- dplyr::bind_rows(combined_data, trans_data)
  }
  
  if (has_lfq && exists("lfq_data") && nrow(lfq_data) > 0) {
    combined_datasets <- c(combined_datasets, "LFQ Proteome")
    combined_data <- dplyr::bind_rows(combined_data, lfq_data)
  }
  
  if (has_tmt && exists("tmt_data") && nrow(tmt_data) > 0) {
    combined_datasets <- c(combined_datasets, "TMT Proteome")
    combined_data <- dplyr::bind_rows(combined_data, tmt_data)
  }
  
  # If we have data from multiple datasets, create a combined visualization
  if (length(combined_datasets) > 1 && nrow(combined_data) > 0) {
    cat("Creating combined visualization for datasets:", paste(combined_datasets, collapse=", "), "\n")
    
    # Create combined plot
    p_combined <- ggplot2::ggplot(combined_data, 
                                  ggplot2::aes(x = reorder(interaction(dataset, immuno_tmt_pattern), -n), 
                                               y = n, fill = dataset)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_text(ggplot2::aes(label = paste0(n, " (", round(percentage, 1), "%)")), 
                         vjust = -0.5, size = 3) +
      ggplot2::scale_fill_brewer(palette = "Set1") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Correlation Patterns Across All Datasets",
        x = "",
        y = "Number of Peptides",
        fill = "Dataset"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "right"
      )
    
    # Save combined plot
    for (format in output_formats) {
      file_path <- file.path(viz_dir, paste0("all_datasets_correlation.", format))
      cat("Saving combined correlation plot to:", file_path, "\n")
      
      if (format == "pdf") {
        ggplot2::ggsave(file_path, p_combined, width = plot_width * 1.5, height = plot_height)
      } else if (format == "png") {
        ggplot2::ggsave(file_path, p_combined, width = plot_width * 1.5, height = plot_height, dpi = png_dpi)
      }
      
      file_paths[[paste0("combined_", format)]] <- file_path
    }
  }
  
  return(file_paths)
}

#' Create scatter plots comparing fold changes between omics datasets
#' 
#' @param integrated_data Data frame with fold changes from multiple omics datasets
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization files
create_omics_scatter_plots <- function(integrated_data, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Initialize file paths list
  file_paths <- list()
  
  # Check which omics data we have
  has_transcriptome <- "log2_fold_change_transcriptome" %in% colnames(integrated_data)
  has_lfq <- "log2_fold_change_lfq" %in% colnames(integrated_data)
  has_tmt <- "log2_fold_change_tmt" %in% colnames(integrated_data)
  
  # Need to have at least two datasets to compare
  datasets_available <- sum(has_transcriptome, has_lfq, has_tmt)
  if (datasets_available < 2) {
    warning("At least 2 omics datasets with fold changes are needed for scatter plots")
    return(NULL)
  }
  
  # Create transcriptome vs LFQ scatter plot
  if (has_transcriptome && has_lfq) {
    # Filter to rows with data in both datasets
    scatter_data <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_transcriptome) & !is.na(log2_fold_change_lfq))
    
    if (nrow(scatter_data) > 0) {
      # Calculate correlation
      correlation <- cor(scatter_data$log2_fold_change_transcriptome, 
                         scatter_data$log2_fold_change_lfq, 
                         method = "spearman", 
                         use = "pairwise.complete.obs")
      
      # Create the plot
      p <- ggplot2::ggplot(scatter_data, 
                           ggplot2::aes(x = log2_fold_change_transcriptome, 
                                        y = log2_fold_change_lfq)) +
        ggplot2::geom_point(alpha = 0.7, color = "steelblue", size = 2) +
        ggplot2::geom_smooth(method = "lm", color = "darkred", alpha = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Gene Expression Correlation: Transcriptome vs LFQ Proteome",
          subtitle = paste0("Spearman correlation: ", round(correlation, 3)),
          x = "Transcriptome Log2 Fold Change",
          y = "LFQ Proteome Log2 Fold Change"
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold")
        )
      
      # Save in each requested format
      for (format in output_formats) {
        file_path <- file.path(viz_dir, paste0("transcriptome_vs_lfq_scatter.", format))
        
        if (format == "pdf") {
          ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height)
        } else if (format == "png") {
          ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, dpi = png_dpi)
        }
        
        file_paths[[paste0("transcriptome_lfq_", format)]] <- file_path
      }
    }
  }
  
  # Create transcriptome vs TMT scatter plot
  if (has_transcriptome && has_tmt) {
    # Filter to rows with data in both datasets
    scatter_data <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_transcriptome) & !is.na(log2_fold_change_tmt))
    
    if (nrow(scatter_data) > 0) {
      # Calculate correlation
      correlation <- cor(scatter_data$log2_fold_change_transcriptome, 
                         scatter_data$log2_fold_change_tmt, 
                         method = "spearman", 
                         use = "pairwise.complete.obs")
      
      # Create the plot
      p <- ggplot2::ggplot(scatter_data, 
                           ggplot2::aes(x = log2_fold_change_transcriptome, 
                                        y = log2_fold_change_tmt)) +
        ggplot2::geom_point(alpha = 0.7, color = "steelblue", size = 2) +
        ggplot2::geom_smooth(method = "lm", color = "darkred", alpha = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Gene Expression Correlation: Transcriptome vs TMT Proteome",
          subtitle = paste0("Spearman correlation: ", round(correlation, 3)),
          x = "Transcriptome Log2 Fold Change",
          y = "TMT Proteome Log2 Fold Change"
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold")
        )
      
      # Save in each requested format
      for (format in output_formats) {
        file_path <- file.path(viz_dir, paste0("transcriptome_vs_tmt_scatter.", format))
        
        if (format == "pdf") {
          ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height)
        } else if (format == "png") {
          ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, dpi = png_dpi)
        }
        
        file_paths[[paste0("transcriptome_tmt_", format)]] <- file_path
      }
    }
  }
  
  # Create LFQ vs TMT scatter plot
  if (has_lfq && has_tmt) {
    # Filter to rows with data in both datasets
    scatter_data <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_lfq) & !is.na(log2_fold_change_tmt))
    
    if (nrow(scatter_data) > 0) {
      # Calculate correlation
      correlation <- cor(scatter_data$log2_fold_change_lfq, 
                         scatter_data$log2_fold_change_tmt, 
                         method = "spearman", 
                         use = "pairwise.complete.obs")
      
      # Create the plot
      p <- ggplot2::ggplot(scatter_data, 
                           ggplot2::aes(x = log2_fold_change_lfq, 
                                        y = log2_fold_change_tmt)) +
        ggplot2::geom_point(alpha = 0.7, color = "steelblue", size = 2) +
        ggplot2::geom_smooth(method = "lm", color = "darkred", alpha = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Gene Expression Correlation: LFQ vs TMT Proteome",
          subtitle = paste0("Spearman correlation: ", round(correlation, 3)),
          x = "LFQ Proteome Log2 Fold Change",
          y = "TMT Proteome Log2 Fold Change"
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold")
        )
      
      # Save in each requested format
      for (format in output_formats) {
        file_path <- file.path(viz_dir, paste0("lfq_vs_tmt_scatter.", format))
        
        if (format == "pdf") {
          ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height)
        } else if (format == "png") {
          ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, dpi = png_dpi)
        }
        
        file_paths[[paste0("lfq_tmt_", format)]] <- file_path
      }
    }
  }
  
  return(file_paths)
}

#' Create a Venn diagram of upregulated genes across omics datasets
#' 
#' @param integrated_data Data frame with fold changes from multiple omics datasets
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization file
#' 
# Check if we can create the Venn diagram
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  warning("VennDiagram package not available. Cannot create Venn diagram.")
  return(NULL)
}

create_upregulated_venn_diagram <- function(integrated_data, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Initialize file paths list
  file_paths <- list()
  
  # Extract genes upregulated in each omics dataset
  genes_immuno_up <- NULL
  if ("log2_fold_change_immuno" %in% colnames(integrated_data)) {
    genes_immuno_up <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_immuno) & log2_fold_change_immuno > 1) %>%
      dplyr::pull(primary_gene) %>%
      unique()
  } else {
    # Extract any upregulated genes using immunopeptide data if log2_fold_change_immuno doesn't exist
    genes_immuno_up <- integrated_data %>%
      dplyr::pull(primary_gene) %>%
      unique()
  }
  
  has_transcriptome <- "log2_fold_change_transcriptome" %in% colnames(integrated_data)
  has_lfq <- "log2_fold_change_lfq" %in% colnames(integrated_data)
  has_tmt <- "log2_fold_change_tmt" %in% colnames(integrated_data)
  
  if (has_transcriptome) {
    genes_trans_up <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > 1) %>%
      dplyr::pull(primary_gene) %>%
      unique()
  } else {
    genes_trans_up <- c()
  }
  
  if (has_lfq) {
    genes_lfq_up <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_lfq) & log2_fold_change_lfq > 1) %>%
      dplyr::pull(primary_gene) %>%
      unique()
  } else {
    genes_lfq_up <- c()
  }
  
  if (has_tmt) {
    genes_tmt_up <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_tmt) & log2_fold_change_tmt > 1) %>%
      dplyr::pull(primary_gene) %>%
      unique()
  } else {
    genes_tmt_up <- c()
  }
  
  # Count datasets with upregulated genes
  omics_count <- sum(length(genes_immuno_up) > 0, 
                     length(genes_trans_up) > 0, 
                     length(genes_lfq_up) > 0, 
                     length(genes_tmt_up) > 0)
  
  if (omics_count < 2) {
    warning("At least 2 omics datasets with upregulated genes are needed for Venn diagram")
    return(NULL)
  }
  
  # Create a base file path for the visualizations
  base_file_path <- file.path(viz_dir, "upregulated_genes_venn")
  for (format in output_formats) {
    file_paths[[format]] <- paste0(base_file_path, ".", format)
  }
  
  # Define colors for the Venn diagram
  venn_colors <- c("steelblue", "firebrick", "forestgreen", "darkorange")
  
  # Create the appropriate Venn diagram based on available datasets
  if (omics_count == 2) {
    # Determine which two datasets to use
    if (length(genes_immuno_up) > 0 && length(genes_trans_up) > 0) {
      category_names <- c("Immunopeptidome", "Transcriptome")
      set1 <- genes_immuno_up
      set2 <- genes_trans_up
    } else if (length(genes_immuno_up) > 0 && length(genes_lfq_up) > 0) {
      category_names <- c("Immunopeptidome", "LFQ Proteome")
      set1 <- genes_immuno_up
      set2 <- genes_lfq_up
    } else if (length(genes_immuno_up) > 0 && length(genes_tmt_up) > 0) {
      category_names <- c("Immunopeptidome", "TMT Proteome")
      set1 <- genes_immuno_up
      set2 <- genes_tmt_up
    } else if (length(genes_trans_up) > 0 && length(genes_lfq_up) > 0) {
      category_names <- c("Transcriptome", "LFQ Proteome")
      set1 <- genes_trans_up
      set2 <- genes_lfq_up
    } else if (length(genes_trans_up) > 0 && length(genes_tmt_up) > 0) {
      category_names <- c("Transcriptome", "TMT Proteome")
      set1 <- genes_trans_up
      set2 <- genes_tmt_up
    } else if (length(genes_lfq_up) > 0 && length(genes_tmt_up) > 0) {
      category_names <- c("LFQ Proteome", "TMT Proteome")
      set1 <- genes_lfq_up
      set2 <- genes_tmt_up
    } else {
      # Fallback - this should not happen given earlier checks
      warning("Could not determine which datasets to use for Venn diagram")
      return(NULL)
    }
    
    # Create 2-way Venn diagram
    venn_file <- file_paths$png
    VennDiagram::venn.diagram(
      x = list(set1, set2),
      category.names = category_names,
      filename = venn_file,
      output = TRUE,
      imagetype = "png",
      height = plot_height * png_dpi,
      width = plot_width * png_dpi,
      resolution = png_dpi,
      compression = "lzw",
      lwd = 2,
      lty = "solid",
      fill = venn_colors[1:2],
      cex = 2,
      fontfamily = "sans",
      cat.cex = 1.5,
      cat.fontfamily = "sans"
    )
    
    # Also create PDF if needed
    if ("pdf" %in% output_formats) {
      venn_file_pdf <- file_paths$pdf
      VennDiagram::venn.diagram(
        x = list(set1, set2),
        category.names = category_names,
        filename = venn_file_pdf,
        output = TRUE,
        imagetype = "pdf",
        height = plot_height,
        width = plot_width,
        lwd = 2,
        lty = "solid",
        fill = venn_colors[1:2],
        cex = 2,
        fontfamily = "sans",
        cat.cex = 1.5,
        cat.fontfamily = "sans"
      )
    }
  } else if (omics_count == 3) {
    # Determine which three datasets to use
    if (length(genes_immuno_up) > 0 && length(genes_trans_up) > 0 && length(genes_lfq_up) > 0) {
      category_names <- c("Immunopeptidome", "Transcriptome", "LFQ Proteome")
      set1 <- genes_immuno_up
      set2 <- genes_trans_up
      set3 <- genes_lfq_up
    } else if (length(genes_immuno_up) > 0 && length(genes_trans_up) > 0 && length(genes_tmt_up) > 0) {
      category_names <- c("Immunopeptidome", "Transcriptome", "TMT Proteome")
      set1 <- genes_immuno_up
      set2 <- genes_trans_up
      set3 <- genes_tmt_up
    } else if (length(genes_immuno_up) > 0 && length(genes_lfq_up) > 0 && length(genes_tmt_up) > 0) {
      category_names <- c("Immunopeptidome", "LFQ Proteome", "TMT Proteome")
      set1 <- genes_immuno_up
      set2 <- genes_lfq_up
      set3 <- genes_tmt_up
    } else if (length(genes_trans_up) > 0 && length(genes_lfq_up) > 0 && length(genes_tmt_up) > 0) {
      category_names <- c("Transcriptome", "LFQ Proteome", "TMT Proteome")
      set1 <- genes_trans_up
      set2 <- genes_lfq_up
      set3 <- genes_tmt_up
    } else {
      # This should not happen given earlier checks
      warning("Could not determine which three datasets to use for Venn diagram")
      return(NULL)
    }
    
    # Create 3-way Venn diagram
    venn_file <- file_paths$png
    VennDiagram::venn.diagram(
      x = list(set1, set2, set3),
      category.names = category_names,
      filename = venn_file,
      output = TRUE,
      imagetype = "png",
      height = plot_height * png_dpi,
      width = plot_width * png_dpi,
      resolution = png_dpi,
      compression = "lzw",
      lwd = 2,
      lty = "solid",
      fill = venn_colors[1:3],
      cex = 2,
      fontfamily = "sans",
      cat.cex = 1.5,
      cat.fontfamily = "sans"
    )
    
    # Also create PDF if needed
    if ("pdf" %in% output_formats) {
      venn_file_pdf <- file_paths$pdf
      VennDiagram::venn.diagram(
        x = list(set1, set2, set3),
        category.names = category_names,
        filename = venn_file_pdf,
        output = TRUE,
        imagetype = "pdf",
        height = plot_height,
        width = plot_width,
        lwd = 2,
        lty = "solid",
        fill = venn_colors[1:3],
        cex = 2,
        fontfamily = "sans",
        cat.cex = 1.5,
        cat.fontfamily = "sans"
      )
    }
  } else if (omics_count == 4) {
    # Use all four datasets
    # Create 4-way Venn diagram
    venn_file <- file_paths$png
    VennDiagram::venn.diagram(
      x = list(genes_immuno_up, genes_trans_up, genes_lfq_up, genes_tmt_up),
      category.names = c("Immunopeptidome", "Transcriptome", "LFQ Proteome", "TMT Proteome"),
      filename = venn_file,
      output = TRUE,
      imagetype = "png",
      height = plot_height * png_dpi,
      width = plot_width * png_dpi,
      resolution = png_dpi,
      compression = "lzw",
      lwd = 2,
      lty = "solid",
      fill = venn_colors,
      cex = 1.5,
      fontfamily = "sans",
      cat.cex = 1.2,
      cat.fontfamily = "sans"
    )
    
    # Also create PDF if needed
    if ("pdf" %in% output_formats) {
      venn_file_pdf <- file_paths$pdf
      VennDiagram::venn.diagram(
        x = list(genes_immuno_up, genes_trans_up, genes_lfq_up, genes_tmt_up),
        category.names = c("Immunopeptidome", "Transcriptome", "LFQ Proteome", "TMT Proteome"),
        filename = venn_file_pdf,
        output = TRUE,
        imagetype = "pdf",
        height = plot_height,
        width = plot_width,
        lwd = 2,
        lty = "solid",
        fill = venn_colors,
        cex = 1.5,
        fontfamily = "sans",
        cat.cex = 1.2,
        cat.fontfamily = "sans"
      )
    }
  }
  
  return(file_paths)
}

#' Create bar chart showing overlap of upregulated genes
#' 
#' @param integrated_data Data frame with fold changes from multiple omics datasets
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization files
create_overlap_bar_chart <- function(integrated_data, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Initialize file paths list
  file_paths <- list()
  
  # Extract genes upregulated in each omics dataset
  genes_immuno_up <- NULL
  if ("log2_fold_change_immuno" %in% colnames(integrated_data)) {
    genes_immuno_up <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_immuno) & log2_fold_change_immuno > 1) %>%
      dplyr::pull(primary_gene) %>%
      unique()
  } else {
    # Extract any peptide genes as immuno genes
    genes_immuno_up <- integrated_data %>%
      dplyr::pull(primary_gene) %>%
      unique()
  }
  
  has_transcriptome <- "log2_fold_change_transcriptome" %in% colnames(integrated_data)
  has_lfq <- "log2_fold_change_lfq" %in% colnames(integrated_data)
  has_tmt <- "log2_fold_change_tmt" %in% colnames(integrated_data)
  
  if (has_transcriptome) {
    genes_trans_up <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_transcriptome) & log2_fold_change_transcriptome > 1) %>%
      dplyr::pull(primary_gene) %>%
      unique()
  } else {
    genes_trans_up <- c()
  }
  
  if (has_lfq) {
    genes_lfq_up <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_lfq) & log2_fold_change_lfq > 1) %>%
      dplyr::pull(primary_gene) %>%
      unique()
  } else {
    genes_lfq_up <- c()
  }
  
  if (has_tmt) {
    genes_tmt_up <- integrated_data %>%
      dplyr::filter(!is.na(log2_fold_change_tmt) & log2_fold_change_tmt > 1) %>%
      dplyr::pull(primary_gene) %>%
      unique()
  } else {
    genes_tmt_up <- c()
  }
  
  # Count datasets
  omics_count <- sum(length(genes_immuno_up) > 0, 
                     length(genes_trans_up) > 0, 
                     length(genes_lfq_up) > 0, 
                     length(genes_tmt_up) > 0)
  
  if (omics_count < 2) {
    warning("At least 2 omics datasets with upregulated genes are needed for overlap chart")
    return(NULL)
  }
  
  # Create overlap summary
  overlap_summary <- data.frame(
    Category = character(),
    Count = integer(),
    Type = character(),
    stringsAsFactors = FALSE
  )
  
  # Helper function to safely calculate overlaps
  safe_intersection <- function(set1, set2) {
    if (length(set1) == 0 || length(set2) == 0) return(character(0))
    return(intersect(set1, set2))
  }
  
  safe_setdiff <- function(set1, set2) {
    if (length(set1) == 0) return(character(0))
    if (length(set2) == 0) return(set1)
    return(setdiff(set1, set2))
  }
  
  # Calculate all possible overlaps based on available datasets
  if (length(genes_immuno_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "Immunopeptidome only",
        Count = length(safe_setdiff(genes_immuno_up, 
                                    union(union(genes_trans_up, genes_lfq_up), genes_tmt_up))),
        Type = "Single dataset",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (length(genes_trans_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "Transcriptome only",
        Count = length(safe_setdiff(genes_trans_up, 
                                    union(union(genes_immuno_up, genes_lfq_up), genes_tmt_up))),
        Type = "Single dataset",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (length(genes_lfq_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "LFQ only",
        Count = length(safe_setdiff(genes_lfq_up, 
                                    union(union(genes_immuno_up, genes_trans_up), genes_tmt_up))),
        Type = "Single dataset",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (length(genes_tmt_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "TMT only",
        Count = length(safe_setdiff(genes_tmt_up, 
                                    union(union(genes_immuno_up, genes_trans_up), genes_lfq_up))),
        Type = "Single dataset",
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Two-way overlaps
  if (length(genes_immuno_up) > 0 && length(genes_trans_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "Immuno + Trans",
        Count = length(safe_setdiff(
          safe_intersection(genes_immuno_up, genes_trans_up),
          union(genes_lfq_up, genes_tmt_up)
        )),
        Type = "Two datasets",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (length(genes_immuno_up) > 0 && length(genes_lfq_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "Immuno + LFQ",
        Count = length(safe_setdiff(
          safe_intersection(genes_immuno_up, genes_lfq_up),
          union(genes_trans_up, genes_tmt_up)
        )),
        Type = "Two datasets",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (length(genes_immuno_up) > 0 && length(genes_tmt_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "Immuno + TMT",
        Count = length(safe_setdiff(
          safe_intersection(genes_immuno_up, genes_tmt_up),
          union(genes_trans_up, genes_lfq_up)
        )),
        Type = "Two datasets",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (length(genes_trans_up) > 0 && length(genes_lfq_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "Trans + LFQ",
        Count = length(safe_setdiff(
          safe_intersection(genes_trans_up, genes_lfq_up),
          union(genes_immuno_up, genes_tmt_up)
        )),
        Type = "Two datasets",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (length(genes_trans_up) > 0 && length(genes_tmt_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "Trans + TMT",
        Count = length(safe_setdiff(
          safe_intersection(genes_trans_up, genes_tmt_up),
          union(genes_immuno_up, genes_lfq_up)
        )),
        Type = "Two datasets",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (length(genes_lfq_up) > 0 && length(genes_tmt_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "LFQ + TMT",
        Count = length(safe_setdiff(
          safe_intersection(genes_lfq_up, genes_tmt_up),
          union(genes_immuno_up, genes_trans_up)
        )),
        Type = "Two datasets",
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Three-way overlaps
  if (length(genes_immuno_up) > 0 && length(genes_trans_up) > 0 && length(genes_lfq_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "Immuno + Trans + LFQ",
        Count = length(safe_setdiff(
          safe_intersection(safe_intersection(genes_immuno_up, genes_trans_up), genes_lfq_up),
          genes_tmt_up
        )),
        Type = "Three datasets",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (length(genes_immuno_up) > 0 && length(genes_trans_up) > 0 && length(genes_tmt_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "Immuno + Trans + TMT",
        Count = length(safe_setdiff(
          safe_intersection(safe_intersection(genes_immuno_up, genes_trans_up), genes_tmt_up),
          genes_lfq_up
        )),
        Type = "Three datasets",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (length(genes_immuno_up) > 0 && length(genes_lfq_up) > 0 && length(genes_tmt_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "Immuno + LFQ + TMT",
        Count = length(safe_setdiff(
          safe_intersection(safe_intersection(genes_immuno_up, genes_lfq_up), genes_tmt_up),
          genes_trans_up
        )),
        Type = "Three datasets",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (length(genes_trans_up) > 0 && length(genes_lfq_up) > 0 && length(genes_tmt_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "Trans + LFQ + TMT",
        Count = length(safe_setdiff(
          safe_intersection(safe_intersection(genes_trans_up, genes_lfq_up), genes_tmt_up),
          genes_immuno_up
        )),
        Type = "Three datasets",
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Four-way overlap
  if (length(genes_immuno_up) > 0 && length(genes_trans_up) > 0 && 
      length(genes_lfq_up) > 0 && length(genes_tmt_up) > 0) {
    overlap_summary <- rbind(
      overlap_summary,
      data.frame(
        Category = "All datasets",
        Count = length(safe_intersection(
          safe_intersection(genes_immuno_up, genes_trans_up),
          safe_intersection(genes_lfq_up, genes_tmt_up)
        )),
        Type = "All datasets",
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Filter out zero counts and sort
  overlap_summary <- overlap_summary %>%
    dplyr::filter(Count > 0) %>%
    dplyr::arrange(dplyr::desc(Count))
  
  # If we have no rows with counts, return NULL
  if(nrow(overlap_summary) == 0) {
    warning("No overlaps found with counts > 0")
    return(NULL)
  }
  
  # Order categories by count
  overlap_summary$Category <- factor(overlap_summary$Category, 
                                     levels = overlap_summary$Category[order(overlap_summary$Count, 
                                                                             decreasing = TRUE)])
  
  # Create bar chart
  p <- ggplot2::ggplot(overlap_summary, 
                       ggplot2::aes(x = Category, y = Count, fill = Type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = Count), vjust = -0.5) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      title = "Overlap of Upregulated Genes Across Omics Datasets",
      x = "",
      y = "Number of Genes",
      fill = "Overlap Type"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(size = 14, face = "bold")
    )
  
  # Save the plots in each format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0("upregulated_genes_overlap.", format))
    
    if (format == "pdf") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height)
    } else if (format == "png") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, dpi = png_dpi)
    }
    
    file_paths[[format]] <- file_path
  }
  
  return(file_paths)
}

#' from here new
#' Create a detection status barplot
#' 
#' @param integrated_data Data frame with integrated multi-omics data
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @param dataset Dataset to create barplot for ("transcriptome", "lfq", or "tmt")
#' @return Path to saved visualization file
create_detection_barplot <- function(integrated_data, viz_dir, config, dataset = "transcriptome") {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Check which dataset to use
  if (dataset == "transcriptome") {
    status_col <- "transcriptome_status"
    title <- "Detected Peptides by Transcriptome Status"
    file_prefix <- "detection_transcriptome"
  } else if (dataset == "lfq") {
    status_col <- "lfq_status"
    title <- "Detected Peptides by LFQ Proteome Status"
    file_prefix <- "detection_lfq"
  } else if (dataset == "tmt") {
    status_col <- "tmt_status"
    title <- "Detected Peptides by TMT Proteome Status"
    file_prefix <- "detection_tmt"
  } else {
    stop("Invalid dataset specified. Must be 'transcriptome', 'lfq', or 'tmt'")
  }
  
  # Check if we have the required column
  if (!status_col %in% colnames(integrated_data)) {
    warning("Required column '", status_col, "' not found in integrated data")
    return(NULL)
  }
  
  # Count by status category
  status_counts <- integrated_data %>%
    dplyr::filter(!is.na(!!dplyr::sym(status_col))) %>%
    dplyr::count(!!dplyr::sym(status_col)) %>%
    dplyr::mutate(
      percentage = n / sum(n) * 100,
      label = paste0(n, " (", round(percentage, 1), "%)")
    )
  
  # Set colors based on dataset
  if (dataset == "transcriptome") {
    status_colors <- c(
      "Up in transcriptome" = "red3",
      "Unchanged in transcriptome" = "grey50",
      "Down in transcriptome" = "blue3"
    )
  } else if (dataset == "lfq") {
    status_colors <- c(
      "Up in LFQ proteome" = "red3",
      "Unchanged in LFQ proteome" = "grey50",
      "Down in LFQ proteome" = "blue3"
    )
  } else if (dataset == "tmt") {
    status_colors <- c(
      "Up in TMT proteome" = "red3",
      "Unchanged in TMT proteome" = "grey50",
      "Down in TMT proteome" = "blue3"
    )
  }
  
  # Create the plot
  status_counts$status <- factor(
    status_counts[[status_col]],
    levels = names(status_colors)
  )
  
  p <- ggplot2::ggplot(status_counts, 
                       ggplot2::aes(x = status, y = n, fill = status)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.5) +
    ggplot2::scale_fill_manual(values = status_colors) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = "",
      y = "Number of Peptides",
      fill = "Status"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(size = 14, face = "bold")
    )
  
  # File paths for each format
  file_paths <- list()
  
  # Save in each requested format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0(file_prefix, ".", format))
    
    if (format == "pdf") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height)
    } else if (format == "png") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, dpi = png_dpi)
    }
    
    file_paths[[format]] <- file_path
  }
  
  return(file_paths)
}

#' Create a multi-omics heatmap for shared peptides
#' 
#' @param integrated_data Data frame with integrated multi-omics data
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization files
create_multi_omics_heatmap <- function(integrated_data, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Determine which omics data sources we have
  has_transcriptome <- "log2_fold_change_transcriptome" %in% colnames(integrated_data)
  has_lfq <- "log2_fold_change_lfq" %in% colnames(integrated_data)
  has_tmt <- "log2_fold_change_tmt" %in% colnames(integrated_data)
  
  # Need at least one omics source
  if (!has_transcriptome && !has_lfq && !has_tmt) {
    warning("No fold change data available for multi-omics heatmap")
    return(NULL)
  }
  
  # Create a matrix for the heatmap
  heatmap_cols <- c()
  if (has_transcriptome) heatmap_cols <- c(heatmap_cols, "log2_fold_change_transcriptome")
  if (has_lfq) heatmap_cols <- c(heatmap_cols, "log2_fold_change_lfq")
  if (has_tmt) heatmap_cols <- c(heatmap_cols, "log2_fold_change_tmt")
  
  # Filter to rows that have at least one non-NA value in the fold change columns
  filtered_data <- integrated_data %>%
    dplyr::filter_at(
      dplyr::vars(dplyr::any_of(heatmap_cols)), 
      dplyr::any_vars(!is.na(.))
    )
  
  if (nrow(filtered_data) == 0) {
    warning("No data available for multi-omics heatmap after filtering")
    return(NULL)
  }
  
  # Sort by number of samples peptide is found in
  filtered_data <- filtered_data %>%
    dplyr::arrange(dplyr::desc(sample_count), primary_gene)
  
  # Create labels for each peptide-gene pair
  row_labels <- paste0(filtered_data$primary_gene, " (", filtered_data$Peptide, ")")
  
  # Limit to top peptides if there are too many
  max_peptides <- 50
  if (length(row_labels) > max_peptides) {
    cat("Limiting heatmap to top", max_peptides, "peptides\n")
    filtered_data <- filtered_data[1:max_peptides, ]
    row_labels <- row_labels[1:max_peptides]
  }
  
  # Create the matrix for visualization
  heatmap_matrix <- filtered_data %>%
    dplyr::select(dplyr::all_of(heatmap_cols)) %>%
    as.matrix()
  
  # Replace NA with 0 for visualization (this is just for display)
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  
  # Set row names
  rownames(heatmap_matrix) <- row_labels
  
  # Set column names
  colnames(heatmap_matrix) <- c(
    if (has_transcriptome) "Transcriptome" else NULL,
    if (has_lfq) "LFQ Proteome" else NULL,
    if (has_tmt) "TMT Proteome" else NULL
  )
  
  # Create row annotations if we have neoantigen classifications
  row_annotation <- NULL
  if ("public_neoantigen_classification" %in% colnames(filtered_data)) {
    row_annotation <- data.frame(
      Classification = filtered_data$public_neoantigen_classification,
      row.names = row_labels
    )
  }
  
  # File paths for each format
  file_paths <- list()
  
  # Create heatmap for each format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0("multi_omics_heatmap.", format))
    
    if (format == "pdf") {
      pdf(file_path, width = plot_width, height = max(8, nrow(heatmap_matrix)/3))
    } else if (format == "png") {
      png(file_path, width = plot_width * png_dpi, height = max(800, nrow(heatmap_matrix)*40), 
          res = png_dpi)
    }
    
    # Create the heatmap
    pheatmap::pheatmap(
      heatmap_matrix,
      main = "Multi-Omics Integration: Log2 Fold Changes",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-3, 3, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_row = row_annotation,
      display_numbers = TRUE,
      number_format = "%.1f",
      fontsize_row = 8,
      fontsize_col = 10
    )
    
    dev.off()
    
    file_paths[[format]] <- file_path
  }
  
  return(file_paths)
}

#' Create a Sankey diagram showing flow between omics platforms
#' 
#' @param integrated_data Data frame with integrated multi-omics data
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization files
create_multi_omics_sankey <- function(integrated_data, viz_dir, config) {
  # Check if we have the required columns for a Sankey diagram
  has_transcriptome <- "transcriptome_status" %in% colnames(integrated_data)
  has_lfq <- "lfq_status" %in% colnames(integrated_data)
  has_tmt <- "tmt_status" %in% colnames(integrated_data)
  
  # Need at least one additional omics platform
  if (!has_transcriptome && !has_lfq && !has_tmt) {
    warning("No omics status columns available for Sankey diagram")
    return(NULL)
  }
  
  # Prepare data for Sankey diagram
  # We'll use the networkD3 package
  # First, ensure we have it installed
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    install.packages("networkD3")
  }
  
  # Determine which datasets to include
  datasets <- c("Immunopeptidome")
  if (has_transcriptome) datasets <- c(datasets, "Transcriptome")
  if (has_lfq) datasets <- c(datasets, "LFQ")
  if (has_tmt) datasets <- c(datasets, "TMT")
  
  # Create nodes data frame
  nodes <- data.frame(
    name = c(
      # Start with detection status
      "Detected",
      
      # Add transcriptome statuses if available
      if (has_transcriptome) c(
        "Up in transcriptome",
        "Unchanged in transcriptome",
        "Down in transcriptome"
      ) else c(),
      
      # Add LFQ statuses if available
      if (has_lfq) c(
        "Up in LFQ proteome",
        "Unchanged in LFQ proteome",
        "Down in LFQ proteome"
      ) else c(),
      
      # Add TMT statuses if available
      if (has_tmt) c(
        "Up in TMT proteome",
        "Unchanged in TMT proteome",
        "Down in TMT proteome"
      ) else c()
    )
  )
  
  # Create links data frame
  # We'll count the flows between different statuses
  links <- data.frame(
    source = integer(),
    target = integer(),
    value = integer()
  )
  
  # Add links from immunopeptidome to transcriptome
  if (has_transcriptome) {
    # Count peptides in each transcriptome status category
    trans_counts <- integrated_data %>%
      dplyr::filter(!is.na(transcriptome_status)) %>%
      dplyr::count(transcriptome_status)
    
    # Add links
    for (i in 1:nrow(trans_counts)) {
      status <- trans_counts$transcriptome_status[i]
      count <- trans_counts$n[i]
      
      # Find index of this status in nodes
      target_idx <- match(status, nodes$name) - 1  # 0-based indexing
      
      # Add link from Detected (0) to this status
      links <- rbind(links, data.frame(
        source = 0,  # Detected
        target = target_idx,
        value = count
      ))
    }
  }
  
  # Add links from transcriptome to LFQ if both are available
  if (has_transcriptome && has_lfq) {
    # Count flows between transcriptome and LFQ statuses
    trans_lfq_counts <- integrated_data %>%
      dplyr::filter(!is.na(transcriptome_status) & !is.na(lfq_status)) %>%
      dplyr::count(transcriptome_status, lfq_status)
    
    # Add links
    for (i in 1:nrow(trans_lfq_counts)) {
      trans_status <- trans_lfq_counts$transcriptome_status[i]
      lfq_status <- trans_lfq_counts$lfq_status[i]
      count <- trans_lfq_counts$n[i]
      
      # Find indices in nodes
      source_idx <- match(trans_status, nodes$name) - 1  # 0-based indexing
      target_idx <- match(lfq_status, nodes$name) - 1  # 0-based indexing
      
      # Add link
      links <- rbind(links, data.frame(
        source = source_idx,
        target = target_idx,
        value = count
      ))
    }
  } else if (has_lfq) {
    # If we have LFQ but not transcriptome, link directly from immunopeptidome
    lfq_counts <- integrated_data %>%
      dplyr::filter(!is.na(lfq_status)) %>%
      dplyr::count(lfq_status)
    
    # Add links
    for (i in 1:nrow(lfq_counts)) {
      status <- lfq_counts$lfq_status[i]
      count <- lfq_counts$n[i]
      
      # Find index of this status in nodes
      target_idx <- match(status, nodes$name) - 1  # 0-based indexing
      
      # Add link from Detected (0) to this status
      links <- rbind(links, data.frame(
        source = 0,  # Detected
        target = target_idx,
        value = count
      ))
    }
  }
  
  # Add links to TMT
  if (has_lfq && has_tmt) {
    # Link from LFQ to TMT
    lfq_tmt_counts <- integrated_data %>%
      dplyr::filter(!is.na(lfq_status) & !is.na(tmt_status)) %>%
      dplyr::count(lfq_status, tmt_status)
    
    # Add links
    for (i in 1:nrow(lfq_tmt_counts)) {
      lfq_status <- lfq_tmt_counts$lfq_status[i]
      tmt_status <- lfq_tmt_counts$tmt_status[i]
      count <- lfq_tmt_counts$n[i]
      
      # Find indices in nodes
      source_idx <- match(lfq_status, nodes$name) - 1  # 0-based indexing
      target_idx <- match(tmt_status, nodes$name) - 1  # 0-based indexing
      
      # Add link
      links <- rbind(links, data.frame(
        source = source_idx,
        target = target_idx,
        value = count
      ))
    }
  } else if (has_transcriptome && has_tmt) {
    # Link from transcriptome to TMT
    trans_tmt_counts <- integrated_data %>%
      dplyr::filter(!is.na(transcriptome_status) & !is.na(tmt_status)) %>%
      dplyr::count(transcriptome_status, tmt_status)
    
    # Add links
    for (i in 1:nrow(trans_tmt_counts)) {
      trans_status <- trans_tmt_counts$transcriptome_status[i]
      tmt_status <- trans_tmt_counts$tmt_status[i]
      count <- trans_tmt_counts$n[i]
      
      # Find indices in nodes
      source_idx <- match(trans_status, nodes$name) - 1  # 0-based indexing
      target_idx <- match(tmt_status, nodes$name) - 1  # 0-based indexing
      
      # Add link
      links <- rbind(links, data.frame(
        source = source_idx,
        target = target_idx,
        value = count
      ))
    }
  } else if (has_tmt) {
    # Link directly from immunopeptidome to TMT
    tmt_counts <- integrated_data %>%
      dplyr::filter(!is.na(tmt_status)) %>%
      dplyr::count(tmt_status)
    
    # Add links
    for (i in 1:nrow(tmt_counts)) {
      status <- tmt_counts$tmt_status[i]
      count <- tmt_counts$n[i]
      
      # Find index of this status in nodes
      target_idx <- match(status, nodes$name) - 1  # 0-based indexing
      
      # Add link from Detected (0) to this status
      links <- rbind(links, data.frame(
        source = 0,  # Detected
        target = target_idx,
        value = count
      ))
    }
  }
  
  # Extract visualization settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # File paths
  file_paths <- list()
  
  # Create Sankey diagram using networkD3
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0("multi_omics_sankey.", format))
    
    if (format == "html") {
      # Create interactive HTML Sankey
      sankey <- networkD3::sankeyNetwork(
        Links = links,
        Nodes = nodes,
        Source = "source",
        Target = "target",
        Value = "value",
        NodeID = "name",
        units = "peptides",
        fontSize = 12,
        nodeWidth = 30
      )
      
      # Save as HTML
      htmlwidgets::saveWidget(sankey, file_path, selfcontained = TRUE)
      
      file_paths[["html"]] <- file_path
    } else if (format == "png") {
      # For static image formats, we'll use plotly to render the Sankey
      # then capture it as an image
      if (requireNamespace("plotly", quietly = TRUE)) {
        # Create Sankey with plotly
        p <- plotly::plot_ly(
          type = "sankey",
          orientation = "h",
          node = list(
            label = nodes$name,
            pad = 15,
            thickness = 20,
            line = list(
              color = "black",
              width = 0.5
            )
          ),
          link = list(
            source = links$source,
            target = links$target,
            value = links$value
          )
        )
        
        # Add title
        p <- plotly::layout(p, title = "Multi-Omics Integration Flow")
        
        # Save as image
        plotly::orca(p, file_path, width = plot_width * 100, height = plot_height * 100)
        
        file_paths[["png"]] <- file_path
      } else {
        warning("plotly package not available for static Sankey diagram rendering")
      }
    }
  }
  
  return(file_paths)
}

#' Create a table of candidate public neoantigens
#' 
#' @param integrated_data Data frame with integrated multi-omics data
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization files
create_neoantigen_table <- function(integrated_data, viz_dir, config) {
  # Check if we have neoantigen classifications
  if (!"public_neoantigen_classification" %in% colnames(integrated_data)) {
    warning("No neoantigen classification column found")
    return(NULL)
  }
  
  # Filter to potential neoantigens
  neoantigens <- integrated_data %>%
    dplyr::filter(public_neoantigen_classification != "Not a public neoantigen") %>%
    dplyr::arrange(dplyr::desc(neoantigen_score), dplyr::desc(sample_count))
  
  if (nrow(neoantigens) == 0) {
    warning("No potential neoantigens found")
    return(NULL)
  }
  
  # Select relevant columns
  table_data <- neoantigens %>%
    dplyr::select(
      Peptide,
      primary_gene,
      sample_count,
      dplyr::contains("log2_fold_change"),
      public_neoantigen_classification
    )
  
  # Write to CSV file
  file_path <- file.path(viz_dir, "candidate_neoantigens.csv")
  write.csv(table_data, file_path, row.names = FALSE)
  
  # If we have knitr/kableExtra, create a formatted HTML table
  if (requireNamespace("knitr", quietly = TRUE) && 
      requireNamespace("kableExtra", quietly = TRUE)) {
    
    # Create a nicely formatted table
    html_table <- knitr::kable(table_data, format = "html", caption = "Candidate Public Neoantigens") %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      kableExtra::column_spec(1, bold = TRUE) %>%  # Peptide column
      kableExtra::column_spec(2, italic = TRUE)    # Gene column
    
    # Color code the fold changes
    for (col in grep("log2_fold_change", colnames(table_data), value = TRUE)) {
      col_idx <- which(colnames(table_data) == col)
      html_table <- html_table %>%
        kableExtra::column_spec(
          col_idx,
          color = "white",
          background = kableExtra::spec_color(table_data[[col]], end = 0.7, option = "D", 
                                              direction = 1, scale_from = c(-3, 3))
        )
    }
    
    # Write to HTML file
    html_file <- file.path(viz_dir, "candidate_neoantigens.html")
    cat(as.character(html_table), file = html_file)
    
    return(list(csv = file_path, html = html_file))
  }
  
  return(list(csv = file_path))
}

#' Generate multi-omics visualizations
#' 
#' @param multi_omics_results Results from multi-omics integration
#' @param viz_dir Directory to save visualizations
#' @param config Configuration settings
#' @return List of paths to visualizations
generate_multi_omics_visualizations <- function(multi_omics_results, viz_dir, config) {
  # Make sure output directory exists
  if (!dir.exists(viz_dir)) {
    dir.create(viz_dir, recursive = TRUE)
  }
  
  # Initialize results
  viz_paths <- list()
  
  # Generate visualizations for each threshold
  if (!is.null(multi_omics_results$threshold_results)) {
    for (threshold_name in names(multi_omics_results$threshold_results)) {
      threshold_data <- multi_omics_results$threshold_results[[threshold_name]]
      threshold_value <- as.numeric(gsub("shared_", "", threshold_name))
      
      # Only proceed if we have potential neoantigens
      if ("public_neoantigen_classification" %in% colnames(threshold_data)) {
        neoantigens <- threshold_data %>%
          dplyr::filter(public_neoantigen_classification != "Not a public neoantigen")
        
        if (nrow(neoantigens) > 0) {
          # Create heatmap
          tryCatch({
            viz_paths[[paste0(threshold_name, "_heatmap")]] <- create_public_neoantigens_heatmap(
              neoantigens, viz_dir, config
            )
          }, error = function(e) {
            cat("Error creating neoantigen heatmap for threshold", threshold_value, ":", e$message, "\n")
          })
        }
      }
    }
  }
  
  return(viz_paths)
}