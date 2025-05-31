#' HLA-I Analysis Visualization Functions
#' Functions for creating standardized visualizations for peptide data
#' @author Your Name
#' @version 1.0

source("peptide_core_utils.R")

# Create a heatmap of peptide presence/absence or intensity
create_peptide_heatmap <- function(peptide_matrix, value_cols, 
                                   is_intensity = TRUE, log_transform = TRUE,
                                   peptide_col = "Peptide", 
                                   annotation_cols = NULL,
                                   title = "Peptide Heatmap",
                                   cluster_rows = FALSE,
                                   cluster_cols = FALSE) {
  
  # Extract matrix for heatmap
  matrix_data <- peptide_matrix %>%
    select(all_of(c(peptide_col, value_cols))) %>%
    column_to_rownames(peptide_col)
  
  # Convert to matrix
  heat_matrix <- as.matrix(matrix_data)
  
  # Log transform intensity if needed
  if (is_intensity && log_transform) {
    heat_matrix <- log10(heat_matrix + 1)
    title <- paste0(title, " (log10)")
  }
  
  # Create presence/absence matrix if needed
  if (!is_intensity) {
    heat_matrix <- (heat_matrix > 0) * 1
  }
  
  # Create row annotations if provided
  row_annotation <- NULL
  if (!is.null(annotation_cols)) {
    row_annotation <- peptide_matrix %>%
      select(all_of(c(peptide_col, annotation_cols))) %>%
      column_to_rownames(peptide_col)
  }
  
  # Create color palette
  if (is_intensity) {
    colors <- get_color_palette("sequential", 100)
  } else {
    colors <- get_color_palette("binary", 2)
  }
  
  # Create the heatmap
  heatmap <- pheatmap(
    heat_matrix,
    main = title,
    color = colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    annotation_row = row_annotation,
    display_numbers = nrow(heat_matrix) <= 50, # Only show numbers for smaller matrices
    number_format = if(is_intensity) "%.1f" else "%d",
    fontsize_row = max(4, min(10, 300/nrow(heat_matrix))), # Adjust font size
    fontsize_col = 10
  )
  
  return(heatmap)
}

# Create scatter plot for comparing two datasets (e.g., immunopeptidome vs transcriptome)
create_comparison_scatter <- function(data, x_col, y_col, 
                                      color_col = NULL, highlight_col = NULL,
                                      x_label = NULL, y_label = NULL,
                                      title = "Data Comparison", 
                                      subtitle = NULL) {
  
  # Filter out NAs
  plot_data <- data %>%
    filter(!is.na(!!sym(x_col))) %>%
    filter(!is.na(!!sym(y_col)))
  
  if (nrow(plot_data) == 0) {
    warning("No data available for scatter plot:", title)
    return(NULL)
  }
  
  # Set default labels if not provided
  if (is.null(x_label)) x_label <- x_col
  if (is.null(y_label)) y_label <- y_col
  
  # Base plot
  p <- ggplot(plot_data, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    # Add lines at +/- 1 log2FC
    geom_hline(yintercept = 1, linetype = "dotted", color = "darkgray") +
    geom_hline(yintercept = -1, linetype = "dotted", color = "darkgray") +
    geom_vline(xintercept = 1, linetype = "dotted", color = "darkgray") +
    geom_vline(xintercept = -1, linetype = "dotted", color = "darkgray") +
    theme_minimal() +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  # Add color if specified
  if(!is.null(color_col)) {
    p <- p + geom_point(aes(color = !!sym(color_col)), alpha = 0.7) +
      labs(color = color_col)
  } else {
    p <- p + geom_point(alpha = 0.7, color = "steelblue")
  }
  
  # Highlight specific points if specified
  if(!is.null(highlight_col)) {
    # Filter only highlighted points
    highlighted_data <- plot_data %>%
      filter(!!sym(highlight_col) == TRUE)
    
    if(nrow(highlighted_data) > 0) {
      p <- p + 
        geom_point(data = highlighted_data, 
                   aes(x = !!sym(x_col), y = !!sym(y_col)), 
                   color = "red", size = 3, shape = 17)
      
      # Add labels if fewer than 30 highlighted points
      if(nrow(highlighted_data) <= 30 && "Peptide" %in% colnames(highlighted_data)) {
        p <- p + geom_text_repel(data = highlighted_data,
                                 aes(x = !!sym(x_col), y = !!sym(y_col), 
                                     label = Peptide),
                                 size = 3, max.overlaps = 20)
      } else if(nrow(highlighted_data) <= 30 && "primary_gene" %in% colnames(highlighted_data)) {
        p <- p + geom_text_repel(data = highlighted_data,
                                 aes(x = !!sym(x_col), y = !!sym(y_col), 
                                     label = primary_gene),
                                 size = 3, max.overlaps = 20)
      }
    }
  }
  
  return(p)
}

# Enhanced volcano plot function that properly handles statistical significance
create_volcano_plot <- function(data, fc_col = "log2_fold_change", p_val_col = NULL, 
                                label_col = NULL, highlight_col = NULL, sig_threshold = 0.05,
                                fc_threshold = 1, title = NULL, max_labels = 30) {
  
  # Check that fold change column exists
  if (!fc_col %in% colnames(data)) {
    stop("Fold change column '", fc_col, "' not found in data")
  }
  
  # If p-value column exists, use it; otherwise, use fold change and intensity as a "stand-in" 
  if (!is.null(p_val_col) && p_val_col %in% colnames(data)) {
    # Traditional volcano plot with p-values
    plot_data <- data %>%
      filter(!is.na(.data[[fc_col]]), !is.na(.data[[p_val_col]])) %>%
      mutate(
        neg_log10_p = -log10(.data[[p_val_col]]),
        significance = case_when(
          .data[[p_val_col]] < sig_threshold & .data[[fc_col]] > fc_threshold ~ "Up",
          .data[[p_val_col]] < sig_threshold & .data[[fc_col]] < -fc_threshold ~ "Down",
          TRUE ~ "Not significant"
        )
      )
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = .data[[fc_col]], y = neg_log10_p)) +
      geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "darkgray") +
      geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "darkgray")
    
    # Set y-axis label
    y_label <- paste0("-log10(", p_val_col, ")")
    
  } else {
    # For fold change-only data, use fold change on x and intensity on y axes
    # Look for intensity column 
    intensity_cols <- grep("intensity|abundance|expr", colnames(data), ignore.case = TRUE, value = TRUE)
    
    # If no intensity columns found, create a simulated significance
    if (length(intensity_cols) == 0) {
      # Create a simulated "significance" based on fold change magnitude
      plot_data <- data %>%
        filter(!is.na(.data[[fc_col]])) %>%
        mutate(
          simulated_significance = abs(.data[[fc_col]]),
          significance = case_when(
            .data[[fc_col]] > fc_threshold ~ "Up",
            .data[[fc_col]] < -fc_threshold ~ "Down",
            TRUE ~ "Not significant"
          )
        )
      
      # Create the plot with fold change vs fold change magnitude
      p <- ggplot(plot_data, aes(x = .data[[fc_col]], y = simulated_significance)) +
        geom_hline(yintercept = fc_threshold, linetype = "dashed", color = "darkgray") +
        geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "darkgray")
      
      # Set y-axis label
      y_label <- "Fold Change Magnitude (absolute value)"
      
    } else {
      # Use the first intensity column found
      intensity_col <- intensity_cols[1]
      
      # Log-transform the intensity to make the plot more readable
      plot_data <- data %>%
        filter(!is.na(.data[[fc_col]]), !is.na(.data[[intensity_col]])) %>%
        mutate(
          log_intensity = log10(.data[[intensity_col]] + 1),  # Add 1 to avoid log(0)
          significance = case_when(
            .data[[fc_col]] > fc_threshold ~ "Up",
            .data[[fc_col]] < -fc_threshold ~ "Down",
            TRUE ~ "Not significant"
          )
        )
      
      # Create the plot with fold change vs intensity
      p <- ggplot(plot_data, aes(x = .data[[fc_col]], y = log_intensity)) +
        geom_hline(yintercept = log10(10 + 1), linetype = "dashed", color = "darkgray") +
        geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "darkgray")
      
      # Set y-axis label
      y_label <- paste0("Log10(", intensity_col, " + 1)")
    }
  }
  
  # Add points with coloring by significance
  p <- p + geom_point(aes(color = significance), alpha = 0.7) +
    scale_color_manual(values = c(
      "Up" = "red",
      "Down" = "blue",
      "Not significant" = "gray70"
    ))
  
  # Add highlight for specific points if highlight column is provided
  if (!is.null(highlight_col) && highlight_col %in% colnames(data)) {
    # Add highlight information to plot data
    plot_data$highlighted <- plot_data[[highlight_col]]
    
    # If highlight column is logical, add special highlighting
    if (is.logical(plot_data$highlighted)) {
      # Add highlighted points with bigger size
      p <- p + geom_point(
        data = filter(plot_data, highlighted), 
        aes(shape = "Highlighted"), 
        size = 3, 
        stroke = 1.5
      ) +
        scale_shape_manual(values = c("Highlighted" = 1))  # Use hollow circle
    }
  }
  
  # Add labels for top differentially expressed points if label column is provided
  if (!is.null(label_col) && label_col %in% colnames(data)) {
    # Make sure label values are available
    plot_data$label <- as.character(plot_data[[label_col]])
    
    # Select most significant up and down points for labeling
    if (exists("neg_log10_p", plot_data)) {
      # Use p-value significance for selection
      to_label <- plot_data %>%
        filter(significance %in% c("Up", "Down")) %>%
        arrange(desc(neg_log10_p)) %>%
        head(max_labels)
    } else if (exists("log_intensity", plot_data)) {
      # Use fold change and intensity for selection
      up_labels <- plot_data %>%
        filter(significance == "Up") %>%
        arrange(desc(.data[[fc_col]]), desc(log_intensity)) %>%
        head(max_labels/2)
      
      down_labels <- plot_data %>%
        filter(significance == "Down") %>%
        arrange(.data[[fc_col]], desc(log_intensity)) %>%
        head(max_labels/2)
      
      to_label <- bind_rows(up_labels, down_labels)
    } else {
      # Use fold change magnitude for selection
      to_label <- plot_data %>%
        filter(significance %in% c("Up", "Down")) %>%
        arrange(desc(abs(.data[[fc_col]]))) %>%
        head(max_labels)
    }
    
    # Try to use ggrepel for better label placement if available
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = to_label,
        aes(label = label),
        box.padding = 0.5,
        max.overlaps = 30,
        size = 3
      )
    } else {
      # Fallback to standard text labels
      p <- p + geom_text(
        data = to_label,
        aes(label = label),
        hjust = 0, vjust = 0,
        nudge_x = 0.1, nudge_y = 0.1,
        size = 3
      )
    }
  }
  
  # Add appropriate labels
  p <- p + labs(
    title = title %||% "Volcano Plot",
    x = paste("Log2 Fold Change (", gsub("log2_fold_change_*", "", fc_col), ")", sep=""),
    y = y_label,
    color = "Significance",
    shape = NULL
  ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10)
    )
  
  return(p)
}

# Function to create better visualizations for tumor-specific peptides
create_tumor_specific_plot <- function(data, tumor_specific_col = "detection_status", 
                                       intensity_col = NULL, tumor_value = "Tumor-specific",
                                       label_col = NULL, title = NULL, max_peptides = 50) {
  
  # Check that required columns exist
  if (!tumor_specific_col %in% colnames(data)) {
    stop("Tumor-specific column '", tumor_specific_col, "' not found in data")
  }
  
  # Find tumor-specific peptides
  tumor_specific <- data %>%
    filter(.data[[tumor_specific_col]] == tumor_value)
  
  if (nrow(tumor_specific) == 0) {
    stop("No tumor-specific peptides found")
  }
  
  # Find intensity column if not specified
  if (is.null(intensity_col)) {
    # Look for tumor intensity column
    intensity_cols <- grep("intensity.*tumor|tumor.*intensity", 
                           colnames(data), ignore.case = TRUE, value = TRUE)
    
    if (length(intensity_cols) > 0) {
      intensity_col <- intensity_cols[1]
    } else {
      # Try any intensity column
      intensity_cols <- grep("intensity", colnames(data), ignore.case = TRUE, value = TRUE)
      if (length(intensity_cols) > 0) {
        intensity_col <- intensity_cols[1]
      } else {
        stop("No intensity column found. Please specify the intensity_col parameter.")
      }
    }
  }
  
  # Determine label column if not specified
  if (is.null(label_col)) {
    # Look for peptide column first
    if ("Peptide" %in% colnames(data)) {
      label_col <- "Peptide"
    } else if ("peptide" %in% colnames(data)) {
      label_col <- "peptide"
    } else {
      # Try gene column
      if ("primary_gene" %in% colnames(data)) {
        label_col <- "primary_gene"
      } else if ("gene" %in% colnames(data)) {
        label_col <- "gene"
      } else {
        # Use first column as fallback
        label_col <- colnames(data)[1]
      }
    }
  }
  
  # Sort tumor-specific peptides by intensity
  tumor_specific <- tumor_specific %>%
    filter(!is.na(.data[[intensity_col]])) %>%
    arrange(desc(.data[[intensity_col]])) %>%
    # Take top peptides for visualization
    head(max_peptides)
  
  # Add label information
  tumor_specific$label <- as.character(tumor_specific[[label_col]])
  
  # Check if we have fusion information
  has_fusion_info <- any(grepl("fusion|spans_junction", colnames(tumor_specific)))
  
  # Create enhanced bar plot
  if (has_fusion_info) {
    # Identify fusion-related columns for coloring
    fusion_col <- grep("from_fusion", colnames(tumor_specific), value = TRUE)[1]
    junction_col <- grep("spans_junction", colnames(tumor_specific), value = TRUE)[1]
    
    # Determine color mapping based on available columns
    if (!is.null(fusion_col) && !is.null(junction_col)) {
      tumor_specific <- tumor_specific %>%
        mutate(fusion_category = case_when(
          .data[[junction_col]] ~ "Junction-spanning",
          .data[[fusion_col]] ~ "Fusion-derived",
          TRUE ~ "Non-fusion"
        ))
      color_col <- "fusion_category"
    } else if (!is.null(fusion_col)) {
      color_col <- fusion_col
    } else {
      color_col <- NULL
    }
    
    # Create plot with fusion coloring
    p <- ggplot(tumor_specific, aes(x = reorder(label, .data[[intensity_col]]), 
                                    y = .data[[intensity_col]])) +
      geom_bar(stat = "identity", aes(fill = .data[[color_col]])) +
      scale_fill_manual(values = c(
        "Junction-spanning" = "purple",
        "Fusion-derived" = "orange",
        "Non-fusion" = "blue",
        "TRUE" = "orange",  # For boolean columns
        "FALSE" = "blue"    # For boolean columns
      ))
  } else {
    # Check for other interesting categorical columns for coloring
    categorical_cols <- sapply(tumor_specific, function(x) is.factor(x) || is.character(x))
    categorical_cols <- names(categorical_cols)[categorical_cols]
    
    # Exclude label and status columns
    categorical_cols <- setdiff(categorical_cols, c(label_col, tumor_specific_col))
    
    if (length(categorical_cols) > 0) {
      # Use first categorical column for coloring
      color_col <- categorical_cols[1]
      p <- ggplot(tumor_specific, aes(x = reorder(label, .data[[intensity_col]]), 
                                      y = .data[[intensity_col]])) +
        geom_bar(stat = "identity", aes(fill = .data[[color_col]]))
    } else {
      # Basic plot without special coloring
      p <- ggplot(tumor_specific, aes(x = reorder(label, .data[[intensity_col]]), 
                                      y = .data[[intensity_col]])) +
        geom_bar(stat = "identity", fill = "steelblue")
    }
  }
  
  # Complete the plot with appropriate styling
  p <- p + 
    coord_flip() +  # Horizontal bars for better label readability
    labs(
      title = title %||% "Tumor-Specific Peptides",
      x = paste(label_col),
      y = paste("Intensity (", gsub(".*intensity_", "", intensity_col), ")", sep=""),
      fill = "Category"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  return(p)
}

# Create ranked bar plot for tumor-specific peptides
tumor_specific_plot <- create_tumor_specific_plot(
  data = multi_omics_data,
  tumor_specific_col = "detection_status",
  intensity_col = paste0("total_intensity_", config$tumor_id),
  tumor_value = "Tumor-specific",
  label_col = "Peptide",
  title = paste0("Peptides Exclusive to ", config$tumor_id)
)
save_plot(tumor_specific_plot, file.path(viz_dir, "tumor_specific_peptides"))

# Function to create better heatmap for visualizing peptide patterns
create_peptide_heatmap <- function(data, value_cols, is_intensity = TRUE, 
                                   log_transform = TRUE, peptide_col = "Peptide",
                                   annotation_cols = NULL, title = NULL,
                                   cluster_rows = TRUE, cluster_cols = TRUE) {
  
  # Check required packages
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required for heatmap creation")
  }
  
  # Verify that required columns exist
  if (!peptide_col %in% colnames(data)) {
    stop("Peptide column '", peptide_col, "' not found in data")
  }
  
  # Check that all value columns exist
  missing_cols <- value_cols[!value_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    stop("Value columns not found: ", paste(missing_cols, collapse = ", "))
  }
  
  # Extract data for heatmap
  heatmap_data <- as.matrix(data[, value_cols, drop = FALSE])
  rownames(heatmap_data) <- data[[peptide_col]]
  
  # Process values based on parameters
  if (is_intensity) {
    # For intensity data, add a small value to avoid log(0) issues
    heatmap_data[is.na(heatmap_data)] <- 0
    
    if (log_transform) {
      heatmap_data <- log2(heatmap_data + 1)
    }
  } else {
    # For presence/absence data, convert to binary
    heatmap_data[is.na(heatmap_data)] <- 0
    heatmap_data <- ifelse(heatmap_data > 0, 1, 0)
  }
  
  # Clean column names for display
  colnames(heatmap_data) <- gsub("total_intensity_|intensity_", "", colnames(heatmap_data))
  
  # Set up annotation if specified
  row_annotation <- NULL
  if (!is.null(annotation_cols)) {
    valid_anno_cols <- annotation_cols[annotation_cols %in% colnames(data)]
    
    if (length(valid_anno_cols) > 0) {
      row_annotation <- data[, valid_anno_cols, drop = FALSE]
      rownames(row_annotation) <- data[[peptide_col]]
    }
  }
  
  # Set up color palette based on data type
  if (is_intensity) {
    if (log_transform) {
      # For log-transformed intensity data, use a blue-white-red gradient
      color_palette <- colorRampPalette(c("navy", "blue", "white", "red", "darkred"))(100)
      # Center around median
      breaks_mid <- median(heatmap_data, na.rm = TRUE)
      max_abs <- max(abs(range(heatmap_data, na.rm = TRUE) - breaks_mid))
      breaks_range <- c(breaks_mid - max_abs, breaks_mid + max_abs)
    } else {
      # For raw intensity, use a white to blue gradient
      color_palette <- colorRampPalette(c("white", "steelblue", "navy"))(100)
      breaks_range <- range(heatmap_data, na.rm = TRUE)
    }
  } else {
    # For binary data, use a simple two-color palette
    color_palette <- c("white", "darkblue")
    breaks_range <- c(0, 1)
  }
  
  # Create heatmap
  pheatmap::pheatmap(
    heatmap_data,
    color = color_palette,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    annotation_row = row_annotation,
    main = title %||% "Peptide Heatmap",
    fontsize_row = ifelse(nrow(heatmap_data) <= 50, 9, 6),
    fontsize_col = 10,
    scale = "none",  # We've already done any necessary transformations
    show_rownames = nrow(heatmap_data) <= 100,  # Only show row names for smaller heatmaps
    border_color = NA,  # No cell borders for cleaner look
    treeheight_row = ifelse(cluster_rows, 50, 0),
    treeheight_col = ifelse(cluster_cols, 50, 0)
  )
}

# Function to create an enrichment analysis visualization for tumor-specific peptides
create_enrichment_visualization <- function(data, gene_col = "primary_gene", 
                                            group_col = "detection_status", 
                                            group_value = "Tumor-specific",
                                            title = NULL) {
  # This is a simplified enrichment visualization - a proper enrichment analysis
  # would require additional packages like clusterProfiler or gprofiler2
  
  # Extract genes from the specified group
  group_genes <- data %>%
    filter(.data[[group_col]] == group_value) %>%
    pull(gene_col) %>%
    unique()
  
  if (length(group_genes) == 0) {
    stop("No genes found for the specified group")
  }
  
  # Simple word frequency analysis for gene names
  gene_prefixes <- substr(group_genes, 1, 3)
  prefix_counts <- table(gene_prefixes)
  prefix_df <- data.frame(
    prefix = names(prefix_counts),
    count = as.numeric(prefix_counts)
  ) %>%
    arrange(desc(count)) %>%
    filter(count >= 2)  # Only show prefixes with at least 2 occurrences
  
  # Create visualization
  if (nrow(prefix_df) > 0) {
    p <- ggplot(prefix_df, aes(x = reorder(prefix, count), y = count)) +
      geom_bar(stat = "identity", fill = "darkgreen") +
      coord_flip() +
      labs(
        title = title %||% paste("Gene Family Distribution in", group_value, "Genes"),
        subtitle = paste(length(group_genes), "unique genes analyzed"),
        x = "Gene Prefix",
        y = "Count"
      ) +
      theme_minimal()
    
    return(p)
  } else {
    message("Not enough gene prefix patterns found for visualization")
    return(NULL)
  }
}

# Helper function for NULL coalescing
'%||%' <- function(x, y) if(is.null(x)) y else x

# # Create volcano plot for differential expression analysis
# create_volcano_plot <- function(data, fc_col, pval_col = NULL, 
#                                 label_col = NULL, highlight_col = NULL,
#                                 fc_cutoff = 1, pval_cutoff = 0.05,
#                                 title = "Volcano Plot", 
#                                 point_size = 2, alpha = 0.7) {
#   
#   # Filter out NAs for fold change
#   plot_data <- data %>%
#     filter(!is.na(!!sym(fc_col)))
#   
#   # Base plot without p-value
#   if (is.null(pval_col) || !pval_col %in% colnames(plot_data)) {
#     # Use fixed y-value for p-value 
#     p <- ggplot(plot_data, aes(x = !!sym(fc_col), y = -log10(0.05))) +
#       geom_point(alpha = alpha, size = point_size, aes(color = abs(!!sym(fc_col)) > fc_cutoff)) +
#       scale_color_manual(values = c("gray", "red"), 
#                          labels = c(paste("< ", fc_cutoff), paste(">= ", fc_cutoff)),
#                          name = "Log2 FC")
#     
#     # Set y-label for placeholder p-value
#     y_label <- "-log10(p-value) [placeholder]"
#   } else {
#     # Filter out NAs for p-value
#     plot_data <- plot_data %>%
#       filter(!is.na(!!sym(pval_col)))
#     
#     # Convert p-values of 0 to a small value to avoid -log10(0) = Inf
#     plot_data <- plot_data %>%
#       mutate(!!pval_col := ifelse(!!sym(pval_col) == 0, 1e-300, !!sym(pval_col)))
#     
#     # Create the plot with p-value
#     p <- ggplot(plot_data, aes(x = !!sym(fc_col), y = -log10(!!sym(pval_col)))) +
#       geom_point(alpha = alpha, size = point_size, 
#                  aes(color = (abs(!!sym(fc_col)) > fc_cutoff & !!sym(pval_col) < pval_cutoff))) +
#       scale_color_manual(values = c("gray", "red"), 
#                          labels = c("Not significant", "Significant"),
#                          name = "Significance")
#     
#     # Set y-label for actual p-value
#     y_label <- "-log10(p-value)"
#   }
#   
#   # Common plot elements
#   p <- p +
#     geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
#     geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "darkgray") +
#     theme_minimal() +
#     labs(
#       title = title,
#       x = "Log2 Fold Change",
#       y = y_label
#     )
#   
#   # Add labels if specified and column exists
#   if(!is.null(label_col) && label_col %in% colnames(plot_data)) {
#     # Only label significant points
#     if(!is.null(pval_col) && pval_col %in% colnames(plot_data)) {
#       label_data <- plot_data %>%
#         filter(abs(!!sym(fc_col)) > fc_cutoff & !!sym(pval_col) < pval_cutoff)
#     } else {
#       label_data <- plot_data %>%
#         filter(abs(!!sym(fc_col)) > fc_cutoff)
#     }
#     
#     # Limit to top most extreme points
#     if(nrow(label_data) > 0) {
#       # Calculate extreme score safely
#       if(!is.null(pval_col) && pval_col %in% colnames(label_data)) {
#         label_data <- label_data %>%
#           mutate(extreme_score = abs(!!sym(fc_col)) * -log10(!!sym(pval_col))) %>%
#           arrange(desc(extreme_score)) %>%
#           head(25)
#       } else {
#         label_data <- label_data %>%
#           mutate(extreme_score = abs(!!sym(fc_col))) %>%
#           arrange(desc(extreme_score)) %>%
#           head(25)
#       }
#       
#       # Add text labels avoiding the sym() call on NULL values
#       if(!is.null(pval_col) && pval_col %in% colnames(label_data)) {
#         # Use p-value column for y-position
#         p <- p + geom_text_repel(data = label_data,
#                                  aes(x = !!sym(fc_col), 
#                                      y = -log10(!!sym(pval_col)),
#                                      label = !!sym(label_col)),
#                                  size = 3, max.overlaps = 20)
#       } else {
#         # Use fixed y-position
#         p <- p + geom_text_repel(data = label_data,
#                                  aes(x = !!sym(fc_col), 
#                                      y = -log10(0.05),
#                                      label = !!sym(label_col)),
#                                  size = 3, max.overlaps = 20)
#       }
#     }
#   }
#   
#   # Highlight specific points if specified
#   if(!is.null(highlight_col) && highlight_col %in% colnames(plot_data)) {
#     # Filter only highlighted points
#     highlighted_data <- plot_data %>%
#       filter(!!sym(highlight_col) == TRUE)
#     
#     if(nrow(highlighted_data) > 0) {
#       # Add highlighted points avoiding the sym() call on NULL values
#       if(!is.null(pval_col) && pval_col %in% colnames(highlighted_data)) {
#         # Use p-value column for y-position
#         p <- p + geom_point(data = highlighted_data, 
#                             aes(x = !!sym(fc_col), 
#                                 y = -log10(!!sym(pval_col))),
#                             color = "purple", size = 3, shape = 17)
#       } else {
#         # Use fixed y-position
#         p <- p + geom_point(data = highlighted_data, 
#                             aes(x = !!sym(fc_col), 
#                                 y = -log10(0.05)),
#                             color = "purple", size = 3, shape = 17)
#       }
#     }
#   }
#   
#   return(p)
# }

# Create Venn diagram for sample overlaps
create_venn_diagram <- function(data, sample_cols, labels = NULL, 
                                title = "Sample Overlaps", 
                                colors = NULL, output_file = NULL) {
  
  if(length(sample_cols) < 2 || length(sample_cols) > 5) {
    stop("Venn diagrams support between 2 and 5 sets")
  }
  
  # Create list of items in each set
  sets <- list()
  for(i in 1:length(sample_cols)) {
    sets[[i]] <- data %>%
      filter(!!sym(sample_cols[i]) > 0) %>%
      pull(Peptide)
  }
  
  # Set names if provided
  if(!is.null(labels) && length(labels) == length(sample_cols)) {
    names(sets) <- labels
  } else {
    names(sets) <- sample_cols
  }
  
  # Set colors if provided
  if(is.null(colors)) {
    colors <- brewer.pal(length(sample_cols), "Set1")
  }
  
  # Create Venn diagram
  if(!is.null(output_file)) {
    # Save to file
    pdf(output_file, width = 10, height = 8)
    venn_obj <- venn.diagram(
      sets,
      filename = NULL,  # Don't save separately
      fill = colors,
      alpha = 0.5,
      main = title,
      main.cex = 1.5
    )
    grid.draw(venn_obj)
    dev.off()
  } else {
    # Return the Venn diagram object for rendering
    venn_obj <- venn.diagram(
      sets,
      filename = NULL,
      fill = colors,
      alpha = 0.5,
      main = title,
      main.cex = 1.5
    )
    return(venn_obj)
  }
}

# Create UpSet plot for complex set intersections
create_upset_plot <- function(data, sample_cols, labels = NULL,
                              min_size = 1, max_sets = 20,
                              title = "Sample Intersections") {
  
  # Create binary matrix for UpSet
  binary_data <- matrix(0, nrow = nrow(data), ncol = length(sample_cols))
  for(i in 1:length(sample_cols)) {
    binary_data[, i] <- as.numeric(data[[sample_cols[i]]] > 0)
  }
  
  # Set column names
  if(!is.null(labels) && length(labels) == length(sample_cols)) {
    colnames(binary_data) <- labels
  } else {
    colnames(binary_data) <- sample_cols
  }
  
  # Create UpSet plot
  upset_plot <- upset(
    as.data.frame(binary_data),
    nsets = length(sample_cols),
    nintersects = max_sets,
    mb.ratio = c(0.4, 0.6),
    order.by = "freq",
    main.bar.color = "steelblue",
    sets.bar.color = "darkred",
    keep.order = FALSE,
    set_size.show = TRUE,
    text.scale = 1.2,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size"
  )
  
  return(upset_plot)
}

# Create bar plot for categories
create_bar_plot <- function(data, x_col, y_col, fill_col = NULL,
                            stack = FALSE, coord_flip = FALSE,
                            title = NULL, x_label = NULL, y_label = NULL,
                            text_labels = TRUE) {
  
  # Set default labels if not provided
  if(is.null(x_label)) x_label <- x_col
  if(is.null(y_label)) y_label <- y_col
  
  # Basic plot
  p <- ggplot(data, aes_string(x = x_col, y = y_col))
  
  # Add fill if specified
  if(!is.null(fill_col)) {
    if(stack) {
      p <- p + geom_bar(aes_string(fill = fill_col), stat = "identity", position = "stack")
    } else {
      p <- p + geom_bar(aes_string(fill = fill_col), stat = "identity", position = "dodge")
    }
  } else {
    p <- p + geom_bar(fill = "steelblue", stat = "identity")
  }
  
  # Add text labels if requested
  if(text_labels) {
    if(!is.null(fill_col) && !stack) {
      p <- p + geom_text(aes_string(label = y_col), position = position_dodge(width = 0.9), vjust = -0.5)
    } else {
      p <- p + geom_text(aes_string(label = y_col), vjust = -0.5)
    }
  }
  
  # Flip coordinates if requested
  if(coord_flip) {
    p <- p + coord_flip()
  }
  
  # Formatting
  p <- p + 
    theme_minimal() +
    labs(
      title = title,
      x = x_label,
      y = y_label
    )
  
  return(p)
}

# Create histogram of fold changes
create_fold_change_histogram <- function(data, fc_col = "log2_fold_change", 
                                         category_col = "peptide_category",
                                         title = "Distribution of Peptide Fold Changes", 
                                         bin_count = 50) {
  
  # Filter out NAs
  plot_data <- data %>%
    filter(!is.na(!!sym(fc_col)))
  
  # Check if category column exists and create it if not
  if (!category_col %in% colnames(plot_data)) {
    cat("Column", category_col, "not found. Creating it based on fold change values.\n")
    plot_data <- plot_data %>%
      mutate(!!category_col := case_when(
        !!sym(fc_col) > 1 ~ "Up in Tumor (FC > 2)",
        !!sym(fc_col) < -1 ~ "Down in Tumor (FC < 0.5)",
        TRUE ~ "Similar (-1 < log2FC < 1)"
      ))
  }
  
  # Create the plot
  if (category_col %in% colnames(plot_data)) {
    # With category coloring
    p <- ggplot(plot_data, aes(x = !!sym(fc_col), fill = !!sym(category_col))) +
      geom_histogram(bins = bin_count, color = "black", alpha = 0.7) +
      scale_fill_manual(values = c("Up in Tumor (FC > 2)" = "red", 
                                   "Down in Tumor (FC < 0.5)" = "blue", 
                                   "Similar (-1 < log2FC < 1)" = "gray")) +
      theme_minimal() +
      labs(
        title = title,
        x = "Log2 Fold Change",
        y = "Count",
        fill = "Category"
      )
  } else {
    # Without category coloring (fallback)
    p <- ggplot(plot_data, aes(x = !!sym(fc_col))) +
      geom_histogram(bins = bin_count, color = "black", fill = "steelblue", alpha = 0.7) +
      theme_minimal() +
      labs(
        title = title,
        x = "Log2 Fold Change",
        y = "Count"
      )
  }
  
  return(p)
}

# Create barplot of detection status
create_detection_status_plot <- function(data, status_col = "detection_status",
                                         title = "Peptide Detection Status",
                                         coord_flip = FALSE) {
  
  # Create summary
  status_summary <- data %>%
    group_by(!!sym(status_col)) %>%
    summarise(
      count = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(count))
  
  # Create the plot
  p <- ggplot(status_summary, aes(x = !!sym(status_col), y = count, fill = !!sym(status_col))) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), vjust = -0.5) +
    theme_minimal() +
    labs(
      title = title,
      x = "Detection Status",
      y = "Count",
      fill = "Status"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Flip coordinates if requested
  if(coord_flip) {
    p <- p + coord_flip()
  }
  
  return(p)
}

# Create analysis for exclusive peptides
create_exclusive_peptide_analysis <- function(data, tumor_id, normal_id, 
                                              tumor_intensity_col = NULL,
                                              normal_intensity_col = NULL) {
  # Create default column names if not provided
  if(is.null(tumor_intensity_col)) tumor_intensity_col <- paste0("total_intensity_", tumor_id)
  if(is.null(normal_intensity_col)) normal_intensity_col <- paste0("total_intensity_", normal_id)
  
  # Find peptides that are exclusively in tumor (not detected in normal)
  tumor_exclusive_peptides <- data %>%
    filter(!!sym(tumor_intensity_col) > 0 & !!sym(normal_intensity_col) == 0) %>%
    arrange(desc(!!sym(tumor_intensity_col)))
  
  # Find peptides that are exclusively in normal (not detected in tumor)
  normal_exclusive_peptides <- data %>%
    filter(!!sym(normal_intensity_col) > 0 & !!sym(tumor_intensity_col) == 0) %>%
    arrange(desc(!!sym(normal_intensity_col)))
  
  # Count of exclusive peptides
  cat("\nExclusive peptide counts:\n")
  cat("Peptides found only in", tumor_id, ":", nrow(tumor_exclusive_peptides), "\n")
  cat("Peptides found only in", normal_id, ":", nrow(normal_exclusive_peptides), "\n")
  
  # Extract fusion peptides if from_fusion column exists
  if("from_fusion" %in% colnames(data)) {
    tumor_exclusive_fusion <- tumor_exclusive_peptides %>%
      filter(from_fusion) %>%
      arrange(desc(if("spans_junction" %in% colnames(.)) spans_junction else TRUE), 
              desc(!!sym(tumor_intensity_col)))
    
    normal_exclusive_fusion <- normal_exclusive_peptides %>%
      filter(from_fusion) %>%
      arrange(desc(if("spans_junction" %in% colnames(.)) spans_junction else TRUE), 
              desc(!!sym(normal_intensity_col)))
    
    if(nrow(tumor_exclusive_fusion) > 0 || nrow(normal_exclusive_fusion) > 0) {
      cat("\nFusion peptide exclusive counts:\n")
      if(nrow(tumor_exclusive_fusion) > 0) {
        cat("Fusion peptides found only in", tumor_id, ":", nrow(tumor_exclusive_fusion), "\n")
      }
      if(nrow(normal_exclusive_fusion) > 0) {
        cat("Fusion peptides found only in", normal_id, ":", nrow(normal_exclusive_fusion), "\n")
      }
      
      return(list(
        tumor_exclusive = tumor_exclusive_peptides,
        normal_exclusive = normal_exclusive_peptides,
        tumor_exclusive_fusion = if(nrow(tumor_exclusive_fusion) > 0) tumor_exclusive_fusion else NULL,
        normal_exclusive_fusion = if(nrow(normal_exclusive_fusion) > 0) normal_exclusive_fusion else NULL
      ))
    }
  }
  
  return(list(
    tumor_exclusive = tumor_exclusive_peptides,
    normal_exclusive = normal_exclusive_peptides
  ))
}

# Add sanity checks visualization
create_sanity_check_report <- function(data, tumor_id, normal_id,
                                       tumor_intensity_col = NULL,
                                       normal_intensity_col = NULL,
                                       fc_col = "log2_fold_change",
                                       category_col = "peptide_category") {
  
  # Create default column names if not provided
  if(is.null(tumor_intensity_col)) tumor_intensity_col <- paste0("total_intensity_", tumor_id)
  if(is.null(normal_intensity_col)) normal_intensity_col <- paste0("total_intensity_", normal_id)
  
  # Store results
  sanity_results <- list(
    passed = character(),
    failed = list()
  )
  
  # 1. Check consistency in detection status vs intensity values
  detection_sanity_check <- data %>%
    mutate(
      status_check = case_when(
        detection_status == paste0(tumor_id, "-specific") & !!sym(normal_intensity_col) > 0 ~ "FAIL",
        detection_status == paste0(normal_id, "-specific") & !!sym(tumor_intensity_col) > 0 ~ "FAIL",
        detection_status == "Detected in both" & (!!sym(tumor_intensity_col) == 0 | !!sym(normal_intensity_col) == 0) ~ "FAIL",
        TRUE ~ "PASS"
      )
    )
  
  failed_detection <- detection_sanity_check %>% filter(status_check == "FAIL")
  if(nrow(failed_detection) > 0) {
    sanity_results$failed$detection <- failed_detection
  } else {
    sanity_results$passed <- c(sanity_results$passed, "Detection status is consistent with intensity values")
  }
  
  # 2. Check peptide category assignment
  category_sanity_check <- data %>%
    mutate(
      expected_category = case_when(
        !!sym(fc_col) > 1 ~ "Up in Tumor (FC > 2)",
        !!sym(fc_col) < -1 ~ "Down in Tumor (FC < 0.5)",
        TRUE ~ "Similar (-1 < log2FC < 1)"
      ),
      category_check = ifelse(expected_category == !!sym(category_col), "PASS", "FAIL")
    ) %>%
    filter(category_check == "FAIL")
  
  if(nrow(category_sanity_check) > 0) {
    sanity_results$failed$category <- category_sanity_check
  } else {
    sanity_results$passed <- c(sanity_results$passed, "Peptide category assignments are consistent")
  }
  
  # 3. Check for NAs in key analytical columns
  na_check <- data %>%
    summarise(
      NA_in_fold_change = sum(is.na(!!sym(fc_col))),
      NA_in_detection = sum(is.na(detection_status)),
      NA_in_category = sum(is.na(!!sym(category_col)))
    )
  
  if(any(unlist(na_check) > 0)) {
    sanity_results$failed$na_values <- na_check
  } else {
    sanity_results$passed <- c(sanity_results$passed, "No missing values in key analytical columns")
  }
  
  # Basic statistics
  basic_stats <- data %>%
    summarise(
      min_fold_change = min(!!sym(fc_col), na.rm = TRUE),
      max_fold_change = max(!!sym(fc_col), na.rm = TRUE),
      mean_fold_change = mean(!!sym(fc_col), na.rm = TRUE),
      median_fold_change = median(!!sym(fc_col), na.rm = TRUE),
      min_intensity_tumor = min(!!sym(tumor_intensity_col), na.rm = TRUE),
      max_intensity_tumor = max(!!sym(tumor_intensity_col), na.rm = TRUE),
      min_intensity_normal = min(!!sym(normal_intensity_col), na.rm = TRUE),
      max_intensity_normal = max(!!sym(normal_intensity_col), na.rm = TRUE)
    )
  
  sanity_results$stats <- basic_stats
  
  return(sanity_results)
}

# Create interactive dashboard with multiple visualizations
create_interactive_dashboard <- function(plot_list, title = "Peptide Analysis Dashboard",
                                         output_file = "dashboard.html",
                                         summary_text = NULL) {
  
  # Create HTML header with CSS styling
  html_header <- '
  <!DOCTYPE html>
  <html>
  <head>
    <meta charset="UTF-8">
    <title>%s</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
      body {
        font-family: Arial, sans-serif;
        margin: 0;
        padding: 20px;
        background-color: #f5f5f5;
      }
      .dashboard-container {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
        gap: 20px;
        margin-bottom: 30px;
      }
      .dashboard-item {
        background-color: white;
        border-radius: 5px;
        box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        padding: 15px;
        height: 450px;
      }
      iframe {
        width: 100%%;
        height: 100%%;
        border: none;
      }
      h1 {
        color: #333;
        text-align: center;
        margin-bottom: 30px;
      }
      h3 {
        color: #555;
        margin-top: 0;
        margin-bottom: 15px;
        border-bottom: 1px solid #eee;
        padding-bottom: 8px;
      }
      .summary {
        background-color: white;
        border-radius: 5px;
        box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        padding: 20px;
        margin-bottom: 20px;
      }
      .footer {
        text-align: center;
        margin-top: 30px;
        color: #777;
        font-size: 0.9em;
      }
    </style>
  </head>
  <body>
    <h1>%s</h1>
  '
  
  # Format header with title
  html_header <- sprintf(html_header, title, title)
  
  # Add summary section if provided
  if(!is.null(summary_text)) {
    summary_html <- '
    <div class="summary">
      <h2>Analysis Summary</h2>
      <p>%s</p>
    </div>
    '
    html_header <- paste0(html_header, sprintf(summary_html, summary_text))
  }
  
  # Function to create an iframe element for each plot
  create_plot_iframe <- function(file_path, title) {
    iframe_template <- '
    <div class="dashboard-item">
      <h3>%s</h3>
      <iframe src="%s"></iframe>
    </div>
    '
    
    # Format the iframe div
    sprintf(iframe_template, title, basename(file_path))
  }
  
  # Start building the dashboard HTML
  dashboard_html <- html_header
  
  # Add plots in groups of 2
  for(i in seq(1, length(plot_list), by = 2)) {
    # Open dashboard container
    dashboard_html <- paste0(dashboard_html, '\n<div class="dashboard-container">\n')
    
    # Add first plot in this group
    dashboard_html <- paste0(dashboard_html, 
                             create_plot_iframe(plot_list[[i]]$path, plot_list[[i]]$title))
    
    # Add second plot if available
    if(i + 1 <= length(plot_list)) {
      dashboard_html <- paste0(dashboard_html, 
                               create_plot_iframe(plot_list[[i+1]]$path, plot_list[[i+1]]$title))
    }
    
    # Close dashboard container
    dashboard_html <- paste0(dashboard_html, '\n</div>\n')
  }
  
  # Add footer with date
  footer_html <- '
  <div class="footer">
    <p>%s | Generated: %s</p>
  </div>
  </body>
  </html>
  '
  dashboard_html <- paste0(dashboard_html, 
                           sprintf(footer_html, title, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  # Write the complete dashboard HTML to a file
  writeLines(dashboard_html, output_file)
  cat("Dashboard created successfully:", output_file, "\n")
}