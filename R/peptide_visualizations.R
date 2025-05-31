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

# Create volcano plot for differential expression analysis
create_volcano_plot <- function(data, fc_col, pval_col = NULL, 
                                label_col = NULL, highlight_col = NULL,
                                fc_cutoff = 1, pval_cutoff = 0.05,
                                title = "Volcano Plot", 
                                point_size = 2, alpha = 0.7) {
  
  # Filter out NAs for fold change
  plot_data <- data %>%
    filter(!is.na(!!sym(fc_col)))
  
  # Base plot without p-value
  if (is.null(pval_col) || !pval_col %in% colnames(plot_data)) {
    p <- ggplot(plot_data, aes(x = !!sym(fc_col), y = -log10(0.05))) +
      geom_point(alpha = alpha, size = point_size, aes(color = abs(!!sym(fc_col)) > fc_cutoff)) +
      scale_color_manual(values = c("gray", "red"), 
                         labels = c(paste("< ", fc_cutoff), paste(">= ", fc_cutoff)),
                         name = "Log2 FC")
  } else {
    # Filter out NAs for p-value
    plot_data <- plot_data %>%
      filter(!is.na(!!sym(pval_col)))
    
    # Convert p-values of 0 to a small value to avoid -log10(0) = Inf
    plot_data <- plot_data %>%
      mutate(!!pval_col := ifelse(!!sym(pval_col) == 0, 1e-300, !!sym(pval_col)))
    
    # Create the plot with p-value
    p <- ggplot(plot_data, aes(x = !!sym(fc_col), y = -log10(!!sym(pval_col)))) +
      geom_point(alpha = alpha, size = point_size, 
                 aes(color = (abs(!!sym(fc_col)) > fc_cutoff & !!sym(pval_col) < pval_cutoff))) +
      scale_color_manual(values = c("gray", "red"), 
                         labels = c("Not significant", "Significant"),
                         name = "Significance")
  }
  
  # Common plot elements
  p <- p +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "darkgray") +
    theme_minimal() +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-log10(p-value)"
    )
  
  # Add labels if specified
  if(!is.null(label_col) && label_col %in% colnames(plot_data)) {
    # Only label significant points
    if(!is.null(pval_col) && pval_col %in% colnames(plot_data)) {
      label_data <- plot_data %>%
        filter(abs(!!sym(fc_col)) > fc_cutoff & !!sym(pval_col) < pval_cutoff)
    } else {
      label_data <- plot_data %>%
        filter(abs(!!sym(fc_col)) > fc_cutoff)
    }
    
    # Limit to top most extreme points
    label_data <- label_data %>%
      mutate(extreme_score = abs(!!sym(fc_col)) * -log10(ifelse(is.null(pval_col), 0.05, !!sym(pval_col)))) %>%
      arrange(desc(extreme_score)) %>%
      head(25)
    
    if(nrow(label_data) > 0) {
      p <- p + geom_text_repel(data = label_data,
                               aes(x = !!sym(fc_col), 
                                   y = -log10(ifelse(is.null(pval_col), 0.05, !!sym(pval_col))), 
                                   label = !!sym(label_col)),
                               size = 3, max.overlaps = 20)
    }
  }
  
  # Highlight specific points if specified
  if(!is.null(highlight_col)) {
    # Filter only highlighted points
    highlighted_data <- plot_data %>%
      filter(!!sym(highlight_col) == TRUE)
    
    if(nrow(highlighted_data) > 0) {
      p <- p + 
        geom_point(data = highlighted_data, 
                   aes(x = !!sym(fc_col), 
                       y = -log10(ifelse(is.null(pval_col), 0.05, !!sym(pval_col)))), 
                   color = "purple", size = 3, shape = 17)
    }
  }
  
  return(p)
}

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