# Visualization utilities for peptide analysis

# Helper function for safe column selection
# Add this at the top of your visualization_utils.R file

#' Safely select columns from a dataframe
#'
#' This function attempts to select columns from a dataframe even if they have 
#' slight variations in naming. It handles spaces vs underscores and skips columns
#' that don't exist without causing errors.
#'
#' @param df Data frame to select from
#' @param cols Character vector of column names to select
#' @return Data frame with selected columns
safe_select <- function(df, cols) {
  available_cols <- intersect(cols, colnames(df))
  if (length(available_cols) == 0) {
    warning("None of the requested columns exist in the data frame")
    return(df) # Return original dataframe if no columns match
  }
  return(df[, available_cols, drop = FALSE])
}

#' Create a bar chart of peptide distribution by sample count
#' 
#' @param sharing_summary Data frame with sample count distribution
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization file
create_sample_count_distribution <- function(sharing_summary, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Create the plot
  p <- ggplot2::ggplot(sharing_summary, 
                       ggplot2::aes(x = factor(sample_count), y = n,
                                    text = paste0("Sample count: ", sample_count,
                                                  "<br>Peptides: ", n,
                                                  "<br>Percentage: ", round(percentage, 1), "%"))) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(
      title = "Peptide Distribution by Sample Count",
      x = "Number of Samples",
      y = "Number of Peptides"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    )
  
  # Create file paths for each format
  file_paths <- list()
  
  # Save in each requested format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0("sample_count_distribution.", format))
    
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

#' Create a pie chart of unique vs shared peptides
#' 
#' @param unique_shared Data frame with unique vs shared counts
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Path to saved visualization file
create_unique_shared_pie <- function(unique_shared, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Create the plot
  p <- ggplot2::ggplot(unique_shared, ggplot2::aes(x = "", y = count, fill = category)) +
    ggplot2::geom_bar(width = 1, stat = "identity") +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::labs(
      title = "Unique vs Shared Peptides",
      fill = "Category"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(category, "\n", count, " (", round(percentage, 1), "%)")),
      position = ggplot2::position_stack(vjust = 0.5)
    ) +
    ggplot2::scale_fill_manual(values = c("lightblue", "lightgreen"))
  
  # Create file paths for each format
  file_paths <- list()
  
  # Save in each requested format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0("unique_vs_shared_pie.", format))
    
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

#' Create a bar chart of peptide length distribution
#' 
#' @param length_stats Data frame with length distribution
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Paths to saved visualization files
create_length_distribution_plot <- function(length_stats, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Create the plot
  p <- ggplot2::ggplot(length_stats, 
                       ggplot2::aes(x = factor(`Peptide Length`), y = count)) +
    ggplot2::geom_bar(stat = "identity", fill = "darkgreen") +
    ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.5, size = 3) +
    ggplot2::labs(
      title = "Peptide Length Distribution",
      subtitle = paste0("Total peptides: ", sum(length_stats$count)),
      x = "Peptide Length (amino acids)",
      y = "Number of Peptides"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    )
  
  # Create file paths for each format
  file_paths <- list()
  
  # Save in each requested format
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0("peptide_length_distribution.", format))
    
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

#' Create a heatmap of peptide length distribution by sample
#' 
#' @param length_by_sample Data frame with length by sample
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Paths to saved visualization files
create_length_by_sample_heatmap <- function(length_by_sample, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Create a matrix for the heatmap
  heatmap_data <- length_by_sample %>%
    tidyr::pivot_wider(
      names_from = `Peptide Length`,
      values_from = peptide_count,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("SampleID")
  
  # Save heatmap in each format
  file_paths <- list()
  
  for (format in output_formats) {
    file_path <- file.path(viz_dir, paste0("peptide_length_by_sample.", format))
    
    if (format == "pdf") {
      pdf(file_path, width = plot_width, height = plot_height)
    } else if (format == "png") {
      png(file_path, width = plot_width * png_dpi, height = plot_height * png_dpi, res = png_dpi)
    }
    
    # Create the heatmap
    pheatmap::pheatmap(
      as.matrix(heatmap_data),
      main = "Peptide Length Distribution by Sample",
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = TRUE,
      number_format = "%d",
      fontsize_number = 10,
      color = colorRampPalette(c("white", "steelblue", "darkblue"))(100)
    )
    
    dev.off()
    
    file_paths[[format]] <- file_path
  }
  
  return(file_paths)
}

#' Create bar chart comparing original vs spiked samples
#' 
#' @param spike_matrix_data Data frame with spike-in data
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Paths to saved visualization files
create_spike_comparison_plot <- function(spike_matrix_data, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Extract sample IDs from data
  original_sample <- unique(spike_matrix_data$SampleID)[1]
  spiked_sample <- unique(spike_matrix_data$SampleID)[2]
  
  # Make sure total_intensity is numeric
  spike_matrix_data$total_intensity <- as.numeric(spike_matrix_data$total_intensity)
  
  # Create the plot
  p <- ggplot2::ggplot(spike_matrix_data, 
                       ggplot2::aes(x = Peptide, y = total_intensity, 
                                    fill = SampleID)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::scale_fill_manual(values = c("steelblue", "orange"),
                               labels = c(paste0(original_sample, " (Original)"), 
                                          paste0(spiked_sample, " (Spiked)"))) +
    ggplot2::labs(
      title = paste("Spike-in Peptide Comparison:", original_sample, "vs", spiked_sample),
      subtitle = "Comparison of peptide intensities between original and spiked samples",
      x = "Peptide Sequence",
      y = "Total Intensity",
      fill = "Sample"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12),
      axis.title = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 12)
    )
  
  # Check if we have valid data for log scale
  if(all(spike_matrix_data$total_intensity >= 0) && sum(spike_matrix_data$total_intensity) > 0) {
    # Also create a log scale version
    p_log <- p + 
      ggplot2::scale_y_log10(labels = scales::scientific) +
      ggplot2::labs(
        subtitle = "Comparison of peptide intensities (log scale)",
        y = "Total Intensity (log scale)"
      )
  } else {
    # Create a version without log scale if data is zero or negative
    p_log <- p + 
      ggplot2::labs(
        subtitle = "No valid data for log scale plot",
        y = "Total Intensity"
      )
  }
  
  # Save regular scale plots
  file_paths <- list()
  
  for (format in output_formats) {
    # Regular scale
    file_path <- file.path(viz_dir, paste0("spike_comparison.", format))
    
    if (format == "pdf") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height)
    } else if (format == "png") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, dpi = png_dpi)
    }
    
    file_paths[[paste0("regular_", format)]] <- file_path
    
    # Log scale
    file_path_log <- file.path(viz_dir, paste0("spike_comparison_log.", format))
    
    if (format == "pdf") {
      ggplot2::ggsave(file_path_log, p_log, width = plot_width, height = plot_height)
    } else if (format == "png") {
      ggplot2::ggsave(file_path_log, p_log, width = plot_width, height = plot_height, dpi = png_dpi)
    }
    
    file_paths[[paste0("log_", format)]] <- file_path_log
  }
  
  return(file_paths)
}

#' Create a heatmap of spike-in peptide detection
#' 
#' @param spike_matrix_data Data frame with spike-in data
#' @param viz_dir Directory to save visualization
#' @param config Configuration settings
#' @return Paths to saved visualization files
create_spike_detection_heatmap <- function(spike_matrix_data, viz_dir, config) {
  # Extract viz settings
  plot_width <- config$visualization$plot_width
  plot_height <- config$visualization$plot_height
  output_formats <- config$visualization$output_formats
  png_dpi <- config$visualization$png_dpi
  
  # Create a matrix for presence/absence
  detection_matrix <- spike_matrix_data %>%
    dplyr::select(Peptide, SampleID, detected) %>%
    tidyr::pivot_wider(
      names_from = SampleID,
      values_from = detected
    ) %>%
    tibble::column_to_rownames("Peptide")
  
  # Create a matrix for intensity
  intensity_matrix <- spike_matrix_data %>%
    dplyr::select(Peptide, SampleID, total_intensity) %>%
    tidyr::pivot_wider(
      names_from = SampleID,
      values_from = total_intensity
    ) %>%
    tibble::column_to_rownames("Peptide")
  
  # Make sure data is numeric
  detection_matrix <- as.matrix(detection_matrix)
  mode(detection_matrix) <- "numeric"
  
  intensity_matrix <- as.matrix(intensity_matrix)
  intensity_matrix <- matrix(as.numeric(intensity_matrix), 
                             nrow = nrow(intensity_matrix),
                             dimnames = dimnames(intensity_matrix))
  
  # Log transform intensity for better visualization (adding 1 to avoid log(0))
  log_intensity_matrix <- log10(intensity_matrix + 1)
  
  # Save heatmaps in each format
  file_paths <- list()
  
  for (format in output_formats) {
    # Presence/absence heatmap
    p_file <- file.path(viz_dir, paste0("spike_detection_heatmap.", format))
    
    if (format == "pdf") {
      pdf(p_file, width = plot_width, height = plot_height)
    } else if (format == "png") {
      png(p_file, width = plot_width * png_dpi, height = plot_height * png_dpi, res = png_dpi)
    }
    
    # Check if matrix is valid before plotting
    if(nrow(detection_matrix) > 0 && ncol(detection_matrix) > 0) {
      pheatmap::pheatmap(
        detection_matrix,
        main = "Spike-in Peptide Detection",
        color = c("white", "darkblue"),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        display_numbers = TRUE,
        fontsize_number = 10
      )
    } else {
      # Create a simple error plot if matrix is invalid
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(1, 1, "Not enough data for heatmap", cex = 1.5)
    }
    
    dev.off()
    
    file_paths[[paste0("detection_", format)]] <- p_file
    
    # Intensity heatmap
    i_file <- file.path(viz_dir, paste0("spike_intensity_heatmap.", format))
    
    if (format == "pdf") {
      pdf(i_file, width = plot_width, height = plot_height)
    } else if (format == "png") {
      png(i_file, width = plot_width * png_dpi, height = plot_height * png_dpi, res = png_dpi)
    }
    
    # Check if matrix is valid before plotting
    if(nrow(log_intensity_matrix) > 0 && ncol(log_intensity_matrix) > 0) {
      pheatmap::pheatmap(
        log_intensity_matrix,
        main = "Spike-in Peptide Intensity (log10)",
        color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        display_numbers = TRUE,
        number_format = "%.1f",
        fontsize_number = 10
      )
    } else {
      # Create a simple error plot if matrix is invalid
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(1, 1, "Not enough data for heatmap", cex = 1.5)
    }
    
    dev.off()
    
    file_paths[[paste0("intensity_", format)]] <- i_file
  }
  
  return(file_paths)
}

# Functions to add to src/core/visualization_utils.R

#' Create fusion peptide detection heatmap
#'
#' @param fusion_matrix Matrix of fusion peptide presence/absence
#' @param output_dir Output directory for saving the plot
#' @param config Configuration list
#' @return Path to saved visualization
create_fusion_detection_heatmap <- function(fusion_matrices, output_dir, config) {
  cat("Creating fusion peptide detection heatmap...\n")
  
  # Extract matrices
  presence_matrix <- fusion_matrices$presence_matrix
  intensity_matrix <- fusion_matrices$intensity_matrix
  
  # If no fusion peptides found, return NULL
  if (nrow(presence_matrix) == 0) {
    cat("No fusion peptides to visualize\n")
    return(NULL)
  }
  
  # Prepare file names
  file_base <- file.path(output_dir, "fusion_peptide_heatmap")
  file_paths <- list()
  
  # Create heatmaps for each output format
  for (format in config$visualization$output_formats) {
    file_name <- paste0(file_base, ".", format)
    
    # Open device
    if (format == "pdf") {
      pdf(file_name, width = config$visualization$plot_width, height = config$visualization$plot_height)
    } else if (format == "png") {
      png(file_name, width = config$visualization$plot_width, height = config$visualization$plot_height, 
          units = "in", res = config$visualization$png_dpi)
    }
    
    # Create presence/absence heatmap with pheatmap
    heatmap_colors <- c("white", "darkblue")
    heatmap_title <- "DNAJB1:PRKACA Fusion Peptide Detection Across Samples"
    
    pheatmap(
      presence_matrix,
      color = heatmap_colors,
      main = heatmap_title,
      cluster_rows = FALSE,  # Don't cluster rows to preserve peptide order
      cluster_cols = TRUE,   # Cluster columns to group similar samples
      legend = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 10
    )
    
    # Close device
    dev.off()
    
    file_paths[["presence"]] <- file_name
    cat(sprintf("Saved %s presence heatmap to %s\n", format, file_name))
    
    # Create intensity heatmap if intensities are available
    if (any(intensity_matrix > 0)) {
      # Log-transform intensities for better visualization
      intensity_matrix_log <- log10(intensity_matrix + 1)
      
      # Prepare file name for intensity heatmap
      intensity_file_name <- paste0(file.path(output_dir, "fusion_peptide_intensity"), ".", format)
      
      # Open device
      if (format == "pdf") {
        pdf(intensity_file_name, width = config$visualization$plot_width, height = config$visualization$plot_height)
      } else if (format == "png") {
        png(intensity_file_name, width = config$visualization$plot_width, height = config$visualization$plot_height, 
            units = "in", res = config$visualization$png_dpi)
      }
      
      # Create intensity heatmap
      intensity_title <- "DNAJB1:PRKACA Fusion Peptide Intensity Across Samples (log10 scale)"
      
      pheatmap(
        intensity_matrix_log,
        main = intensity_title,
        cluster_rows = FALSE, 
        cluster_cols = FALSE, 
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize_row = 8,
        fontsize_col = 10
      )
      
      # Close device
      dev.off()
      
      file_paths[["intensity"]] <- intensity_file_name
      cat(sprintf("Saved %s intensity heatmap to %s\n", format, intensity_file_name))
    }
  }
  
  return(file_paths)
}

#' Create fusion peptide sequence visualization
#'
#' @param fusion_results Dataframe of fusion peptides
#' @param output_dir Output directory for saving the plot
#' @param config Configuration list
#' @return Path to saved visualization
create_fusion_sequence_plot <- function(fusion_results, fusion_peptides, output_dir, config) {
  cat("Creating fusion peptide sequence visualization...\n")
  
  # If no fusion peptides found, return NULL
  if (nrow(fusion_results) == 0) {
    cat("No fusion peptides to visualize\n")
    return(NULL)
  }
  
  # Define the fusion sequences
  dnajb1_seq <- "RKREIFDRYGE"
  prkaca_seq <- "VKEFLAKAKED"
  
  # Full fusion sequence for reference
  full_fusion_seq <- paste0(dnajb1_seq, "|", prkaca_seq)
  junction_point <- nchar(dnajb1_seq) + 1  # Position of the junction
  
  # Get unique peptides with their details
  unique_peptides <- fusion_results %>%
    select(Peptide, seq1_part, seq2_part, seq1_contribution, seq2_contribution) %>%
    distinct()
  
  # Sort by position and length
  unique_peptides <- unique_peptides %>%
    arrange(desc(seq1_contribution), seq2_contribution, desc(nchar(Peptide)))
  
  # Count samples for each peptide
  peptide_sample_counts <- fusion_results %>%
    group_by(Peptide) %>%
    summarize(sample_count = n_distinct(SampleID)) %>%
    arrange(desc(sample_count))
  
  # Join with unique_peptides
  unique_peptides <- unique_peptides %>%
    left_join(peptide_sample_counts, by = "Peptide")
  
  # Prepare file names
  file_base <- file.path(output_dir, "fusion_peptide_sequence")
  file_paths <- list()
  
  # Create plot for each output format
  for (format in config$visualization$output_formats) {
    file_name <- paste0(file_base, ".", format)
    
    # Open device
    if (format == "pdf") {
      pdf(file_name, width = config$visualization$plot_width, height = max(8, nrow(unique_peptides) * 0.3))
    } else if (format == "png") {
      png(file_name, width = config$visualization$plot_width, height = max(8, nrow(unique_peptides) * 0.3), 
          units = "in", res = config$visualization$png_dpi)
    }
    
    # Set up plotting area with larger margins
    par(mar = c(4, 10, 6, 10))
    
    # Calculate total width for plot
    total_width <- nchar(full_fusion_seq)
    
    # Create empty plot
    plot(1, 1, type = "n", 
         xlim = c(0, total_width), 
         ylim = c(0, nrow(unique_peptides) + 1),
         xlab = "Amino Acid Position", ylab = "", 
         main = "DNAJB1:PRKACA Fusion Peptides",
         yaxt = "n")
    
    # Add full sequence at the top
    text(total_width/2, nrow(unique_peptides) + 0.8, 
         labels = paste("Fusion Sequence:", full_fusion_seq), 
         cex = 0.9, font = 2)
    
    # Add vertical line at junction
    abline(v = junction_point - 0.5, col = "red", lty = 2)
    text(junction_point - 0.5, nrow(unique_peptides) + 0.4, 
         labels = "Fusion Point", col = "red", cex = 0.8)
    
    # Add peptide labels on y-axis
    axis(2, at = 1:nrow(unique_peptides), 
         labels = unique_peptides$Peptide, 
         las = 1, cex.axis = 0.8)
    
    # Add sample count on right y-axis
    axis(4, at = 1:nrow(unique_peptides), 
         labels = paste(unique_peptides$sample_count, "samples"), 
         las = 1, cex.axis = 0.7, col.axis = "darkgreen")
    
    # Draw peptides
    for (i in 1:nrow(unique_peptides)) {
      # Get peptide parts and position info
      peptide <- unique_peptides$Peptide[i]
      seq1_part <- unique_peptides$seq1_part[i]
      seq2_part <- unique_peptides$seq2_part[i]
      seq1_contribution <- unique_peptides$seq1_contribution[i]
      
      # Calculate starting position (right-aligned to junction)
      start_pos <- junction_point - seq1_contribution
      
      # Draw line for DNAJB1 part (blue)
      segments(start_pos, i, junction_point - 0.5, i, col = "blue", lwd = 3)
      
      # Draw line for PRKACA part (red)
      segments(junction_point - 0.5, i, start_pos + nchar(peptide), i, col = "red", lwd = 3)
      
      # Add peptide text
      text(start_pos + nchar(peptide)/2, i, peptide, cex = 0.7, pos = 3)
    }
    
    # Add legend
    legend("bottomright", 
           legend = c("DNAJB1", "PRKACA", "Fusion Point"), 
           col = c("blue", "red", "red"),
           lwd = c(3, 3, 1), 
           lty = c(1, 1, 2),
           bg = "white",
           cex = 0.8)
    
    # Close device
    dev.off()
    
    file_paths[[format]] <- file_name
    cat(sprintf("Saved %s sequence plot to %s\n", format, file_name))
  }
  
  return(file_paths)
}

#' Create fusion peptide distribution bar plot
#'
#' @param fusion_metrics List of fusion peptide metrics
#' @param output_dir Output directory for saving the plot
#' @param config Configuration list
#' @return Path to saved visualization
create_fusion_distribution_plot <- function(fusion_metrics, output_dir, config) {
  cat("Creating fusion peptide distribution plot...\n")
  
  # If no fusion peptides found, return NULL
  if (fusion_metrics$total_fusion_peptides == 0) {
    cat("No fusion peptides to visualize\n")
    return(NULL)
  }
  
  # Prepare file names
  file_base <- file.path(output_dir, "fusion_peptide_distribution")
  file_paths <- list()
  
  # Create plot for each output format
  for (format in config$visualization$output_formats) {
    file_name <- paste0(file_base, ".", format)
    
    # Open device
    if (format == "pdf") {
      pdf(file_name, width = config$visualization$plot_width, height = config$visualization$plot_height)
    } else if (format == "png") {
      png(file_name, width = config$visualization$plot_width, height = config$visualization$plot_height, 
          units = "in", res = config$visualization$png_dpi)
    }
    
    # Create bar plot of sample counts
    p <- fusion_metrics$sample_counts %>%
      ggplot(aes(x = reorder(SampleID, -unique_peptides), y = unique_peptides)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      labs(
        title = "DNAJB1:PRKACA Fusion Peptides by Sample",
        x = "Sample",
        y = "Number of Unique Fusion Peptides"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p)
    
    # Close device
    dev.off()
    
    file_paths[[format]] <- file_name
    cat(sprintf("Saved %s distribution plot to %s\n", format, file_name))
  }
  
  return(file_paths)
}