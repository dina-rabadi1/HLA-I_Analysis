#!/usr/bin/env Rscript

#' HLA Analysis in R
#' 
#' Combined script for HLA frequency analysis, visualization, and validation
#' Converted from Python modules
#' 
#' # Setting directory
setwd("~/Documents/Github/HLA-I_Analysis/")

# Define paths
RAW_DATA_DIR <- "rawdata"
OUTPUT_DIR <- "HLA_frequencies"

library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ComplexHeatmap)
library(viridis)
library(readxl) # For reading Excel files

#' ----------------------------------------
#' Data loading and preparation
#' ----------------------------------------

#' Load and prepare HLA data
#' @param input_file Path to input file
#' @return Prepared dataframe
load_hla_data <- function(input_file) {
  # Check file extension
  file_ext <- tools::file_ext(input_file)
  
  if (file_ext == "csv") {
    # Read CSV with stringsAsFactors=FALSE to avoid factor conversion issues
    df <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (file_ext %in% c("xlsx", "xls")) {
    # Read Excel file
    df <- as.data.frame(readxl::read_excel(input_file))
  } else {
    stop("Unsupported file format. Please provide a CSV or Excel file.")
  }
  
  # Print the column names for debugging
  cat("Columns in the dataset:", paste(colnames(df), collapse=", "), "\n")
  
  # Return the prepared dataframe
  return(df)
}

#' ----------------------------------------
#' Frequency calculation functions
#' ----------------------------------------

#' Calculate frequencies for HLA alleles
#' @param df Data frame with HLA typing data
#' @param locus HLA locus (A, B, or C)
#' @return Data frame with frequency calculations
calculate_frequencies <- function(df, locus) {
  # Count total number of patients
  total_patients <- nrow(df)
  cat("Total patients:", total_patients, "\n")
  
  # Extract all alleles for this locus (each patient has 2)
  col1 <- paste0(locus, "1")
  col2 <- paste0(locus, "2")
  
  # Check if columns exist
  if (!col1 %in% colnames(df) || !col2 %in% colnames(df)) {
    stop(paste0("Columns ", col1, " and/or ", col2, " not found in the dataset."))
  }
  
  # Extract alleles
  all_alleles <- c(df[[col1]], df[[col2]])
  all_alleles <- all_alleles[!is.na(all_alleles)]
  
  cat("Number of alleles for locus", locus, ":", length(all_alleles), "\n")
  
  # If no alleles, return empty dataframe with proper structure
  if (length(all_alleles) == 0) {
    return(data.frame(
      Allele_Frequency = numeric(0),
      Patient_Frequency = numeric(0),
      Allele_Percentage = numeric(0),
      Patient_Percentage = numeric(0)
    ))
  }
  
  # Calculate allele frequencies (proportion of all alleles)
  allele_counts <- table(all_alleles)
  total_alleles <- length(all_alleles)
  
  # Convert to data frame properly
  allele_freq <- data.frame(
    Allele = names(allele_counts),
    Allele_Frequency = as.numeric(allele_counts) / total_alleles,
    stringsAsFactors = FALSE
  )
  
  # Calculate patient frequencies (proportion of patients with allele)
  patient_freq <- data.frame(
    Allele = names(allele_counts),
    Patient_Frequency = 0,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(patient_freq)) {
    allele <- patient_freq$Allele[i]
    # Count patients with this allele in either position
    patients_with_allele <- sum(
      (df[[col1]] == allele) | (df[[col2]] == allele),
      na.rm = TRUE
    )
    patient_freq$Patient_Frequency[i] <- patients_with_allele / total_patients
  }
  
  # Merge the two frequency data frames
  result <- merge(allele_freq, patient_freq, by = "Allele")
  
  # Convert to percentages
  result$Allele_Percentage <- result$Allele_Frequency * 100
  result$Patient_Percentage <- result$Patient_Frequency * 100
  
  # Set row names to alleles and remove the Allele column
  rownames(result) <- result$Allele
  result$Allele <- NULL
  
  # Sort by patient percentage
  result <- result[order(result$Patient_Percentage, decreasing = TRUE),]
  
  # Print the frequency counts for debugging
  cat("Number of unique alleles for locus", locus, ":", nrow(result), "\n")
  if (nrow(result) > 0) {
    cat("Top allele:", rownames(result)[1], "with frequency:", result$Patient_Percentage[1], "%\n")
  }
  
  return(result)
}

#' Run frequency analysis for all loci
#' @param df Data frame with HLA typing data
#' @return List with frequency results and metadata
run_frequency_analysis <- function(df) {
  # Calculate frequencies for each locus
  results <- list()
  frequencies <- list()
  
  for (locus in c("A", "B", "C")) {
    cat("\nProcessing locus:", locus, "\n")
    frequencies[[locus]] <- calculate_frequencies(df, locus)
  }
  
  results$frequencies <- frequencies
  results$metadata <- list(
    sample_size = nrow(df),
    timestamp = Sys.time()
  )
  
  return(results)
}

#' ----------------------------------------
#' Visualization functions
#' ----------------------------------------

#' Calculate Shannon diversity index
#' @param freq_series Frequency series
#' @return Shannon diversity index
calculate_shannon_diversity <- function(freq_series) {
  # Convert to proportions if values are percentages
  if (sum(freq_series, na.rm = TRUE) > 1) {
    proportions <- freq_series / 100
  } else {
    proportions <- freq_series
  }
  
  # Ensure no zeros which would cause -Inf in log calculation
  proportions <- proportions[proportions > 0]
  
  # If no valid proportions, return 0
  if (length(proportions) == 0) {
    return(0)
  }
  
  # Calculate Shannon diversity index
  shannon <- -sum(proportions * log(proportions))
  return(shannon)
}

#' Create individual frequency plot with optional y-axis scaling
#' @param freq_df Frequency data frame
#' @param locus HLA locus
#' @param output_path Output file path
#' @param y_max Optional maximum for y-axis
#' @return List with summary statistics
create_frequency_plot <- function(freq_df, locus, output_path, y_max = NULL) {
  # Check if the dataframe is empty
  if (nrow(freq_df) == 0) {
    cat("No data available for locus", locus, "- skipping plot\n")
    return(list(
      unique_alleles = 0,
      shannon_diversity_allele = 0,
      shannon_diversity_patient = 0,
      top_allele = NA,
      top_frequency = 0
    ))
  }
  
  # Reshape data for ggplot
  plot_data <- data.frame(
    Allele = rownames(freq_df),
    Allele_Percentage = freq_df$Allele_Percentage,
    Patient_Percentage = freq_df$Patient_Percentage
  )
  
  plot_data_long <- reshape2::melt(plot_data, 
                                   id.vars = "Allele", 
                                   variable.name = "Type",
                                   value.name = "Percentage")
  
  # Calculate Shannon diversity
  shannon_allele <- calculate_shannon_diversity(freq_df$Allele_Frequency)
  shannon_patient <- calculate_shannon_diversity(freq_df$Patient_Frequency)
  
  # Create plot
  p <- ggplot(plot_data_long, aes(x = Allele, y = Percentage, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("Allele_Percentage" = "skyblue", 
                                 "Patient_Percentage" = "lightgreen")) +
    labs(title = paste0("HLA-", locus, " Frequencies"),
         subtitle = paste0("n=", nrow(freq_df), " alleles, Shannon Diversity: ",
                           sprintf("%.2f", shannon_allele), " (allele), ",
                           sprintf("%.2f", shannon_patient), " (patient)"),
         x = "Alleles", y = "Frequency (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Set y-axis limit if specified
  if (!is.null(y_max)) {
    p <- p + ylim(0, y_max)
  }
  
  # Save plot
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  ggsave(output_path, p, width = 12, height = 6)
  
  # Return summary stats
  return(list(
    unique_alleles = nrow(freq_df),
    shannon_diversity_allele = shannon_allele,
    shannon_diversity_patient = shannon_patient,
    top_allele = if(nrow(freq_df) > 0) rownames(freq_df)[which.max(freq_df$Patient_Percentage)] else NA,
    top_frequency = if(nrow(freq_df) > 0) max(freq_df$Patient_Percentage) else 0
  ))
}

#' Create sorted frequency plot with alleles ordered by frequency
#' @param freq_df Frequency data frame
#' @param locus HLA locus
#' @param output_path Output file path
#' @param y_max Optional maximum for y-axis
#' @return List with summary statistics
create_sorted_frequency_plot <- function(freq_df, locus, output_path, y_max = NULL) {
  # Check if the dataframe is empty
  if (nrow(freq_df) == 0) {
    cat("No data available for locus", locus, "- skipping sorted plot\n")
    return(list(
      unique_alleles = 0,
      shannon_diversity_allele = 0,
      shannon_diversity_patient = 0,
      top_allele = NA,
      top_frequency = 0
    ))
  }
  
  # Reshape data for ggplot
  plot_data <- data.frame(
    Allele = rownames(freq_df),
    Allele_Percentage = freq_df$Allele_Percentage,
    Patient_Percentage = freq_df$Patient_Percentage
  )
  
  # Sort by patient percentage (highest to lowest)
  plot_data <- plot_data[order(plot_data$Patient_Percentage, decreasing = TRUE),]
  
  # Set factor levels to maintain sort order in the plot
  plot_data$Allele <- factor(plot_data$Allele, levels = plot_data$Allele)
  
  plot_data_long <- reshape2::melt(plot_data, 
                                   id.vars = "Allele", 
                                   variable.name = "Type",
                                   value.name = "Percentage")
  
  # Calculate Shannon diversity
  shannon_allele <- calculate_shannon_diversity(freq_df$Allele_Frequency)
  shannon_patient <- calculate_shannon_diversity(freq_df$Patient_Frequency)
  
  # Create plot with sorted data
  p <- ggplot(plot_data_long, aes(x = Allele, y = Percentage, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("Allele_Percentage" = "skyblue", 
                                 "Patient_Percentage" = "lightgreen")) +
    labs(title = paste0("HLA-", locus, " Frequencies (Sorted by Frequency)"),
         subtitle = paste0("n=", nrow(freq_df), " alleles, Shannon Diversity: ",
                           sprintf("%.2f", shannon_allele), " (allele), ",
                           sprintf("%.2f", shannon_patient), " (patient)"),
         x = "Alleles", y = "Frequency (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Set y-axis limit if specified
  if (!is.null(y_max)) {
    p <- p + ylim(0, y_max)
  }
  
  # Save plot
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  ggsave(output_path, p, width = 12, height = 6)
  
  # Return summary stats
  return(list(
    unique_alleles = nrow(freq_df),
    shannon_diversity_allele = shannon_allele,
    shannon_diversity_patient = shannon_patient,
    top_allele = if(nrow(freq_df) > 0) rownames(freq_df)[which.max(freq_df$Patient_Percentage)] else NA,
    top_frequency = if(nrow(freq_df) > 0) max(freq_df$Patient_Percentage) else 0
  ))
}

#' Update create_visualizations function to include sorted frequency plots
#' @param results_dict Results dictionary from frequency analysis
#' @param output_dir Output directory
#' @param df Original data frame
create_visualizations <- function(results_dict, output_dir, df) {
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Find max frequency across all loci for consistent y-axis
  max_freq <- 0
  for (locus in names(results_dict$frequencies)) {
    freq_df <- results_dict$frequencies[[locus]]
    if (nrow(freq_df) > 0) {
      max_freq <- max(max_freq,
                      max(freq_df$Allele_Percentage),
                      max(freq_df$Patient_Percentage))
    }
  }
  y_max <- ceiling(max_freq / 5) * 5
  
  # Create both scaled and unscaled versions
  locus_summaries <- list()
  for (locus in names(results_dict$frequencies)) {
    freq_df <- results_dict$frequencies[[locus]]
    
    # Skip if no data
    if (nrow(freq_df) == 0) {
      cat("No data for locus", locus, "- skipping visualization\n")
      next
    }
    
    # Original unscaled version
    unscaled_path <- file.path(output_dir, paste0("HLA-", locus, "_frequencies_unscaled.png"))
    locus_summaries[[locus]] <- create_frequency_plot(freq_df, locus, unscaled_path)
    
    # Original scaled version
    scaled_path <- file.path(output_dir, paste0("HLA-", locus, "_frequencies_scaled.png"))
    create_frequency_plot(freq_df, locus, scaled_path, y_max)
    
    # NEW: Sorted unscaled version
    sorted_unscaled_path <- file.path(output_dir, paste0("HLA-", locus, "_frequencies_sorted_unscaled.png"))
    create_sorted_frequency_plot(freq_df, locus, sorted_unscaled_path)
    
    # NEW: Sorted scaled version
    sorted_scaled_path <- file.path(output_dir, paste0("HLA-", locus, "_frequencies_sorted_scaled.png"))
    create_sorted_frequency_plot(freq_df, locus, sorted_scaled_path, y_max)
  }
  
  # Create all other visualizations
  create_combined_visualizations(results_dict, df, output_dir)
  
  # Print summary statistics
  cat("\nHLA Allele Diversity Summary:\n")
  for (locus in names(locus_summaries)) {
    stats <- locus_summaries[[locus]]
    cat(paste0("\nHLA-", locus, ":\n"))
    cat(paste0("Number of unique alleles: ", stats$unique_alleles, "\n"))
    if (!is.na(stats$top_allele)) {
      cat(paste0("Most common allele: ", stats$top_allele, " (", 
                 sprintf("%.1f", stats$top_frequency), "% of patients)\n"))
    } else {
      cat("No alleles found for this locus\n")
    }
    cat(paste0("Shannon diversity index (allele): ", 
               sprintf("%.3f", stats$shannon_diversity_allele), "\n"))
    cat(paste0("Shannon diversity index (patient): ", 
               sprintf("%.3f", stats$shannon_diversity_patient), "\n"))
  }
}

#' Create top 20 frequency plot
#' @param data Frequency data
#' @param column Column to plot
#' @param freq_type Frequency type (Patient or Allele)
#' @param output_dir Output directory
create_top_20_plot <- function(data, column, freq_type, output_dir) {
  # Skip if no data
  if (nrow(data) == 0) {
    cat("No data for top 20", freq_type, "frequencies - skipping plot\n")
    return()
  }
  
  # Limit to top 20 if more rows
  if (nrow(data) > 20) {
    data <- head(data, 20)
  }
  
  # Create label with allele and locus
  data$Label <- paste(data$Locus, data$Allele, sep="-")
  
  p <- ggplot(data, aes_string(x = "reorder(Label, -get(column))", y = column)) +
    geom_bar(stat = "identity", fill = ifelse(freq_type == "Patient", 
                                              "lightgreen", "skyblue")) +
    geom_text(aes(label = sprintf("%.1f%%", get(column))), 
              vjust = -0.5, size = 3) +
    labs(title = paste0("Top ", nrow(data), " HLA Frequencies by ", freq_type, " Percentage"),
         x = "HLA Alleles", y = paste0(freq_type, " Frequency (%)")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plot
  output_path <- file.path(output_dir, paste0("top_", nrow(data), "_", 
                                              tolower(freq_type), "_frequencies.png"))
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  ggsave(output_path, p, width = 15, height = 8)
}

#' Create locus comparison heatmap
#' @param locus1 First HLA locus
#' @param locus2 Second HLA locus
#' @param freq_dict Frequency dictionary
#' @param df Original data frame
#' @param output_dir Output directory
create_locus_comparison_heatmap <- function(locus1, locus2, freq_dict, df, output_dir) {
  # Check if we have data for both loci
  if (nrow(freq_dict[[locus1]]) == 0 || nrow(freq_dict[[locus2]]) == 0) {
    cat("Insufficient data for locus comparison between", locus1, "and", locus2, 
        "- skipping heatmap\n")
    return()
  }
  
  # Get top 10 alleles from each locus (or all if fewer than 10)
  top_alleles1 <- head(freq_dict[[locus1]][order(freq_dict[[locus1]]$Patient_Percentage, 
                                                 decreasing = TRUE),], 
                       min(10, nrow(freq_dict[[locus1]])))
  
  top_alleles2 <- head(freq_dict[[locus2]][order(freq_dict[[locus2]]$Patient_Percentage, 
                                                 decreasing = TRUE),], 
                       min(10, nrow(freq_dict[[locus2]])))
  
  # Create matrix of co-occurrences
  matrix <- matrix(0, nrow = nrow(top_alleles1), ncol = nrow(top_alleles2))
  rownames(matrix) <- rownames(top_alleles1)
  colnames(matrix) <- rownames(top_alleles2)
  
  # Fill matrix
  for (i in 1:nrow(df)) {
    alleles1 <- c(df[[paste0(locus1, "1")]][i], df[[paste0(locus1, "2")]][i])
    alleles2 <- c(df[[paste0(locus2, "1")]][i], df[[paste0(locus2, "2")]][i])
    
    for (a1 in alleles1) {
      for (a2 in alleles2) {
        if (!is.na(a1) && !is.na(a2) && 
            a1 %in% rownames(matrix) && a2 %in% colnames(matrix)) {
          matrix[a1, a2] <- matrix[a1, a2] + 1
        }
      }
    }
  }
  
  # Convert to percentages
  matrix <- (matrix / nrow(df) * 100) |> round(0)
  
  # Create directory for output
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create heatmap with ComplexHeatmap
  tryCatch({
    # Create PDF version
    pdf_file <- file.path(output_dir, paste0("HLA_", locus1, "_vs_", locus2, "_heatmap.pdf"))
    pdf(pdf_file, width = 12, height = 8)
    
    # Set up colors for heatmap
    color_func <- colorRampPalette(c("white", "yellow", "orange", "red"))
    colors <- color_func(100)
    
    # Create heatmap with text labels
    heatmap_result <- ComplexHeatmap::Heatmap(
      matrix,
      name = "Co-occurrence (%)",
      col = colors,
      rect_gp = grid::gpar(col = "white", lwd = 1),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.text(matrix[i, j], x, y, gp = grid::gpar(fontsize = 10))
      },
      column_title = paste0("HLA-", locus2, " Alleles"),
      row_title = paste0("HLA-", locus1, " Alleles"),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_side = "left",
      column_names_side = "top"
    )
    
    # Draw the heatmap with title included in the heatmap
    ComplexHeatmap::draw(
      heatmap_result, 
      padding = unit(c(2, 2, 2, 2), "cm"),
      column_title = paste0("HLA-", locus1, " vs HLA-", locus2, " Co-occurrence (%)")
    )
    
    # Close the PDF device
    dev.off()
    
    # Create PNG version
    png_file <- file.path(output_dir, paste0("HLA_", locus1, "_vs_", locus2, "_heatmap.png"))
    png(png_file, width = 1200, height = 800, res = 100)
    
    # Draw the same heatmap for PNG
    ComplexHeatmap::draw(
      heatmap_result, 
      padding = unit(c(2, 2, 2, 2), "cm"),
      column_title = paste0("HLA-", locus1, " vs HLA-", locus2, " Co-occurrence (%)")
    )
    
    # Close the PNG device
    dev.off()
    
    cat("Created heatmap for HLA-", locus1, " vs HLA-", locus2, "\n")
  }, error = function(e) {
    cat("Error creating heatmap for HLA-", locus1, " vs HLA-", locus2, ":", conditionMessage(e), "\n")
  })
}

#' Create clustered heatmap of HLA correlations
#' @param freq_dict Frequency dictionary
#' @param df Original data frame
#' @param output_dir Output directory
create_clustered_heatmap <- function(freq_dict, df, output_dir) {
  # Combine frequencies from all loci
  all_freqs <- data.frame()
  
  for (locus in c("A", "B", "C")) {
    # Skip loci with no data
    if (nrow(freq_dict[[locus]]) == 0) {
      next
    }
    
    freq_df <- freq_dict[[locus]]
    temp_df <- data.frame(
      Allele = paste0(locus, "-", rownames(freq_df)),
      Patient_Percentage = freq_df$Patient_Percentage,
      stringsAsFactors = FALSE
    )
    all_freqs <- rbind(all_freqs, temp_df)
  }
  
  # Check if we have enough data
  if (nrow(all_freqs) < 2) {
    cat("Insufficient data for clustered heatmap - skipping\n")
    return()
  }
  
  # Get top 30 overall frequencies (or all if fewer than 30)
  top_freqs <- head(all_freqs[order(all_freqs$Patient_Percentage, decreasing = TRUE),], 
                    min(30, nrow(all_freqs)))
  
  # Create binary matrix for presence/absence of each allele
  binary_matrix <- matrix(0, nrow = nrow(df), ncol = nrow(top_freqs))
  colnames(binary_matrix) <- top_freqs$Allele
  
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(top_freqs)) {
      allele_info <- strsplit(top_freqs$Allele[j], "-")[[1]]
      locus <- allele_info[1]
      allele_name <- paste(allele_info[-1], collapse = "-") # In case allele has dash in name
      
      if (!is.na(df[[paste0(locus, "1")]][i]) && 
          df[[paste0(locus, "1")]][i] == allele_name) {
        binary_matrix[i, j] <- 1
      } else if (!is.na(df[[paste0(locus, "2")]][i]) && 
                 df[[paste0(locus, "2")]][i] == allele_name) {
        binary_matrix[i, j] <- 1
      }
    }
  }
  
  # Calculate correlation matrix
  corr_matrix <- cor(binary_matrix) |> round(2)
  
  # Create directory for output
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create clustered heatmap using tryCatch to handle errors
  tryCatch({
    # PDF version
    pdf_file <- file.path(output_dir, "HLA_clustered_correlations.pdf")
    pdf(pdf_file, width = 15, height = 12)
    
    # Set up colors for heatmap
    colors <- colorRampPalette(c("blue", "white", "red"))(100)
    
    # Create clustered heatmap
    heatmap_result <- ComplexHeatmap::Heatmap(
      corr_matrix,
      name = "Correlation",
      col = colors,
      rect_gp = grid::gpar(col = "white", lwd = 1),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.text(corr_matrix[i, j], x, y, gp = grid::gpar(fontsize = 8))
      },
      column_title = "HLA Alleles",
      row_title = "HLA Alleles",
      cluster_rows = TRUE,
      cluster_columns = TRUE
    )
    
    # Draw the heatmap with title included
    ComplexHeatmap::draw(
      heatmap_result, 
      padding = unit(c(2, 2, 2, 2), "cm"),
      column_title = "HLA Clustered Correlations"
    )
    
    # Close the PDF device
    dev.off()
    
    # PNG version
    png_file <- file.path(output_dir, "HLA_clustered_correlations.png")
    png(png_file, width = 1500, height = 1200, res = 100)
    
    # Draw again for PNG
    ComplexHeatmap::draw(
      heatmap_result, 
      padding = unit(c(2, 2, 2, 2), "cm"),
      column_title = "HLA Clustered Correlations"
    )
    
    # Close the PNG device
    dev.off()
    
    cat("Created clustered correlation heatmap\n")
  }, error = function(e) {
    cat("Error creating clustered heatmap:", conditionMessage(e), "\n")
  })
}

#' Create all combined visualizations
#' @param results_dict Results dictionary
#' @param df Original data frame
#' @param output_dir Output directory
create_combined_visualizations <- function(results_dict, df, output_dir) {
  # Create overall top 20 frequencies across all loci
  all_frequencies <- data.frame(
    Locus = character(),
    Allele = character(),
    Allele_Percentage = numeric(),
    Patient_Percentage = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (locus in names(results_dict$frequencies)) {
    # Skip loci with no data
    if (nrow(results_dict$frequencies[[locus]]) == 0) {
      next
    }
    
    freq_df <- results_dict$frequencies[[locus]]
    temp_df <- data.frame(
      Locus = paste0("HLA-", locus),
      Allele = rownames(freq_df),
      Allele_Percentage = freq_df$Allele_Percentage,
      Patient_Percentage = freq_df$Patient_Percentage,
      stringsAsFactors = FALSE
    )
    all_frequencies <- rbind(all_frequencies, temp_df)
  }
  
  # Skip if no data at all
  if (nrow(all_frequencies) == 0) {
    cat("No frequency data available - skipping combined visualizations\n")
    return()
  }
  
  # Top 20 Patient Frequencies (or all if fewer than 20)
  top_patient <- all_frequencies[order(all_frequencies$Patient_Percentage, 
                                       decreasing = TRUE),]
  create_top_20_plot(top_patient, "Patient_Percentage", "Patient", output_dir)
  
  # Top 20 Allele Frequencies (or all if fewer than 20)
  top_allele <- all_frequencies[order(all_frequencies$Allele_Percentage, 
                                      decreasing = TRUE),]
  create_top_20_plot(top_allele, "Allele_Percentage", "Allele", output_dir)
  
  # Create locus comparisons
  create_locus_comparison_heatmap("A", "B", results_dict$frequencies, df, output_dir)
  create_locus_comparison_heatmap("B", "C", results_dict$frequencies, df, output_dir)
  create_locus_comparison_heatmap("A", "C", results_dict$frequencies, df, output_dir)
  create_clustered_heatmap(results_dict$frequencies, df, output_dir)
}

#' ----------------------------------------
#' Validation functions
#' ----------------------------------------

#' Validate HLA frequency calculations and data integrity
#' @param df Data frame with HLA typing data
#' @param results_dict Results dictionary from frequency analysis
#' @return List with validation results
validate_calculations <- function(df, results_dict) {
  validation_results <- list()
  
  # 1. Check total patient counts
  total_patients <- nrow(df)
  validation_results$total_patients <- total_patients
  
  # 2. Validate allele counts for each locus
  for (locus in c("A", "B", "C")) {
    # Check if necessary columns exist
    col1 <- paste0(locus, "1")
    col2 <- paste0(locus, "2")
    
    if (!col1 %in% colnames(df) || !col2 %in% colnames(df)) {
      cat("Warning: Columns", col1, "and/or", col2, "not found. Skipping validation for locus", locus, "\n")
      next
    }
    
    # Each patient should have exactly 2 alleles per locus
    expected_allele_count <- total_patients * 2
    actual_allele_count <- sum(!is.na(df[[col1]])) + sum(!is.na(df[[col2]]))
    
    validation_results[[paste0(locus, "_allele_count_match")]] <- list(
      expected = expected_allele_count,
      actual = actual_allele_count,
      valid = expected_allele_count == actual_allele_count
    )
    
    # 3. Validate frequency calculations if data exists
    if (locus %in% names(results_dict$frequencies) && nrow(results_dict$frequencies[[locus]]) > 0) {
      freq_df <- results_dict$frequencies[[locus]]
      allele_freq_sum <- sum(freq_df$Allele_Percentage)
      validation_results[[paste0(locus, "_allele_freq_sum")]] <- list(
        sum = allele_freq_sum,
        valid = abs(100 - allele_freq_sum) < 0.01 # Should sum to 100%
      )
      
      # Patient frequencies can sum to >100% because patients can have multiple alleles
      patient_freq_sum <- sum(freq_df$Patient_Percentage)
      validation_results[[paste0(locus, "_patient_freq_sum")]] <- patient_freq_sum
    }
  }
  
  # 4. Validate haplotype counts (if available)
  if ("haplotypes" %in% names(results_dict) && "Count" %in% names(results_dict$haplotypes)) {
    total_haplotypes <- sum(results_dict$haplotypes$Count)
    validation_results$haplotype_count_match <- list(
      expected = total_patients,
      actual = total_haplotypes,
      valid = total_patients == total_haplotypes
    )
  }
  
  # 5. Validate group comparisons
  if ("Group" %in% colnames(df)) {
    for (group in unique(df$Group[!is.na(df$Group)])) {
      group_size <- sum(df$Group == group, na.rm = TRUE)
      validation_results[[paste0("group_", group, "_size")]] <- group_size
    }
  }
  
  return(validation_results)
}

#' Print a formatted validation report
#' @param validation_results Validation results list
print_validation_report <- function(validation_results) {
  cat("\n=== HLA Analysis Validation Report ===\n")
  cat(paste0("Total Patients: ", validation_results$total_patients, "\n"))
  
  cat("\n--- Allele Count Validation ---\n")
  for (locus in c("A", "B", "C")) {
    locus_key <- paste0(locus, "_allele_count_match")
    
    if (!locus_key %in% names(validation_results)) {
      cat(paste0("\nHLA-", locus, ": No validation data available\n"))
      next
    }
    
    result <- validation_results[[locus_key]]
    cat(paste0("\nHLA-", locus, ":\n"))
    cat(paste0("Expected alleles: ", result$expected, "\n"))
    cat(paste0("Actual alleles: ", result$actual, "\n"))
    cat(paste0("Valid: ", result$valid, "\n"))
    
    freq_sum_key <- paste0(locus, "_allele_freq_sum")
    if (freq_sum_key %in% names(validation_results)) {
      freq_sum <- validation_results[[freq_sum_key]]
      cat(paste0("Allele frequency sum: ", sprintf("%.2f", freq_sum$sum), "%\n"))
      cat(paste0("Valid frequency sum: ", freq_sum$valid, "\n"))
      
      patient_freq_key <- paste0(locus, "_patient_freq_sum")
      if (patient_freq_key %in% names(validation_results)) {
        cat(paste0("Patient frequency sum: ", 
                   sprintf("%.2f", validation_results[[patient_freq_key]]), "%\n"))
      }
    }
  }
  
  if ("haplotype_count_match" %in% names(validation_results)) {
    cat("\n--- Haplotype Validation ---\n")
    hap_result <- validation_results$haplotype_count_match
    cat(paste0("Expected haplotypes: ", hap_result$expected, "\n"))
    cat(paste0("Actual haplotypes: ", hap_result$actual, "\n"))
    cat(paste0("Valid: ", hap_result$valid, "\n"))
  }
  
  cat("\n--- Group Sizes ---\n")
  has_groups <- FALSE
  for (key in names(validation_results)) {
    if (grepl("^group_", key)) {
      has_groups <- TRUE
      group_name <- gsub("^group_", "", key)
      cat(paste0(group_name, ": ", validation_results[[key]], " patients\n"))
    }
  }
  
  # Fix the typo in line 721 (change TRUEs to TRUE)
  # Then add the rest of the code:
  
  if (!has_groups) {
    cat("No group information available\n")
  }
}

#' Main function to run HLA analysis
#' @param input_file Input file path (CSV or Excel)
#' @param output_dir Output directory
#' @return Results list
run_hla_analysis <- function(input_file, output_dir) {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load and prepare data
  cat("Loading data from:", input_file, "\n")
  df <- load_hla_data(input_file)
  
  # Print data summary
  cat("Data summary:\n")
  cat("Total patients:", nrow(df), "\n")
  cat("Columns:", paste(colnames(df), collapse=", "), "\n\n")
  
  # Run frequency analysis
  cat("Running frequency analysis...\n")
  results <- run_frequency_analysis(df)
  
  # Validate calculations
  cat("Validating calculations...\n")
  validation <- validate_calculations(df, results)
  print_validation_report(validation)
  
  # Create visualizations
  cat("Creating visualizations...\n")
  create_visualizations(results, output_dir, df)
  
  cat(paste0("\nAnalysis complete! Results saved to: ", output_dir, "\n"))
  
  # Return results for further analysis if needed
  return(list(
    results = results,
    validation = validation
  ))
}

# ----------------------------------------
# Run the analysis
# ----------------------------------------

# Try both file formats - first CSV, then Excel if CSV fails
csv_path <- file.path(RAW_DATA_DIR, "HLA-I_all.csv")
excel_path <- file.path(RAW_DATA_DIR, "HLA-I_all.xlsx")

cat("Starting HLA analysis...\n")

# Try to read the CSV file first
if (file.exists(csv_path)) {
  cat("Using CSV file:", csv_path, "\n")
  tryCatch({
    results <- run_hla_analysis(csv_path, OUTPUT_DIR)
    cat("CSV analysis complete!\n")
  }, error = function(e) {
    cat("Error with CSV file:", conditionMessage(e), "\n")
    cat("Trying Excel file instead...\n")
    
    if (file.exists(excel_path)) {
      results <- run_hla_analysis(excel_path, OUTPUT_DIR)
      cat("Excel analysis complete!\n")
    } else {
      stop("Excel file not found either.")
    }
  })
} else if (file.exists(excel_path)) {
  # If CSV doesn't exist, try Excel
  cat("CSV file not found. Using Excel file:", excel_path, "\n")
  results <- run_hla_analysis(excel_path, OUTPUT_DIR)
  cat("Excel analysis complete!\n")
} else {
  stop("Could not find input data files in the rawdata directory.")
}