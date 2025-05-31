# Script to correlate peptide intensity peptide data from 2cv and 3cv files
# This uses the original imp_014 data files

# Setting directory
setwd("~/Documents/Github/HLA-I_Analysis/")

# Load required packages
library(tidyverse)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(gridExtra)

directory <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata/imp_014_rawdata"

# Function to list and categorize TSV files
list_tsv_files <- function(directory) {
  # List all files in directory
  all_files <- list.files(path = directory, pattern = ".*_peptides\\.tsv$", full.names = TRUE)
  
  if (length(all_files) == 0) {
    stop("No peptide TSV files found in ", directory)
  }
  
  # Categorize by CV type
  files_2cv <- all_files[grepl("_2CV_", all_files)]
  files_3cv <- all_files[grepl("_3CV_", all_files)]
  
  cat("Found", length(files_2cv), "2cv files and", length(files_3cv), "3cv files\n")
  
  return(list(files_2cv = files_2cv, files_3cv = files_3cv))
}

# Function to extract sample ID from filename
extract_sample_id <- function(filename) {
  # Extract sample ID using regex
  # Pattern like: DDA_2CV_117_peptides.tsv
  filename <- basename(filename)
  match <- str_extract(filename, "DDA_[23]CV_([^_]+)")
  sample_id <- gsub("DDA_[23]CV_", "", match)
  
  if (is.na(sample_id)) {
    return("unknown")
  }
  
  return(sample_id)
}

# Function to read peptide files
read_peptide_files <- function(file_list, cv_type) {
  if (length(file_list) == 0) {
    return(NULL)
  }
  
  all_data <- list()
  
  for (i in seq_along(file_list)) {
    file <- file_list[i]
    cat("Reading", cv_type, "file:", basename(file), "\n")
    
    # Read the TSV file
    data <- read_tsv(file, show_col_types = FALSE)
    
    # Add metadata columns
    data$SourceFile <- basename(file)
    data$SampleID <- extract_sample_id(file)
    data$CVType <- cv_type
    
    all_data[[i]] <- data
  }
  
  # Combine all data frames
  combined_data <- bind_rows(all_data)
  
  return(combined_data)
}

# Function to perform correlation analysis for matched peptides
perform_correlation_analysis <- function(data_2cv, data_3cv) {
  # Get unique sample IDs
  samples_2cv <- unique(data_2cv$SampleID)
  samples_3cv <- unique(data_3cv$SampleID)
  
  # Find common sample IDs
  common_samples <- intersect(samples_2cv, samples_3cv)
  
  if (length(common_samples) == 0) {
    stop("No common sample IDs found between 2CV and 3CV datasets")
  }
  
  cat("Found", length(common_samples), "common sample IDs\n")
  
  # Initialize results list
  results <- list(
    by_sample = list(),
    overall = NULL
  )
  
  # Initialize dataframe to store all matched peptides for overall analysis
  all_matched_peptides <- data.frame()
  
  # Process each sample
  for (sample_id in common_samples) {
    cat("Processing sample:", sample_id, "\n")
    
    # Filter data for this sample
    sample_2cv <- data_2cv %>% filter(SampleID == sample_id)
    sample_3cv <- data_3cv %>% filter(SampleID == sample_id)
    
    # Find common peptides
    peptides_2cv <- sample_2cv$Peptide
    peptides_3cv <- sample_3cv$Peptide
    common_peptides <- intersect(peptides_2cv, peptides_3cv)
    
    cat("  Found", length(common_peptides), "peptides in both 2CV and 3CV\n")
    
    if (length(common_peptides) < 5) {
      warning("Too few common peptides for sample ", sample_id, ". Skipping.")
      next
    }
    
    # Create matched data
    matched_data <- data.frame(
      Peptide = common_peptides,
      SampleID = sample_id,
      Intensity_2CV = sample_2cv %>% 
        filter(Peptide %in% common_peptides) %>% 
        select(Peptide, Intensity) %>% 
        rename(Intensity_2CV = Intensity) %>% 
        pull(Intensity_2CV),
      Intensity_3CV = sample_3cv %>% 
        filter(Peptide %in% common_peptides) %>% 
        select(Peptide, Intensity) %>% 
        rename(Intensity_3CV = Intensity) %>% 
        pull(Intensity_3CV)
    )
    
    # Add this after reading in your data but before correlation analysis
    # Check for negative or zero intensity values in 2CV
    negative_2cv <- sum(data_2cv$Intensity < 0)
    zero_2cv <- sum(data_2cv$Intensity == 0)
    cat("Found", negative_2cv, "negative intensity values and", zero_2cv, "zero values in 2CV data\n")
    
    # Check for negative or zero intensity values in 3CV
    negative_3cv <- sum(data_3cv$Intensity < 0)
    zero_3cv <- sum(data_3cv$Intensity == 0)
    cat("Found", negative_3cv, "negative intensity values and", zero_3cv, "zero values in 3CV data\n")
    
    # If there are negative values, let's see their distribution
    if (negative_2cv > 0) {
      cat("Range of negative 2CV intensities:", range(data_2cv$Intensity[data_2cv$Intensity < 0]), "\n")
    }
    if (negative_3cv > 0) {
      cat("Range of negative 3CV intensities:", range(data_3cv$Intensity[data_3cv$Intensity < 0]), "\n")
    }
    
    # Calculate correlations
    pearson_cor <- cor(matched_data$Intensity_2CV, matched_data$Intensity_3CV, method = "pearson")
    spearman_cor <- cor(matched_data$Intensity_2CV, matched_data$Intensity_3CV, method = "spearman")
    
    # Calculate R-squared (coefficient of determination)
    lm_model <- lm(Intensity_3CV ~ Intensity_2CV, data = matched_data)
    r_squared <- summary(lm_model)$r.squared
    
    # Store results for this sample
    results$by_sample[[sample_id]] <- list(
      matched_data = matched_data,
      pearson = pearson_cor,
      spearman = spearman_cor,
      r_squared = r_squared
    )
    
    # Add to overall dataset
    all_matched_peptides <- bind_rows(all_matched_peptides, matched_data)
  }
  
  # Overall analysis if we have data
  if (nrow(all_matched_peptides) > 0) {
    # Calculate overall correlations
    pearson_cor_overall <- cor(all_matched_peptides$Intensity_2CV, all_matched_peptides$Intensity_3CV, method = "pearson")
    spearman_cor_overall <- cor(all_matched_peptides$Intensity_2CV, all_matched_peptides$Intensity_3CV, method = "spearman")
    
    # Calculate overall R-squared
    lm_model_overall <- lm(Intensity_3CV ~ Intensity_2CV, data = all_matched_peptides)
    r_squared_overall <- summary(lm_model_overall)$r.squared
    
    results$overall <- list(
      matched_data = all_matched_peptides,
      pearson = pearson_cor_overall,
      spearman = spearman_cor_overall,
      r_squared = r_squared_overall
    )
  }
  
  return(results)
}

# Function to create visualization
create_visualizations <- function(correlation_results) {
  # Create directory for plots if it doesn't exist
  dir.create("plots", showWarnings = FALSE)
  
  create_scatter_plot <- function(data, title, filename) {
    # Filter out zeros and negative values for log transformation
    data_filtered <- data %>% filter(Intensity_2CV > 0, Intensity_3CV > 0)
    
    # Create scatter plot
    p <- ggplot(data_filtered, aes(x = Intensity_2CV, y = Intensity_3CV)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", color = "red") +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        title = title,
        subtitle = sprintf(
          "Pearson: %.3f, Spearman: %.3f, R²: %.3f (n=%d, filtered out %d zeros/negatives)",
          cor(data_filtered$Intensity_2CV, data_filtered$Intensity_3CV, method = "pearson"),
          cor(data_filtered$Intensity_2CV, data_filtered$Intensity_3CV, method = "spearman"),
          summary(lm(Intensity_3CV ~ Intensity_2CV, data = data_filtered))$r.squared,
          nrow(data_filtered),
          nrow(data) - nrow(data_filtered)
        ),
        x = "2CV Intensity (log10)",
        y = "3CV Intensity (log10)"
      ) +
      theme_minimal()
    
    ggsave(filename, plot = p, width = 10, height = 8)
    return(p)
  }
  
  # Create plots for each sample
  sample_plots <- list()
  for (sample_id in names(correlation_results$by_sample)) {
    data <- correlation_results$by_sample[[sample_id]]$matched_data
    title <- paste("Peptide Intensity Correlation for Sample", sample_id)
    filename <- paste0("plots/scatter_", sample_id, ".png")
    
    sample_plots[[sample_id]] <- create_scatter_plot(data, title, filename)
    cat("Created scatter plot for sample", sample_id, "\n")
  }
  
  # Create overall plot
  if (!is.null(correlation_results$overall)) {
    overall_data <- correlation_results$overall$matched_data
    title <- "Overall Peptide Intensity Correlation (All Samples)"
    filename <- "plots/scatter_overall.png"
    
    overall_plot <- create_scatter_plot(overall_data, title, filename)
    cat("Created overall scatter plot\n")
    
    # Create heatmap
    # For the heatmap, we'll bin the data into a 2D histogram
    # First, filter out zeros and negative values before log transformation
    heatmap_data <- overall_data %>%
      filter(Intensity_2CV > 0, Intensity_3CV > 0) %>%  # Filter out non-positive values
      mutate(
        log_2CV = log10(Intensity_2CV),
        log_3CV = log10(Intensity_3CV)
      ) %>%
      mutate(
        bin_x = cut(log_2CV, breaks = 50),
        bin_y = cut(log_3CV, breaks = 50)
      ) %>%
      group_by(bin_x, bin_y) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(
        bin_x_numeric = as.numeric(as.character(bin_x)),  # More robust conversion
        bin_y_numeric = as.numeric(as.character(bin_y))
      )
    
    heatmap_plot <- ggplot(heatmap_data, aes(x = bin_x_numeric, y = bin_y_numeric, fill = count)) +
      geom_tile() +
      scale_fill_viridis_c(trans = "log1p") +
      labs(
        title = "Heatmap of Peptide Intensity Correlation",
        x = "2CV Intensity (log10, binned)",
        y = "3CV Intensity (log10, binned)",
        fill = "Count"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank())
    
    ggsave("plots/heatmap_overall.png", plot = heatmap_plot, width = 10, height = 8)
    cat("Created overall heatmap\n")
  }
  
  # Create summary table
  summary_table <- data.frame(
    SampleID = character(),
    PeptideCount = integer(),
    Pearson = numeric(),
    Spearman = numeric(),
    R_Squared = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (sample_id in names(correlation_results$by_sample)) {
    result <- correlation_results$by_sample[[sample_id]]
    summary_table <- bind_rows(summary_table, data.frame(
      SampleID = sample_id,
      PeptideCount = nrow(result$matched_data),
      Pearson = result$pearson,
      Spearman = result$spearman,
      R_Squared = result$r_squared
    ))
  }
  
  if (!is.null(correlation_results$overall)) {
    summary_table <- bind_rows(summary_table, data.frame(
      SampleID = "OVERALL",
      PeptideCount = nrow(correlation_results$overall$matched_data),
      Pearson = correlation_results$overall$pearson,
      Spearman = correlation_results$overall$spearman,
      R_Squared = correlation_results$overall$r_squared
    ))
  }
  
  # Save summary table
  write.csv(summary_table, "correlation_summary.csv", row.names = FALSE)
  cat("Created correlation summary table\n")
  
  return(list(
    sample_plots = sample_plots,
    overall_plot = if (!is.null(correlation_results$overall)) overall_plot else NULL,
    heatmap = if (!is.null(correlation_results$overall)) heatmap_plot else NULL,
    summary_table = summary_table
  ))
}

# Main function to run the analysis
run_peptide_correlation_analysis <- function(directory = "data") {
  # List TSV files
  files <- list_tsv_files(directory)
  
  # Read peptide files
  data_2cv <- read_peptide_files(files$files_2cv, "2CV")
  data_3cv <- read_peptide_files(files$files_3cv, "3CV")
  
  if (is.null(data_2cv) || is.null(data_3cv)) {
    stop("Missing data for one or both CV types")
  }
  
  cat("Read", nrow(data_2cv), "rows from 2CV files and", nrow(data_3cv), "rows from 3CV files\n")
  
  # Perform correlation analysis
  correlation_results <- perform_correlation_analysis(data_2cv, data_3cv)
  
  # Create visualizations
  visualization_results <- create_visualizations(correlation_results)
  
  # Return results
  return(list(
    correlation = correlation_results,
    visualization = visualization_results
  ))
}

# Run the analysis
results <- run_peptide_correlation_analysis(directory)

# Improved PDF creation function with better heatmap handling
create_summary_pdf_robust <- function(correlation_results, visualization_results) {
  # Load required packages
  library(ggplot2)
  library(gridExtra)
  library(grid)
  
  # Create plots directory if it doesn't exist
  dir.create("plots", showWarnings = FALSE)
  
  # Set the PDF path in the plots directory
  pdf_path <- file.path("plots", "peptide_correlation_summary.pdf")
  cat("Creating PDF at:", pdf_path, "\n")
  
  # Create PDF with try-catch to handle errors
  tryCatch({
    # Start PDF device
    pdf(pdf_path, width = 11, height = 11)
    
    # Title page
    plot.new()
    text(x = 0.5, y = 0.7, "Peptide Intensity Correlation Analysis Summary", 
         cex = 2, font = 2)
    text(x = 0.5, y = 0.5, paste("Analysis Date:", format(Sys.Date(), "%B %d, %Y")),
         cex = 1.2)
    
    # Summary statistics table
    plot.new()
    grid.text("Correlation Statistics Summary", x = 0.5, y = 0.95, 
              gp = gpar(fontsize = 18, fontface = "bold"))
    
    # Format numeric columns to 3 decimal places
    summary_table <- visualization_results$summary_table
    for(col in c("Pearson", "Spearman", "R_Squared")) {
      summary_table[[col]] <- sprintf("%.3f", summary_table[[col]])
    }
    
    # Draw table using grid.table for simplicity
    grid.table(summary_table, rows = NULL)
    
    # Overall correlation plot
    if (!is.null(visualization_results$overall_plot)) {
      print(visualization_results$overall_plot + 
              ggtitle("Overall Peptide Intensity Correlation (All Samples)"))
    }
    
    # Sample plots - 4 per page
    sample_ids <- names(visualization_results$sample_plots)
    num_samples <- length(sample_ids)
    
    # Process in batches of 4
    for (i in seq(1, num_samples, by = 4)) {
      # Get current batch
      end_idx <- min(i + 3, num_samples)
      current_batch <- sample_ids[i:end_idx]
      
      # Create a list of plots for this page
      current_plots <- list()
      for (j in 1:length(current_batch)) {
        current_plots[[j]] <- visualization_results$sample_plots[[current_batch[j]]]
      }
      
      # Add dummy plots if needed to fill the grid
      while (length(current_plots) < 4) {
        current_plots[[length(current_plots) + 1]] <- ggplot() + theme_void()
      }
      
      # Print the grid of plots
      grid.arrange(
        grobs = current_plots,
        ncol = 2,
        top = "Sample Correlation Plots"
      )
    }
    
    # Create a much simpler heatmap that won't cause issues
    if (!is.null(correlation_results$overall$matched_data)) {
      # Safe heatmap creation
      tryCatch({
        cat("Creating simplified heatmap...\n")
        
        # Get data
        overall_data <- correlation_results$overall$matched_data
        
        # Filter out zeros and log transform safely
        heat_data <- overall_data %>%
          filter(Intensity_2CV > 0, Intensity_3CV > 0) %>%
          mutate(
            log10_2CV = log10(Intensity_2CV),
            log10_3CV = log10(Intensity_3CV)
          )
        
        # Use hexbin instead of tile for better performance
        simple_heatmap <- ggplot(heat_data, aes(x = log10_2CV, y = log10_3CV)) +
          stat_bin2d(bins = 30) +  # Use stat_bin2d instead of manually creating bins
          scale_fill_viridis_c(trans = "log1p", name = "Count") +
          labs(
            title = "Peptide Intensity Correlation Heatmap",
            x = "2CV Intensity (log10)",
            y = "3CV Intensity (log10)"
          ) +
          theme_minimal() +
          theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
        
        print(simple_heatmap)
        cat("Heatmap created successfully\n")
      }, error = function(e) {
        # If heatmap creation fails, just show a message and continue
        cat("Error creating heatmap:", e$message, "\n")
        cat("Skipping heatmap and continuing with PDF creation\n")
        
        # Create an empty plot with an error message
        plot.new()
        text(0.5, 0.5, "Heatmap could not be created due to an error", cex = 1.5)
        text(0.5, 0.4, "See console for details", cex = 1.2)
      })
    }
    
    # Interpretation page
    plot.new()
    text(x = 0.5, y = 0.95, "Interpretation of Results", cex = 1.5, font = 2)
    
    # Get values for interpretation safely
    pearson_overall <- as.numeric(summary_table$Pearson[nrow(summary_table)])
    r_squared_overall <- as.numeric(summary_table$R_Squared[nrow(summary_table)])
    pearson_values <- as.numeric(summary_table$Pearson[-nrow(summary_table)])
    pearson_min <- min(pearson_values)
    pearson_max <- max(pearson_values)
    
    # Create interpretation text
    interpretation <- paste(
      "Summary of findings:",
      "",
      "1. Overall correlation: Strong correlation between 2CV and 3CV samples",
      paste("   Pearson correlation:", pearson_overall),
      paste("   R-squared:", r_squared_overall),
      "",
      "2. Sample variability: Sample correlations range from",
      paste("   ", pearson_min, "to", pearson_max, "(Pearson)"),
      "",
      "3. Data filtering: Zero intensity values were excluded from log-scale visualizations",
      "",
      "4. Heatmap: Diagonal pattern confirms strong correlation between 2CV and 3CV",
      sep = "\n"
    )
    
    # Add interpretation text
    text(x = 0.1, y = 0.8, interpretation, adj = c(0, 1), cex = 0.9, family = "mono")
    
    # Close PDF device
    dev.off()
    
    cat("Summary PDF successfully created at:", pdf_path, "\n")
  }, 
  error = function(e) {
    # Make sure PDF device is closed if an error occurs
    if (names(dev.cur()) != "null device") dev.off()
    cat("Error creating PDF:", e$message, "\n")
  })
}

create_summary_pdf_robust(results$correlation, results$visualization)

# Function to create a polished PDF with no cut-off elements
create_summary_pdf_polished <- function(correlation_results, visualization_results) {
  library(ggplot2)
  library(gridExtra)
  library(grid)
  
  # Create plots directory if it doesn't exist
  dir.create("plots", showWarnings = FALSE)
  
  # Set the PDF path in the plots directory
  pdf_path <- "plots/peptide_correlation_summary.pdf"
  cat("Creating polished PDF at:", pdf_path, "\n")
  
  # Create PDF with larger dimensions and margins
  pdf(pdf_path, width = 11, height = 11)
  
  # Title page
  grid.newpage()
  grid.text("Peptide Intensity Correlation Analysis Summary", 
            x = 0.5, y = 0.7, gp = gpar(fontsize = 24, fontface = "bold"))
  grid.text(paste("Analysis Date:", format(Sys.Date(), "%B %d, %Y")),
            x = 0.5, y = 0.5, gp = gpar(fontsize = 16))
  
  # Summary statistics table
  grid.newpage()
  grid.text("Correlation Statistics Summary", x = 0.5, y = 0.95, 
            gp = gpar(fontsize = 18, fontface = "bold"))
  
  # Format numeric columns to 3 decimal places
  summary_table <- visualization_results$summary_table
  summary_table$Pearson <- sprintf("%.3f", summary_table$Pearson)
  summary_table$Spearman <- sprintf("%.3f", summary_table$Spearman)
  summary_table$R_Squared <- sprintf("%.3f", summary_table$R_Squared)
  
  # Create table with better spacing
  table_theme <- ttheme_minimal(
    core = list(fg_params = list(fontsize = 10, hjust = 0.5, x = 0.5),
                bg_params = list(fill = c("white", "grey95"), col = NA)),
    colhead = list(fg_params = list(fontsize = 11, fontface = "bold", hjust = 0.5, x = 0.5),
                   bg_params = list(fill = "grey90", col = NA))
  )
  
  grid.table(summary_table, rows = NULL, theme = table_theme)
  
  # Overall correlation plot
  if (!is.null(visualization_results$overall_plot)) {
    # Create a new version of the plot with better formatting
    overall_data <- correlation_results$overall$matched_data %>%
      filter(Intensity_2CV > 0, Intensity_3CV > 0)
    
    # Create improved plot with better margins
    overall_plot <- ggplot(overall_data, aes(x = Intensity_2CV, y = Intensity_3CV)) +
      geom_point(alpha = 0.5, size = 0.7) +  # Smaller points
      geom_smooth(method = "lm", color = "red") +
      scale_x_log10(labels = scales::scientific) +
      scale_y_log10(labels = scales::scientific) +
      labs(
        title = "Overall Peptide Intensity Correlation (All Samples)",
        subtitle = sprintf(
          "Pearson: %.3f, Spearman: %.3f, R²: %.3f (n=%d, filtered out %d zeros/negatives)",
          cor(overall_data$Intensity_2CV, overall_data$Intensity_3CV, method = "pearson"),
          cor(overall_data$Intensity_2CV, overall_data$Intensity_3CV, method = "spearman"),
          summary(lm(Intensity_3CV ~ Intensity_2CV, data = overall_data))$r.squared,
          nrow(overall_data),
          nrow(correlation_results$overall$matched_data) - nrow(overall_data)
        ),
        x = "2CV Intensity (log10)",
        y = "3CV Intensity (log10)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.margin = margin(30, 30, 30, 30, "pt"),  # Larger margins
        panel.grid.minor = element_blank()  # Remove minor grid lines
      )
    
    print(overall_plot)
  }
  
  # Sample plots - 4 per page with smaller points and better margins
  sample_ids <- names(correlation_results$by_sample)
  num_samples <- length(sample_ids)
  
  # Process in batches of 4
  for (i in seq(1, num_samples, by = 4)) {
    grid.newpage()
    
    # Set up the layout
    pushViewport(viewport(layout = grid.layout(3, 2, heights = c(0.1, 1, 1), widths = c(1, 1))))
    
    # Add the page title
    grid.text("Sample Correlation Plots", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
              gp = gpar(fontsize = 16, fontface = "bold"))
    
    # Get indices for this page
    end_idx <- min(i + 3, num_samples)
    batch_ids <- sample_ids[i:end_idx]
    
    # Calculate positions for each plot
    positions <- list(
      c(2, 1), c(2, 2),
      c(3, 1), c(3, 2)
    )
    
    # Create and place each plot
    for (j in 1:length(batch_ids)) {
      sample_id <- batch_ids[j]
      pos <- positions[[j]]
      
      # Get data for this sample
      sample_data <- correlation_results$by_sample[[sample_id]]$matched_data %>%
        filter(Intensity_2CV > 0, Intensity_3CV > 0)
      
      # Create improved plot
      p <- ggplot(sample_data, aes(x = Intensity_2CV, y = Intensity_3CV)) +
        geom_point(alpha = 0.5, size = 0.7) +
        geom_smooth(method = "lm", color = "red") +
        scale_x_log10(labels = scales::scientific) +
        scale_y_log10(labels = scales::scientific) +
        labs(
          title = paste("Sample", sample_id),
          subtitle = sprintf(
            "Pearson: %.3f, Spearman: %.3f, R²: %.3f (n=%d, filtered out %d zeros)",
            cor(sample_data$Intensity_2CV, sample_data$Intensity_3CV, method = "pearson"),
            cor(sample_data$Intensity_2CV, sample_data$Intensity_3CV, method = "spearman"),
            summary(lm(Intensity_3CV ~ Intensity_2CV, data = sample_data))$r.squared,
            nrow(sample_data),
            nrow(correlation_results$by_sample[[sample_id]]$matched_data) - nrow(sample_data)
          ),
          x = "2CV Intensity (log10)",
          y = "3CV Intensity (log10)"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 8, hjust = 0.5),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          plot.margin = margin(15, 15, 15, 15, "pt"),
          panel.grid.minor = element_blank()
        )
      
      # Print the plot in the right position
      print(p, vp = viewport(layout.pos.row = pos[1], layout.pos.col = pos[2]))
    }
    
    popViewport()
  }
  
  # Simplified heatmap that's guaranteed to work
  if (!is.null(correlation_results$overall$matched_data)) {
    # Get the data
    heat_data <- correlation_results$overall$matched_data %>%
      filter(Intensity_2CV > 0, Intensity_3CV > 0) %>%
      mutate(
        log10_2CV = log10(Intensity_2CV),
        log10_3CV = log10(Intensity_3CV)
      )
    
    # Create a simplified 2D histogram with hexbins
    heat_plot <- ggplot(heat_data, aes(x = log10_2CV, y = log10_3CV)) +
      stat_bin2d(bins = 30) +
      scale_fill_viridis_c(option = "viridis", trans = "log1p") +
      labs(
        title = "Peptide Intensity Correlation Heatmap",
        x = "2CV Intensity (log10)",
        y = "3CV Intensity (log10)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = margin(20, 20, 20, 20, "pt")
      )
    
    print(heat_plot)
  }
  
  # Interpretation page
  grid.newpage()
  grid.text("Interpretation of Results", x = 0.5, y = 0.95, 
            gp = gpar(fontsize = 18, fontface = "bold"))
  
  # Extract values from the original summary table to ensure we use the correct values
  pearson_overall <- visualization_results$summary_table$Pearson[nrow(visualization_results$summary_table)]
  r_squared_overall <- visualization_results$summary_table$R_Squared[nrow(visualization_results$summary_table)]
  pearson_min <- min(visualization_results$summary_table$Pearson[-nrow(visualization_results$summary_table)])
  pearson_max <- max(visualization_results$summary_table$Pearson[-nrow(visualization_results$summary_table)])
  
  # Format them
  pearson_overall_str <- sprintf("%.3f", pearson_overall)
  r_squared_overall_str <- sprintf("%.3f", r_squared_overall)
  pearson_min_str <- sprintf("%.3f", pearson_min)
  pearson_max_str <- sprintf("%.3f", pearson_max)
  
  # Create text
  interpretation_text <- paste(
    "Summary of findings:",
    "",
    "1. Overall correlation: Strong correlation between 2CV and 3CV samples",
    paste("   Pearson correlation:", pearson_overall_str),
    paste("   R-squared:", r_squared_overall_str),
    "",
    "2. Sample variability: Sample correlations range from",
    paste("   ", pearson_min_str, "to", pearson_max_str, "(Pearson)"),
    "",
    "3. Data filtering: Zero intensity values were excluded from log-scale visualizations.",
    "   These represent peptides detected in both samples but with intensity below",
    "   the quantification threshold in at least one of the samples.",
    "",
    "4. Heatmap interpretation: The heatmap shows the density of peptide intensities,",
    "   with the diagonal pattern confirming the strong correlation between 2CV and 3CV measurements.",
    sep = "\n"
  )
  
  # Add text to the page with better formatting
  grid.text(interpretation_text, x = 0.05, y = 0.8, just = c("left", "top"),
            gp = gpar(fontsize = 12))
  
  # Close the PDF device
  dev.off()
  cat("Polished summary PDF created at:", pdf_path, "\n")
}

create_summary_pdf_polished(results$correlation, results$visualization)

# Print summary of results
cat("\n=== CORRELATION ANALYSIS SUMMARY ===\n")
print(results$visualization$summary_table)
cat("\nVisualization files have been saved in the 'plots' directory\n")