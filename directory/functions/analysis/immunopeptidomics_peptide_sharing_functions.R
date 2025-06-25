# ===============================
# PEPTIDE SHARING ANALYSIS FUNCTIONS
# ===============================
# Functions for comprehensive peptide sharing and distribution analysis
# Author: Generated for immunopeptidomics analysis

# Function to prepare peptide sharing data
prepare_peptide_sharing_data <- function(combined_data, min_length = 8, max_length = 12, exclude_samples = NULL, include_samples = NULL) {
  
  cat("\n=== PEPTIDE SHARING DATA PREPARATION ===\n")
  
  # Apply sample filtering first
  data <- combined_data
  original_count <- nrow(data)
  
  # Apply inclusion filter
  if (!is.null(include_samples)) {
    missing_samples <- setdiff(include_samples, unique(data$SampleID))
    if (length(missing_samples) > 0) {
      warning("WARNING: Requested samples not found in data: ", paste(missing_samples, collapse = ", "))
    }
    data <- data %>% filter(SampleID %in% include_samples)
    cat("After including samples", paste(include_samples, collapse = ", "), ":", nrow(data), "entries\n")
  }
  
  # Apply exclusion filter
  if (!is.null(exclude_samples)) {
    found_samples <- intersect(exclude_samples, unique(data$SampleID))
    if (length(found_samples) > 0) {
      data <- data %>% filter(!SampleID %in% exclude_samples)
      cat("After excluding samples", paste(found_samples, collapse = ", "), ":", nrow(data), "entries\n")
    }
  }
  
  # Check for required columns and map if needed
  required_mapping <- list(
    "Peptide" = c("Peptide", "Sequence", "peptide"),
    "SampleID" = c("SampleID", "Sample", "sample_id"),
    "Peptide Length" = c("Peptide Length", "Length", "peptide_length", "PeptideLength"),
    "Gene" = c("Gene", "gene", "Gene Name", "GeneName"),
    "Protein" = c("Protein", "protein", "Protein Name", "ProteinName")
  )
  
  # Map column names
  for (target_col in names(required_mapping)) {
    possible_cols <- required_mapping[[target_col]]
    found_col <- intersect(possible_cols, colnames(data))
    
    if (length(found_col) > 0) {
      if (found_col[1] != target_col) {
        data[[target_col]] <- data[[found_col[1]]]
        cat("Mapped column:", found_col[1], "->", target_col, "\n")
      }
    } else if (target_col == "Peptide Length") {
      # Calculate peptide length if not available
      data$`Peptide Length` <- nchar(data$Peptide)
      cat("Calculated peptide length from sequences\n")
    } else if (target_col %in% c("Gene", "Protein")) {
      # These are optional, set to NA if missing
      data[[target_col]] <- NA
      cat("Warning: Column", target_col, "not found - setting to NA\n")
    }
  }
  
  # Filter by peptide length
  data <- data %>% 
    filter(`Peptide Length` >= min_length & `Peptide Length` <= max_length)
  
  cat("After length filtering (", min_length, "-", max_length, "mers):", nrow(data), "entries\n")
  cat("Unique peptides:", n_distinct(data$Peptide), "\n")
  cat("Unique samples:", n_distinct(data$SampleID), "\n")
  
  return(data)
}

# Function to create peptide report with sample distribution
create_peptide_report <- function(peptide_data) {
  
  cat("Creating peptide sample distribution report...\n")
  
  peptide_report <- peptide_data %>%
    select(Peptide, SampleID, `Peptide Length`, Gene, Protein) %>%
    distinct() %>%
    group_by(Peptide, `Peptide Length`) %>%
    summarize(
      sample_list = paste(sort(SampleID), collapse = ", "),
      sample_count = n_distinct(SampleID),
      genes = paste(unique(na.omit(Gene)), collapse = "; "),
      proteins = paste(unique(na.omit(Protein)), collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(desc(sample_count), Peptide)
  
  # Clean up empty gene/protein fields
  peptide_report$genes[peptide_report$genes == ""] <- NA
  peptide_report$proteins[peptide_report$proteins == ""] <- NA
  
  return(peptide_report)
}

# Function to perform peptide sharing analysis
perform_peptide_sharing_analysis <- function(peptide_data, output_dir, dataset_name, timestamp) {
  
  cat("\n=== PEPTIDE SHARING ANALYSIS ===\n")
  
  # Create analysis subdirectory
  sharing_dir <- file.path(output_dir, "peptide_sharing")
  plots_dir <- file.path(sharing_dir, "plots")
  tables_dir <- file.path(sharing_dir, "tables")
  
  dir.create(sharing_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create peptide report
  peptide_report <- create_peptide_report(peptide_data)
  
  # Get unique samples
  unique_samples <- sort(unique(peptide_data$SampleID))
  num_samples <- length(unique_samples)
  
  cat("Analysis includes", nrow(peptide_report), "unique peptides across", num_samples, "samples\n")
  
  # 1. Sample count distribution
  sample_count_dist <- peptide_report %>%
    count(sample_count) %>%
    mutate(percentage = n / sum(n) * 100)
  
  # 2. Unique vs shared peptides
  unique_shared <- data.frame(
    category = c("Unique", "Shared"),
    count = c(
      sum(sample_count_dist$n[sample_count_dist$sample_count == 1]),
      sum(sample_count_dist$n[sample_count_dist$sample_count > 1])
    )
  ) %>%
    mutate(percentage = count / sum(count) * 100)
  
  # 3. Length distribution
  length_dist <- peptide_report %>%
    count(`Peptide Length`) %>%
    mutate(percentage = n / sum(n) * 100)
  
  # 4. Average sharing by length
  avg_sharing_by_length <- peptide_report %>%
    group_by(`Peptide Length`) %>%
    summarise(
      avg_samples = mean(sample_count),
      median_samples = median(sample_count),
      total_peptides = n(),
      unique_count = sum(sample_count == 1),
      shared_count = sum(sample_count > 1),
      percent_shared = round(shared_count / total_peptides * 100, 1),
      .groups = "drop"
    )
  
  # 5. Robust sharing analysis
  sharing_thresholds <- c(0.25, 0.5, 0.75, 0.9)
  samples_needed <- ceiling(num_samples * sharing_thresholds)
  threshold_labels <- paste0(sharing_thresholds * 100, "% of samples (", samples_needed, "/", num_samples, ")")
  names(samples_needed) <- threshold_labels
  
  robustly_shared_peptides <- list()
  for (i in 1:length(samples_needed)) {
    threshold <- names(samples_needed)[i]
    min_samples <- samples_needed[i]
    
    robust_peptides <- peptide_report %>%
      filter(sample_count >= min_samples) %>%
      arrange(desc(sample_count), Peptide)
    
    robustly_shared_peptides[[threshold]] <- robust_peptides
    cat("Peptides present in at least", threshold, ":", nrow(robust_peptides), "\n")
  }
  
  robust_summary <- data.frame(
    threshold = names(robustly_shared_peptides),
    peptide_count = sapply(robustly_shared_peptides, nrow),
    min_samples = samples_needed
  )
  
  # 6. Gene analysis (if gene data available)
  gene_analysis <- NULL
  if (!all(is.na(peptide_report$genes))) {
    # Get top genes by peptide count
    top_genes <- peptide_report %>%
      filter(!is.na(genes) & genes != "") %>%
      mutate(gene_list = strsplit(genes, "; ")) %>%
      unnest(gene_list) %>%
      count(gene_list, sort = TRUE) %>%
      head(15) %>%
      pull(gene_list)
    
    if (length(top_genes) > 0) {
      gene_analysis <- list(top_genes = top_genes)
    }
  }
  
  # 7. Sample-sample Jaccard similarity
  jaccard_results <- calculate_jaccard_similarity(peptide_data, unique_samples)
  
  # Save all tables
  cat("Saving peptide sharing tables...\n")
  write_csv(peptide_report, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_01_peptide_sample_distribution.csv")))
  write_csv(sample_count_dist, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_02_sample_count_distribution.csv")))
  write_csv(unique_shared, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_03_unique_vs_shared_summary.csv")))
  write_csv(length_dist, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_04_length_distribution.csv")))
  write_csv(avg_sharing_by_length, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_05_sharing_by_length.csv")))
  write_csv(robust_summary, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_06_robust_sharing_summary.csv")))
  
  # Save robust peptide sets
  for (i in 1:length(robustly_shared_peptides)) {
    threshold <- names(robustly_shared_peptides)[i]
    robust_set <- robustly_shared_peptides[[threshold]]
    
    filename <- paste0(timestamp, "_", dataset_name, "_07_robust_peptides_", sprintf("%02d", i), ".csv")
    write_csv(robust_set, file.path(tables_dir, filename))
  }
  
  # Save Jaccard similarity results
  write_csv(jaccard_results$similarity_df, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_08_sample_jaccard_similarity.csv")))
  write_csv(jaccard_results$similarity_matrix, file.path(tables_dir, paste0(timestamp, "_", dataset_name, "_09_jaccard_similarity_matrix.csv")))
  
  # Return results for plotting
  results <- list(
    peptide_report = peptide_report,
    sample_count_dist = sample_count_dist,
    unique_shared = unique_shared,
    length_dist = length_dist,
    avg_sharing_by_length = avg_sharing_by_length,
    robustly_shared_peptides = robustly_shared_peptides,
    robust_summary = robust_summary,
    gene_analysis = gene_analysis,
    jaccard_results = jaccard_results,
    unique_samples = unique_samples,
    num_samples = num_samples,
    plots_dir = plots_dir,
    tables_dir = tables_dir
  )
  
  return(results)
}

# Function to calculate Jaccard similarity between samples
calculate_jaccard_similarity <- function(peptide_data, unique_samples) {
  
  cat("Calculating Jaccard similarity between samples...\n")
  
  # Function to count shared peptides between two samples
  count_shared_peptides <- function(sample1, sample2) {
    peptides_sample1 <- peptide_data %>%
      filter(SampleID == sample1) %>%
      pull(Peptide) %>%
      unique()
    
    peptides_sample2 <- peptide_data %>%
      filter(SampleID == sample2) %>%
      pull(Peptide) %>%
      unique()
    
    shared_count <- length(intersect(peptides_sample1, peptides_sample2))
    return(shared_count)
  }
  
  # Create sample-sample sharing matrix
  sample_sharing_matrix <- matrix(0, nrow = length(unique_samples), ncol = length(unique_samples))
  rownames(sample_sharing_matrix) <- unique_samples
  colnames(sample_sharing_matrix) <- unique_samples
  
  # Fill the matrix
  for (i in 1:length(unique_samples)) {
    for (j in 1:length(unique_samples)) {
      sample1 <- unique_samples[i]
      sample2 <- unique_samples[j]
      
      if (i == j) {
        # Diagonal: total peptides for that sample
        total_peptides <- peptide_data %>%
          filter(SampleID == sample1) %>%
          pull(Peptide) %>%
          unique() %>%
          length()
        sample_sharing_matrix[i, j] <- total_peptides
      } else {
        # Off diagonal: shared peptides
        sample_sharing_matrix[i, j] <- count_shared_peptides(sample1, sample2)
      }
    }
  }
  
  # Convert to data frame for analysis
  similarity_df <- as.data.frame(sample_sharing_matrix)
  similarity_df$Sample1 <- rownames(sample_sharing_matrix)
  similarity_df <- similarity_df %>%
    pivot_longer(cols = unique_samples,
                 names_to = "Sample2",
                 values_to = "shared_peptides") %>%
    rowwise() %>%
    mutate(
      sample1_total = sample_sharing_matrix[Sample1, Sample1],
      sample2_total = sample_sharing_matrix[Sample2, Sample2],
      jaccard_index = shared_peptides / (sample1_total + sample2_total - shared_peptides)
    )
  
  # Create matrix format for easier viewing
  similarity_matrix <- similarity_df %>%
    select(Sample1, Sample2, jaccard_index) %>%
    pivot_wider(names_from = Sample2, values_from = jaccard_index)
  
  return(list(
    similarity_df = similarity_df,
    similarity_matrix = similarity_matrix,
    sharing_matrix = sample_sharing_matrix
  ))
}

# Fixed function to create UpSet plot data and plot
create_upset_plot_fixed <- function(peptide_data, unique_samples, plots_dir, dataset_name, timestamp) {
  
  cat("Creating UpSet plot for sample intersections...\n")
  
  tryCatch({
    # Create a proper binary matrix for UpSetR
    # Each row is a peptide, each column is a sample (1 = present, 0 = absent)
    
    # Get all unique peptides
    all_peptides <- unique(peptide_data$Peptide)
    
    # Create the binary matrix
    upset_matrix <- data.frame(row.names = all_peptides)
    
    # Add a column for each sample
    for (sample in unique_samples) {
      # Get peptides for this sample
      sample_peptides <- peptide_data %>%
        filter(SampleID == sample) %>%
        pull(Peptide) %>%
        unique()
      
      # Create binary vector (1 if peptide present in sample, 0 if not)
      upset_matrix[[sample]] <- as.integer(all_peptides %in% sample_peptides)
    }
    
    cat("UpSet matrix dimensions:", nrow(upset_matrix), "x", ncol(upset_matrix), "\n")
    cat("Sample columns:", paste(colnames(upset_matrix), collapse = ", "), "\n")
    cat("Total peptides with at least one sample:", sum(rowSums(upset_matrix) > 0), "\n")
    
    # Check if we have valid data
    if (nrow(upset_matrix) == 0 || ncol(upset_matrix) == 0) {
      cat("⚠ No data available for UpSet plot\n")
      return(FALSE)
    }
    
    # Remove peptides that don't appear in any sample (shouldn't happen but just in case)
    upset_matrix <- upset_matrix[rowSums(upset_matrix) > 0, ]
    
    if (nrow(upset_matrix) == 0) {
      cat("⚠ No peptides found in any sample for UpSet plot\n")
      return(FALSE)
    }
    
    cat("Final UpSet matrix:", nrow(upset_matrix), "peptides x", ncol(upset_matrix), "samples\n")
    
    # Create the UpSet plot
    png_file <- file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_07_upset_sample_intersections.png"))
    
    png(png_file, width = 14, height = 10, units = "in", res = 300)
    
    # Create UpSet plot with better parameters
    UpSetR::upset(
      upset_matrix,
      sets = unique_samples,
      sets.bar.color = "steelblue",
      main.bar.color = "darkgreen",
      matrix.color = "red",
      order.by = "freq",
      decreasing = TRUE,
      keep.order = FALSE,
      number.angles = 30,
      point.size = 3,
      line.size = 1.2,
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size",
      text.scale = c(1.3, 1.3, 1.2, 1.2, 1.5, 1.2),
      set_size.angles = 45,
      set_size.show = TRUE,
      nintersects = min(20, 2^length(unique_samples) - 1)  # Limit number of intersections shown
    )
    
    dev.off()
    
    cat("✓ UpSet plot saved successfully to:", png_file, "\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("✗ Error creating UpSet plot:", e$message, "\n")
    cat("This might be due to data format issues or package conflicts\n")
    
    # Create an alternative visualization as fallback
    cat("Creating alternative intersection plot...\n")
    
    tryCatch({
      # Create a simple intersection size bar plot as alternative
      intersection_data <- peptide_data %>%
        select(Peptide, SampleID) %>%
        distinct() %>%
        group_by(Peptide) %>%
        summarise(
          sample_count = n_distinct(SampleID),
          samples = paste(sort(SampleID), collapse = " & "),
          .groups = "drop"
        ) %>%
        count(sample_count, samples) %>%
        arrange(desc(sample_count), desc(n))
      
      # Limit to top 15 intersections for readability
      if (nrow(intersection_data) > 15) {
        intersection_data <- intersection_data %>%
          slice_head(n = 15) %>%
          mutate(samples = ifelse(nchar(samples) > 30, 
                                  paste0(substr(samples, 1, 27), "..."), 
                                  samples))
      }
      
      alt_plot <- ggplot(intersection_data, aes(x = reorder(samples, n), y = n)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
        geom_text(aes(label = n), hjust = -0.1, size = 3) +
        coord_flip() +
        labs(
          title = "Sample Intersections (Alternative View)",
          subtitle = paste("Top intersections -", nrow(intersection_data), "shown"),
          x = "Sample Combination",
          y = "Number of Shared Peptides"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 10)
        )
      
      alt_file <- file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_07_alternative_sample_intersections.png"))
      ggsave(alt_file, alt_plot, width = 12, height = 8, dpi = 300)
      
      cat("✓ Alternative intersection plot saved to:", alt_file, "\n")
      return(TRUE)
      
    }, error = function(e2) {
      cat("✗ Failed to create alternative plot:", e2$message, "\n")
      return(FALSE)
    })
  })
}

# Function to create static peptide sharing plots
create_peptide_sharing_plots <- function(sharing_results, dataset_name, timestamp) {
  
  cat("Creating peptide sharing visualizations...\n")
  
  plots_dir <- sharing_results$plots_dir
  plots <- list()
  
  # 1. Sample count distribution
  plots$sample_count_dist <- ggplot(sharing_results$sample_count_dist, 
                                    aes(x = factor(sample_count), y = n)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")), 
              vjust = -0.3, size = 3) +
    labs(title = "Peptide Distribution by Sample Count",
         subtitle = paste("Total peptides:", sum(sharing_results$sample_count_dist$n)),
         x = "Number of Samples",
         y = "Number of Peptides") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # 2. Unique vs shared pie chart
  plots$unique_shared_pie <- ggplot(sharing_results$unique_shared, 
                                    aes(x = "", y = count, fill = category)) +
    geom_col(width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("Unique" = "#4CAF50", "Shared" = "#FF9800")) +
    labs(title = "Unique vs Shared Peptides",
         subtitle = paste("Total peptides:", sum(sharing_results$unique_shared$count)),
         fill = "Category") +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom") +
    geom_text(aes(label = paste0(category, "\n", count, "\n(", round(percentage, 1), "%)")), 
              position = position_stack(vjust = 0.5), size = 4)
  
  # 3. Peptide length distribution
  plots$length_dist <- ggplot(sharing_results$length_dist, 
                              aes(x = factor(`Peptide Length`), y = n)) +
    geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
    geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")), 
              vjust = -0.3, size = 3) +
    labs(title = "Peptide Length Distribution",
         subtitle = paste("Peptides ranging from", min(sharing_results$length_dist$`Peptide Length`), 
                          "to", max(sharing_results$length_dist$`Peptide Length`), "amino acids"),
         x = "Peptide Length",
         y = "Number of Peptides") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # 4. Average sharing by length
  plots$avg_sharing_length <- ggplot(sharing_results$avg_sharing_by_length, 
                                     aes(x = factor(`Peptide Length`), y = avg_samples)) +
    geom_bar(stat = "identity", fill = "orange", alpha = 0.7) +
    geom_text(aes(label = round(avg_samples, 2)), vjust = -0.3, size = 3) +
    labs(title = "Average Sample Sharing by Peptide Length",
         subtitle = "Higher values indicate more sharing across samples",
         x = "Peptide Length",
         y = "Average Number of Samples") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # 5. Robust sharing summary
  plots$robust_sharing <- ggplot(sharing_results$robust_summary, 
                                 aes(x = reorder(threshold, -peptide_count), y = peptide_count)) +
    geom_bar(stat = "identity", fill = "darkblue", alpha = 0.7) +
    geom_text(aes(label = peptide_count), vjust = -0.3, size = 3) +
    labs(title = "Robustly Shared Peptides by Threshold",
         subtitle = "Number of peptides found in different proportions of samples",
         x = "Sharing Threshold",
         y = "Number of Peptides") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 6. Jaccard similarity heatmap
  jaccard_matrix <- sharing_results$jaccard_results$jaccard_index
  jaccard_long <- sharing_results$jaccard_results$similarity_df %>%
    select(Sample1, Sample2, jaccard_index)
  
  plots$jaccard_heatmap <- ggplot(jaccard_long, 
                                  aes(x = Sample2, y = Sample1, fill = jaccard_index)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, name = "Jaccard\nIndex") +
    geom_text(aes(label = round(jaccard_index, 3)), size = 3) +
    labs(title = "Sample Similarity (Jaccard Index)",
         subtitle = "1.0 = identical peptide sets, 0.0 = no shared peptides",
         x = "Sample", y = "Sample") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 7. UpSet plot for sample intersections (FIXED VERSION)
  if (sharing_results$num_samples >= 2) {
    # Use the fixed UpSet plot function
    upset_success <- create_upset_plot_fixed(
      peptide_data = sharing_results$peptide_report %>%
        select(Peptide, sample_list) %>%
        separate_rows(sample_list, sep = ", ") %>%
        rename(SampleID = sample_list),
      unique_samples = sharing_results$unique_samples,
      plots_dir = plots_dir,
      dataset_name = dataset_name,
      timestamp = timestamp
    )
    
    if (!upset_success) {
      cat("⚠ UpSet plot creation failed, but alternative may be available\n")
    }
  } else {
    cat("⚠ Need at least 2 samples for intersection analysis\n")
  }
  
  # Save individual plots
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_01_peptide_count_distribution.png")), 
         plots$sample_count_dist, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_02_unique_vs_shared_pie.png")), 
         plots$unique_shared_pie, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_03_length_distribution.png")), 
         plots$length_dist, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_04_avg_sharing_by_length.png")), 
         plots$avg_sharing_length, width = 10, height = 8, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_05_robust_sharing_summary.png")), 
         plots$robust_sharing, width = 12, height = 8, dpi = 300)
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_06_jaccard_similarity_heatmap.png")), 
         plots$jaccard_heatmap, width = 10, height = 8, dpi = 300)
  
  # Create combined plots
  combined_plot1 <- grid.arrange(
    plots$sample_count_dist, plots$unique_shared_pie,
    plots$length_dist, plots$avg_sharing_length,
    ncol = 2, nrow = 2
  )
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_08_combined_basic_analysis.png")), 
         combined_plot1, width = 16, height = 12, dpi = 300)
  
  combined_plot2 <- grid.arrange(
    plots$robust_sharing, plots$jaccard_heatmap,
    ncol = 2, nrow = 1
  )
  
  ggsave(file.path(plots_dir, paste0(timestamp, "_", dataset_name, "_09_combined_advanced_analysis.png")), 
         combined_plot2, width = 16, height = 8, dpi = 300)
  
  cat("✓ Peptide sharing plots saved!\n")
  
  return(plots)
}