#!/usr/bin/env Rscript
#' HLA-I Analysis - Fusion Protein Analysis
#' analyze_fusion.R
#' Specialized script for analyzing fusion protein-derived peptides
#' @author Your Name
#' @version 1.0

#' Main function for fusion protein analysis
#' @param config Configuration object
#' @param dirs Directory structure created by create_output_directories
#' @return Invisibly returns a list of results
analyze_fusion <- function(config, dirs = NULL) {
  # Create output directories if not provided
  if (is.null(dirs)) {
    dirs <- create_output_directories(config$data_path, config$output_name)
  }
  
  # Validate required fusion parameters
  if (is.null(config$fusion_sequence) || is.null(config$fusion_parts) || is.null(config$junction_position)) {
    stop("Fusion analysis requires fusion_sequence, fusion_parts, and junction_position in config")
  }
  
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
    experiment_type = "multi_sample",  # Treat as multi-sample for flexibility
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
  
  # Perform multi-sample analysis first
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
  
  # 3. Identify fusion peptides (core functionality of this script)
  cat("\n## 3. Identifying fusion peptides...\n")
  
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
  
  # Calculate fusion peptide metrics
  if (nrow(fusion_peptides) > 0) {
    # Add sequence position information
    fusion_peptides <- fusion_peptides %>%
      mutate(
        fusion_peptide_length = nchar(Peptide),
        
        # Try to find peptide in fusion sequence
        fusion_pos = if_else(
          from_fusion,
          as.integer(stringr::str_locate(config$fusion_sequence, Peptide)[1]),
          NA_integer_
        ),
        
        # Determine relative position to junction
        junction_relative_position = if_else(
          !is.na(fusion_pos),
          fusion_pos - config$junction_position,
          NA_integer_
        ),
        
        # Classify as upstream, downstream, or spanning
        position_class = case_when(
          spans_junction ~ "Junction-spanning",
          !is.na(junction_relative_position) & junction_relative_position < 0 & 
            junction_relative_position + fusion_peptide_length <= 0 ~ "Upstream (Part 1)",
          !is.na(junction_relative_position) & junction_relative_position >= 0 ~ "Downstream (Part 2)",
          TRUE ~ "Other"
        )
      )
  }
  
  # 4. Perform detailed analysis of junction-spanning peptides
  cat("\n## 4. Analyzing junction-spanning peptides...\n")
  
  junction_peptides <- fusion_peptides %>%
    filter(spans_junction)
  
  if (nrow(junction_peptides) > 0) {
    # Calculate positions relative to junction
    junction_peptides <- junction_peptides %>%
      mutate(
        amino_acids_before_junction = config$junction_position - fusion_pos,
        amino_acids_after_junction = fusion_peptide_length - amino_acids_before_junction,
        amino_acids_ratio = amino_acids_before_junction / amino_acids_after_junction,
        
        # Extract sequences from each part
        part1_seq = substring(Peptide, 1, amino_acids_before_junction),
        part2_seq = substring(Peptide, amino_acids_before_junction + 1, fusion_peptide_length)
      )
    
    cat("Detailed analysis of", nrow(junction_peptides), "junction-spanning peptides:\n")
    
    # Create a summary of junction-spanning peptides
    junction_summary <- junction_peptides %>%
      select(
        Peptide, fusion_peptide_length, 
        amino_acids_before_junction, amino_acids_after_junction, 
        part1_seq, part2_seq
      ) %>%
      arrange(desc(amino_acids_before_junction))
    
    print(junction_summary)
  } else {
    cat("No junction-spanning peptides found\n")
  }
  
  # 5. Generate visualizations
  cat("\n## 5. Generating visualizations...\n")
  
  # Prepare visualization directory
  viz_dir <- dirs$viz_dir
  
  # a) Heatmap of fusion peptides across samples
  if (nrow(fusion_peptides) > 0 && length(sample_ids) > 0) {
    # Get intensity columns for all samples
    intensity_cols <- paste0("total_intensity_", sample_ids)
    
    # Create heatmap of fusion peptides
    fusion_heatmap <- create_peptide_heatmap(
      fusion_peptides,
      value_cols = intensity_cols,
      is_intensity = TRUE,
      log_transform = TRUE,
      peptide_col = "Peptide",
      annotation_cols = c("fusion_peptide_type", "spans_junction", "position_class"),
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
  }
  
  # b) Bar plot of fusion peptide types
  if (nrow(fusion_peptides) > 0 && "fusion_peptide_type" %in% colnames(fusion_peptides)) {
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
  
  # c) Position distribution plot relative to the junction
  if (nrow(fusion_peptides) > 0 && "junction_relative_position" %in% colnames(fusion_peptides)) {
    # Filter out NA values
    position_data <- fusion_peptides %>%
      filter(!is.na(junction_relative_position)) %>%
      mutate(
        # Add an artificial stacking variable for peptides at same position
        stack_pos = dense_rank(interaction(Peptide, junction_relative_position))
      )
    
    if (nrow(position_data) > 0) {
      position_plot <- ggplot(position_data, 
                              aes(x = junction_relative_position, 
                                  y = stack_pos, 
                                  color = position_class, 
                                  label = Peptide)) +
        geom_point(size = 3) +
        geom_segment(aes(x = junction_relative_position, 
                         xend = junction_relative_position + fusion_peptide_length - 1,
                         y = stack_pos, 
                         yend = stack_pos),
                     size = 2, alpha = 0.6) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
        annotate("text", x = 0, y = 0, label = "Junction", color = "red", angle = 90, vjust = -0.5) +
        scale_color_brewer(palette = "Set1") +
        labs(
          title = "Fusion Peptide Positions Relative to Junction",
          subtitle = "Each line represents a peptide, vertical red line is the junction",
          x = "Position Relative to Junction",
          y = "Peptide Index",
          color = "Position Class"
        ) +
        theme_minimal() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      
      save_plot(position_plot, file.path(viz_dir, "fusion_peptide_positions"), 
                width = 12, height = max(8, nrow(position_data)/3))
      
      # Add text labels for junction-spanning peptides in a separate plot
      junction_peptides_data <- position_data %>%
        filter(spans_junction)
      
      if (nrow(junction_peptides_data) > 0) {
        junction_peptide_plot <- position_plot +
          geom_text_repel(
            data = junction_peptides_data,
            aes(label = paste0(Peptide, " (", part1_seq, "|", part2_seq, ")")),
            size = 3, max.overlaps = 30,
            nudge_y = 0.5
          ) +
          labs(title = "Junction-Spanning Peptides with Sequences",
               subtitle = "Vertical bar | marks the junction position within each peptide")
        
        save_plot(junction_peptide_plot, file.path(viz_dir, "junction_peptide_positions"), 
                  width = 12, height = max(8, nrow(junction_peptides_data)*1.5))
      }
    }
  }
  
  # d) Sample distribution of fusion peptides (presence heatmap)
  if (nrow(fusion_peptides) > 0 && length(sample_ids) > 1) {
    # Get presence/absence matrix for samples
    presence_cols <- paste0("total_intensity_", sample_ids)
    
    # Create binary presence/absence heatmap
    presence_heatmap <- create_peptide_heatmap(
      fusion_peptides,
      value_cols = presence_cols,
      is_intensity = FALSE,  # Binary presence/absence
      peptide_col = "Peptide",
      annotation_cols = c("fusion_peptide_type", "spans_junction", "position_class"),
      title = "Fusion Peptide Presence Across Samples",
      cluster_rows = FALSE,
      cluster_cols = TRUE
    )
    
    pdf(file.path(viz_dir, "fusion_peptide_presence.pdf"), width = 10, height = max(8, nrow(fusion_peptides)/3))
    print(presence_heatmap)
    dev.off()
    
    png(file.path(viz_dir, "fusion_peptide_presence.png"), width = 800, height = max(600, nrow(fusion_peptides)*40), res = 100)
    print(presence_heatmap)
    dev.off()
  }
  
  # e) Spiked peptide visualizations (if applicable)
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
        title = paste0("Detection of Spiked Fusion Junction Peptides Across ", length(sample_ids), " Samples"),
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
          title = "Intensity of Detected Spiked Fusion Peptides",
          x = "Peptide-Sample",
          y = "Total Intensity",
          fill = "Type"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      
      save_plot(intensity_plot, file.path(viz_dir, "spiked_peptide_intensity"))
    }
  }
  
  # 6. Generate interactive visualizations if requested
  if (config$generate_interactive) {
    cat("\n## 6. Creating interactive visualizations...\n")
    
    # Convert appropriate plots to interactive versions
    interactive_viz_list <- list()
    
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
    
    # Position distribution plot
    if (exists("position_plot")) {
      position_interactive <- ggplotly(position_plot)
      htmlwidgets::saveWidget(position_interactive, 
                              file.path(viz_dir, "interactive_fusion_positions.html"), 
                              selfcontained = TRUE)
      
      interactive_viz_list[[length(interactive_viz_list) + 1]] <- list(
        path = file.path(viz_dir, "interactive_fusion_positions.html"),
        title = "Fusion Peptide Positions"
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
        "Analysis of fusion peptides across ", length(sample_ids), " samples. ",
        "A total of ", nrow(fusion_peptides), " fusion-derived peptides were identified, ",
        if("spans_junction" %in% colnames(fusion_peptides)) 
          paste0("of which ", sum(fusion_peptides$spans_junction), " span the fusion junction. ") else "",
        if(exists("spiked_matrix_data")) 
          paste0("Analysis included ", nrow(spiked_matrix_data), " spiked fusion peptides.") else ""
      )
      
      # Create the dashboard
      create_interactive_dashboard(
        interactive_viz_list,
        title = "Fusion Peptide Analysis Dashboard",
        output_file = file.path(viz_dir, "interactive_dashboard.html"),
        summary_text = summary_text
      )
    }
  }
  
  # 7. Generate Excel reports
  cat("\n## 7. Generating Excel reports...\n")
  
  # Create Excel sheets
  excel_sheets <- list()
  
  # Add fusion peptides sheet
  excel_sheets[["All_Fusion_Peptides"]] <- fusion_peptides
  
  # Add junction-spanning peptides sheet if available
  if (exists("junction_peptides") && nrow(junction_peptides) > 0) {
    excel_sheets[["Junction_Spanning_Peptides"]] <- junction_peptides
  }
  
  # Add fusion position details
  if ("position_class" %in% colnames(fusion_peptides)) {
    excel_sheets[["Fusion_Position_Details"]] <- fusion_peptides %>%
      select(Peptide, fusion_pos, junction_relative_position, position_class, 
             any_of(c("amino_acids_before_junction", "amino_acids_after_junction", 
                      "part1_seq", "part2_seq")))
  }
  
  # Add spiked peptides sheet if available
  if (exists("spiked_matrix_data")) {
    excel_sheets[["Spiked_Peptides"]] <- spiked_matrix_data
  }
  
  # Add sample matrix showing where each fusion peptide was found
  if (length(sample_ids) > 1) {
    presence_summary <- fusion_peptides %>%
      select(Peptide, fusion_peptide_type, spans_junction, position_class, 
             all_of(paste0("total_intensity_", sample_ids)))
    
    excel_sheets[["Sample_Distribution"]] <- presence_summary
  }
  
  # Generate Excel report
  excel_output_file <- file.path(dirs$excel_dir, 
                                 paste0(config$output_name, "_results.xlsx"))
  generate_excel_report(excel_sheets, excel_output_file)
  
  # 8. Save results for future use
  save(fusion_peptides, file = file.path(dirs$data_dir, "fusion_peptides.RData"))
  
  if (exists("junction_peptides") && nrow(junction_peptides) > 0) {
    save(junction_peptides, file = file.path(dirs$data_dir, "junction_peptides.RData"))
  }
  
  if (exists("spiked_matrix_data")) {
    save(spiked_matrix_data, file = file.path(dirs$data_dir, "spiked_peptides.RData"))
  }
  
  # Save multi-sample analysis results
  save(multi_sample_analysis, file = file.path(dirs$data_dir, "multi_sample_analysis.RData"))
  
  # 9. Print summary information
  cat("\n## 9. Analysis summary:\n")
  cat("\nFusion protein analysis complete!\n")
  cat("Results saved to:", dirs$main_dir, "\n")
  
  # Fusion protein details
  cat("\nFusion protein details:\n")
  cat("- Part 1 (upstream) length:", config$junction_position, "amino acids\n")
  cat("- Part 2 (downstream) length:", nchar(config$fusion_sequence) - config$junction_position, "amino acids\n")
  cat("- Total fusion protein length:", nchar(config$fusion_sequence), "amino acids\n")
  
  # Fusion peptide summary
  cat("\nFusion peptide summary:\n")
  cat("Total fusion-derived peptides:", nrow(fusion_peptides), "\n")
  
  if ("spans_junction" %in% colnames(fusion_peptides)) {
    cat("Junction-spanning peptides:", sum(fusion_peptides$spans_junction), "\n")
  }
  
  if ("position_class" %in% colnames(fusion_peptides)) {
    position_counts <- fusion_peptides %>%
      count(position_class) %>%
      mutate(percentage = round(n / sum(n) * 100, 1))
    
    print(position_counts)
  }
  
  # Spiked peptide summary if available
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
  
  # Return results invisibly
  invisible(list(
    fusion_peptides = fusion_peptides,
    junction_peptides = if(exists("junction_peptides")) junction_peptides else NULL,
    spiked_peptides = if(exists("spiked_matrix_data")) spiked_matrix_data else NULL,
    directories = dirs
  ))
}