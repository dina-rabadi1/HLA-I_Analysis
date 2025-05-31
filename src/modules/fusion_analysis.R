# src/modules/fusion_analysis.R
# DNAJB1:PRKACA fusion protein analysis module

#' Run DNAJB1:PRKACA fusion protein analysis
#'
#' This function analyzes peptides that span the fusion junction between
#' DNAJB1 and PRKACA proteins across all samples.
#'
#' @param peptide_data Dataframe containing peptide data
#' @param config Configuration list
#' @return List containing analysis results
run_fusion_analysis <- function(peptide_data, config) {
  cat("Starting DNAJB1:PRKACA fusion analysis...\n")
  
  # Define the fusion sequences
  dnajb1_seq <- "KREIFDRYGEE"
  prkaca_seq <- "VKEFLAKAKED"
  
  # Generate all possible fusion-spanning peptides (8-12mers)
  fusion_peptides <- generate_fusion_peptides(dnajb1_seq, prkaca_seq, 
                                              min_length = config$analysis$min_length,
                                              max_length = config$analysis$max_length)
  
  # Identify fusion peptides in the data
  fusion_results <- identify_fusion_peptides(peptide_data, fusion_peptides, config)
  
  # Calculate metrics
  fusion_metrics <- calculate_fusion_metrics(fusion_results, fusion_peptides, config)
  
  # Prepare matrix for visualization
  fusion_matrix <- prepare_fusion_matrix(fusion_results, config)
  
  # Return results
  return(list(
    fusion_results = fusion_results,
    fusion_metrics = fusion_metrics,
    fusion_matrix = fusion_matrix,
    fusion_peptides = fusion_peptides
  ))
}

#' Generate all possible fusion-spanning peptides
#'
#' @param seq1 First protein sequence
#' @param seq2 Second protein sequence
#' @param min_length Minimum peptide length
#' @param max_length Maximum peptide length
#' @return Dataframe of fusion peptides
generate_fusion_peptides <- function(seq1, seq2, min_length, max_length) {
  cat("Generating all possible fusion-spanning peptides...\n")
  
  # Create empty list to store peptides
  fusion_peptides <- list()
  
  # For each possible length
  for (len in min_length:max_length) {
    # For each possible split point
    for (i in 1:(len-1)) {
      # Need at least 1 amino acid from each protein
      if (i <= nchar(seq1) && (len - i) <= nchar(seq2)) {
        # Extract segments from each protein
        seq1_part <- substr(seq1, nchar(seq1) - i + 1, nchar(seq1))
        seq2_part <- substr(seq2, 1, len - i)
        
        # Combine to create fusion peptide
        fusion_peptide <- paste0(seq1_part, seq2_part)
        
        # Add to list with details
        fusion_peptides[[length(fusion_peptides) + 1]] <- list(
          peptide = fusion_peptide,
          length = len,
          seq1_contribution = i,
          seq2_contribution = len - i,
          seq1_part = seq1_part,
          seq2_part = seq2_part
        )
      }
    }
  }
  
  # Convert list to dataframe
  fusion_peptides_df <- do.call(rbind, lapply(fusion_peptides, function(x) {
    data.frame(
      peptide = x$peptide,
      length = x$length,
      seq1_contribution = x$seq1_contribution,
      seq2_contribution = x$seq2_contribution,
      seq1_part = x$seq1_part,
      seq2_part = x$seq2_part,
      stringsAsFactors = FALSE
    )
  }))
  
  cat(sprintf("Generated %d possible fusion-spanning peptides\n", nrow(fusion_peptides_df)))
  return(fusion_peptides_df)
}

#' Identify fusion peptides in the peptide data
#'
#' @param peptide_data Dataframe containing peptide data
#' @param fusion_peptides Dataframe of fusion peptides
#' @param config Configuration list
#' @return Dataframe of identified fusion peptides
identify_fusion_peptides <- function(peptide_data, fusion_peptides, config) {
  cat("Identifying fusion peptides in the dataset...\n")
  
  # Get the list of samples to analyze
  samples <- config$samples$include
  
  # Create a list to store results
  results <- list()
  
  # Check each sample
  for (sample in samples) {
    cat(sprintf("Processing sample %s...\n", sample))
    
    # Filter peptide data for current sample
    sample_data <- peptide_data %>% 
      filter(SampleID == sample)
    
    # Match fusion peptides
    matches <- sample_data %>%
      filter(Peptide %in% fusion_peptides$peptide)
    
    # Skip if no matches
    if (nrow(matches) == 0) {
      cat(sprintf("No fusion peptides found in sample %s\n", sample))
      next
    }
    
    # Add fusion peptide details
    matches <- matches %>%
      left_join(fusion_peptides, by = c("Peptide" = "peptide"))
    
    # Add to results
    results[[sample]] <- matches
  }
  
  # Combine results
  if (length(results) > 0) {
    all_results <- bind_rows(results)
    cat(sprintf("Found %d fusion peptides across all samples\n", nrow(all_results)))
    return(all_results)
  } else {
    cat("No fusion peptides found in any sample\n")
    return(data.frame())
  }
}

#' Calculate fusion peptide metrics
#'
#' @param fusion_results Dataframe of fusion peptides
#' @param fusion_peptides Dataframe of all possible fusion peptides
#' @param config Configuration list
#' @return List of metrics
calculate_fusion_metrics <- function(fusion_results, fusion_peptides, config) {
  cat("Calculating fusion peptide metrics...\n")
  
  # If no fusion peptides found, return empty metrics
  if (nrow(fusion_results) == 0) {
    return(list(
      total_fusion_peptides = 0,
      unique_fusion_peptides = 0,
      samples_with_fusion = 0,
      samples_analyzed = length(config$samples$include),
      peptide_length_counts = data.frame(
        length = config$analysis$min_length:config$analysis$max_length,
        count = 0
      ),
      sample_counts = data.frame(
        SampleID = config$samples$include,
        unique_peptides = 0,
        total_peptides = 0
      )
    ))
  }
  
  # Calculate metrics
  total_fusion_peptides <- nrow(fusion_results)
  unique_fusion_peptides <- length(unique(fusion_results$Peptide))
  samples_with_fusion <- length(unique(fusion_results$SampleID))
  
  # Count by peptide length
  peptide_length_counts <- fusion_results %>%
    group_by(length) %>%
    summarize(count = n()) %>%
    arrange(length)
  
  # Ensure all lengths in the range are included
  all_lengths <- data.frame(
    length = config$analysis$min_length:config$analysis$max_length
  )
  peptide_length_counts <- all_lengths %>%
    left_join(peptide_length_counts, by = "length") %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  # Count by sample
  sample_counts <- fusion_results %>%
    group_by(SampleID) %>%
    summarize(
      unique_peptides = n_distinct(Peptide),
      total_peptides = n(),
      avg_intensity = mean(Intensity, na.rm = TRUE)
    )
  
  # Ensure all samples are included
  all_samples <- data.frame(
    SampleID = config$samples$include
  )
  sample_counts <- all_samples %>%
    left_join(sample_counts, by = "SampleID") %>%
    mutate(
      unique_peptides = ifelse(is.na(unique_peptides), 0, unique_peptides),
      total_peptides = ifelse(is.na(total_peptides), 0, total_peptides),
      avg_intensity = ifelse(is.na(avg_intensity), 0, avg_intensity)
    )
  
  # Calculate detection statistics
  total_possible_peptides <- nrow(fusion_peptides)
  detected_peptides <- unique_fusion_peptides
  detection_rate <- detected_peptides / total_possible_peptides * 100
  
  # Return metrics
  return(list(
    total_fusion_peptides = total_fusion_peptides,
    unique_fusion_peptides = unique_fusion_peptides,
    samples_with_fusion = samples_with_fusion,
    samples_analyzed = length(config$samples$include),
    total_possible_peptides = total_possible_peptides,
    detection_rate = detection_rate,
    peptide_length_counts = peptide_length_counts,
    sample_counts = sample_counts
  ))
}

#' Prepare fusion matrix for visualization
#'
#' @param fusion_results Dataframe of fusion peptides
#' @param config Configuration list
#' @return Matrix for visualization
prepare_fusion_matrix <- function(fusion_results, config) {
  cat("Preparing fusion peptide matrix for visualization...\n")
  
  # If no fusion peptides found, return empty matrix
  if (nrow(fusion_results) == 0) {
    # Create a matrix with sample columns but no rows
    empty_matrix <- matrix(0, 
                           nrow = 0, 
                           ncol = length(config$samples$include))
    colnames(empty_matrix) <- config$samples$include
    return(list(
      presence_matrix = empty_matrix,
      intensity_matrix = empty_matrix
    ))
  }
  
  # Get all samples and peptides
  all_samples <- config$samples$include
  all_peptides <- unique(fusion_results$Peptide)
  
  # Create a matrix of presence/absence
  presence_matrix <- matrix(0, 
                            nrow = length(all_peptides), 
                            ncol = length(all_samples))
  
  rownames(presence_matrix) <- all_peptides
  colnames(presence_matrix) <- all_samples
  
  # Create a matrix of intensities
  intensity_matrix <- matrix(0, 
                             nrow = length(all_peptides), 
                             ncol = length(all_samples))
  
  rownames(intensity_matrix) <- all_peptides
  colnames(intensity_matrix) <- all_samples
  
  # Fill the matrices
  for (peptide in all_peptides) {
    peptide_data <- fusion_results %>%
      filter(Peptide == peptide)
    
    for (sample in unique(peptide_data$SampleID)) {
      sample_data <- peptide_data %>%
        filter(SampleID == sample)
      
      peptide_idx <- which(all_peptides == peptide)
      sample_idx <- which(all_samples == sample)
      
      if (length(peptide_idx) > 0 && length(sample_idx) > 0) {
        # Mark as present
        presence_matrix[peptide_idx, sample_idx] <- 1
        
        # Average intensity if multiple matches
        avg_intensity <- mean(sample_data$Intensity, na.rm = TRUE)
        intensity_matrix[peptide_idx, sample_idx] <- avg_intensity
      }
    }
  }
  
  # Order peptides by their position in the fusion junction
  peptide_metadata <- fusion_results %>%
    select(Peptide, seq1_contribution, seq2_contribution) %>%
    distinct() %>%
    arrange(desc(seq1_contribution), seq2_contribution)
  
  # Reorder matrices
  ordered_peptides <- peptide_metadata$Peptide
  ordered_idx <- match(ordered_peptides, rownames(presence_matrix))
  ordered_idx <- ordered_idx[!is.na(ordered_idx)]
  
  if (length(ordered_idx) > 0) {
    presence_matrix <- presence_matrix[ordered_idx, , drop = FALSE]
    intensity_matrix <- intensity_matrix[ordered_idx, , drop = FALSE]
  }
  
  return(list(
    presence_matrix = presence_matrix,
    intensity_matrix = intensity_matrix
  ))
}