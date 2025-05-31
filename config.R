#' HLA-I Analysis Pipeline Configuration
#' config.R
#' Central configuration file for the HLA-I peptide analysis pipeline
#' @author Your Name
#' @version 1.0

# Load required base packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, ggplot2)

#' Function to create a standard configuration object
#' @param analysis_type Type of analysis ("tumor_normal" or "multi_sample")
#' @param data_path Path to the data directory
#' @param output_name Name for output directories and files
#' @param ... Additional parameters specific to each analysis type
#' @return Configuration list object
create_config <- function(analysis_type = c("tumor_normal", "multi_sample", "fusion_analysis"),
                          data_path = NULL,
                          output_name = NULL,
                          ...) {
  
  # Check required parameters
  analysis_type <- match.arg(analysis_type)
  
  if (is.null(data_path)) {
    stop("data_path is required")
  }
  
  # Set default output name if not provided
  if (is.null(output_name)) {
    output_name <- paste0(analysis_type, "_", format(Sys.time(), "%Y%m%d_%H%M"))
  }
  
  # Base configuration common to all analysis types
  config <- list(
    analysis_type = analysis_type,
    data_path = data_path,
    output_name = output_name,
    timestamp = Sys.time(),
    peptide_length_filter = c(8, 12),  # Default peptide length filter
    generate_interactive = TRUE,       # Default to generating interactive visualizations
    packages = c(                      # Common required packages
      "tidyverse", "openxlsx", "ggplot2", "pheatmap", 
      "plotly", "VennDiagram", "RColorBrewer", "ggrepel"
    )
  )
  
  # Add additional parameters based on analysis type
  extra_params <- list(...)
  
  if (analysis_type == "tumor_normal") {
    # Default tumor-normal specific parameters
    tumor_normal_defaults <- list(
      tumor_id = "T",
      normal_id = "N",
      transcriptome_path = NULL,  # Set to NULL if not available
      lfq_path = NULL,            # Set to NULL if not available
      tmt_path = NULL,            # Set to NULL if not available
      sample_pattern = ".*_([TN])_.*"  # Default pattern to extract T or N samples
    )
    
    # Merge user-provided parameters with defaults
    for (param in names(tumor_normal_defaults)) {
      if (is.null(extra_params[[param]])) {
        config[[param]] <- tumor_normal_defaults[[param]]
      } else {
        config[[param]] <- extra_params[[param]]
      }
    }
    
    # Add additional packages needed for this analysis
    config$packages <- unique(c(config$packages, "UpSetR"))
    
  } else if (analysis_type == "multi_sample") {
    # Default multi-sample specific parameters
    multi_sample_defaults <- list(
      sample_pattern = ".*_([0-9]+[A-Za-z]?)_.*",  # Pattern to extract sample IDs
      spiked_peptides = NULL,     # Set to a vector of spiked peptide sequences if applicable
      spiked_sample = NULL,       # Set to the sample ID that was spiked
      min_samples_shared = 2      # Minimum number of samples for a peptide to be considered "shared"
    )
    
    # Merge user-provided parameters with defaults
    for (param in names(multi_sample_defaults)) {
      if (is.null(extra_params[[param]])) {
        config[[param]] <- multi_sample_defaults[[param]]
      } else {
        config[[param]] <- extra_params[[param]]
      }
    }
    
    # Add additional packages needed for this analysis
    config$packages <- unique(c(config$packages, "UpSetR"))
    
  } else if (analysis_type == "fusion_analysis") {
    # Default fusion-specific parameters
    fusion_defaults <- list(
      fusion_parts = list(
        part1 = NULL,  # Must be specified by user
        part2 = NULL   # Must be specified by user
      ),
      fusion_sequence = NULL,  # Must be specified by user
      junction_position = NULL,  # Must be specified by user
      sample_pattern = ".*_([0-9]+[A-Za-z]?)_.*"  # Pattern to extract sample IDs (same as multi-sample)
    )
    
    # Check if fusion parameters are provided
    if (is.null(extra_params$fusion_parts) || 
        is.null(extra_params$fusion_sequence) ||
        is.null(extra_params$junction_position)) {
      warning("Fusion analysis requires 'fusion_parts', 'fusion_sequence', and 'junction_position'\n",
              "Using default DNAJB1-PRKACA fusion parameters")
      
      # Set default fusion protein (DNAJB1-PRKACA)
      config$fusion_parts <- list(
        part1 = "MDQAIKCYQFSSSAEPDLFRGGGMPSSEDRAEDGGSQPPASGNGPAEPTEEGGSPAPGPGR",
        part2 = "DFGFAKIVDGVAFYAKNLDPPLMAFIKIMLGKGGSGEIKELRGEVNILEIPTLQIKLCINGV"
      )
      
      config$fusion_sequence <- paste0(config$fusion_parts$part1, config$fusion_parts$part2)
      config$junction_position <- nchar(config$fusion_parts$part1)
    } else {
      for (param in names(fusion_defaults)) {
        if (!is.null(extra_params[[param]])) {
          config[[param]] <- extra_params[[param]]
        }
      }
    }
  }
  
  # Merge any remaining custom parameters
  for (param in names(extra_params)) {
    if (!(param %in% names(config))) {
      config[[param]] <- extra_params[[param]]
    }
  }
  
  return(config)
}

#' Function to save a configuration
#' @param config Configuration object
#' @param file Path to save the configuration (RDS format)
save_config <- function(config, file = NULL) {
  if (is.null(file)) {
    file <- file.path(config$data_path, paste0(config$output_name, "_config.rds"))
  }
  
  saveRDS(config, file)
  
  # Also save a text version for easier inspection
  text_file <- sub("\\.rds$", ".txt", file)
  
  # Convert config to text representation
  config_text <- capture.output({
    cat("HLA-I Analysis Pipeline Configuration\n")
    cat("====================================\n")
    cat("Analysis type:", config$analysis_type, "\n")
    cat("Timestamp:", format(config$timestamp), "\n\n")
    
    for (param in names(config)) {
      if (param != "timestamp") {
        if (is.list(config[[param]]) && !is.data.frame(config[[param]])) {
          cat(param, ":\n")
          if (param == "fusion_parts") {
            cat("  part1:", substr(config[[param]]$part1, 1, 40), "...\n")
            cat("  part2:", substr(config[[param]]$part2, 1, 40), "...\n")
          } else {
            for (subparam in names(config[[param]])) {
              cat("  ", subparam, ":", config[[param]][[subparam]], "\n")
            }
          }
        } else if (is.character(config[[param]]) && length(config[[param]]) > 1) {
          cat(param, ":", paste(head(config[[param]], 5), collapse = ", "))
          if (length(config[[param]]) > 5) cat("...")
          cat("\n")
        } else {
          cat(param, ":", config[[param]], "\n")
        }
      }
    }
  })
  
  writeLines(config_text, text_file)
  
  cat("Configuration saved to:", file, "and", text_file, "\n")
  
  return(invisible(config))
}

#' Function to load a configuration
#' @param file Path to the saved configuration file
#' @return Configuration object
load_config <- function(file) {
  if (!file.exists(file)) {
    stop("Configuration file does not exist:", file)
  }
  
  config <- readRDS(file)
  
  # Ensure config has required elements
  required_elements <- c("analysis_type", "data_path", "output_name", "timestamp")
  missing_elements <- required_elements[!required_elements %in% names(config)]
  
  if (length(missing_elements) > 0) {
    stop("Configuration file is missing required elements:", 
         paste(missing_elements, collapse = ", "))
  }
  
  cat("Loaded configuration from:", file, "\n")
  cat("Analysis type:", config$analysis_type, "\n")
  cat("Timestamp:", format(config$timestamp), "\n")
  
  return(config)
}

#' Function to create example configurations
#' @param output_dir Directory to save example configurations
create_example_configs <- function(output_dir = ".") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Example 1: Tumor-Normal Analysis
  tumor_normal_config <- create_config(
    analysis_type = "tumor_normal",
    data_path = "data/patient_148",
    output_name = "patient_148_TN_analysis",
    tumor_id = "148T",
    normal_id = "148N",
    transcriptome_path = "data/patient_148/148_transcriptome.xlsx",
    lfq_path = "data/patient_148/148_lfq_proteome.xlsx",
    tmt_path = "data/patient_148/148_tmt_proteome.xlsx",
    fusion_parts = list(
      part1 = "MDQAIKCYQFSSSAEPDLFRGGGMPSSEDRAEDGGSQPPASGNGPAEPTEEGGSPAPGPGR",
      part2 = "DFGFAKIVDGVAFYAKNLDPPLMAFIKIMLGKGGSGEIKELRGEVNILEIPTLQIKLCINGV"
    ),
    fusion_sequence = "MDQAIKCYQFSSSAEPDLFRGGGMPSSEDRAEDGGSQPPASGNGPAEPTEEGGSPAPGPGRDFGFAKIVDGVAFYAKNLDPPLMAFIKIMLGKGGSGEIKELRGEVNILEIPTLQIKLCINGV",
    junction_position = 60
  )
  
  save_config(tumor_normal_config, file.path(output_dir, "example_tumor_normal_config.rds"))
  
  # Example 2: Multi-Sample Analysis
  multi_sample_config <- create_config(
    analysis_type = "multi_sample",
    data_path = "data/pdx_models",
    output_name = "pdx_models_comparison",
    sample_pattern = ".*_PDX([0-9]+)_.*",
    fusion_parts = list(
      part1 = "MDQAIKCYQFSSSAEPDLFRGGGMPSSEDRAEDGGSQPPASGNGPAEPTEEGGSPAPGPGR",
      part2 = "DFGFAKIVDGVAFYAKNLDPPLMAFIKIMLGKGGSGEIKELRGEVNILEIPTLQIKLCINGV"
    ),
    fusion_sequence = "MDQAIKCYQFSSSAEPDLFRGGGMPSSEDRAEDGGSQPPASGNGPAEPTEEGGSPAPGPGRDFGFAKIVDGVAFYAKNLDPPLMAFIKIMLGKGGSGEIKELRGEVNILEIPTLQIKLCINGV",
    junction_position = 60,
    min_samples_shared = 3
  )
  
  save_config(multi_sample_config, file.path(output_dir, "example_multi_sample_config.rds"))
  
  # Example 3: Fusion Analysis with Spiked Peptides
  fusion_config <- create_config(
    analysis_type = "fusion_analysis",
    data_path = "data/spiked_experiment",
    output_name = "fusion_spiked_analysis",
    fusion_parts = list(
      part1 = "MDQAIKCYQFSSSAEPDLFRGGGMPSSEDRAEDGGSQPPASGNGPAEPTEEGGSPAPGPGR",
      part2 = "DFGFAKIVDGVAFYAKNLDPPLMAFIKIMLGKGGSGEIKELRGEVNILEIPTLQIKLCINGV"
    ),
    fusion_sequence = "MDQAIKCYQFSSSAEPDLFRGGGMPSSEDRAEDGGSQPPASGNGPAEPTEEGGSPAPGPGRDFGFAKIVDGVAFYAKNLDPPLMAFIKIMLGKGGSGEIKELRGEVNILEIPTLQIKLCINGV",
    junction_position = 60,
    spiked_peptides = c(
      "EGGSPAPGP",
      "PGPGRDFGF",
      "GSPAPGPGRD",
      "PAPGPGRDFG",
      "SPAPGPGRDF",
      "APGPGRDFGF"
    ),
    spiked_sample = "51S"
  )
  
  save_config(fusion_config, file.path(output_dir, "example_fusion_spiked_config.rds"))
  
  cat("Example configurations created in", output_dir, "\n")
}