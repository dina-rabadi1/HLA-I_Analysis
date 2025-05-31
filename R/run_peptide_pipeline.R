#!/usr/bin/env Rscript
#' HLA-I Analysis Pipeline
#' run_peptide_pipeline.R
#' Master script for running the modular HLA-I peptide analysis pipeline
#' @author Your Name
#' @version 1.0

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Load configuration system
source("config.R")

# Show usage information if no arguments provided
if (length(args) == 0) {
  cat("Usage: Rscript run_peptide_pipeline.R <config_file> [options]\n")
  cat("OR:    Rscript run_peptide_pipeline.R --create-config <analysis_type> <data_path> <output_name> [options]\n\n")
  
  cat("Available analysis types:\n")
  cat("  tumor_normal   - Compare tumor and normal samples\n")
  cat("  multi_sample   - Analyze peptides across multiple samples\n")
  cat("  fusion_analysis - Focus on fusion protein analysis\n\n")
  
  cat("Examples:\n")
  cat("  Rscript run_peptide_pipeline.R configs/my_analysis_config.rds\n")
  cat("  Rscript run_peptide_pipeline.R --create-config tumor_normal data/patient_148 p148_analysis --tumor_id 148T --normal_id 148N\n\n")
  
  # Create example configs if requested
  if (length(args) > 0 && args[1] == "--examples") {
    cat("Creating example configuration files...\n")
    create_example_configs("example_configs")
  }
  
  quit(status = 0)
}

# Check if creating a new config
if (args[1] == "--create-config") {
  if (length(args) < 4) {
    stop("Creating a config requires at least: analysis_type, data_path, and output_name")
  }
  
  analysis_type <- args[2]
  data_path <- args[3]
  output_name <- args[4]
  
  # Process additional parameters
  extra_params <- list()
  if (length(args) > 4) {
    param_args <- args[5:length(args)]
    for (i in seq(1, length(param_args), by = 2)) {
      if (i + 1 <= length(param_args)) {
        param_name <- sub("^--", "", param_args[i])
        param_value <- param_args[i + 1]
        
        # Convert strings to appropriate types when possible
        if (param_value == "TRUE" || param_value == "true") {
          param_value <- TRUE
        } else if (param_value == "FALSE" || param_value == "false") {
          param_value <- FALSE
        } else if (param_value == "NULL" || param_value == "null") {
          param_value <- NULL
        } else if (grepl("^[0-9]+$", param_value)) {
          param_value <- as.integer(param_value)
        } else if (grepl("^[0-9]*\\.[0-9]+$", param_value)) {
          param_value <- as.numeric(param_value)
        }
        
        extra_params[[param_name]] <- param_value
      }
    }
  }
  
  # Create the configuration
  config <- do.call(create_config, c(
    list(analysis_type = analysis_type, data_path = data_path, output_name = output_name),
    extra_params
  ))
  
  # Save the configuration
  config_file <- file.path(data_path, paste0(output_name, "_config.rds"))
  save_config(config, config_file)
  
  # Ask if the user wants to run the analysis now
  cat("\nConfiguration created and saved to:", config_file, "\n")
  cat("Run analysis now? (y/n): ")
  run_now <- readline()
  
  if (tolower(run_now) != "y") {
    cat("Exiting. Run the analysis later with: Rscript run_peptide_pipeline.R", config_file, "\n")
    quit(status = 0)
  }
} else {
  # Load existing configuration
  config_file <- args[1]
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file)
  }
  
  config <- load_config(config_file)
}

# Load required modules based on analysis type
source("peptide_core_utils.R")
source("peptide_data_processing.R")
source("peptide_visualizations.R")
source("peptide_integration.R")

# Load required packages
cat("\nLoading required packages...\n")
load_required_packages(config$packages)

# Create output directories
cat("\nCreating output directories...\n")
dirs <- create_output_directories(config$data_path, config$output_name)

# Run the appropriate analysis based on configuration
if (config$analysis_type == "tumor_normal") {
  cat("\nRunning tumor-normal analysis...\n")
  source("analyze_tumor_normal.R")
  
  # Execute tumor-normal analysis function
  analyze_tumor_normal(config, dirs)
  
} else if (config$analysis_type == "multi_sample") {
  cat("\nRunning multi-sample analysis...\n")
  source("analyze_multi_sample.R")
  
  # Execute multi-sample analysis function
  analyze_multi_sample(config, dirs)
  
} else if (config$analysis_type == "fusion_analysis") {
  cat("\nRunning fusion protein analysis...\n")
  source("analyze_fusion.R")
  
  # Execute fusion analysis function
  analyze_fusion(config, dirs)
  
} else {
  stop("Unknown analysis type: ", config$analysis_type)
}

# Save the final configuration including any updates
save_config(config, file.path(dirs$main_dir, paste0(config$output_name, "_final_config.rds")))

# Create a readme file
readme_file <- file.path(dirs$main_dir, "README.txt")
cat(
  "HLA-I Peptide Analysis", "\n",
  "=====================", "\n\n",
  "Analysis type: ", config$analysis_type, "\n",
  "Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
  "This directory contains the results of the HLA-I peptide analysis.", "\n",
  "The analysis was performed using the HLA-I Analysis Pipeline.", "\n\n",
  "Directories:", "\n",
  "- ", basename(dirs$viz_dir), ": Visualizations (PDF, PNG, and interactive HTML)", "\n",
  "- ", basename(dirs$excel_dir), ": Excel reports with detailed results", "\n",
  "- ", basename(dirs$data_dir), ": Processed data files (RData format)", "\n\n",
  "Configuration:", "\n",
  "- A copy of the configuration file can be found at:", "\n  ",
  file.path(dirs$main_dir, paste0(config$output_name, "_final_config.rds")), "\n\n",
  "For more information, please refer to the pipeline documentation.", "\n",
  file = readme_file
)

cat("\nAnalysis complete!\n")
cat("Results saved to:", dirs$main_dir, "\n")
cat("See", readme_file, "for more information.\n")