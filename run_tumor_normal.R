# Load required libraries
library(tidyverse)
library(ggplot2)
library(openxlsx)
library(pheatmap)
library(ggrepel)  # For text label repulsion in plots
library(RColorBrewer)  # For color palettes

# Create a tumor-normal configuration
config <- list(
  data_path = "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/data/imp_014_rawdata",  # Your data path
  output_name = "tumor_normal_analysis",
  tumor_id = "148T",  # Your tumor ID pattern
  normal_id = "148N",  # Your normal ID pattern
  peptide_length_filter = c(8, 12),
  generate_interactive = TRUE
)

# Create output directories first
dirs <- create_output_directories(config$data_path, config$output_name)

# Source required utility scripts
source("peptide_core_utils.R")

# Source the processing scripts (make sure they don't reference viz_dir directly)
source("peptide_data_processing.R")
source("peptide_integration.R")
source("peptide_visualizations.R")

# Run tumor-normal analysis
source("analyze_tumor_normal.R")  # Using the fixed version you showed
analyze_tumor_normal(config, dirs)