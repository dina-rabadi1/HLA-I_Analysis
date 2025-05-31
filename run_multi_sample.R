# run_multi_sample.R
# Load required libraries
library(tidyverse)

# Source required scripts
source("peptide_core_utils.R")
source("peptide_data_processing.R")
source("peptide_visualizations.R")
source("peptide_integration.R")
source("config.R")

# Load the configuration
config <- load_config("configs/multi_sample_config.rds")

# Run the analysis
source("analyze_multi_sample.R")
# The main function is called inside analyze_multi_sample.R