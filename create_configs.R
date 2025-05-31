# Source the config file
# Run config.R then create_configs.R
source("config.R")

# Create a directory for configs if it doesn't exist
dir.create("configs", showWarnings = FALSE)

# Example 1: Tumor-Normal Analysis
tumor_normal_config <- create_config(
  analysis_type = "tumor_normal",
  data_path = "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/data/imp_014_rawdata",  # Replace with your actual data path
  output_name = "tumor_normal_analysis",
  tumor_id = "148T",  # Replace with your actual tumor ID
  normal_id = "148N",  # Replace with your actual normal ID
  
  # Optional paths - set to NULL if not available
  transcriptome_path = "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/data/Normalized_Gene_counts_FLCdb_Panel_1.xlsx",  # Path to transcriptome data file
  lfq_path = "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/data/Levin2023/adg7038_Table_S2_LFQ.xlsx",  # Path to LFQ proteomics data
  tmt_path = "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/data/Levin2023/adg7038_Table_S1_TMT.xlsx",  # Path to TMT proteomics data
  
  # Fusion analysis parameters (optional)
  fusion_parts = list(
    part1 = "SEQUENCE_PART1",  # Replace with actual sequence
    part2 = "SEQUENCE_PART2"   # Replace with actual sequence
  ),
  fusion_sequence = "SEQUENCE_PART1SEQUENCE_PART2",  # Combined sequence
  junction_position = 12  # Position where the fusion occurs
)

# Save tumor-normal config
save_config(tumor_normal_config, "configs/tumor_normal_config.rds")

# Example 2: Multi-Sample Analysis
multi_sample_config <- create_config(
  analysis_type = "multi_sample",
  data_path = "path/to/your/data",  # Replace with your actual data path
  output_name = "multi_sample_analysis",
  sample_pattern = ".*_([0-9]+[A-Za-z]?)_.*",  # Pattern to extract sample IDs
  
  # Optional spiked peptide analysis
  spiked_peptides = NULL,  # Vector of spiked peptide sequences
  spiked_sample = NULL     # Sample ID that was spiked
)

# Save multi-sample config
save_config(multi_sample_config, "configs/multi_sample_config.rds")