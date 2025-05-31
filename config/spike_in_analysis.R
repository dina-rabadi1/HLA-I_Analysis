# Configuration for spike-in peptide analysis
# This focuses on 51 vs 51S spike-in comparison

config <- list(
  # Input data
  input = list(
    # Raw data directory
    raw_data_dir = "rawdata/imp_014_rawdata",
    
    # Peptide file (will be created if concatenation is enabled)
    peptide_file = "results/peptide_analysis/data/combined_peptides.tsv",
    
    # Other files
    transcriptome_file = "rawdata/Normalized_Gene_counts_FLCdb_Panel_1.xlsx",
    lfq_file = "rawdata/Levin2023/adg7038_Table_S2_LFQ.xlsx",
    tmt_file = "rawdata/Levin2023/adg7038_Table_S1_TMT.xlsx"
  ),
  
  # Preprocessing settings
  preprocessing = list(
    # Whether to perform 2CV/3CV concatenation
    concatenate_2cv_3cv = TRUE,
    
    # Patterns for file matching
    pattern_2cv = ".*_2CV_.*_peptides\\.tsv$",
    pattern_3cv = ".*_3CV_.*_peptides\\.tsv$"
  ),
  
  # Samples to include/exclude
  samples = list(
    # Include both 51 and 51S for the spike-in comparison
    include = c("51", "51S"),
    
    # No samples to exclude from this focused analysis
    exclude_from_shared = c()
  ),
  
  # Output settings
  output = list(
    # Base directory for results
    base_dir = "results/spike_in_analysis",
    
    # Create Excel report
    create_excel = TRUE
  ),
  
  # Analysis settings
  analysis = list(
    # Peptide length filtering
    min_length = 8,
    max_length = 12,
    
    # Which analyses to run
    do_shared_peptide = TRUE,
    do_amino_acid = FALSE,
    do_fusion = FALSE,
    do_multi_omics = FALSE,
    do_length_distribution = TRUE,
    
    # Spike-in analysis settings
    sample_comparisons = list(
      # Enable 51 vs 51S comparison
      do_51_vs_51S = TRUE,
      
      # Sample IDs
      original_sample = "51",
      spiked_sample = "51S",
      
      # Known spike-in peptides
      spiked_peptides = c(
        "EEVKEFLAK", 
        "YGEEVKEFL", 
        "RYGEEVKEF", 
        "RYGEEVKEFL", 
        "EIFDRYGEEV", 
        "IFDRYGEEV"
      )
    )
  ),
  
  # Visualization settings
  visualization = list(
    # Generate both PDF and PNG
    output_formats = c("pdf", "png"),
    
    # Plot dimensions
    plot_width = 10,
    plot_height = 8,
    
    # PNG resolution (dots per inch)
    png_dpi = 300
  )
)