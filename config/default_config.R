# Default configuration for HLA-I peptide analysis
# This file defines settings for the analysis pipeline
# config/default_config.R
# This is the default config but also the all-sample (excluding 148N and 51S)

config <- list(
  # Input data
  input = list(
    # 2CV and 3CV data directory
    raw_data_dir = "rawdata/imp_014_rawdata",
    
    # Combined peptide file (will be created if concatenation is enabled)
    peptide_file = "results/peptide_analysis/data/combined_peptides.tsv",
    
    # Transcriptome data (if needed)
    transcriptome_file = "rawdata/Normalized_Gene_counts_FLCdb_Panel_1.xlsx",
    
    # Proteome data (if needed)
    lfq_file = "rawdata/Levin2023/adg7038_Table_S2_LFQ.xlsx",
    tmt_file = "rawdata/Levin2023/adg7038_Table_S1_TMT.xlsx"
  ),
  
  # Data preprocessing
  preprocessing = list(
    # Whether to perform 2CV/3CV concatenation
    concatenate_2cv_3cv = TRUE,
    
    # Pattern to match 2CV files
    pattern_2cv = ".*_2CV_.*_peptides\\.tsv$",
    
    # Pattern to match 3CV files
    pattern_3cv = ".*_3CV_.*_peptides\\.tsv$"
  ),
  
  # Samples to include/exclude
  samples = list(
    # All sample IDs to analyze
    include = c("51", "57", "59", "62", "63", "88", "117", "123", "148T"),
    
    # Samples to exclude from shared peptide analysis
    exclude_from_shared = c("51S", "148N")
  ),
  
  # Output settings
  output = list(
    # Base directory for results (will create subfolders)
    base_dir = "results/all_sample_analysis",
    
    # Create Excel reports
    create_excel = TRUE
  ),

  # Analysis settings
  analysis = list(
    # Peptide length filtering
    min_length = 8,
    max_length = 12,
    
    # Whether to perform specific analyses
    do_shared_peptide = TRUE,
    do_amino_acid = TRUE,
    do_fusion = TRUE,
    do_multi_omics = TRUE,
    do_length_distribution = TRUE,
    
    # Empty sample_comparisons section (used by spike-in analysis)
    sample_comparisons = list(
      do_51_vs_51S = FALSE
    )
  ),
  
  # Visualization settings
  visualization = list(
    # Color schemes
    color_scheme = "viridis",
    
    # Generate both PDF and PNG
    output_formats = c("pdf", "png"),
    
    # Plot dimensions
    plot_width = 10,
    plot_height = 8,
    
    # PNG resolution (dots per inch)
    png_dpi = 300
  )
)