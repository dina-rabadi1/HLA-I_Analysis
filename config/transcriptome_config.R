# Configuration for transcriptome analysis with all samples
# config/transcriptome_config.R

config <- list(
  # Input data
  input = list(
    # 2CV and 3CV data directory
    raw_data_dir = "rawdata/imp_014_rawdata",
    
    # Combined peptide file (will be created if concatenation is enabled)
    peptide_file = "results/peptide_analysis/data/combined_peptides.tsv",
    
    # Transcriptome data
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
    include = c("51", "59", "63", "88", "117", "123", "148T"),
    
    # Samples to exclude from shared peptide analysis
    exclude_from_shared = c("51S", "57", "62", "148N")
  ),
  
  # Output settings
  output = list(
    # Base directory for results (will create subfolders)
    base_dir = "results/transcriptome_analysis",
    
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
    do_length_distribution = TRUE,
    
    # Enable transcriptome analysis
    do_transcriptome = TRUE,
    
    # Enable multi-omics analysis
    do_multi_omics = TRUE,
    
    # Tumor-normal analysis (for RU148 pair)
    tumor_normal = list(
      do_tumor_normal = TRUE,
      pairs = list(
        list(tumor = "148T", normal = "148N")
      )
    ),
    
    # Transcriptome analysis settings
    transcriptome = list(
      # Define sharing thresholds for shared peptide analysis
      # Will analyze peptides shared in at least these many samples
      sharing_thresholds = c(5, 6, 7)
    ),
    
    # Sample comparisons (disable for now)
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