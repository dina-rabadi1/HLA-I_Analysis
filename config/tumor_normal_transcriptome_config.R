# Configuration for tumor-normal transcriptome analysis (RU148 only)
# config/tumor_normal_transcriptome_config.R

config <- list(
  # Input data
  input = list(
    # 2CV and 3CV data directory
    raw_data_dir = "rawdata/imp_014_rawdata",
    
    # Generated combined peptide file
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
    # Only RU148 tumor and normal samples
    include = c("148T", "148N"),
    
    # No samples to exclude from shared peptide analysis
    exclude_from_shared = c()
  ),
  
  # Output settings
  output = list(
    # Base directory for results (will create subfolders)
    base_dir = "results/tumor_normal_transcriptome",
    
    # Create Excel reports
    create_excel = TRUE
  ),
  
  # Analysis settings
  analysis = list(
    # Peptide length filtering
    min_length = 8,
    max_length = 12,
    
    # Whether to perform specific analyses
    do_shared_peptide = FALSE,  # Not needed for paired analysis
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