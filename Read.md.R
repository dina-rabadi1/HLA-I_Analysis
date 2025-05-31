# HLA-I Peptide Analysis Pipeline

A flexible, modular R pipeline for analyzing immunopeptidome data, with special focus on tumor/normal comparisons, multi-sample analysis, and fusion peptide identification.

## Overview

This pipeline provides a comprehensive framework for analyzing HLA-I peptide data from mass spectrometry experiments. It supports:
  
  - **Tumor/Normal Comparisons**: Compare peptide presentation between tumor and matched normal samples
- **Multi-Sample Analysis**: Analyze peptide sharing across multiple samples (PDX models, patient cohorts)
- **Fusion Peptide Analysis**: Identify and characterize fusion protein-derived peptides
- **Multi-Omics Integration**: Combine immunopeptidome data with transcriptome and proteome (LFQ and TMT) data

## Directory Structure

```
HLA-I_Analysis/
  ├── R/
  │   ├── analyze_tumor_normal.R      # Tumor-normal analysis script
│   ├── analyze_multi_sample.R      # Multi-sample analysis script
│   ├── analyze_fusion.R            # Fusion peptide analysis script
│   ├── peptide_core_utils.R        # Core utility functions
│   ├── peptide_data_processing.R   # Data processing functions
│   ├── peptide_visualizations.R    # Visualization functions
│   ├── peptide_integration.R       # Multi-omics integration functions
│   ├── config.R                    # Configuration management
│   └── run_peptide_pipeline.R      # Master script to run analyses
├── data/                           # Directory for input data
  │   ├── patient_148/                # Example patient data folder
  │   ├── pdx_models/                 # Example PDX models folder
  │   └── ...
├── configs/                        # Saved configuration files
  │   ├── example_tumor_normal_config.rds
│   ├── example_multi_sample_config.rds
│   └── ...
└── results/                        # Analysis results will be saved here
  ├── patient_148_TN_analysis/    # Example output folder
  │   ├── visualizations/         # Generated plots (PDF, PNG, HTML)
  │   ├── excel_reports/          # Excel result files
  │   └── processed_data/         # Saved R data objects
  └── ...
```

## Installation

### Prerequisites

- R 4.1.0 or higher
- Required R packages:
  - tidyverse
- ggplot2
- pheatmap
- openxlsx
- plotly
- VennDiagram
- RColorBrewer
- ggrepel
- UpSetR
- htmlwidgets

### Setup

1. Clone this repository:
  ```
git clone https://github.com/yourusername/HLA-I_Analysis.git
cd HLA-I_Analysis
```

2. Install required R packages:
  ```R
# From R console
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, ggplot2, pheatmap, openxlsx, plotly, 
               VennDiagram, RColorBrewer, ggrepel, UpSetR, htmlwidgets)
```

## Quick Start

### Using the Command Line

1. Create a new configuration:
  ```
Rscript run_peptide_pipeline.R --create-config tumor_normal data/patient_148 p148_analysis --tumor_id 148T --normal_id 148N
```

2. Or run with an existing configuration:
  ```
Rscript run_peptide_pipeline.R configs/my_analysis_config.rds
```

3. To see example configurations:
  ```
Rscript run_peptide_pipeline.R --examples
```

### Using R

```R
# Load configuration system
source("config.R")

# Create a configuration
my_config <- create_config(
  analysis_type = "tumor_normal",
  data_path = "data/patient_148",
  output_name = "p148_analysis",
  tumor_id = "148T",
  normal_id = "148N",
  transcriptome_path = "data/patient_148/transcriptome.xlsx",
  lfq_path = "data/patient_148/lfq_proteome.xlsx"
)

# Save the configuration
save_config(my_config, "configs/my_config.rds")

# Run the pipeline
source("run_peptide_pipeline.R")
config <- load_config("configs/my_config.rds")
source("analyze_tumor_normal.R")
analyze_tumor_normal(config)
```

## Module Descriptions

### Core Modules

- **peptide_core_utils.R**: Contains shared utility functions like loading packages, creating directories, and defining standard color palettes
- **peptide_data_processing.R**: Functions for processing peptide data, creating matrices, and handling tumor/normal and multi-sample analyses
- **peptide_visualizations.R**: Standard visualization functions for heatmaps, scatter plots, volcano plots, etc.
- **peptide_integration.R**: Functions for integrating multiple omics data types and identifying public neoantigens

### Analysis Modules

- **analyze_tumor_normal.R**: Complete tumor vs. normal analysis workflow
- **analyze_multi_sample.R**: Workflow for analyzing peptides across multiple samples
- **analyze_fusion.R**: Specialized workflow for fusion peptide analysis

### Configuration and Pipeline Control

- **config.R**: System for creating, saving, and loading analysis configurations
- **run_peptide_pipeline.R**: Master script that ties everything together

## Analysis Types

### Tumor/Normal Analysis

Compares peptide presentation between tumor and matched normal samples, with options for:
  - Differential expression analysis
- Multi-omics integration (immunopeptidome, transcriptome, proteome)
- Identification of tumor-specific peptides
- Public neoantigen identification
- Fusion peptide detection

### Multi-Sample Analysis

Analyzes peptide sharing across multiple samples, with:
  - Shared vs. private peptide identification
- Sample clustering based on peptide profiles
- Overlap visualization using UpSet plots and Venn diagrams
- Fusion peptide detection across samples

### Fusion Peptide Analysis

Specialized analysis for fusion protein-derived peptides:
  - Detection of peptides spanning fusion junctions
- Positional analysis relative to the junction
- Visualization of peptide distribution along the fusion sequence
- Analysis of synthetic spiked fusion peptides

## Output

Each analysis produces:
  
  1. **Visualizations**:
  - Static plots (PDF, PNG)
- Interactive HTML visualizations
- Interactive dashboard

2. **Excel Reports**:
  - Complete analysis results in spreadsheet format
- Multiple sheets organized by analysis component

3. **Processed Data**:
  - R data objects for further analysis

## Example Usage

### Tumor/Normal Analysis

```R
# Create configuration
tumor_normal_config <- create_config(
  analysis_type = "tumor_normal",
  data_path = "data/patient_148",
  output_name = "patient_148_TN_analysis",
  tumor_id = "148T",
  normal_id = "148N",
  transcriptome_path = "data/patient_148/transcriptome.xlsx",
  lfq_path = "data/patient_148/lfq_proteome.xlsx",
  fusion_parts = list(
    part1 = "MDQAIKCYQFSSSAEPDLFRGGGMPSSEDRAEDGGSQPPASGNGPAEPTEEGGSPAPGPGR",
    part2 = "DFGFAKIVDGVAFYAKNLDPPLMAFIKIMLGKGGSGEIKELRGEVNILEIPTLQIKLCINGV"
  ),
  fusion_sequence = "MDQAIKCYQFSSSAEPDLFRGGGMPSSEDRAEDGGSQPPASGNGPAEPTEEGGSPAPGPGRDFGFAKIVDGVAFYAKNLDPPLMAFIKIMLGKGGSGEIKELRGEVNILEIPTLQIKLCINGV",
  junction_position = 60
)

# Run analysis
source("analyze_tumor_normal.R")
results <- analyze_tumor_normal(tumor_normal_config)
```

### Multi-Sample Analysis

```R
# Create configuration
multi_sample_config <- create_config(
  analysis_type = "multi_sample",
  data_path = "data/pdx_models",
  output_name = "pdx_models_comparison",
  sample_pattern = ".*_PDX([0-9]+)_.*"
)

# Run analysis
source("analyze_multi_sample.R")
results <- analyze_multi_sample(multi_sample_config)
```

### Fusion Peptide Analysis

```R
# Create configuration
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

# Run analysis
source("analyze_fusion.R")
results <- analyze_fusion(fusion_config)
```

## Customization

The pipeline is designed to be modular and customizable. You can:
  
  1. Modify visualization parameters in `peptide_visualizations.R`
2. Add new data processing functions in `peptide_data_processing.R`
3. Define custom analysis workflows by creating new analysis scripts
4. Extend the configuration system in `config.R` with new parameters

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- This pipeline was developed for analyzing HLA-I peptide data from mass spectrometry experiments
- Special thanks to [Your Name/Lab] for providing test data and use cases