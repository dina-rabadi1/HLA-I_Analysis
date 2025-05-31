# HLA-I Analysis Pipeline Structure Visualization (Final Fixed Version)
# This script creates a visualization with no cut-offs anywhere

# Install required packages if not already installed
if(!require(igraph)) install.packages("igraph")

# Load required libraries
library(igraph)

# Create a parameter consistency diagram with no cut-offs
create_final_diagram <- function() {
  # Create a plot with very large dimensions and more compact design
  png("hla_pipeline_parameters_compact.png", width = 1600, height = 900, res = 100)
  
  # Set up empty plot with wider dimensions and more margin space
  # Increase right and left margins substantially
  par(mar = c(2, 2, 3, 2))
  plot(0, type = "n", xlim = c(0, 16), ylim = c(0, 10), 
       xlab = "", ylab = "", main = "",
       axes = FALSE)
  
  # Draw config box - moved further left
  rect(0.5, 7.5, 4, 9, col = "#FFFFE0", border = "black")
  text(2.25, 8.5, "Configuration Files", font = 2, cex = 1)
  text(0.7, 8.2, "default_config.R", adj = 0, cex = 0.8)
  text(0.7, 7.9, "transcriptome_config.R", adj = 0, cex = 0.8)
  text(0.7, 7.6, "tumor_normal_only.R", adj = 0, cex = 0.8)
  
  # Draw source box - centered better
  rect(7, 7.5, 10.5, 9, col = "#ADD8E6", border = "black")
  text(8.75, 8.5, "Source Code", font = 2, cex = 1)
  text(7.2, 8.2, "src/core/peptide_core_utils.R", adj = 0, cex = 0.8)
  text(7.2, 7.9, "src/modules/peptide_integration.R", adj = 0, cex = 0.8)
  text(7.2, 7.6, "src/modules/peptide_search", adj = 0, cex = 0.8)
  
  # Draw legend in a clearly visible box to the right
  rect(12, 6, 15.5, 9, col = "#F5F5F5", border = "gray50", lty = 2)
  
  # Add parameter consistency legend inside the box
  text(12.2, 8.5, "Parameter Consistency Flow:", adj = 0, font = 2, cex = 0.9)
  
  # Color dots for legend
  points(12.5, 8, pch = 19, col = "red", cex = 1)
  text(13, 8, "Configuration parameters", adj = 0, cex = 0.8)
  
  points(12.5, 7.5, pch = 19, col = "blue", cex = 1)
  text(13, 7.5, "Shared code utilities", adj = 0, cex = 0.8)
  
  # Example parameters - inside legend box
  text(12.2, 7, "Example Parameters:", adj = 0, font = 2, cex = 0.9)
  text(12.5, 6.6, "- HLA typing parameters", adj = 0, cex = 0.75)
  text(12.5, 6.3, "- Peptide search thresholds", adj = 0, cex = 0.75)
  
  # Benefits box under the legend box
  rect(12, 3.5, 15.5, 5.5, col = "#E8F8E8", border = "darkgreen", lty = 2)
  text(12.2, 5.1, "Benefits of Parameter Consistency:", adj = 0, font = 2, cex = 0.9, col = "darkgreen")
  text(12.5, 4.7, "- Reproducible analyses", adj = 0, cex = 0.75, col = "darkgreen")
  text(12.5, 4.3, "- Comparable results", adj = 0, cex = 0.75, col = "darkgreen")
  text(12.5, 3.9, "- Consistent file paths", adj = 0, cex = 0.75, col = "darkgreen")
  
  # Draw pipeline scripts with proper spacing
  rect(1.5, 5.5, 4.5, 6.5, col = "#F08080", border = "black")
  text(3, 6, "run_peptide_pipeline.R", cex = 0.9)
  
  rect(7, 5.5, 10, 6.5, col = "#F08080", border = "black")
  text(8.5, 6, "run_peptide_transcriptome_pipeline.R", cex = 0.8)
  
  # Draw primary results directories with proper spacing
  rect(1, 3.5, 3.5, 4.5, col = "#90EE90", border = "black")
  text(2.25, 4, "results/peptide_analysis", cex = 0.8)
  
  rect(4.5, 3.5, 7, 4.5, col = "#90EE90", border = "black")
  text(5.75, 4, "results/transcriptome_analysis", cex = 0.8)
  
  rect(8, 3.5, 10.5, 4.5, col = "#90EE90", border = "black")
  text(9.25, 4, "results/tumor_normal_transcriptome", cex = 0.7)
  
  # Secondary results directories with proper spacing
  rect(1, 1.5, 3.5, 2.5, col = "#90EE90", border = "black")
  text(2.25, 2, "results/148tumor_normal", cex = 0.8)
  
  rect(4.5, 1.5, 7, 2.5, col = "#90EE90", border = "black")
  text(5.75, 2, "results/all_sample_analysis", cex = 0.8)
  
  rect(8, 1.5, 10.5, 2.5, col = "#90EE90", border = "black")
  text(9.25, 2, "results/spike_in_analysis", cex = 0.8)
  
  # Integration box
  rect(4.5, 0, 7, 1, col = "#B19CD9", border = "black")
  text(5.75, 0.5, "Integrated Analysis", font = 2, cex = 0.9)
  
  # Draw parameter flow arrows with proper spacing
  # From config to pipelines
  arrows(2.25, 7.5, 3, 6.5, length = 0.1, lwd = 2, col = "red")
  arrows(2.25, 7.5, 8.5, 6.5, length = 0.1, lwd = 2, col = "red")
  
  # From source to pipelines
  arrows(8.75, 7.5, 3, 6.5, length = 0.1, lwd = 2, col = "blue")
  arrows(8.75, 7.5, 8.5, 6.5, length = 0.1, lwd = 2, col = "blue")
  
  # From pipelines to primary results
  arrows(3, 5.5, 2.25, 4.5, length = 0.1, lwd = 1.5)
  arrows(3, 5.5, 5.75, 4.5, length = 0.1, lwd = 1.5)
  arrows(8.5, 5.5, 5.75, 4.5, length = 0.1, lwd = 1.5)
  arrows(8.5, 5.5, 9.25, 4.5, length = 0.1, lwd = 1.5)
  
  # From primary results to secondary results
  arrows(2.25, 3.5, 2.25, 2.5, length = 0.1, lwd = 1.5)
  arrows(5.75, 3.5, 5.75, 2.5, length = 0.1, lwd = 1.5)
  arrows(9.25, 3.5, 9.25, 2.5, length = 0.1, lwd = 1.5)
  
  # From secondary results to final integration
  arrows(2.25, 1.5, 5, 1, length = 0.1, lwd = 1.5)
  arrows(5.75, 1.5, 5.75, 1, length = 0.1, lwd = 1.5)
  arrows(9.25, 1.5, 6.5, 1, length = 0.1, lwd = 1.5)
  
  # Add main title with better positioning
  text(8, 9.5, "HLA-I Analysis Pipeline - Parameter Consistency Flow", font = 2, cex = 1.3)
  
  dev.off()
  
  message("Final parameter consistency diagram created with more compact layout: hla_pipeline_parameters_compact.png")
  return(TRUE)
}

# Run the final diagram function
create_final_diagram()

# Also create an even safer version with very compact layout
create_safe_diagram <- function() {
  # Create a plot with standard dimensions but very compact layout
  png("hla_pipeline_parameters_safe.png", width = 1200, height = 800, res = 100)
  
  # Set up empty plot with compact dimensions
  par(mar = c(1, 1, 2, 1))
  plot(0, type = "n", xlim = c(0, 10), ylim = c(0, 10), 
       xlab = "", ylab = "", main = "",
       axes = FALSE)
  
  # Draw config box - top left
  rect(0.5, 8, 3, 9.5, col = "#FFFFE0", border = "black")
  text(1.75, 9, "Configuration Files", font = 2, cex = 0.8)
  text(0.7, 8.7, "default_config.R", adj = 0, cex = 0.7)
  text(0.7, 8.4, "transcriptome_config.R", adj = 0, cex = 0.7)
  text(0.7, 8.1, "tumor_normal_only.R", adj = 0, cex = 0.7)
  
  # Draw source box - top right
  rect(7, 8, 9.5, 9.5, col = "#ADD8E6", border = "black")
  text(8.25, 9, "Source Code", font = 2, cex = 0.8)
  text(7.2, 8.7, "src/core/peptide_core_utils.R", adj = 0, cex = 0.7)
  text(7.2, 8.4, "src/modules/peptide_integration.R", adj = 0, cex = 0.7)
  text(7.2, 8.1, "src/modules/peptide_search", adj = 0, cex = 0.7)
  
  # Draw legend in center top
  text(4.5, 9, "Parameter Flow:", adj = 0.5, font = 2, cex = 0.8)
  points(4, 8.5, pch = 19, col = "red", cex = 1)
  text(4.3, 8.5, "Config params", adj = 0, cex = 0.7)
  points(4, 8, pch = 19, col = "blue", cex = 1)
  text(4.3, 8, "Shared utilities", adj = 0, cex = 0.7)
  
  # Draw pipeline scripts
  rect(1.5, 6, 3.5, 7, col = "#F08080", border = "black")
  text(2.5, 6.5, "run_peptide_pipeline.R", cex = 0.7)
  
  rect(6.5, 6, 8.5, 7, col = "#F08080", border = "black")
  text(7.5, 6.5, "run_peptide_transcriptome_pipeline.R", cex = 0.6)
  
  # Draw primary results directories
  rect(1, 4, 3, 5, col = "#90EE90", border = "black")
  text(2, 4.5, "results/peptide_analysis", cex = 0.7)
  
  rect(4, 4, 6, 5, col = "#90EE90", border = "black")
  text(5, 4.5, "results/transcriptome_analysis", cex = 0.7)
  
  rect(7, 4, 9, 5, col = "#90EE90", border = "black")
  text(8, 4.5, "results/tumor_normal_transcriptome", cex = 0.6)
  
  # Secondary results directories
  rect(1, 2, 3, 3, col = "#90EE90", border = "black")
  text(2, 2.5, "results/148tumor_normal", cex = 0.7)
  
  rect(4, 2, 6, 3, col = "#90EE90", border = "black")
  text(5, 2.5, "results/all_sample_analysis", cex = 0.7)
  
  rect(7, 2, 9, 3, col = "#90EE90", border = "black")
  text(8, 2.5, "results/spike_in_analysis", cex = 0.7)
  
  # Benefits text - small and compact
  text(4, 1.5, "Benefits: Reproducible analyses, Comparable results", cex = 0.6, col = "darkgreen")
  
  # Integration box
  rect(4, 0.5, 6, 1.5, col = "#B19CD9", border = "black")
  text(5, 1, "Integrated Analysis", font = 2, cex = 0.7)
  
  # Draw parameter flow arrows
  # From config to pipelines
  arrows(1.75, 8, 2.5, 7, length = 0.1, lwd = 1.5, col = "red")
  arrows(1.75, 8, 7.5, 7, length = 0.1, lwd = 1.5, col = "red")
  
  # From source to pipelines
  arrows(8.25, 8, 2.5, 7, length = 0.1, lwd = 1.5, col = "blue")
  arrows(8.25, 8, 7.5, 7, length = 0.1, lwd = 1.5, col = "blue")
  
  # From pipelines to primary results
  arrows(2.5, 6, 2, 5, length = 0.1, lwd = 1)
  arrows(2.5, 6, 5, 5, length = 0.1, lwd = 1)
  arrows(7.5, 6, 5, 5, length = 0.1, lwd = 1)
  arrows(7.5, 6, 8, 5, length = 0.1, lwd = 1)
  
  # From primary results to secondary results
  arrows(2, 4, 2, 3, length = 0.1, lwd = 1)
  arrows(5, 4, 5, 3, length = 0.1, lwd = 1)
  arrows(8, 4, 8, 3, length = 0.1, lwd = 1)
  
  # From secondary results to final integration
  arrows(2, 2, 4.5, 1.5, length = 0.1, lwd = 1)
  arrows(5, 2, 5, 1.5, length = 0.1, lwd = 1)
  arrows(8, 2, 5.5, 1.5, length = 0.1, lwd = 1)
  
  # Add main title
  text(5, 9.8, "HLA-I Analysis Pipeline - Parameter Consistency Flow", font = 2, cex = 1)
  
  dev.off()
  
  message("Safe compact version created: hla_pipeline_parameters_safe.png")
  return(TRUE)
}

# Run both diagram functions
create_final_diagram()
create_safe_diagram()

# Display a message in the RStudio console
cat("\nTwo versions of the HLA-I Analysis Pipeline visualization created:\n")
cat("1. hla_pipeline_parameters_compact.png - More detailed version\n")
cat("2. hla_pipeline_parameters_safe.png - Extra-compact version to ensure no cutoffs\n")
cat("Please use whichever version displays better in your environment.\n")# HLA-I Analysis Pipeline Structure Visualization (Final Version)
# This script creates a visualization with completely fixed spacing and legends

# Install required packages if not already installed
if(!require(igraph)) install.packages("igraph")

# Load required libraries
library(igraph)

# Create a parameter consistency diagram with completely fixed spacing
create_final_diagram <- function() {
  # Create a plot with even larger dimensions to fix cut-off issues
  png("hla_pipeline_parameters_final.png", width = 1800, height = 1100, res = 120)
  
  # Set up empty plot with wider dimensions and more margin space
  par(mar = c(2, 2, 3, 10))  # Increase right margin for legends
  plot(0, type = "n", xlim = c(0, 16), ylim = c(0, 10), 
       xlab = "", ylab = "", main = "",
       axes = FALSE)
  
  # Draw config box - with proper spacing
  rect(1, 7.5, 5, 9, col = "#FFFFE0", border = "black")
  text(3, 8.5, "Configuration Files", font = 2, cex = 1.1)
  text(1.5, 8.2, "default_config.R", adj = 0, cex = 0.9)
  text(1.5, 7.9, "transcriptome_config.R", adj = 0, cex = 0.9)
  text(1.5, 7.6, "tumor_normal_only.R", adj = 0, cex = 0.9)
  
  # Draw source box - with proper spacing
  rect(10, 7.5, 14, 9, col = "#ADD8E6", border = "black")
  text(12, 8.5, "Source Code", font = 2, cex = 1.1)
  text(10.5, 8.2, "src/core/peptide_core_utils.R", adj = 0, cex = 0.9)
  text(10.5, 7.9, "src/modules/peptide_integration.R", adj = 0, cex = 0.9)
  text(10.5, 7.6, "src/modules/peptide_search", adj = 0, cex = 0.9)
  
  # Draw pipeline scripts with proper spacing
  rect(2, 5.5, 6, 6.5, col = "#F08080", border = "black")
  text(4, 6, "run_peptide_pipeline.R", cex = 1)
  
  rect(9, 5.5, 13, 6.5, col = "#F08080", border = "black")
  text(11, 6, "run_peptide_transcriptome_pipeline.R", cex = 1)
  
  # Draw primary results directories with proper spacing
  rect(1, 3.5, 4, 4.5, col = "#90EE90", border = "black")
  text(2.5, 4, "results/peptide_analysis", cex = 0.9)
  
  rect(6, 3.5, 9, 4.5, col = "#90EE90", border = "black")
  text(7.5, 4, "results/transcriptome_analysis", cex = 0.9)
  
  rect(11, 3.5, 14, 4.5, col = "#90EE90", border = "black")
  text(12.5, 4, "results/tumor_normal_transcriptome", cex = 0.9)
  
  # Secondary results directories with proper spacing
  rect(1, 1.5, 4, 2.5, col = "#90EE90", border = "black")
  text(2.5, 2, "results/148tumor_normal", cex = 0.9)
  
  rect(6, 1.5, 9, 2.5, col = "#90EE90", border = "black")
  text(7.5, 2, "results/all_sample_analysis", cex = 0.9)
  
  rect(11, 1.5, 14, 2.5, col = "#90EE90", border = "black")
  text(12.5, 2, "results/spike_in_analysis", cex = 0.9)
  
  # Integration box
  rect(6, 0, 9, 1, col = "#B19CD9", border = "black")
  text(7.5, 0.5, "Integrated Analysis", font = 2, cex = 1)
  
  # Draw parameter flow arrows with proper spacing
  # From config to pipelines
  arrows(3, 7.5, 4, 6.5, length = 0.1, lwd = 2, col = "red")
  arrows(3, 7.5, 11, 6.5, length = 0.1, lwd = 2, col = "red")
  
  # From source to pipelines
  arrows(12, 7.5, 4, 6.5, length = 0.1, lwd = 2, col = "blue")
  arrows(12, 7.5, 11, 6.5, length = 0.1, lwd = 2, col = "blue")
  
  # From pipelines to primary results
  arrows(4, 5.5, 2.5, 4.5, length = 0.1, lwd = 1.5)
  arrows(4, 5.5, 7.5, 4.5, length = 0.1, lwd = 1.5)
  arrows(11, 5.5, 7.5, 4.5, length = 0.1, lwd = 1.5)
  arrows(11, 5.5, 12.5, 4.5, length = 0.1, lwd = 1.5)
  
  # From primary results to secondary results
  arrows(2.5, 3.5, 2.5, 2.5, length = 0.1, lwd = 1.5)
  arrows(7.5, 3.5, 7.5, 2.5, length = 0.1, lwd = 1.5)
  arrows(12.5, 3.5, 12.5, 2.5, length = 0.1, lwd = 1.5)
  
  # From secondary results to final integration
  arrows(2.5, 1.5, 6.5, 1, length = 0.1, lwd = 1.5)
  arrows(7.5, 1.5, 7.5, 1, length = 0.1, lwd = 1.5)
  arrows(12.5, 1.5, 8.5, 1, length = 0.1, lwd = 1.5)
  
  # Add clean parameter consistency legend
  text(14.7, 8.5, "Parameter Consistency Flow:", adj = 0, font = 2, cex = 1)
  
  # Color dots for legend
  points(14.7, 8, pch = 19, col = "red", cex = 1.2)
  text(15.2, 8, "Configuration parameters", adj = 0, cex = 0.9)
  
  points(14.7, 7.5, pch = 19, col = "blue", cex = 1.2)
  text(15.2, 7.5, "Shared code utilities", adj = 0, cex = 0.9)
  
  # Example parameters - properly placed
  text(14.7, 6.5, "Example Parameters:", adj = 0, font = 2, cex = 1)
  text(15, 6, "- HLA typing parameters", adj = 0, cex = 0.9)
  text(15, 5.5, "- Peptide search thresholds", adj = 0, cex = 0.9)
  text(15, 5, "- File paths across analyses", adj = 0, cex = 0.9)
  text(15, 4.5, "- Common statistical methods", adj = 0, cex = 0.9)
  
  # Add benefits of consistency - properly placed
  text(14.7, 3.5, "Benefits of Parameter Consistency:", adj = 0, font = 2, cex = 1, col = "darkgreen")
  text(15, 3, "- Reproducible analyses across pipelines", adj = 0, cex = 0.9, col = "darkgreen")
  text(15, 2.5, "- Comparable results between data types", adj = 0, cex = 0.9, col = "darkgreen")
  
  # Add main title with better positioning
  text(8, 9.5, "HLA-I Analysis Pipeline - Parameter Consistency Flow", font = 2, cex = 1.5)
  
  dev.off()
  
  message("Final parameter consistency diagram created with all issues fixed: hla_pipeline_parameters_final.png")
  return(TRUE)
}

# Run the final diagram function
create_final_diagram()

# Display a message in the RStudio console
cat("\nHLA-I Analysis Pipeline visualization completed with all issues fixed.\n")
cat("- All text and legends are now fully visible\n")
cat("- Proper spacing between all components\n")
cat("- Consistent labeling throughout\n")
cat("The visualization file is in your working directory: hla_pipeline_parameters_final.png\n")
# # HLA-I Analysis Pipeline Structure Visualization (Final Version)
# # This script creates a visualization with improved spacing to avoid text overlap
# 
# # Install required packages if not already installed
# if(!require(igraph)) install.packages("igraph")
# 
# # Load required libraries
# library(igraph)
# 
# # Create a parameter consistency diagram with greatly improved spacing
# create_parameter_consistency_diagram <- function() {
#   # Create a simple plot showing how parameters flow through the pipeline
#   png("hla_pipeline_parameters.png", width = 1600, height = 1000, res = 120)
#   
#   # Set up empty plot with wider dimensions
#   plot(0, type = "n", xlim = c(0, 16), ylim = c(0, 12), 
#        xlab = "", ylab = "", main = "",
#        axes = FALSE)
#   
#   # Draw config box - moved left
#   rect(1, 9, 5, 11, col = "#FFFFE0", border = "black")
#   text(3, 10.5, "Configuration Files", font = 2, cex = 1.1)
#   text(1.5, 10, "default_config.R", adj = 0, cex = 0.9)
#   text(1.5, 9.5, "transcriptome_config.R", adj = 0, cex = 0.9)
#   text(1.5, 9, "tumor_normal_only.R", adj = 0, cex = 0.9)
#   
#   # Draw source box - moved right
#   rect(10, 9, 14, 11, col = "#ADD8E6", border = "black")
#   text(12, 10.5, "Source Code", font = 2, cex = 1.1)
#   text(10.5, 10, "src/core/peptide_core_utils.R", adj = 0, cex = 0.9)
#   text(10.5, 9.5, "src/modules/peptide_integration.R", adj = 0, cex = 0.9)
#   text(10.5, 9, "src/modules/peptide_search", adj = 0, cex = 0.9)
#   
#   # Draw pipeline scripts (much more spaced out)
#   rect(2, 6, 6, 7.5, col = "#F08080", border = "black")
#   text(4, 6.75, "run_peptide_pipeline.R", cex = 1)
#   
#   rect(9, 6, 13, 7.5, col = "#F08080", border = "black")
#   text(11, 6.75, "run_peptide_transcriptome_pipeline.R", cex = 1)
#   
#   # Draw primary results directories with more spacing
#   rect(1, 3.5, 4, 4.5, col = "#90EE90", border = "black")
#   text(2.5, 4, "results/peptide_analysis", cex = 0.9)
#   
#   rect(6, 3.5, 9, 4.5, col = "#90EE90", border = "black")
#   text(7.5, 4, "results/transcriptome_analysis", cex = 0.9)
#   
#   rect(11, 3.5, 14, 4.5, col = "#90EE90", border = "black")
#   text(12.5, 4, "results/tumor_normal_transcriptome", cex = 0.9)
#   
#   # Secondary results directories with better spacing
#   rect(1, 1.5, 4, 2.5, col = "#90EE90", border = "black")
#   text(2.5, 2, "results/148tumor_normal", cex = 0.9)
#   
#   rect(6, 1.5, 9, 2.5, col = "#90EE90", border = "black")
#   text(7.5, 2, "results/all_sample_analysis", cex = 0.9)
#   
#   rect(11, 1.5, 14, 2.5, col = "#90EE90", border = "black")
#   text(12.5, 2, "results/spike_in_analysis", cex = 0.9)
#   
#   # Integration box
#   rect(6, 0, 9, 1, col = "#B19CD9", border = "black")
#   text(7.5, 0.5, "Integrated Analysis", font = 2, cex = 1)
#   
#   # Draw parameter flow arrows with better spacing
#   # From config to pipelines
#   arrows(3, 9, 4, 7.5, length = 0.1, lwd = 2, col = "red")
#   arrows(3, 9, 11, 7.5, length = 0.1, lwd = 2, col = "red")
#   
#   # From source to pipelines
#   arrows(12, 9, 4, 7.5, length = 0.1, lwd = 2, col = "blue")
#   arrows(12, 9, 11, 7.5, length = 0.1, lwd = 2, col = "blue")
#   
#   # From pipelines to primary results
#   arrows(4, 6, 2.5, 4.5, length = 0.1, lwd = 1.5)
#   arrows(4, 6, 7.5, 4.5, length = 0.1, lwd = 1.5)
#   arrows(11, 6, 7.5, 4.5, length = 0.1, lwd = 1.5)
#   arrows(11, 6, 12.5, 4.5, length = 0.1, lwd = 1.5)
#   
#   # From primary results to secondary results
#   arrows(2.5, 3.5, 2.5, 2.5, length = 0.1, lwd = 1.5)
#   arrows(7.5, 3.5, 7.5, 2.5, length = 0.1, lwd = 1.5)
#   arrows(12.5, 3.5, 12.5, 2.5, length = 0.1, lwd = 1.5)
#   
#   # From secondary results to final integration
#   arrows(2.5, 1.5, 6.5, 1, length = 0.1, lwd = 1.5)
#   arrows(7.5, 1.5, 7.5, 1, length = 0.1, lwd = 1.5)
#   arrows(12.5, 1.5, 8.5, 1, length = 0.1, lwd = 1.5)
#   
#   # Add legend with ample spacing
#   # Parameter consistency key - top right with more spacing
#   text(15, 10.5, "Parameter Consistency Flow:", adj = 0, font = 2, cex = 1)
#   points(15.2, 9.8, pch = 19, col = "red", cex = 1.2)
#   text(15.7, 9.8, "Configuration parameters", adj = 0, cex = 0.9)
#   points(15.2, 9.2, pch = 19, col = "blue", cex = 1.2)
#   text(15.7, 9.2, "Shared code utilities", adj = 0, cex = 0.9)
#   
#   # Example parameters - below the legend with good spacing
#   text(15, 8, "Example Parameters:", adj = 0, font = 2, cex = 1)
#   text(15.3, 7.4, "- HLA typing parameters", adj = 0, cex = 0.9)
#   text(15.3, 6.8, "- Peptide search thresholds", adj = 0, cex = 0.9)
#   text(15.3, 6.2, "- File paths across analyses", adj = 0, cex = 0.9)
#   text(15.3, 5.6, "- Common statistical methods", adj = 0, cex = 0.9)
#   
#   # Add benefits of consistency - moved to more visible location
#   rect(14, 4, 16, 5, border = "darkgreen", lty = 2, lwd = 2)
#   text(15, 4.8, "Benefits of Parameter Consistency:", adj = 0.5, font = 2, cex = 0.9, col = "darkgreen")
#   text(14.2, 4.4, "- Reproducible analyses across pipelines", adj = 0, cex = 0.8, col = "darkgreen")
#   text(14.2, 4.1, "- Comparable results between data types", adj = 0, cex = 0.8, col = "darkgreen")
#   
#   # Add main title with better positioning
#   text(8, 11.5, "HLA-I Analysis Pipeline - Parameter Consistency Flow", font = 2, cex = 1.5)
#   
#   dev.off()
#   
#   message("Parameter consistency diagram created with improved spacing: hla_pipeline_parameters.png")
#   return(TRUE)
# }
# 
# # Run the diagram function
# create_parameter_consistency_diagram()
# 
# # Display a message in the RStudio console
# cat("\nHLA-I Analysis Pipeline visualization completed with improved spacing.\n")
# cat("All text should now be clearly readable without overlap.\n")
# cat("The visualization file is in your working directory: hla_pipeline_parameters.png\n")
# 
# # # HLA-I Analysis Pipeline Structure Visualization (Enhanced Version)
# # # This script creates a detailed visual representation of the HLA-I Analysis pipeline structure,
# # # showing files within src and config directories to highlight parameter consistency.
# # 
# # # Install required packages if not already installed
# # if(!require(igraph)) install.packages("igraph")
# # 
# # # Load required libraries
# # library(igraph)
# # 
# # # Create a parameter consistency diagram with improved spacing and additional results directories
# # create_parameter_consistency_diagram <- function() {
# #   # Create a simple plot showing how parameters flow through the pipeline
# #   png("hla_pipeline_parameters.png", width = 1200, height = 900, res = 120)
# #   
# #   # Set up empty plot
# #   plot(0, type = "n", xlim = c(0, 12), ylim = c(0, 11), 
# #        xlab = "", ylab = "", main = "HLA-I Analysis Pipeline - Parameter Consistency",
# #        axes = FALSE)
# #   
# #   # Draw config box
# #   rect(1, 8.5, 4, 10, col = "#FFFFE0", border = "black")
# #   text(2.5, 9.7, "Configuration Files", font = 2)
# #   text(1.2, 9.4, "default_config.R", adj = 0, cex = 0.8)
# #   text(1.2, 9.1, "transcriptome_config.R", adj = 0, cex = 0.8)
# #   text(1.2, 8.8, "tumor_normal_only.R", adj = 0, cex = 0.8)
# #   
# #   # Draw source box
# #   rect(6, 8.5, 9, 10, col = "#ADD8E6", border = "black")
# #   text(7.5, 9.7, "Source Code", font = 2)
# #   text(6.2, 9.4, "src/core/peptide_core_utils.R", adj = 0, cex = 0.8)
# #   text(6.2, 9.1, "src/modules/peptide_integration.R", adj = 0, cex = 0.8)
# #   text(6.2, 8.8, "src/modules/peptide_search", adj = 0, cex = 0.8)
# #   
# #   # Draw pipeline scripts (more centered)
# #   rect(3, 6.5, 5, 7.5, col = "#F08080", border = "black")
# #   text(4, 7, "run_peptide_pipeline.R", cex = 0.9)
# #   
# #   rect(7, 6.5, 9, 7.5, col = "#F08080", border = "black")
# #   text(8, 7, "run_peptide_transcriptome_pipeline.R", cex = 0.8)
# #   
# #   # Draw ALL results directories from your listing
# #   # Main results directories
# #   rect(2.5, 4, 4.5, 5, col = "#90EE90", border = "black")
# #   text(3.5, 4.5, "results/peptide_analysis", cex = 0.8)
# #   
# #   rect(5, 4, 7, 5, col = "#90EE90", border = "black")
# #   text(6, 4.5, "results/transcriptome_analysis", cex = 0.8)
# #   
# #   rect(7.5, 4, 9.5, 5, col = "#90EE90", border = "black")
# #   text(8.5, 4.5, "results/tumor_normal_transcriptome", cex = 0.7)
# #   
# #   # Additional results directories you mentioned
# #   rect(1, 2.5, 3, 3.5, col = "#90EE90", border = "black")
# #   text(2, 3, "results/148tumor_normal", cex = 0.8)
# #   
# #   rect(3.5, 2.5, 5.5, 3.5, col = "#90EE90", border = "black")
# #   text(4.5, 3, "results/all_sample_analysis", cex = 0.8)
# #   
# #   rect(6, 2.5, 8, 3.5, col = "#90EE90", border = "black")
# #   text(7, 3, "results/spike_in_analysis", cex = 0.8)
# #   
# #   # Integration box
# #   rect(4, 1, 6, 2, col = "#B19CD9", border = "black")
# #   text(5, 1.5, "Integrated Analysis", font = 2, cex = 0.9)
# #   
# #   # Draw parameter flow arrows with better spacing
# #   # From config to pipelines
# #   arrows(2.5, 8.5, 4, 7.5, length = 0.1, lwd = 1.5, col = "red")
# #   arrows(2.5, 8.5, 8, 7.5, length = 0.1, lwd = 1.5, col = "red")
# #   
# #   # From source to pipelines
# #   arrows(7.5, 8.5, 4, 7.5, length = 0.1, lwd = 1.5, col = "blue")
# #   arrows(7.5, 8.5, 8, 7.5, length = 0.1, lwd = 1.5, col = "blue")
# #   
# #   # From pipelines to primary results
# #   arrows(4, 6.5, 3.5, 5, length = 0.1, lwd = 1.5)
# #   arrows(4, 6.5, 6, 5, length = 0.1, lwd = 1.5)
# #   arrows(8, 6.5, 6, 5, length = 0.1, lwd = 1.5)
# #   arrows(8, 6.5, 8.5, 5, length = 0.1, lwd = 1.5)
# #   
# #   # From primary results to secondary results
# #   arrows(3.5, 4, 2, 3.5, length = 0.1, lwd = 1.5)
# #   arrows(6, 4, 4.5, 3.5, length = 0.1, lwd = 1.5)
# #   arrows(8.5, 4, 7, 3.5, length = 0.1, lwd = 1.5)
# #   
# #   # From results to final integration
# #   arrows(2, 2.5, 4.5, 2, length = 0.1, lwd = 1.5)
# #   arrows(4.5, 2.5, 5, 2, length = 0.1, lwd = 1.5)
# #   arrows(7, 2.5, 5.5, 2, length = 0.1, lwd = 1.5)
# #   
# #   # Add explanation - moved to a clearer location
# #   text(10, 9, "Parameter Consistency Flow:", adj = 0, font = 2)
# #   points(10.2, 8.5, pch = 19, col = "red")
# #   text(10.5, 8.5, "Configuration parameters", adj = 0, cex = 0.8)
# #   points(10.2, 8, pch = 19, col = "blue")
# #   text(10.5, 8, "Shared code utilities", adj = 0, cex = 0.8)
# #   
# #   # Add key parameter examples - moved for better spacing
# #   text(10, 7, "Example Parameters:", adj = 0, font = 2)
# #   text(10, 6.6, "- HLA typing parameters", adj = 0, cex = 0.8)
# #   text(10, 6.2, "- Peptide search thresholds", adj = 0, cex = 0.8)
# #   text(10, 5.8, "- File paths across analyses", adj = 0, cex = 0.8)
# #   text(10, 5.4, "- Common statistical methods", adj = 0, cex = 0.8)
# #   
# #   # Add benefits of consistency
# #   rect(9, 4, 11.5, 5, border = "darkgreen", lty = 2)
# #   text(10.25, 4.8, "Benefits of Parameter Consistency:", adj = 0.5, font = 2, col = "darkgreen")
# #   text(9.2, 4.5, "- Reproducible analyses across pipelines", adj = 0, cex = 0.8, col = "darkgreen")
# #   text(9.2, 4.2, "- Comparable results between data types", adj = 0, cex = 0.8, col = "darkgreen")
# #   
# #   # Add title showing purpose of diagram
# #   title(main = "HLA-I Analysis Pipeline - Parameter Consistency Flow", cex.main = 1.3)
# #   
# #   dev.off()
# #   
# #   message("Parameter consistency diagram created: hla_pipeline_parameters.png")
# #   return(TRUE)
# # }
# # 
# # # Run the diagram function
# # create_parameter_consistency_diagram()
# # 
# # # Display a message in the RStudio console
# # cat("\nHLA-I Analysis Pipeline Structure has been visualized.\n")
# # cat("This visualization highlights parameter consistency across pipeline components.\n")
# # cat("The diagram includes all results directories and improved spacing.\n")
# # cat("The visualization file is in your working directory: hla_pipeline_parameters.png\n")
# # 
# # # # HLA-I Analysis Pipeline Structure Visualization (Enhanced Version)
# # # # This script creates a detailed visual representation of the HLA-I Analysis pipeline structure,
# # # # showing files within src and config directories to highlight parameter consistency.
# # # 
# # # # Install required packages if not already installed
# # # if(!require(igraph)) install.packages("igraph")
# # # if(!require(ggplot2)) install.packages("ggplot2")
# # # if(!require(dplyr)) install.packages("dplyr")
# # # 
# # # # Load required libraries
# # # library(igraph)
# # # library(ggplot2)
# # # 
# # # # Create a detailed diagram showing the pipeline components and files
# # # create_detailed_diagram <- function() {
# # #   # Define the directory structure with files
# # #   structure <- list(
# # #     # Source code directory
# # #     "src" = list(
# # #       "core" = c("peptide_core_utils.R", "peptide_data_processing.R"),
# # #       "modules" = c("peptide_integration.R", "peptide_search"),
# # #       "reports" = c("visualization_utils.R")
# # #     ),
# # #     
# # #     # Configuration files
# # #     "config" = c(
# # #       "default_config.R", 
# # #       "transcriptome_config.R", 
# # #       "tumor_normal_only.R", 
# # #       "tumor_normal_transcriptome_config.R",
# # #       "spike_in_analysis.R"
# # #     ),
# # #     
# # #     # Main pipeline scripts
# # #     "pipeline_scripts" = c(
# # #       "run_peptide_pipeline.R",
# # #       "run_peptide_transcriptome_pipeline.R"
# # #     ),
# # #     
# # #     # Results directories
# # #     "results" = list(
# # #       "peptide_analysis" = c("peptide_data.RData", "peptide_results.csv"),
# # #       "transcriptome_analysis" = c("transcript_data.RData", "expression.csv"),
# # #       "tumor_normal_transcriptome" = c("differential_expression.csv", "fusion_candidates.csv")
# # #     ),
# # #     
# # #     # Data directories
# # #     "data" = c("processed_data.RData"),
# # #     "rawdata" = c("sample_files.fastq")
# # #   )
# # #   
# # #   # Create a dataframe for all nodes
# # #   nodes <- data.frame(
# # #     name = character(),
# # #     type = character(),
# # #     parent = character(),
# # #     stringsAsFactors = FALSE
# # #   )
# # #   
# # #   # Create a dataframe for all edges
# # #   edges <- data.frame(
# # #     from = character(),
# # #     to = character(),
# # #     type = character(),
# # #     stringsAsFactors = FALSE
# # #   )
# # #   
# # #   # Add top-level directories as nodes
# # #   for (dir_name in names(structure)) {
# # #     if (dir_name != "pipeline_scripts") {
# # #       nodes <- rbind(nodes, data.frame(
# # #         name = dir_name,
# # #         type = ifelse(dir_name == "config", "config", 
# # #                       ifelse(dir_name == "src", "source",
# # #                              ifelse(dir_name %in% c("data", "rawdata"), "data", "results"))),
# # #         parent = NA,
# # #         stringsAsFactors = FALSE
# # #       ))
# # #     }
# # #   }
# # #   
# # #   # Add subdirectories and files as nodes
# # #   for (dir_name in names(structure)) {
# # #     dir_content <- structure[[dir_name]]
# # #     
# # #     if (dir_name == "pipeline_scripts") {
# # #       # Add pipeline scripts as top-level nodes
# # #       for (script in dir_content) {
# # #         nodes <- rbind(nodes, data.frame(
# # #           name = script,
# # #           type = "pipeline",
# # #           parent = NA,
# # #           stringsAsFactors = FALSE
# # #         ))
# # #       }
# # #     } else if (is.list(dir_content)) {
# # #       # For directories with subdirectories
# # #       for (subdir_name in names(dir_content)) {
# # #         # Add subdirectory
# # #         full_subdir_name <- paste0(dir_name, "/", subdir_name)
# # #         nodes <- rbind(nodes, data.frame(
# # #           name = full_subdir_name,
# # #           type = ifelse(dir_name == "src", "source", 
# # #                         ifelse(dir_name == "config", "config", 
# # #                                ifelse(dir_name %in% c("data", "rawdata"), "data", "results"))),
# # #           parent = dir_name,
# # #           stringsAsFactors = FALSE
# # #         ))
# # #         
# # #         # Add files in subdirectory
# # #         for (file in dir_content[[subdir_name]]) {
# # #           full_file_name <- paste0(full_subdir_name, "/", file)
# # #           nodes <- rbind(nodes, data.frame(
# # #             name = full_file_name,
# # #             type = ifelse(dir_name == "src", "source_file", 
# # #                           ifelse(dir_name == "config", "config_file", 
# # #                                  ifelse(dir_name %in% c("data", "rawdata"), "data_file", "results_file"))),
# # #             parent = full_subdir_name,
# # #             stringsAsFactors = FALSE
# # #           ))
# # #         }
# # #       }
# # #     } else {
# # #       # For directories with just files
# # #       for (file in dir_content) {
# # #         full_file_name <- paste0(dir_name, "/", file)
# # #         nodes <- rbind(nodes, data.frame(
# # #           name = full_file_name,
# # #           type = ifelse(dir_name == "src", "source_file", 
# # #                         ifelse(dir_name == "config", "config_file", 
# # #                                ifelse(dir_name %in% c("data", "rawdata"), "data_file", "results_file"))),
# # #           parent = dir_name,
# # #           stringsAsFactors = FALSE
# # #         ))
# # #       }
# # #     }
# # #   }
# # #   
# # #   # Create parent-child edges
# # #   for (i in 1:nrow(nodes)) {
# # #     if (!is.na(nodes$parent[i])) {
# # #       edges <- rbind(edges, data.frame(
# # #         from = nodes$parent[i],
# # #         to = nodes$name[i],
# # #         type = "contains",
# # #         stringsAsFactors = FALSE
# # #       ))
# # #     }
# # #   }
# # #   
# # #   # Add pipeline data flow edges
# # #   pipeline_edges <- data.frame(
# # #     from = c(
# # #       "rawdata", "data", "data",
# # #       "config/default_config.R", "config/transcriptome_config.R", 
# # #       "config/tumor_normal_only.R", "config/tumor_normal_transcriptome_config.R",
# # #       "src/core/peptide_core_utils.R", "src/core/peptide_data_processing.R",
# # #       "src/modules/peptide_integration.R", "src/modules/peptide_search",
# # #       "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #       "results/peptide_analysis", "results/transcriptome_analysis"
# # #     ),
# # #     to = c(
# # #       "data", "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #       "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #       "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #       "run_peptide_pipeline.R", "run_peptide_pipeline.R",
# # #       "run_peptide_pipeline.R", "run_peptide_pipeline.R",
# # #       "results/peptide_analysis", "results/transcriptome_analysis",
# # #       "results/tumor_normal_transcriptome", "results/tumor_normal_transcriptome"
# # #     ),
# # #     type = "flow",
# # #     stringsAsFactors = FALSE
# # #   )
# # #   
# # #   edges <- rbind(edges, pipeline_edges)
# # #   
# # #   # Create the graph
# # #   g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
# # #   
# # #   # Define colors for node types
# # #   colors <- c(
# # #     "pipeline" = "#F08080",        # lightcoral
# # #     "source" = "#ADD8E6",          # lightblue
# # #     "source_file" = "#87CEEB",     # skyblue
# # #     "config" = "#FFFFE0",          # lightyellow
# # #     "config_file" = "#FFFACD",     # lemonchiffon
# # #     "results" = "#90EE90",         # lightgreen
# # #     "results_file" = "#98FB98",    # palegreen
# # #     "data" = "#D3D3D3",            # lightgrey
# # #     "data_file" = "#E8E8E8"        # whitesmoke
# # #   )
# # #   
# # #   # Set node colors and shapes based on type
# # #   V(g)$color <- colors[V(g)$type]
# # #   V(g)$shape <- ifelse(grepl("file", V(g)$type), "rectangle", "circle")
# # #   
# # #   # Set edge types
# # #   E(g)$lty <- ifelse(E(g)$type == "contains", 2, 1)  # dashed for contains, solid for flow
# # #   E(g)$arrow.size <- ifelse(E(g)$type == "contains", 0.3, 0.5)
# # #   
# # #   # Create a layout that emphasizes flow
# # #   # This is a complex layout that tries to position files within their directories
# # #   # while maintaining the flow of the pipeline
# # #   
# # #   # Start with a hierarchical layout
# # #   layout_matrix <- NULL
# # #   
# # #   # Manually position nodes to show the flow clearly
# # #   layout_matrix <- matrix(0, nrow = vcount(g), ncol = 2)
# # #   rownames(layout_matrix) <- V(g)$name
# # #   
# # #   # Set positions for main directories
# # #   layout_matrix["src", ] <- c(2, 7)
# # #   layout_matrix["config", ] <- c(5, 7)
# # #   layout_matrix["data", ] <- c(8, 7)
# # #   layout_matrix["rawdata", ] <- c(8, 8)
# # #   layout_matrix["results", ] <- c(5, 1)
# # #   
# # #   # Set positions for src subdirectories
# # #   layout_matrix["src/core", ] <- c(1, 6)
# # #   layout_matrix["src/modules", ] <- c(2, 6)
# # #   layout_matrix["src/reports", ] <- c(3, 6)
# # #   
# # #   # Set positions for files in src/core
# # #   if ("src/core/peptide_core_utils.R" %in% V(g)$name)
# # #     layout_matrix["src/core/peptide_core_utils.R", ] <- c(0.5, 5.5)
# # #   if ("src/core/peptide_data_processing.R" %in% V(g)$name)
# # #     layout_matrix["src/core/peptide_data_processing.R", ] <- c(1.5, 5.5)
# # #   
# # #   # Set positions for files in src/modules
# # #   if ("src/modules/peptide_integration.R" %in% V(g)$name)
# # #     layout_matrix["src/modules/peptide_integration.R", ] <- c(1.7, 5.5)
# # #   if ("src/modules/peptide_search" %in% V(g)$name)
# # #     layout_matrix["src/modules/peptide_search", ] <- c(2.3, 5.5)
# # #   
# # #   # Set positions for files in src/reports
# # #   if ("src/reports/visualization_utils.R" %in% V(g)$name)
# # #     layout_matrix["src/reports/visualization_utils.R", ] <- c(3, 5.5)
# # #   
# # #   # Set positions for config files
# # #   for (i in seq_along(structure$config)) {
# # #     file_name <- paste0("config/", structure$config[i])
# # #     if (file_name %in% V(g)$name)
# # #       layout_matrix[file_name, ] <- c(5, 8 - 0.3*i)
# # #   }
# # #   
# # #   # Set positions for pipeline scripts
# # #   layout_matrix["run_peptide_pipeline.R", ] <- c(3, 4)
# # #   layout_matrix["run_peptide_transcriptome_pipeline.R", ] <- c(7, 4)
# # #   
# # #   # Set positions for results directories
# # #   layout_matrix["results/peptide_analysis", ] <- c(3, 2)
# # #   layout_matrix["results/transcriptome_analysis", ] <- c(7, 2)
# # #   layout_matrix["results/tumor_normal_transcriptome", ] <- c(5, 0)
# # #   
# # #   # Set positions for results files
# # #   if ("results/peptide_analysis/peptide_data.RData" %in% V(g)$name)
# # #     layout_matrix["results/peptide_analysis/peptide_data.RData", ] <- c(2.5, 1.5)
# # #   if ("results/peptide_analysis/peptide_results.csv" %in% V(g)$name)
# # #     layout_matrix["results/peptide_analysis/peptide_results.csv", ] <- c(3.5, 1.5)
# # #   
# # #   if ("results/transcriptome_analysis/transcript_data.RData" %in% V(g)$name)
# # #     layout_matrix["results/transcriptome_analysis/transcript_data.RData", ] <- c(6.5, 1.5)
# # #   if ("results/transcriptome_analysis/expression.csv" %in% V(g)$name)
# # #     layout_matrix["results/transcriptome_analysis/expression.csv", ] <- c(7.5, 1.5)
# # #   
# # #   if ("results/tumor_normal_transcriptome/differential_expression.csv" %in% V(g)$name)
# # #     layout_matrix["results/tumor_normal_transcriptome/differential_expression.csv", ] <- c(4.5, -0.5)
# # #   if ("results/tumor_normal_transcriptome/fusion_candidates.csv" %in% V(g)$name)
# # #     layout_matrix["results/tumor_normal_transcriptome/fusion_candidates.csv", ] <- c(5.5, -0.5)
# # #   
# # #   # Set positions for data and rawdata files
# # #   if ("data/processed_data.RData" %in% V(g)$name)
# # #     layout_matrix["data/processed_data.RData", ] <- c(8, 6.5)
# # #   if ("rawdata/sample_files.fastq" %in% V(g)$name)
# # #     layout_matrix["rawdata/sample_files.fastq", ] <- c(8, 7.5)
# # #   
# # #   # Create the plot with highlighting parameter consistency
# # #   pdf("hla_pipeline_detailed_structure.pdf", width = 12, height = 10)
# # #   
# # #   # Adjust margins to allow space for labels
# # #   par(mar = c(1, 1, 2, 1))
# # #   
# # #   # Plot the graph
# # #   plot(g, 
# # #        layout = layout_matrix,
# # #        vertex.size = ifelse(grepl("file", V(g)$type), 15, 25),
# # #        vertex.shape = V(g)$shape,
# # #        vertex.color = V(g)$color,
# # #        vertex.frame.color = "gray50",
# # #        vertex.label.cex = 0.6,
# # #        vertex.label.color = "black",
# # #        edge.arrow.size = E(g)$arrow.size,
# # #        edge.lty = E(g)$lty,
# # #        edge.color = ifelse(E(g)$type == "contains", "gray60", "black"),
# # #        main = "HLA-I Analysis Pipeline - Detailed Structure")
# # #   
# # #   # Add a legend
# # #   legend("topright", 
# # #          legend = c(
# # #            "Pipeline Scripts", 
# # #            "Source Directories", "Source Files",
# # #            "Config Directory", "Config Files",
# # #            "Results Directories", "Result Files",
# # #            "Data Directories", "Data Files",
# # #            "Contains Relationship", "Data Flow"
# # #          ), 
# # #          pch = c(
# # #            21, 
# # #            21, 22,
# # #            21, 22,
# # #            21, 22,
# # #            21, 22,
# # #            NA, NA
# # #          ),
# # #          lty = c(
# # #            NA, 
# # #            NA, NA,
# # #            NA, NA,
# # #            NA, NA,
# # #            NA, NA,
# # #            2, 1
# # #          ),
# # #          lwd = c(
# # #            NA, 
# # #            NA, NA,
# # #            NA, NA,
# # #            NA, NA,
# # #            NA, NA,
# # #            1, 1
# # #          ),
# # #          col = c(
# # #            colors["pipeline"], 
# # #            colors["source"], colors["source_file"],
# # #            colors["config"], colors["config_file"],
# # #            colors["results"], colors["results_file"],
# # #            colors["data"], colors["data_file"],
# # #            "gray60", "black"
# # #          ),
# # #          pt.bg = c(
# # #            colors["pipeline"], 
# # #            colors["source"], colors["source_file"],
# # #            colors["config"], colors["config_file"],
# # #            colors["results"], colors["results_file"],
# # #            colors["data"], colors["data_file"],
# # #            NA, NA
# # #          ),
# # #          pt.cex = c(
# # #            1.5, 
# # #            1.5, 1,
# # #            1.5, 1,
# # #            1.5, 1,
# # #            1.5, 1,
# # #            NA, NA
# # #          ),
# # #          title = "Component Types", 
# # #          cex = 0.7,
# # #          box.lty = 2)
# # #   
# # #   # Add annotation explaining parameter consistency
# # #   text(x = 2, y = -1.5, labels = "Parameter Consistency:", 
# # #        adj = 0, font = 2, cex = 0.8)
# # #   text(x = 2, y = -2, 
# # #        labels = "Config files ensure consistent parameters across pipeline components", 
# # #        adj = 0, cex = 0.7)
# # #   text(x = 2, y = -2.5, 
# # #        labels = "Core utilities are shared between pipelines to maintain consistency", 
# # #        adj = 0, cex = 0.7)
# # #   
# # #   dev.off()
# # #   
# # #   # Create a PNG version
# # #   png("hla_pipeline_detailed_structure.png", width = 1200, height = 1000, res = 120)
# # #   
# # #   # Adjust margins
# # #   par(mar = c(1, 1, 2, 1))
# # #   
# # #   # Plot the graph
# # #   plot(g, 
# # #        layout = layout_matrix,
# # #        vertex.size = ifelse(grepl("file", V(g)$type), 15, 25),
# # #        vertex.shape = V(g)$shape,
# # #        vertex.color = V(g)$color,
# # #        vertex.frame.color = "gray50",
# # #        vertex.label.cex = 0.6,
# # #        vertex.label.color = "black",
# # #        edge.arrow.size = E(g)$arrow.size,
# # #        edge.lty = E(g)$lty,
# # #        edge.color = ifelse(E(g)$type == "contains", "gray60", "black"),
# # #        main = "HLA-I Analysis Pipeline - Detailed Structure")
# # #   
# # #   # Add the same legend
# # #   legend("topright", 
# # #          legend = c(
# # #            "Pipeline Scripts", 
# # #            "Source Directories", "Source Files",
# # #            "Config Directory", "Config Files",
# # #            "Results Directories", "Result Files",
# # #            "Data Directories", "Data Files",
# # #            "Contains Relationship", "Data Flow"
# # #          ), 
# # #          pch = c(
# # #            21, 
# # #            21, 22,
# # #            21, 22,
# # #            21, 22,
# # #            21, 22,
# # #            NA, NA
# # #          ),
# # #          lty = c(
# # #            NA, 
# # #            NA, NA,
# # #            NA, NA,
# # #            NA, NA,
# # #            NA, NA,
# # #            2, 1
# # #          ),
# # #          lwd = c(
# # #            NA, 
# # #            NA, NA,
# # #            NA, NA,
# # #            NA, NA,
# # #            NA, NA,
# # #            1, 1
# # #          ),
# # #          col = c(
# # #            colors["pipeline"], 
# # #            colors["source"], colors["source_file"],
# # #            colors["config"], colors["config_file"],
# # #            colors["results"], colors["results_file"],
# # #            colors["data"], colors["data_file"],
# # #            "gray60", "black"
# # #          ),
# # #          pt.bg = c(
# # #            colors["pipeline"], 
# # #            colors["source"], colors["source_file"],
# # #            colors["config"], colors["config_file"],
# # #            colors["results"], colors["results_file"],
# # #            colors["data"], colors["data_file"],
# # #            NA, NA
# # #          ),
# # #          pt.cex = c(
# # #            1.5, 
# # #            1.5, 1,
# # #            1.5, 1,
# # #            1.5, 1,
# # #            1.5, 1,
# # #            NA, NA
# # #          ),
# # #          title = "Component Types", 
# # #          cex = 0.7,
# # #          box.lty = 2)
# # #   
# # #   # Add annotation explaining parameter consistency
# # #   text(x = 2, y = -1.5, labels = "Parameter Consistency:", 
# # #        adj = 0, font = 2, cex = 0.8)
# # #   text(x = 2, y = -2, 
# # #        labels = "Config files ensure consistent parameters across pipeline components", 
# # #        adj = 0, cex = 0.7)
# # #   text(x = 2, y = -2.5, 
# # #        labels = "Core utilities are shared between pipelines to maintain consistency", 
# # #        adj = 0, cex = 0.7)
# # #   
# # #   dev.off()
# # #   
# # #   message("Detailed pipeline visualization complete! Output saved as:")
# # #   message("- hla_pipeline_detailed_structure.pdf")
# # #   message("- hla_pipeline_detailed_structure.png")
# # #   
# # #   return(TRUE)
# # # }
# # # 
# # # # Alternative simpler version that highlights parameter consistency
# # # create_parameter_consistency_diagram <- function() {
# # #   # Create a simple plot showing how parameters flow through the pipeline
# # #   png("hla_pipeline_parameters.png", width = 1200, height = 800, res = 120)
# # #   
# # #   # Set up empty plot
# # #   plot(0, type = "n", xlim = c(0, 10), ylim = c(0, 10), 
# # #        xlab = "", ylab = "", main = "HLA-I Analysis Pipeline - Parameter Consistency",
# # #        axes = FALSE)
# # #   
# # #   # Draw config box
# # #   rect(1, 8, 4, 9.5, col = "#FFFFE0", border = "black")
# # #   text(2.5, 9.2, "Configuration Files", font = 2)
# # #   text(1.2, 8.8, "default_config.R", adj = 0, cex = 0.8)
# # #   text(1.2, 8.5, "transcriptome_config.R", adj = 0, cex = 0.8)
# # #   text(1.2, 8.2, "tumor_normal_only.R", adj = 0, cex = 0.8)
# # #   
# # #   # Draw source box
# # #   rect(1, 6, 4, 7.5, col = "#ADD8E6", border = "black")
# # #   text(2.5, 7.2, "Source Code", font = 2)
# # #   text(1.2, 6.8, "src/core/peptide_core_utils.R", adj = 0, cex = 0.8)
# # #   text(1.2, 6.5, "src/modules/peptide_integration.R", adj = 0, cex = 0.8)
# # #   text(1.2, 6.2, "src/modules/peptide_search", adj = 0, cex = 0.8)
# # #   
# # #   # Draw pipeline scripts
# # #   rect(3, 4, 5, 5, col = "#F08080", border = "black")
# # #   text(4, 4.5, "run_peptide_pipeline.R", cex = 0.9)
# # #   
# # #   rect(6, 4, 8, 5, col = "#F08080", border = "black")
# # #   text(7, 4.5, "run_peptide_transcriptome_pipeline.R", cex = 0.8)
# # #   
# # #   # Draw results
# # #   rect(3, 2, 5, 3, col = "#90EE90", border = "black")
# # #   text(4, 2.5, "results/peptide_analysis", cex = 0.8)
# # #   
# # #   rect(6, 2, 8, 3, col = "#90EE90", border = "black")
# # #   text(7, 2.5, "results/transcriptome_analysis", cex = 0.8)
# # #   
# # #   rect(4.5, 0.5, 6.5, 1.5, col = "#90EE90", border = "black")
# # #   text(5.5, 1, "results/tumor_normal_transcriptome", cex = 0.8)
# # #   
# # #   # Draw parameter flow arrows
# # #   # From config to pipelines
# # #   arrows(2.5, 8, 4, 5, length = 0.1, lwd = 1.5, col = "red")
# # #   arrows(2.5, 8, 7, 5, length = 0.1, lwd = 1.5, col = "red")
# # #   
# # #   # From source to pipelines
# # #   arrows(2.5, 6, 4, 5, length = 0.1, lwd = 1.5, col = "blue")
# # #   arrows(2.5, 6, 7, 5, length = 0.1, lwd = 1.5, col = "blue")
# # #   
# # #   # From pipelines to results
# # #   arrows(4, 4, 4, 3, length = 0.1, lwd = 1.5)
# # #   arrows(7, 4, 7, 3, length = 0.1, lwd = 1.5)
# # #   
# # #   # From results to final integration
# # #   arrows(4, 2, 5.5, 1.5, length = 0.1, lwd = 1.5)
# # #   arrows(7, 2, 5.5, 1.5, length = 0.1, lwd = 1.5)
# # #   
# # #   # Add explanation
# # #   text(8.5, 7.5, "Parameter Consistency Flow:", adj = 0, font = 2)
# # #   points(8.7, 7, pch = 19, col = "red")
# # #   text(9, 7, "Configuration parameters", adj = 0, cex = 0.8)
# # #   points(8.7, 6.5, pch = 19, col = "blue")
# # #   text(9, 6.5, "Shared code utilities", adj = 0, cex = 0.8)
# # #   
# # #   # Add key parameter examples
# # #   text(5.5, 8.5, "Example Parameters:", adj = 0, font = 2)
# # #   text(5.5, 8.1, "- HLA typing parameters", adj = 0, cex = 0.8)
# # #   text(5.5, 7.7, "- Peptide search thresholds", adj = 0, cex = 0.8)
# # #   text(5.5, 7.3, "- File paths across analyses", adj = 0, cex = 0.8)
# # #   text(5.5, 6.9, "- Common statistical methods", adj = 0, cex = 0.8)
# # #   
# # #   # Add benefits of consistency
# # #   rect(6.5, 5.5, 9.5, 6.5, border = "darkgreen", lty = 2)
# # #   text(8, 6.2, "Benefits of Parameter Consistency:", adj = 0.5, font = 2, col = "darkgreen")
# # #   text(6.7, 5.9, "- Reproducible analyses across pipelines", adj = 0, cex = 0.8, col = "darkgreen")
# # #   text(6.7, 5.6, "- Comparable results between data types", adj = 0, cex = 0.8, col = "darkgreen")
# # #   
# # #   dev.off()
# # #   
# # #   message("Parameter consistency diagram created: hla_pipeline_parameters.png")
# # #   return(TRUE)
# # # }
# # # 
# # # # Run both diagram functions
# # # create_detailed_diagram()
# # # create_parameter_consistency_diagram()
# # # 
# # # # Display a message in the RStudio console
# # # cat("\nHLA-I Analysis Pipeline Structure has been visualized with detailed file information.\n")
# # # cat("This visualization highlights parameter consistency across pipeline components.\n")
# # # cat("You can find the visualization files in your working directory:\n")
# # # cat("- hla_pipeline_detailed_structure.pdf/png: Detailed structure with all files\n")
# # # cat("- hla_pipeline_parameters.png: Focused view of parameter consistency\n")
# # 
# # # # HLA-I Analysis Pipeline Structure Visualization
# # # # This script creates a visual representation of the HLA-I Analysis pipeline structure,
# # # # focusing on run_peptide_pipeline.R and run_peptide_transcriptome_pipeline.R.
# # # 
# # # # Install required packages if not already installed
# # # if(!require(igraph)) install.packages("igraph")
# # # if(!require(ggplot2)) install.packages("ggplot2")
# # # if(!require(ggraph)) install.packages("ggraph")
# # # if(!require(dplyr)) install.packages("dplyr") # Added dplyr package
# # # if(!require(tidygraph)) install.packages("tidygraph") # Added tidygraph package
# # # 
# # # # Load required libraries
# # # library(igraph)
# # # library(ggplot2)
# # # library(ggraph)
# # # library(dplyr)      # Added dplyr library
# # # library(tidygraph)  # Added tidygraph library
# # # 
# # # # Create a simplified diagram showing the main pipeline components
# # # create_simplified_diagram <- function() {
# # #   # Create a directed graph for the pipeline
# # #   edges <- data.frame(
# # #     from = c(
# # #       "rawdata", "data", "data", 
# # #       "config", "config", 
# # #       "src/core", "src/core", 
# # #       "src/modules", "src/modules",
# # #       "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #       "results/peptide_analysis", "results/transcriptome_analysis"
# # #     ),
# # #     to = c(
# # #       "data", "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #       "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #       "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #       "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #       "results/peptide_analysis", "results/transcriptome_analysis",
# # #       "results/tumor_normal_transcriptome", "results/tumor_normal_transcriptome"
# # #     )
# # #   )
# # #   
# # #   # Create vertices (nodes) and define their types
# # #   vertices <- data.frame(
# # #     name = unique(c(as.character(edges$from), as.character(edges$to))),
# # #     stringsAsFactors = FALSE
# # #   )
# # #   
# # #   # Add node types
# # #   vertices$type <- "other"
# # #   vertices$type[grep("run_", vertices$name)] <- "pipeline"
# # #   vertices$type[grep("results", vertices$name)] <- "results"
# # #   vertices$type[grep("src", vertices$name)] <- "source"
# # #   vertices$type[vertices$name == "config"] <- "config"
# # #   vertices$type[vertices$name %in% c("data", "rawdata")] <- "data"
# # #   
# # #   # Create the graph
# # #   g <- graph_from_data_frame(edges, directed = TRUE, vertices = vertices)
# # #   
# # #   # Set layout coordinates to control visualization
# # #   layout_coords <- matrix(c(
# # #     2, 9,   # rawdata
# # #     2, 7,   # data
# # #     2, 5,   # config
# # #     1, 3,   # src/core
# # #     3, 3,   # src/modules
# # #     1, 1,   # run_peptide_pipeline.R
# # #     5, 1,   # run_peptide_transcriptome_pipeline.R
# # #     1, -1,  # results/peptide_analysis
# # #     5, -1,  # results/transcriptome_analysis
# # #     3, -3   # results/tumor_normal_transcriptome
# # #   ), ncol = 2, byrow = TRUE)
# # #   
# # #   rownames(layout_coords) <- vertices$name
# # #   
# # #   # Define colors for node types
# # #   colors <- c(
# # #     "pipeline" = "#F08080",  # lightcoral
# # #     "source" = "#ADD8E6",    # lightblue
# # #     "config" = "#FFFFE0",    # lightyellow
# # #     "results" = "#90EE90",   # lightgreen
# # #     "data" = "#D3D3D3",      # lightgrey
# # #     "other" = "#FFFFFF"      # white
# # #   )
# # #   
# # #   # Set node colors based on type
# # #   V(g)$color <- colors[vertices$type]
# # #   
# # #   # Create plot
# # #   pdf("hla_pipeline_structure.pdf", width = 10, height = 8)
# # #   
# # #   par(mar = c(1, 1, 3, 1))  # Adjust margins
# # #   plot(g, 
# # #        layout = layout_coords,
# # #        vertex.size = 30,
# # #        vertex.color = V(g)$color,
# # #        vertex.label.cex = 0.7,
# # #        edge.arrow.size = 0.5,
# # #        main = "HLA-I Analysis Pipeline Structure")
# # #   
# # #   # Add a legend
# # #   legend("topright", 
# # #          legend = names(colors), 
# # #          fill = colors, 
# # #          title = "Component Type", 
# # #          cex = 0.8)
# # #   
# # #   dev.off()
# # #   
# # #   # Also create a PNG version
# # #   png("hla_pipeline_structure.png", width = 1000, height = 800, res = 100)
# # #   
# # #   par(mar = c(1, 1, 3, 1))  # Adjust margins
# # #   plot(g, 
# # #        layout = layout_coords,
# # #        vertex.size = 30,
# # #        vertex.color = V(g)$color,
# # #        vertex.label.cex = 0.7,
# # #        edge.arrow.size = 0.5,
# # #        main = "HLA-I Analysis Pipeline Structure")
# # #   
# # #   # Add a legend
# # #   legend("topright", 
# # #          legend = names(colors), 
# # #          fill = colors, 
# # #          title = "Component Type", 
# # #          cex = 0.8)
# # #   
# # #   dev.off()
# # #   
# # #   # Print a message to let user know the basic plots are complete
# # #   message("Basic pipeline visualizations saved as:")
# # #   message("- hla_pipeline_structure.pdf")
# # #   message("- hla_pipeline_structure.png")
# # #   
# # #   # Return TRUE to indicate that basic visualization was completed successfully
# # #   return(TRUE)
# # # }
# # # 
# # # # Creating a simpler version without the ggraph portion, which was causing the error
# # # create_basic_visualization <- function() {
# # #   # Run the simplified diagram function first
# # #   success <- create_simplified_diagram()
# # #   
# # #   message("Pipeline visualization complete! Basic visualizations were created successfully.")
# # #   message("If you want to create the advanced ggraph visualization, please ensure the following packages are installed:")
# # #   message("install.packages(c('dplyr', 'tidygraph', 'ggraph'))")
# # #   
# # #   return(success)
# # # }
# # # 
# # # # Adding an optional function to try the ggraph visualization with error handling
# # # try_ggraph_visualization <- function() {
# # #   tryCatch({
# # #     # Load the required packages explicitly here
# # #     requireNamespace("dplyr", quietly = TRUE)
# # #     requireNamespace("tidygraph", quietly = TRUE)
# # #     requireNamespace("ggraph", quietly = TRUE)
# # #     
# # #     # Create a directed graph for the pipeline (repeat from above for encapsulation)
# # #     edges <- data.frame(
# # #       from = c(
# # #         "rawdata", "data", "data", 
# # #         "config", "config", 
# # #         "src/core", "src/core", 
# # #         "src/modules", "src/modules",
# # #         "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #         "results/peptide_analysis", "results/transcriptome_analysis"
# # #       ),
# # #       to = c(
# # #         "data", "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #         "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #         "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #         "run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R",
# # #         "results/peptide_analysis", "results/transcriptome_analysis",
# # #         "results/tumor_normal_transcriptome", "results/tumor_normal_transcriptome"
# # #       )
# # #     )
# # #     
# # #     # Create vertices and define types
# # #     vertices <- data.frame(
# # #       name = unique(c(as.character(edges$from), as.character(edges$to))),
# # #       stringsAsFactors = FALSE
# # #     )
# # #     
# # #     # Add node types
# # #     vertices$type <- "other"
# # #     vertices$type[grep("run_", vertices$name)] <- "pipeline"
# # #     vertices$type[grep("results", vertices$name)] <- "results"
# # #     vertices$type[grep("src", vertices$name)] <- "source"
# # #     vertices$type[vertices$name == "config"] <- "config"
# # #     vertices$type[vertices$name %in% c("data", "rawdata")] <- "data"
# # #     
# # #     # Create the graph
# # #     g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = vertices)
# # #     
# # #     # Set layout coordinates
# # #     layout_coords <- matrix(c(
# # #       2, 9,   # rawdata
# # #       2, 7,   # data
# # #       2, 5,   # config
# # #       1, 3,   # src/core
# # #       3, 3,   # src/modules
# # #       1, 1,   # run_peptide_pipeline.R
# # #       5, 1,   # run_peptide_transcriptome_pipeline.R
# # #       1, -1,  # results/peptide_analysis
# # #       5, -1,  # results/transcriptome_analysis
# # #       3, -3   # results/tumor_normal_transcriptome
# # #     ), ncol = 2, byrow = TRUE)
# # #     
# # #     rownames(layout_coords) <- vertices$name
# # #     
# # #     # Define colors
# # #     colors <- c(
# # #       "pipeline" = "#F08080",  # lightcoral
# # #       "source" = "#ADD8E6",    # lightblue
# # #       "config" = "#FFFFE0",    # lightyellow
# # #       "results" = "#90EE90",   # lightgreen
# # #       "data" = "#D3D3D3",      # lightgrey
# # #       "other" = "#FFFFFF"      # white
# # #     )
# # #     
# # #     # Create the tbl_graph and add coordinates
# # #     g_tbl <- tidygraph::as_tbl_graph(g)
# # #     g_tbl <- g_tbl %>%
# # #       tidygraph::activate(nodes) %>%
# # #       dplyr::mutate(x = layout_coords[as.character(name), 1],
# # #                     y = layout_coords[as.character(name), 2],
# # #                     color = colors[type])
# # #     
# # #     # Create the plot with ggraph
# # #     p <- ggraph::ggraph(g_tbl, layout = "manual", x = x, y = y) +
# # #       ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(4, 'mm')), 
# # #                              end_cap = ggraph::circle(3, 'mm'),
# # #                              start_cap = ggraph::circle(3, 'mm'),
# # #                              alpha = 0.6) +
# # #       ggraph::geom_node_point(aes(color = type), size = 20, alpha = 0.8) +
# # #       ggraph::geom_node_text(aes(label = name), size = 3) +
# # #       ggplot2::scale_color_manual(values = colors) +
# # #       ggraph::theme_graph() +
# # #       ggplot2::labs(title = "HLA-I Analysis Pipeline Structure") +
# # #       ggplot2::theme(legend.position = "right")
# # #     
# # #     # Save the plot
# # #     ggplot2::ggsave("hla_pipeline_structure_ggraph.png", p, width = 10, height = 8, dpi = 100)
# # #     
# # #     message("ggraph visualization successful! Saved as:")
# # #     message("- hla_pipeline_structure_ggraph.png")
# # #     
# # #     return(TRUE)
# # #   }, error = function(e) {
# # #     message("Error creating ggraph visualization: ", e$message)
# # #     message("You can still use the basic visualizations that were created earlier.")
# # #     return(FALSE)
# # #   })
# # # }
# # # 
# # # # Alternative function to visualize using just base R plotting (no package dependencies)
# # # create_simple_base_visualization <- function() {
# # #   # Create the edge data
# # #   connections <- list(
# # #     "rawdata" = "data",
# # #     "data" = c("run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R"),
# # #     "config" = c("run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R"),
# # #     "src/core" = c("run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R"),
# # #     "src/modules" = c("run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R"),
# # #     "run_peptide_pipeline.R" = "results/peptide_analysis",
# # #     "run_peptide_transcriptome_pipeline.R" = "results/transcriptome_analysis",
# # #     "results/peptide_analysis" = "results/tumor_normal_transcriptome",
# # #     "results/transcriptome_analysis" = "results/tumor_normal_transcriptome"
# # #   )
# # #   
# # #   # Node positions
# # #   nodes <- list(
# # #     "rawdata" = c(2, 9),
# # #     "data" = c(2, 7),
# # #     "config" = c(2, 5),
# # #     "src/core" = c(1, 3),
# # #     "src/modules" = c(3, 3),
# # #     "run_peptide_pipeline.R" = c(1, 1),
# # #     "run_peptide_transcriptome_pipeline.R" = c(5, 1),
# # #     "results/peptide_analysis" = c(1, -1),
# # #     "results/transcriptome_analysis" = c(5, -1),
# # #     "results/tumor_normal_transcriptome" = c(3, -3)
# # #   )
# # #   
# # #   # Node colors
# # #   node_types <- list(
# # #     "pipeline" = c("run_peptide_pipeline.R", "run_peptide_transcriptome_pipeline.R"),
# # #     "source" = c("src/core", "src/modules"),
# # #     "config" = c("config"),
# # #     "results" = c("results/peptide_analysis", "results/transcriptome_analysis", "results/tumor_normal_transcriptome"),
# # #     "data" = c("data", "rawdata")
# # #   )
# # #   
# # #   colors <- c(
# # #     "pipeline" = "#F08080",  # lightcoral
# # #     "source" = "#ADD8E6",    # lightblue
# # #     "config" = "#FFFFE0",    # lightyellow
# # #     "results" = "#90EE90",   # lightgreen
# # #     "data" = "#D3D3D3"       # lightgrey
# # #   )
# # #   
# # #   # Create a function to get node color
# # #   get_node_color <- function(node_name) {
# # #     for (type in names(node_types)) {
# # #       if (node_name %in% node_types[[type]]) {
# # #         return(colors[type])
# # #       }
# # #     }
# # #     return("#FFFFFF")  # default color
# # #   }
# # #   
# # #   # Create the visualization
# # #   png("hla_pipeline_structure_simple.png", width = 1000, height = 800, res = 100)
# # #   
# # #   # Setup the plot
# # #   plot(1, type = "n", xlim = c(0, 6), ylim = c(-4, 10), 
# # #        xlab = "", ylab = "", main = "HLA-I Analysis Pipeline Structure",
# # #        axes = FALSE)
# # #   
# # #   # Draw nodes
# # #   for (node_name in names(nodes)) {
# # #     x <- nodes[[node_name]][1]
# # #     y <- nodes[[node_name]][2]
# # #     color <- get_node_color(node_name)
# # #     
# # #     # Draw node
# # #     symbols(x, y, circles = 0.4, inches = FALSE, add = TRUE, bg = color)
# # #     
# # #     # Add label
# # #     text(x, y, node_name, cex = 0.7)
# # #   }
# # #   
# # #   # Draw edges
# # #   for (from_node in names(connections)) {
# # #     from_x <- nodes[[from_node]][1]
# # #     from_y <- nodes[[from_node]][2]
# # #     
# # #     to_nodes <- connections[[from_node]]
# # #     if (!is.list(to_nodes) && !is.vector(to_nodes) && !is.null(to_nodes)) {
# # #       to_nodes <- c(to_nodes)  # Make sure it's a vector
# # #     }
# # #     
# # #     for (to_node in to_nodes) {
# # #       to_x <- nodes[[to_node]][1]
# # #       to_y <- nodes[[to_node]][2]
# # #       
# # #       # Calculate arrow position
# # #       angle <- atan2(to_y - from_y, to_x - from_x)
# # #       from_r <- 0.4  # node radius
# # #       to_r <- 0.4
# # #       
# # #       new_from_x <- from_x + from_r * cos(angle)
# # #       new_from_y <- from_y + from_r * sin(angle)
# # #       new_to_x <- to_x - to_r * cos(angle)
# # #       new_to_y <- to_y - to_r * sin(angle)
# # #       
# # #       # Draw arrow
# # #       arrows(new_from_x, new_from_y, new_to_x, new_to_y, 
# # #              length = 0.1, angle = 20, code = 2, lwd = 1)
# # #     }
# # #   }
# # #   
# # #   # Add legend
# # #   legend("topright", 
# # #          legend = names(colors), 
# # #          fill = colors, 
# # #          title = "Component Type", 
# # #          cex = 0.8)
# # #   
# # #   dev.off()
# # #   
# # #   message("Simple base R visualization saved as 'hla_pipeline_structure_simple.png'")
# # #   return(TRUE)
# # # }
# # # 
# # # # Run the base visualization first (most reliable)
# # # create_basic_visualization()
# # # 
# # # # Create the simple base R visualization (as a fallback)
# # # create_simple_base_visualization()
# # # 
# # # # Try the ggraph visualization if you want (might error)
# # # # Uncomment if you want to try this:
# # # # try_ggraph_visualization()
# # # 
# # # # Display a message in the RStudio console
# # # cat("\nHLA-I Analysis Pipeline Structure has been visualized.\n")
# # # cat("You can find the visualization files in your working directory.\n")