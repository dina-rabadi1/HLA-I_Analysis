# Install required packages if not already installed
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

# Load required libraries
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Set working directory to the raw data folder
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis/rawdata")

# Read the Excel file
gene_data <- read_excel("Normalized_Gene_counts_FLCdb_Panel_1.xlsx")

# Check dimensions of the data
dim(gene_data)
colnames(gene_data)

# Quick check of data values - are they what we expect?
summary(gene_data[c("RU148_N", "RU148_T8", "RU148_T11")])

# Check for missing values
missing_values <- colSums(is.na(gene_data))
print(missing_values)

# Check data distribution
hist(gene_data$RU148_N, main="Distribution of Normal Sample Values", breaks=50)
hist(gene_data$RU148_T8, main="Distribution of Tumor Sample T8 Values", breaks=50)
hist(gene_data$RU148_T11, main="Distribution of Tumor Sample T11 Values", breaks=50)

# Calculate some test fold changes to see their distribution
gene_data$log2FC_T8 <- log2((gene_data$RU148_T8 + 0.1) / (gene_data$RU148_N + 0.1))
gene_data$log2FC_T11 <- log2((gene_data$RU148_T11 + 0.1) / (gene_data$RU148_N + 0.1))

# Visualize fold change distribution
hist(gene_data$log2FC_T8, main="Log2 Fold Change Distribution T8", breaks=50)
hist(gene_data$log2FC_T11, main="Log2 Fold Change Distribution T11", breaks=50)

# Check how many genes would pass different fold change thresholds
cat("Genes with log2FC > 1 in T8:", sum(gene_data$log2FC_T8 > 1, na.rm=TRUE), "\n")
cat("Genes with log2FC > 1 in T11:", sum(gene_data$log2FC_T11 > 1, na.rm=TRUE), "\n")
cat("Genes with log2FC > 0.5 in T8:", sum(gene_data$log2FC_T8 > 0.5, na.rm=TRUE), "\n")
cat("Genes with log2FC > 0.5 in T11:", sum(gene_data$log2FC_T11 > 0.5, na.rm=TRUE), "\n")

# Check expression level distribution
cat("Genes with expression > 10 in T8:", sum(gene_data$RU148_T8 > 10, na.rm=TRUE), "\n")
cat("Genes with expression > 10 in T11:", sum(gene_data$RU148_T11 > 10, na.rm=TRUE), "\n")
cat("Genes with expression > 5 in T8:", sum(gene_data$RU148_T8 > 5, na.rm=TRUE), "\n")
cat("Genes with expression > 5 in T11:", sum(gene_data$RU148_T11 > 5, na.rm=TRUE), "\n")

# Check how many genes meet combined criteria with less stringent thresholds
cat("Genes with log2FC > 0.5 AND expression > 5 in T8:", 
    sum(gene_data$log2FC_T8 > 0.5 & gene_data$RU148_T8 > 5, na.rm=TRUE), "\n")
cat("Genes with log2FC > 0.5 AND expression > 5 in T11:", 
    sum(gene_data$log2FC_T11 > 0.5 & gene_data$RU148_T11 > 5, na.rm=TRUE), "\n")