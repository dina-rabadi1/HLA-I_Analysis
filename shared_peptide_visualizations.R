# R script for interactive peptide sample distribution visualization
# shared_peptide_visualizations.R

# Setting directory
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis")

# Path to the combined peptides file
combined_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/20250516/correlation_results_unique_peptides_all_20250520_101444/unique_peptides_all.tsv"

# Create timestamped main output directory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
input_filename <- basename(combined_file)
input_name_no_ext <- tools::file_path_sans_ext(input_filename)
main_output_dir <- file.path(dirname(combined_file), 
                             paste0("shared_peptide_visualizations_", 
                                    input_name_no_ext, "_", timestamp))
dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)

# Create subdirectories within main output directory
viz_dir <- file.path(main_output_dir, "visualizations")
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

excel_dir <- file.path(main_output_dir, "excel_reports")
dir.create(excel_dir, recursive = TRUE, showWarnings = FALSE)

# Load required libraries
required_packages <- c("tidyverse", "plotly", "htmlwidgets", "openxlsx", 
                       "UpSetR", "VennDiagram", "RColorBrewer", "ggdendro")

# Install and load necessary packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Create a function to add input file information to Excel reports
add_input_info_to_excel <- function(wb, input_file, timestamp) {
  # Create a new worksheet for input information
  addWorksheet(wb, "Input Info")
  
  # Create input information
  input_info <- data.frame(
    Parameter = c("Input File", "Analysis Date", "File Path"),
    Value = c(basename(input_file), timestamp, input_file)
  )
  
  # Write information
  writeData(wb, "Input Info", input_info, startRow = 1, startCol = 1)
  addStyle(wb, "Input Info", headerStyle, rows = 1, cols = 1:ncol(input_info))
  
  # Move worksheet to first position
  worksheetOrder(wb) <- c(length(names(wb)), 1:(length(names(wb))-1))
}

# Read the combined peptides file
cat("Reading combined peptides file...\n")
peptides_data <- read_tsv(combined_file, show_col_types = FALSE)

# Filter to include only 8-12mers
cat("Filtering to include only 8-12mers...\n")
peptides_8_12 <- peptides_data %>%
  filter(`Peptide Length` >= 8 & `Peptide Length` <= 12)

cat("Number of 8-12mer peptides:", n_distinct(peptides_8_12$Peptide), "\n")

# Create a report showing which samples each peptide appears in
# Now including gene and protein information
peptide_report <- peptides_8_12 %>%
  # First get basic peptide and sample information
  select(Peptide, SampleID, `Peptide Length`, Gene, Protein) %>%
  distinct() %>%
  # Group by peptide to consolidate sample information
  group_by(Peptide, `Peptide Length`) %>%
  summarize(
    sample_list = paste(sort(SampleID), collapse=", "),
    sample_count = n_distinct(SampleID),
    # Consolidate gene information
    genes = paste(unique(na.omit(Gene)), collapse="; "),
    # Consolidate protein information
    proteins = paste(unique(na.omit(Protein)), collapse="; "),
    .groups = "drop"
  ) %>%
  arrange(desc(sample_count), Peptide)

# NOW write the report to a file
output_report <- file.path(main_output_dir, "peptide_sample_distribution.tsv")
write_tsv(peptide_report, output_report)
cat("Peptide sample distribution report written to:", output_report, "\n")

# Print summary of the report
cat("\nSample count distribution:\n")
print(table(peptide_report$sample_count))

# Get number of unique samples in the dataset
unique_samples <- peptides_8_12 %>%
  select(SampleID) %>%
  distinct() %>%
  pull()

num_samples <- length(unique_samples)
cat("Number of unique samples in dataset:", num_samples, "\n")

# Create style for Excel headers
headerStyle <- createStyle(textDecoration = "bold", fgFill = "#D9D9D9")

################################################################################
# 1. Distribution of peptides by sample count - WITH EXCEL REPORT
################################################################################

# 1.1 Prepare the data
sample_count_dist <- peptide_report %>%
  count(sample_count) %>%
  mutate(percentage = n / sum(n) * 100)

# 1.2 Create the interactive plot with ggplot2 + plotly
p1 <- ggplot(sample_count_dist, aes(x = factor(sample_count), y = n,
                                    text = paste0("Sample count: ", sample_count,
                                                  "<br>Peptides: ", n,
                                                  "<br>Percentage: ", round(percentage, 1), "%"))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Peptide Distribution by Sample Count",
    x = "Number of Samples",
    y = "Number of Peptides"
  ) +
  theme_minimal()

# Convert to interactive plotly plot
p1_interactive <- ggplotly(p1, tooltip = "text")

# Save as interactive HTML
htmlwidgets::saveWidget(p1_interactive,
                        paste0(viz_dir, "/interactive_sample_count_distribution.html"),
                        selfcontained = TRUE)

# 1.3 Create Excel report
# Add detailed data about peptides in each sample count category
peptides_by_sample_count <- peptide_report %>%
  group_by(sample_count) %>%
  summarise(
    num_peptides = n(),
    percentage = round(n() / nrow(peptide_report) * 100, 2),
    example_peptides = paste0(head(Peptide, 10), collapse = ", "),
    most_common_genes = paste0(head(names(sort(table(unlist(strsplit(genes, "; "))), decreasing = TRUE)), 5), collapse = ", "),
    .groups = "drop"
  )

# Create a workbook for sample count distribution
wb_sample_count <- createWorkbook()

# Add summary worksheet
addWorksheet(wb_sample_count, "Sample Count Summary")
writeData(wb_sample_count, "Sample Count Summary", sample_count_dist, startRow = 1, startCol = 1)

# Add style
addStyle(wb_sample_count, "Sample Count Summary", headerStyle, rows = 1, cols = 1:ncol(sample_count_dist))

# Add detailed worksheet
addWorksheet(wb_sample_count, "Sample Count Details")
writeData(wb_sample_count, "Sample Count Details", peptides_by_sample_count, startRow = 1, startCol = 1)
addStyle(wb_sample_count, "Sample Count Details", headerStyle, rows = 1, cols = 1:ncol(peptides_by_sample_count))

# Add raw data of all peptides with sample counts
addWorksheet(wb_sample_count, "All Peptides")
writeData(wb_sample_count, "All Peptides",
          peptide_report %>% select(Peptide, `Peptide Length`, sample_count, genes),
          startRow = 1, startCol = 1)
addStyle(wb_sample_count, "All Peptides", headerStyle, rows = 1, cols = 1:4)

# Save workbook
saveWorkbook(wb_sample_count, paste0(excel_dir, "/peptide_sample_count_distribution.xlsx"), overwrite = TRUE)

################################################################################
# 2. Unique vs Shared Peptides - WITH EXCEL REPORT
################################################################################

# 2.1 Prepare the data
unique_shared <- data.frame(
  category = c("Unique", "Shared"),
  count = c(
    sum(sample_count_dist$n[sample_count_dist$sample_count == 1]),
    sum(sample_count_dist$n[sample_count_dist$sample_count > 1])
  )
)

unique_shared <- unique_shared %>%
  mutate(percentage = count / sum(count) * 100)

# 2.2 Create the interactive plot
p2 <- plot_ly(unique_shared, labels = ~category, values = ~count, type = 'pie',
              text = ~paste0(category, ": ", count, " peptides (", round(percentage, 1), "%)"),
              hoverinfo = 'text', textposition = 'inside',
              marker = list(colors = c("lightblue", "lightgreen"),
                            line = list(color = '#FFFFFF', width = 1))) %>%
  layout(title = "Unique vs Shared Peptides",
         showlegend = TRUE,
         legend = list(orientation = "h"))

# Save as interactive HTML
htmlwidgets::saveWidget(p2,
                        paste0(viz_dir, "/interactive_unique_vs_shared.html"),
                        selfcontained = TRUE)

# 2.3 Create Excel report
# Create detailed data for unique vs shared peptides
unique_peptides <- peptide_report %>%
  filter(sample_count == 1) %>%
  select(Peptide, `Peptide Length`, genes, sample_list)

shared_peptides <- peptide_report %>%
  filter(sample_count > 1) %>%
  select(Peptide, `Peptide Length`, sample_count, genes, sample_list)

# Create summary by gene for unique vs shared
gene_summary <- peptide_report %>%
  mutate(category = ifelse(sample_count == 1, "Unique", "Shared")) %>%
  # Split genes and unnest
  mutate(gene_list = strsplit(genes, "; ")) %>%
  unnest(gene_list) %>%
  count(category, gene_list) %>%
  spread(category, n, fill = 0) %>%
  mutate(total = Unique + Shared,
         pct_shared = round(Shared / total * 100, 1)) %>%
  arrange(desc(total))

# Create workbook
wb_unique_shared <- createWorkbook()

# Add summary worksheet
addWorksheet(wb_unique_shared, "Summary")
writeData(wb_unique_shared, "Summary", unique_shared, startRow = 1, startCol = 1)
addStyle(wb_unique_shared, "Summary", headerStyle, rows = 1, cols = 1:ncol(unique_shared))

# Add gene summary worksheet
addWorksheet(wb_unique_shared, "Gene Summary")
writeData(wb_unique_shared, "Gene Summary", gene_summary, startRow = 1, startCol = 1)
addStyle(wb_unique_shared, "Gene Summary", headerStyle, rows = 1, cols = 1:ncol(gene_summary))

# Add unique peptides worksheet
addWorksheet(wb_unique_shared, "Unique Peptides")
writeData(wb_unique_shared, "Unique Peptides", unique_peptides, startRow = 1, startCol = 1)
addStyle(wb_unique_shared, "Unique Peptides", headerStyle, rows = 1, cols = 1:ncol(unique_peptides))

# Add shared peptides worksheet
addWorksheet(wb_unique_shared, "Shared Peptides")
writeData(wb_unique_shared, "Shared Peptides", shared_peptides, startRow = 1, startCol = 1)
addStyle(wb_unique_shared, "Shared Peptides", headerStyle, rows = 1, cols = 1:ncol(shared_peptides))

# Save workbook
saveWorkbook(wb_unique_shared, paste0(excel_dir, "/unique_vs_shared_peptides.xlsx"), overwrite = TRUE)

################################################################################
# 3. Peptide Length Distribution - WITH EXCEL REPORT
################################################################################

# 3.1 Prepare the data
length_dist <- peptide_report %>%
  count(`Peptide Length`) %>%
  mutate(percentage = n / sum(n) * 100)

# 3.2 Create the interactive plot
p3 <- ggplot(length_dist, aes(x = factor(`Peptide Length`), y = n,
                              text = paste0("Length: ", `Peptide Length`, " amino acids",
                                            "<br>Count: ", n,
                                            "<br>Percentage: ", round(percentage, 1), "%"))) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  labs(
    title = "Peptide Length Distribution",
    x = "Peptide Length",
    y = "Number of Peptides"
  ) +
  theme_minimal()

# Convert to interactive plotly plot
p3_interactive <- ggplotly(p3, tooltip = "text")

# Save as interactive HTML
htmlwidgets::saveWidget(p3_interactive,
                        paste0(viz_dir, "/interactive_length_distribution.html"),
                        selfcontained = TRUE)

# 3.3 Create Excel report
# Create detailed data for peptide length
peptides_by_length <- peptide_report %>%
  group_by(`Peptide Length`) %>%
  summarise(
    num_peptides = n(),
    percentage = round(n() / nrow(peptide_report) * 100, 2),
    unique_peptides = sum(sample_count == 1),
    shared_peptides = sum(sample_count > 1),
    pct_shared = round(shared_peptides / num_peptides * 100, 1),
    avg_sample_count = mean(sample_count),
    example_peptides = paste0(head(Peptide, 10), collapse = ", "),
    .groups = "drop"
  )

# Create a detailed gene breakdown by length
gene_length_summary <- peptide_report %>%
  # Split genes and unnest
  mutate(gene_list = strsplit(genes, "; ")) %>%
  unnest(gene_list) %>%
  count(`Peptide Length`, gene_list) %>%
  group_by(`Peptide Length`) %>%
  arrange(desc(n)) %>%
  summarise(
    top_genes = paste0(head(paste0(gene_list, " (", n, ")"), 10), collapse = ", "),
    .groups = "drop"
  )

# Create workbook
wb_length <- createWorkbook()

# Add summary worksheet
addWorksheet(wb_length, "Length Summary")
writeData(wb_length, "Length Summary", length_dist, startRow = 1, startCol = 1)
addStyle(wb_length, "Length Summary", headerStyle, rows = 1, cols = 1:ncol(length_dist))

# Add detailed summary
addWorksheet(wb_length, "Length Details")
writeData(wb_length, "Length Details", peptides_by_length, startRow = 1, startCol = 1)
addStyle(wb_length, "Length Details", headerStyle, rows = 1, cols = 1:ncol(peptides_by_length))

# Add gene summary by length
addWorksheet(wb_length, "Gene by Length")
writeData(wb_length, "Gene by Length", gene_length_summary, startRow = 1, startCol = 1)
addStyle(wb_length, "Gene by Length", headerStyle, rows = 1, cols = 1:ncol(gene_length_summary))

# Save workbook
saveWorkbook(wb_length, paste0(excel_dir, "/peptide_length_distribution.xlsx"), overwrite = TRUE)

################################################################################
# 4. Average sharing by peptide length - WITH EXCEL REPORT
################################################################################

# 4.1 Prepare the data
avg_sharing_by_length <- peptide_report %>%
  group_by(`Peptide Length`) %>%
  summarise(
    avg_samples = mean(sample_count),
    median_samples = median(sample_count),
    total_peptides = n(),
    unique_count = sum(sample_count == 1),
    shared_count = sum(sample_count > 1),
    percent_shared = round(shared_count / total_peptides * 100, 1)
  )

# 4.2 Create the interactive plot
p4 <- ggplot(avg_sharing_by_length, aes(x = factor(`Peptide Length`), y = avg_samples,
                                        text = paste0("Length: ", `Peptide Length`, " amino acids",
                                                      "<br>Avg samples: ", round(avg_samples, 2),
                                                      "<br>Median samples: ", median_samples,
                                                      "<br>Total peptides: ", total_peptides,
                                                      "<br>% shared: ", percent_shared, "%"))) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(
    title = "Average Sharing by Peptide Length",
    x = "Peptide Length",
    y = "Average Number of Samples"
  ) +
  theme_minimal()

# Convert to interactive plotly plot
p4_interactive <- ggplotly(p4, tooltip = "text")

# Save as interactive HTML
htmlwidgets::saveWidget(p4_interactive,
                        paste0(viz_dir, "/interactive_avg_sharing_by_length.html"),
                        selfcontained = TRUE)

# 4.3 Create Excel report
# Create detailed data for sharing patterns by length
sharing_patterns <- peptide_report %>%
  group_by(`Peptide Length`, sample_count) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = sample_count, values_from = count, values_fill = 0) %>%
  rowwise() %>%
  mutate(total = sum(c_across(-`Peptide Length`)))

# Create sample count distribution by peptide length
length_sharing_distribution <- peptide_report %>%
  group_by(`Peptide Length`) %>%
  count(sample_count) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  arrange(`Peptide Length`, sample_count)

# Create workbook
wb_sharing <- createWorkbook()

# Add summary worksheet
addWorksheet(wb_sharing, "Sharing Summary")
writeData(wb_sharing, "Sharing Summary", avg_sharing_by_length, startRow = 1, startCol = 1)
addStyle(wb_sharing, "Sharing Summary", headerStyle, rows = 1, cols = 1:ncol(avg_sharing_by_length))

# Add sharing patterns
addWorksheet(wb_sharing, "Sharing Patterns")
writeData(wb_sharing, "Sharing Patterns", sharing_patterns, startRow = 1, startCol = 1)
addStyle(wb_sharing, "Sharing Patterns", headerStyle, rows = 1, cols = 1:ncol(sharing_patterns))

# Add detailed distribution
addWorksheet(wb_sharing, "Detailed Distribution")
writeData(wb_sharing, "Detailed Distribution", length_sharing_distribution, startRow = 1, startCol = 1)
addStyle(wb_sharing, "Detailed Distribution", headerStyle, rows = 1, cols = 1:ncol(length_sharing_distribution))

# Save workbook
saveWorkbook(wb_sharing, paste0(excel_dir, "/peptide_sharing_by_length.xlsx"), overwrite = TRUE)

################################################################################
# 5. Gene sharing heatmap - WITH EXCEL REPORT
################################################################################

# Get top genes by peptide count
top_genes <- peptide_report %>%
  # Split multiple genes
  mutate(gene_list = strsplit(genes, "; ")) %>%
  unnest(gene_list) %>%
  count(gene_list, sort = TRUE) %>%
  head(15) %>%
  pull(gene_list)

# Function to create gene-based heatmap data
create_gene_heatmap_data <- function(genes_to_include) {
  # Create a matrix of gene x sample
  gene_sample_matrix <- peptides_8_12 %>%
    filter(!is.na(Gene)) %>%
    # Split genes if a peptide belongs to multiple genes
    mutate(gene_list = strsplit(Gene, "; ")) %>%
    unnest(gene_list) %>%
    filter(gene_list %in% genes_to_include) %>%
    select(gene_list, SampleID, Peptide) %>%
    distinct() %>%
    count(gene_list, SampleID) %>%
    pivot_wider(
      names_from = SampleID,
      values_from = n,
      values_fill = 0
    )
  
  return(gene_sample_matrix)
}

# Create gene heatmap data
gene_matrix_df <- create_gene_heatmap_data(top_genes)

# Convert to matrix format
gene_matrix <- as.matrix(gene_matrix_df[, -1])
rownames(gene_matrix) <- gene_matrix_df$gene_list

# 5.1 Create interactive heatmap
# Convert matrix to long format for plotly
heatmap_data <- as.data.frame(gene_matrix)
heatmap_data$gene <- rownames(gene_matrix)
heatmap_data_long <- pivot_longer(heatmap_data,
                                  cols = -gene,
                                  names_to = "sample",
                                  values_to = "count")

# Create the interactive heatmap
p5 <- plot_ly(heatmap_data_long,
              x = ~sample,
              y = ~gene,
              z = ~count,
              type = "heatmap",
              colorscale = "Viridis",
              hoverinfo = "text",
              text = ~paste0("Gene: ", gene,
                             "<br>Sample: ", sample,
                             "<br>Peptides: ", count))

p5 <- p5 %>% layout(title = "Peptides by Gene Across Samples",
                    xaxis = list(title = "Sample"),
                    yaxis = list(title = "Gene"))

# Save as interactive HTML
htmlwidgets::saveWidget(p5,
                        paste0(viz_dir, "/interactive_gene_heatmap.html"),
                        selfcontained = TRUE)

# 5.2 Create Excel report
# Calculate shared peptides between samples for each gene
calculate_gene_sharing <- function(genes_to_include) {
  gene_sharing_data <- list()
  
  for (gene in genes_to_include) {
    # Get peptides for this gene
    gene_peptides <- peptides_8_12 %>%
      filter(grepl(gene, Gene, fixed = TRUE)) %>%
      select(Peptide, SampleID) %>%
      distinct()
    
    # Get unique peptides for this gene
    unique_peptides <- gene_peptides %>%
      count(Peptide) %>%
      filter(n == 1) %>%
      nrow()
    
    # Get shared peptides for this gene
    shared_peptides <- gene_peptides %>%
      count(Peptide) %>%
      filter(n > 1) %>%
      nrow()
    
    # Total peptides
    total_peptides <- n_distinct(gene_peptides$Peptide)
    
    # Samples containing this gene
    samples_with_gene <- unique(gene_peptides$SampleID)
    
    # Create row for this gene
    gene_sharing_data[[gene]] <- data.frame(
      gene = gene,
      total_peptides = total_peptides,
      unique_peptides = unique_peptides,
      shared_peptides = shared_peptides,
      percent_shared = ifelse(total_peptides > 0, round(shared_peptides / total_peptides * 100, 1), 0),
      sample_count = length(samples_with_gene),
      samples = paste(samples_with_gene, collapse = ", ")
    )
  }
  
  # Combine all genes
  result <- bind_rows(gene_sharing_data)
  return(result)
}

# Get gene sharing data
gene_sharing <- calculate_gene_sharing(top_genes)

# Get sample-specific gene data
sample_gene_data <- gene_matrix_df %>%
  pivot_longer(cols = -gene_list,
               names_to = "sample",
               values_to = "peptide_count") %>%
  arrange(gene_list, desc(peptide_count))

# Create workbook
wb_gene <- createWorkbook()

# Add gene matrix
addWorksheet(wb_gene, "Gene Matrix")
writeData(wb_gene, "Gene Matrix", gene_matrix_df, startRow = 1, startCol = 1)
addStyle(wb_gene, "Gene Matrix", headerStyle, rows = 1, cols = 1:ncol(gene_matrix_df))

# Add gene sharing summary
addWorksheet(wb_gene, "Gene Sharing Summary")
writeData(wb_gene, "Gene Sharing Summary", gene_sharing, startRow = 1, startCol = 1)
addStyle(wb_gene, "Gene Sharing Summary", headerStyle, rows = 1, cols = 1:ncol(gene_sharing))

# Add sample-gene data
addWorksheet(wb_gene, "Sample Gene Data")
writeData(wb_gene, "Sample Gene Data", sample_gene_data, startRow = 1, startCol = 1)
addStyle(wb_gene, "Sample Gene Data", headerStyle, rows = 1, cols = 1:ncol(sample_gene_data))

# Save workbook
saveWorkbook(wb_gene, paste0(excel_dir, "/gene_peptide_distribution.xlsx"), overwrite = TRUE)

################################################################################
# 6. Sample intersections - WITH EXCEL REPORT
################################################################################

# If we have 2-4 samples, create an interactive Venn diagram
if (num_samples >= 2 && num_samples <= 4) {
  # Create lists of peptides for each sample
  sample_peptides <- list()
  
  for (sample in unique_samples) {
    sample_peptides[[sample]] <- peptides_8_12 %>%
      filter(SampleID == sample) %>%
      pull(Peptide) %>%
      unique()
  }
  
  # Calculate all possible intersections
  intersection_data <- list()
  samples_list <- names(sample_peptides)
  
  # Add single sample data
  for (i in 1:length(samples_list)) {
    s1 <- samples_list[i]
    set_name <- s1
    intersection_data[[set_name]] <- length(sample_peptides[[s1]])
  }
  
  # Add two-sample intersections
  if (length(samples_list) >= 2) {
    for (i in 1:(length(samples_list)-1)) {
      for (j in (i+1):length(samples_list)) {
        s1 <- samples_list[i]
        s2 <- samples_list[j]
        set_name <- paste(s1, s2, sep = " & ")
        intersection_data[[set_name]] <- length(intersect(sample_peptides[[s1]], sample_peptides[[s2]]))
      }
    }
  }
  
  # Add three-sample intersections
  if (length(samples_list) >= 3) {
    for (i in 1:(length(samples_list)-2)) {
      for (j in (i+1):(length(samples_list)-1)) {
        for (k in (j+1):length(samples_list)) {
          s1 <- samples_list[i]
          s2 <- samples_list[j]
          s3 <- samples_list[k]
          set_name <- paste(s1, s2, s3, sep = " & ")
          intersection_data[[set_name]] <- length(intersect(intersect(sample_peptides[[s1]], sample_peptides[[s2]]), sample_peptides[[s3]]))
        }
      }
    }
  }
  
  # Add four-sample intersection
  if (length(samples_list) == 4) {
    s1 <- samples_list[1]
    s2 <- samples_list[2]
    s3 <- samples_list[3]
    s4 <- samples_list[4]
    set_name <- paste(s1, s2, s3, s4, sep = " & ")
    intersection_data[[set_name]] <- length(intersect(intersect(intersect(sample_peptides[[s1]], sample_peptides[[s2]]), sample_peptides[[s3]]), sample_peptides[[s4]]))
  }
  
  # Convert to data frame
  intersection_df <- data.frame(
    intersection = names(intersection_data),
    count = unlist(intersection_data)
  ) %>%
    arrange(desc(count))
  
  # 6.1 Create interactive bar chart of intersections
  p6 <- plot_ly(intersection_df,
                x = ~intersection,
                y = ~count,
                type = "bar",
                marker = list(color = "purple"),
                hoverinfo = "text",
                text = ~paste0("Intersection: ", intersection,
                               "<br>Count: ", count, " peptides"))
  
  p6 <- p6 %>% layout(title = "Sample Intersections",
                      xaxis = list(title = "Sample Combination", tickangle = 45),
                      yaxis = list(title = "Number of Peptides"))
  
  # Save as interactive HTML
  htmlwidgets::saveWidget(p6,
                          paste0(viz_dir, "/interactive_sample_intersections.html"),
                          selfcontained = TRUE)
  
  # 6.2 Create Excel report with detailed intersection data
  # Get peptides in each intersection
  intersection_peptides <- list()
  
  # Single samples
  for (i in 1:length(samples_list)) {
    s1 <- samples_list[i]
    set_name <- s1
    intersection_peptides[[set_name]] <- peptides_8_12 %>%
      filter(SampleID == s1) %>%
      select(Peptide, `Peptide Length`, Gene) %>%
      distinct()
  }
  
  # Two-sample intersections
  if (length(samples_list) >= 2) {
    for (i in 1:(length(samples_list)-1)) {
      for (j in (i+1):length(samples_list)) {
        s1 <- samples_list[i]
        s2 <- samples_list[j]
        set_name <- paste(s1, s2, sep = " & ")
        
        p1 <- peptides_8_12 %>% filter(SampleID == s1) %>% pull(Peptide) %>% unique()
        p2 <- peptides_8_12 %>% filter(SampleID == s2) %>% pull(Peptide) %>% unique()
        shared <- intersect(p1, p2)
        
        intersection_peptides[[set_name]] <- peptides_8_12 %>%
          filter(Peptide %in% shared) %>%
          select(Peptide, `Peptide Length`, Gene) %>%
          distinct()
      }
    }
  }
  
  # Create workbook
  wb_intersections <- createWorkbook()
  
  # Add summary worksheet
  addWorksheet(wb_intersections, "Intersection Summary")
  writeData(wb_intersections, "Intersection Summary", intersection_df, startRow = 1, startCol = 1)
  addStyle(wb_intersections, "Intersection Summary", headerStyle, rows = 1, cols = 1:ncol(intersection_df))
  
  # Add worksheets for each intersection
  for (set_name in names(intersection_peptides)) {
    # Excel doesn't like & in sheet names
    sheet_name <- gsub("&", "and", set_name)
    
    # Excel worksheet names limited to 31 chars
    if (nchar(sheet_name) > 31) {
      sheet_name <- substr(sheet_name, 1, 31)
    }
    
    addWorksheet(wb_intersections, sheet_name)
    writeData(wb_intersections, sheet_name, intersection_peptides[[set_name]], startRow = 1, startCol = 1)
    addStyle(wb_intersections, sheet_name, headerStyle, rows = 1, cols = 1:ncol(intersection_peptides[[set_name]]))
  }
  
  # Save workbook
  saveWorkbook(wb_intersections, paste0(excel_dir, "/sample_intersections.xlsx"), overwrite = TRUE)
}

################################################################################
# 8. Core Peptidome Analysis - Identify robustly shared peptides
################################################################################

# Create output directories for new analyses
core_peptidome_dir <- file.path(main_output_dir, "core_peptidome")
dir.create(core_peptidome_dir, recursive = TRUE, showWarnings = FALSE)

# 8.1 Define robust sharing thresholds
# Calculate the number of samples needed for different sharing thresholds
sharing_thresholds <- c(0.25, 0.5, 0.75, 0.9)
samples_needed <- ceiling(num_samples * sharing_thresholds)

# Create descriptive labels for each threshold
threshold_labels <- paste0(sharing_thresholds * 100, "% of samples (", samples_needed, "/", num_samples, ")")
names(samples_needed) <- threshold_labels

# 8.2 Identify robustly shared peptides at each threshold
robustly_shared_peptides <- list()

for (i in 1:length(samples_needed)) {
  threshold <- names(samples_needed)[i]
  min_samples <- samples_needed[i]
  
  # Filter peptides that appear in at least min_samples samples
  robust_peptides <- peptide_report %>%
    filter(sample_count >= min_samples) %>%
    arrange(desc(sample_count), Peptide)
  
  robustly_shared_peptides[[threshold]] <- robust_peptides
  
  cat("Peptides present in at least", threshold, ":", nrow(robust_peptides), "\n")
}

# 8.3 Create summary table of robust peptide counts
robust_summary <- data.frame(
  threshold = names(robustly_shared_peptides),
  peptide_count = sapply(robustly_shared_peptides, nrow),
  min_samples = samples_needed
)

# 8.4 Visualize the number of robustly shared peptides at each threshold
p_robust <- ggplot(robust_summary, aes(x = threshold, y = peptide_count, 
                                       text = paste0("Threshold: ", threshold,
                                                     "<br>Peptides: ", peptide_count))) +
  geom_bar(stat = "identity", fill = "darkblue") +
  labs(
    title = "Robustly Shared Peptides",
    x = "Sharing Threshold",
    y = "Number of Peptides"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Convert to interactive plotly plot
p_robust_interactive <- ggplotly(p_robust, tooltip = "text")

# Save as interactive HTML
htmlwidgets::saveWidget(p_robust_interactive,
                        paste0(viz_dir, "/interactive_robust_peptides.html"),
                        selfcontained = TRUE)

# 8.5 Create an Excel report for each threshold showing the peptides and their distribution
wb_robust <- createWorkbook()

# Add summary worksheet
addWorksheet(wb_robust, "Robust Sharing Summary")
writeData(wb_robust, "Robust Sharing Summary", robust_summary, startRow = 1, startCol = 1)
addStyle(wb_robust, "Robust Sharing Summary", headerStyle, rows = 1, cols = 1:ncol(robust_summary))

# Add worksheets for each threshold
for (i in 1:length(robustly_shared_peptides)) {
  threshold <- names(robustly_shared_peptides)[i]
  robust_set <- robustly_shared_peptides[[threshold]]
  
  # Clean up threshold name for Excel sheet (remove special characters)
  sheet_name <- gsub("[^A-Za-z0-9]", "_", threshold)
  if (nchar(sheet_name) > 31) {
    sheet_name <- substr(sheet_name, 1, 31)
  }
  
  addWorksheet(wb_robust, sheet_name)
  writeData(wb_robust, sheet_name, robust_set, startRow = 1, startCol = 1)
  addStyle(wb_robust, sheet_name, headerStyle, rows = 1, cols = 1:ncol(robust_set))
}

# Save workbook
saveWorkbook(wb_robust, paste0(excel_dir, "/robust_shared_peptides.xlsx"), overwrite = TRUE)

################################################################################
# 9. Sample-specific distribution of core peptides
################################################################################

# 9.1 Create a function to get sample-specific information for peptide sets
get_sample_distribution <- function(peptide_subset, label) {
  # First get the list of peptides
  peptide_list <- peptide_subset$Peptide
  
  # Find which samples have each peptide
  peptide_sample_dist <- peptides_8_12 %>%
    filter(Peptide %in% peptide_list) %>%
    select(Peptide, SampleID) %>%
    distinct() %>%
    # Count peptides per sample
    count(SampleID) %>%
    # Calculate percentage of the peptide subset present in each sample
    mutate(percentage = n / length(peptide_list) * 100) %>%
    arrange(desc(n))
  
  # Label this dataset
  peptide_sample_dist$threshold <- label
  
  return(peptide_sample_dist)
}

# Get sample distribution for each threshold
sample_distributions <- list()

for (i in 1:length(robustly_shared_peptides)) {
  threshold <- names(robustly_shared_peptides)[i]
  robust_set <- robustly_shared_peptides[[threshold]]
  
  sample_distributions[[threshold]] <- get_sample_distribution(robust_set, threshold)
}

# Combine all distributions into one dataframe
all_sample_distributions <- bind_rows(sample_distributions)

# 9.2 Visualize the sample distribution of core peptides
# Create a heatmap of samples vs thresholds
sample_threshold_matrix <- all_sample_distributions %>%
  select(SampleID, threshold, percentage) %>%
  pivot_wider(names_from = threshold, values_from = percentage)

# Convert to long format for plotting
sample_dist_long <- all_sample_distributions %>%
  select(SampleID, threshold, percentage)

# Create interactive heatmap
p_sample_dist <- plot_ly(sample_dist_long,
                         x = ~threshold,
                         y = ~SampleID,
                         z = ~percentage,
                         type = "heatmap",
                         colorscale = "Viridis",
                         hoverinfo = "text",
                         text = ~paste0("Sample: ", SampleID,
                                        "<br>Threshold: ", threshold,
                                        "<br>% of core peptides: ", round(percentage, 1), "%"))

p_sample_dist <- p_sample_dist %>% layout(
  title = "Sample Coverage of Core Peptidome",
  xaxis = list(title = "Sharing Threshold"),
  yaxis = list(title = "Sample"),
  height = 800  # Make it taller to accommodate all samples
)

# Save as interactive HTML
htmlwidgets::saveWidget(p_sample_dist,
                        paste0(viz_dir, "/interactive_sample_core_coverage.html"),
                        selfcontained = TRUE)

# 9.3 Create Excel report
wb_sample_dist <- createWorkbook()

# Add distribution data
addWorksheet(wb_sample_dist, "Sample Distributions")
writeData(wb_sample_dist, "Sample Distributions", all_sample_distributions, startRow = 1, startCol = 1)
addStyle(wb_sample_dist, "Sample Distributions", headerStyle, rows = 1, cols = 1:ncol(all_sample_distributions))

# Add matrix view
addWorksheet(wb_sample_dist, "Sample Matrix")
writeData(wb_sample_dist, "Sample Matrix", sample_threshold_matrix, startRow = 1, startCol = 1)
addStyle(wb_sample_dist, "Sample Matrix", headerStyle, rows = 1, cols = 1:ncol(sample_threshold_matrix))

# For each threshold, list samples by coverage percentage
for (i in 1:length(sharing_thresholds)) {
  threshold <- names(robustly_shared_peptides)[i]
  
  # Clean up name for Excel
  sheet_name <- paste0("Coverage_", i)
  
  sample_coverage <- sample_distributions[[threshold]] %>%
    arrange(desc(percentage))
  
  addWorksheet(wb_sample_dist, sheet_name)
  writeData(wb_sample_dist, sheet_name, sample_coverage, startRow = 1, startCol = 1)
  addStyle(wb_sample_dist, sheet_name, headerStyle, rows = 1, cols = 1:ncol(sample_coverage))
}

# Save workbook
saveWorkbook(wb_sample_dist, paste0(excel_dir, "/sample_core_coverage.xlsx"), overwrite = TRUE)

################################################################################
# 10. Create detailed sample-by-sample sharing matrix
################################################################################

# 10.1 Create a matrix showing how many peptides are shared between each pair of samples
cat("Creating sample-by-sample sharing matrix...\n")

# Function to count shared peptides between two samples
count_shared_peptides <- function(sample1, sample2, peptide_data) {
  peptides_sample1 <- peptide_data %>%
    filter(SampleID == sample1) %>%
    pull(Peptide) %>%
    unique()
  
  peptides_sample2 <- peptide_data %>%
    filter(SampleID == sample2) %>%
    pull(Peptide) %>%
    unique()
  
  shared_count <- length(intersect(peptides_sample1, peptides_sample2))
  return(shared_count)
}

# Create an empty matrix for sample-sample sharing
sample_sharing_matrix <- matrix(0, nrow = length(unique_samples), ncol = length(unique_samples))
rownames(sample_sharing_matrix) <- unique_samples
colnames(sample_sharing_matrix) <- unique_samples

# Fill the matrix with shared peptide counts
for (i in 1:length(unique_samples)) {
  for (j in 1:length(unique_samples)) {
    sample1 <- unique_samples[i]
    sample2 <- unique_samples[j]
    
    if (i == j) {
      # On diagonal, use total peptides for that sample
      total_peptides <- peptides_8_12 %>%
        filter(SampleID == sample1) %>%
        pull(Peptide) %>%
        unique() %>%
        length()
      
      sample_sharing_matrix[i, j] <- total_peptides
    } else {
      # Off diagonal, count shared peptides
      sample_sharing_matrix[i, j] <- count_shared_peptides(sample1, sample2, peptides_8_12)
    }
  }
}

# Convert matrix to a data frame
sample_sharing_df <- as.data.frame(sample_sharing_matrix)
sample_sharing_df$Sample1 <- rownames(sample_sharing_matrix)
sample_sharing_df <- sample_sharing_df %>%
  pivot_longer(cols = unique_samples,
               names_to = "Sample2",
               values_to = "shared_peptides")

# Calculate Jaccard similarity index (intersection/union)
sample_sharing_df <- sample_sharing_df %>%
  rowwise() %>%
  mutate(
    sample1_total = sample_sharing_matrix[Sample1, Sample1],
    sample2_total = sample_sharing_matrix[Sample2, Sample2],
    jaccard_index = shared_peptides / (sample1_total + sample2_total - shared_peptides)
  )

# 10.2 Visualize the sharing matrix using Jaccard index
p_jaccard <- plot_ly(sample_sharing_df,
                       x = ~Sample2,
                       y = ~Sample1,
                       z = ~jaccard_index,
                       type = "heatmap",
                       colorscale = "Hot",
                       reversescale = TRUE,  # This makes high values red/yellow, low values blue
                       hoverinfo = "text",
                       text = ~paste0("Sample1: ", Sample1,
                                      "<br>Sample2: ", Sample2,
                                      "<br>Shared peptides: ", shared_peptides,
                                      "<br>Jaccard index: ", round(jaccard_index, 3)))

p_jaccard <- p_jaccard %>% layout(
  title = "Sample Similarity (Jaccard Index)",
  xaxis = list(title = "Sample"),
  yaxis = list(title = "Sample")
)

# Save as interactive HTML
htmlwidgets::saveWidget(p_jaccard,
                        paste0(viz_dir, "/interactive_sample_similarity.html"),
                        selfcontained = TRUE)

# 10.3 Create an Excel report
wb_sample_sharing <- createWorkbook()

# Add raw sharing counts
sharing_matrix_wide <- pivot_wider(sample_sharing_df %>% select(Sample1, Sample2, shared_peptides),
                                   names_from = Sample2,
                                   values_from = shared_peptides)

# Add Jaccard index
jaccard_matrix_wide <- pivot_wider(sample_sharing_df %>% select(Sample1, Sample2, jaccard_index),
                                   names_from = Sample2,
                                   values_from = jaccard_index)

# Add worksheets
addWorksheet(wb_sample_sharing, "Shared Peptide Counts")
writeData(wb_sample_sharing, "Shared Peptide Counts", sharing_matrix_wide, startRow = 1, startCol = 1)
addStyle(wb_sample_sharing, "Shared Peptide Counts", headerStyle, rows = 1, cols = 1:ncol(sharing_matrix_wide))

addWorksheet(wb_sample_sharing, "Jaccard Index")
writeData(wb_sample_sharing, "Jaccard Index", jaccard_matrix_wide, startRow = 1, startCol = 1)
addStyle(wb_sample_sharing, "Jaccard Index", headerStyle, rows = 1, cols = 1:ncol(jaccard_matrix_wide))

# Also add the long format data
addWorksheet(wb_sample_sharing, "All Pairs")
writeData(wb_sample_sharing, "All Pairs", sample_sharing_df, startRow = 1, startCol = 1)
addStyle(wb_sample_sharing, "All Pairs", headerStyle, rows = 1, cols = 1:ncol(sample_sharing_df))

# Save workbook
saveWorkbook(wb_sample_sharing, paste0(excel_dir, "/sample_sample_sharing.xlsx"), overwrite = TRUE)

# 10.4 Perform hierarchical clustering to identify similar sample groups
# Convert the Jaccard matrix to a distance matrix
jaccard_dist <- 1 - as.matrix(sample_sharing_matrix) / (outer(diag(sample_sharing_matrix), diag(sample_sharing_matrix), "+") - as.matrix(sample_sharing_matrix))
# Set diagonal elements to 0 (same sample = no distance)
diag(jaccard_dist) <- 0

# Create a hierarchical clustering
hc <- hclust(as.dist(jaccard_dist), method = "complete")

# Create a dendrogram plot
p_dendro <- plot_ly(
  type = "scatter",
  mode = "markers",
  marker = list(opacity = 0) # Hide markers
)

# Add dendrogram
dendro <- ggdendro::dendro_data(hc)
segments <- dendro$segments

for (i in seq_len(nrow(segments))) {
  p_dendro <- p_dendro %>% add_trace(
    x = c(segments$x[i], segments$xend[i]),
    y = c(segments$y[i], segments$yend[i]),
    mode = "lines",
    line = list(color = "black"),
    showlegend = FALSE,
    hoverinfo = "none"
  )
}

# Label the samples
labels <- dendro$labels
p_dendro <- p_dendro %>% add_trace(
  x = labels$x,
  y = labels$y - 0.1, # Slightly adjust labels down
  text = hc$labels[labels$label],
  textposition = "bottom center",
  mode = "text",
  showlegend = FALSE
)

p_dendro <- p_dendro %>% layout(
  title = "Sample Clustering Based on Peptide Sharing",
  xaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE),
  yaxis = list(title = "Distance", zeroline = FALSE)
)

# Save as interactive HTML
htmlwidgets::saveWidget(p_dendro,
                        paste0(viz_dir, "/interactive_sample_clustering.html"),
                        selfcontained = TRUE)

################################################################################
# 11. Intensity Normalization using Spike-in Peptide
################################################################################

# 11.1 First, check if the input data contains intensity information
has_intensity <- FALSE
intensity_column <- NULL

# Check for common intensity column names
possible_intensity_cols <- c(
  "Intensity", "MS1Intensity", "PeptideIntensity", "Area", "PeakArea", 
  "Abundance", "MS1Area", "PrecursorArea", "PrecursorIntensity",
  "Peak Area", "MS1 Intensity", "Precursor Area", "Precursor Intensity"
)

for (col in possible_intensity_cols) {
  if (col %in% colnames(peptides_data)) {
    has_intensity <- TRUE
    intensity_column <- col
    break
  }
}

if (!has_intensity) {
  cat("WARNING: No intensity column found in the data. Skipping intensity normalization.\n")
  cat("If you have intensity data, please ensure it's included in the input file with one of these column names:\n")
  cat(paste(possible_intensity_cols, collapse=", "), "\n")
} else {
  cat("Found intensity data in column:", intensity_column, "\n")
  
  # Create output directory for intensity analysis
  intensity_dir <- file.path(main_output_dir, "intensity_analysis")
  dir.create(intensity_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 11.2 Define the spike-in peptide sequence and samples
  # Using the provided spike-in peptide sequence
  spike_in_sequence <- "YGEEVKEFL"
  
  # Define the spike-in sample and its control
  spike_sample <- "51S"  # Sample with 10ng of synthetic peptide
  control_sample <- "51"  # Control sample with endogenous levels
  
  # Check if these samples exist in the data
  if (!(spike_sample %in% unique_samples) || !(control_sample %in% unique_samples)) {
    cat("WARNING: Spike sample", spike_sample, "or control sample", control_sample, "not found in the data.\n")
    cat("Available samples:", paste(unique_samples, collapse=", "), "\n")
    cat("Skipping spike-based normalization and using TIC normalization instead.\n")
    
    # If spike samples not found, use TIC normalization as fallback
    total_intensity_by_sample <- peptides_data %>%
      group_by(SampleID) %>%
      summarize(TotalIntensity = sum(!!sym(intensity_column), na.rm = TRUE)) %>%
      arrange(SampleID)
    
    # Calculate TIC normalization factors
    median_tic <- median(total_intensity_by_sample$TotalIntensity)
    tic_norm_factors <- total_intensity_by_sample %>%
      mutate(norm_factor = median_tic / TotalIntensity)
    
    # Apply normalization to peptides
    peptide_intensities <- peptides_data %>%
      select(SampleID, Peptide, !!sym(intensity_column)) %>%
      left_join(tic_norm_factors %>% select(SampleID, norm_factor), by = "SampleID") %>%
      mutate(normalized_intensity = !!sym(intensity_column) * norm_factor)
    
    cat("Applied TIC normalization. Results saved to Excel report.\n")
    
    # Create visualizations for TIC normalization
    # Plot original vs normalized total intensities
    tic_viz_data <- peptide_intensities %>%
      group_by(SampleID) %>%
      summarize(
        original_total = sum(!!sym(intensity_column), na.rm = TRUE),
        normalized_total = sum(normalized_intensity, na.rm = TRUE)
      ) %>%
      pivot_longer(cols = c(original_total, normalized_total),
                   names_to = "type",
                   values_to = "intensity")
    
    p_tic_norm <- ggplot(tic_viz_data, aes(x = SampleID, y = intensity, fill = type)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("original_total" = "darkred", "normalized_total" = "darkblue"),
                        labels = c("original_total" = "Original", "normalized_total" = "Normalized")) +
      labs(
        title = "TIC Normalization Effect",
        x = "Sample",
        y = "Total Intensity",
        fill = "Type"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Convert to interactive plotly plot  
    p_tic_norm_interactive <- ggplotly(p_tic_norm)
    
    # Save as interactive HTML
    htmlwidgets::saveWidget(p_tic_norm_interactive,
                            paste0(viz_dir, "/interactive_tic_normalization.html"),
                            selfcontained = TRUE)
    
    # Create Excel report for normalized intensities
    wb_norm <- createWorkbook()
    
    # Add normalization factors
    addWorksheet(wb_norm, "Normalization Factors")
    writeData(wb_norm, "Normalization Factors", tic_norm_factors, startRow = 1, startCol = 1)
    addStyle(wb_norm, "Normalization Factors", headerStyle, rows = 1, cols = 1:ncol(tic_norm_factors))
    
    # Add all normalized peptide intensities (top 5000 rows only to avoid Excel limits)
    normalized_subset <- peptide_intensities %>%
      head(5000)
    
    addWorksheet(wb_norm, "Normalized Intensities")
    writeData(wb_norm, "Normalized Intensities", normalized_subset, startRow = 1, startCol = 1)
    addStyle(wb_norm, "Normalized Intensities", headerStyle, rows = 1, cols = 1:ncol(normalized_subset))
    
    # Save workbook
    saveWorkbook(wb_norm, paste0(excel_dir, "/normalized_intensities_tic.xlsx"), overwrite = TRUE)
    
  } else {
    # 11.3 Extract spike-in peptide data
    spike_peptide_data <- peptides_data %>%
      filter(Peptide == spike_in_sequence) %>%
      select(SampleID, Peptide, !!sym(intensity_column))
    
    if (nrow(spike_peptide_data) == 0) {
      cat("WARNING: Spike-in peptide", spike_in_sequence, "not found in the data.\n")
      cat("Check the sequence or try using partial matching.\n")
      
      # Try partial matching if exact match fails
      cat("Attempting partial matching for spike peptide...\n")
      spike_candidates <- peptides_data %>%
        filter(SampleID %in% c(spike_sample, control_sample)) %>%
        group_by(Peptide) %>%
        summarize(count = n()) %>%
        filter(count == 2) %>%  # Found in both spike and control samples
        pull(Peptide)
      
      cat("Possible spike candidates:", length(spike_candidates), "\n")
      if (length(spike_candidates) > 0) {
        cat("First few candidates:", head(spike_candidates), "\n")
        cat("Using first candidate for normalization (update script with correct sequence if needed).\n")
        
        # Use the first candidate as placeholder
        spike_in_sequence <- spike_candidates[1]
        
        spike_peptide_data <- peptides_data %>%
          filter(Peptide == spike_in_sequence) %>%
          select(SampleID, Peptide, !!sym(intensity_column))
      }
    }
    
    if (nrow(spike_peptide_data) > 0) {
      cat("Found spike-in peptide data:", nrow(spike_peptide_data), "rows\n")
      
      # 11.4 Calculate normalization factors based on spike-in
      spike_in_intensities <- spike_peptide_data %>%
        filter(SampleID %in% c(spike_sample, control_sample)) %>%
        select(SampleID, !!sym(intensity_column))
      
      print(spike_in_intensities)
      
      # Calculate amount of synthetic spike (10ng minus endogenous level)
      spike_intensity <- spike_in_intensities %>%
        filter(SampleID == spike_sample) %>%
        pull(!!sym(intensity_column))
      
      control_intensity <- spike_in_intensities %>%
        filter(SampleID == control_sample) %>%
        pull(!!sym(intensity_column))
      
      if (length(spike_intensity) == 0 || length(control_intensity) == 0) {
        cat("WARNING: Spike peptide not found in both spike and control samples.\n")
        cat("Falling back to TIC normalization.\n")
        
        # Fall back to TIC normalization
        total_intensity_by_sample <- peptides_data %>%
          group_by(SampleID) %>%
          summarize(TotalIntensity = sum(!!sym(intensity_column), na.rm = TRUE)) %>%
          arrange(SampleID)
        
        # Calculate TIC normalization factors
        median_tic <- median(total_intensity_by_sample$TotalIntensity)
        tic_norm_factors <- total_intensity_by_sample %>%
          mutate(norm_factor = median_tic / TotalIntensity)
        
        # Apply normalization to peptides
        peptide_intensities <- peptides_data %>%
          select(SampleID, Peptide, !!sym(intensity_column)) %>%
          left_join(tic_norm_factors %>% select(SampleID, norm_factor), by = "SampleID") %>%
          mutate(normalized_intensity = !!sym(intensity_column) * norm_factor)
        
        # Create visualizations for TIC normalization (same as above)
        tic_viz_data <- peptide_intensities %>%
          group_by(SampleID) %>%
          summarize(
            original_total = sum(!!sym(intensity_column), na.rm = TRUE),
            normalized_total = sum(normalized_intensity, na.rm = TRUE)
          ) %>%
          pivot_longer(cols = c(original_total, normalized_total),
                       names_to = "type",
                       values_to = "intensity")
        
        p_tic_norm <- ggplot(tic_viz_data, aes(x = SampleID, y = intensity, fill = type)) +
          geom_bar(stat = "identity", position = "dodge") +
          scale_fill_manual(values = c("original_total" = "darkred", "normalized_total" = "darkblue"),
                            labels = c("original_total" = "Original", "normalized_total" = "Normalized")) +
          labs(
            title = "TIC Normalization Effect (Fallback Method)",
            x = "Sample",
            y = "Total Intensity",
            fill = "Type"
          ) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        # Convert to interactive plotly plot  
        p_tic_norm_interactive <- ggplotly(p_tic_norm)
        
        # Save as interactive HTML
        htmlwidgets::saveWidget(p_tic_norm_interactive,
                                paste0(viz_dir, "/interactive_tic_normalization.html"),
                                selfcontained = TRUE)
        
        # Create Excel report for normalized intensities (same as above)
        wb_norm <- createWorkbook()
        
        # Add normalization factors
        addWorksheet(wb_norm, "Normalization Factors")
        writeData(wb_norm, "Normalization Factors", tic_norm_factors, startRow = 1, startCol = 1)
        addStyle(wb_norm, "Normalization Factors", headerStyle, rows = 1, cols = 1:ncol(tic_norm_factors))
        
        # Add normalized intensities
        normalized_subset <- peptide_intensities %>%
          head(5000)
        
        addWorksheet(wb_norm, "Normalized Intensities")
        writeData(wb_norm, "Normalized Intensities", normalized_subset, startRow = 1, startCol = 1)
        addStyle(wb_norm, "Normalized Intensities", headerStyle, rows = 1, cols = 1:ncol(normalized_subset))
        
        # Save workbook
        saveWorkbook(wb_norm, paste0(excel_dir, "/normalized_intensities_tic.xlsx"), overwrite = TRUE)
        
      } else {
        synthetic_intensity <- spike_intensity - control_intensity
        
        if (synthetic_intensity <= 0) {
          cat("WARNING: Synthetic intensity calculation resulted in non-positive value.\n")
          cat("Spike intensity:", spike_intensity, "Control intensity:", control_intensity, "\n")
          cat("Using absolute spike intensity instead.\n")
          synthetic_intensity <- spike_intensity
        }
        
        # Normalization factor = 10ng / synthetic intensity
        # This gives us ng per unit of intensity
        norm_factor <- 10 / synthetic_intensity
        
        cat("Spike-in intensity:", spike_intensity, "\n")
        cat("Control intensity:", control_intensity, "\n")
        cat("Synthetic peptide intensity:", synthetic_intensity, "\n")
        cat("Normalization factor (ng/intensity):", norm_factor, "\n")
        
        # 11.5 Apply normalization to all peptides
        # Get all peptide intensities
        peptide_intensities <- peptides_data %>%
          select(SampleID, Peptide, !!sym(intensity_column)) %>%
          # Convert intensity to normalized amount (ng)
          mutate(normalized_amount = !!sym(intensity_column) * norm_factor)
        
        # 11.6 Create visualizations and reports for normalized data
        # Summary of peptide amounts by sample
        sample_amounts <- peptide_intensities %>%
          group_by(SampleID) %>%
          summarize(
            total_intensity = sum(!!sym(intensity_column), na.rm = TRUE),
            total_amount_ng = sum(normalized_amount, na.rm = TRUE),
            peptide_count = n_distinct(Peptide),
            avg_amount_per_peptide = total_amount_ng / peptide_count
          )
        
        # Visualize normalized amounts by sample
        p_norm_amount <- ggplot(sample_amounts, aes(x = SampleID, y = total_amount_ng,
                                                    text = paste0("Sample: ", SampleID,
                                                                  "<br>Total peptide amount: ", round(total_amount_ng, 2), " ng",
                                                                  "<br>Peptide count: ", peptide_count,
                                                                  "<br>Avg peptide amount: ", round(avg_amount_per_peptide, 4), " ng"))) +
          geom_bar(stat = "identity", fill = "darkgreen") +
          labs(
            title = "Normalized Peptide Amounts by Sample",
            x = "Sample",
            y = "Total Peptide Amount (ng)"
          ) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        # Convert to interactive plotly plot
        p_norm_amount_interactive <- ggplotly(p_norm_amount, tooltip = "text")
        
        # Save as interactive HTML
        htmlwidgets::saveWidget(p_norm_amount_interactive,
                                paste0(viz_dir, "/interactive_normalized_amounts.html"),
                                selfcontained = TRUE)
        
        # Create Excel reports for normalized data
        wb_norm <- createWorkbook()
        
        # Add sample summary
        addWorksheet(wb_norm, "Sample Summary")
        writeData(wb_norm, "Sample Summary", sample_amounts, startRow = 1, startCol = 1)
        addStyle(wb_norm, "Sample Summary", headerStyle, rows = 1, cols = 1:ncol(sample_amounts))
        
        # Add normalization parameters
        norm_info <- data.frame(
          Parameter = c("Spike Peptide", "Spike Sample", "Control Sample", 
                        "Spike Intensity", "Control Intensity", "Synthetic Intensity",
                        "Normalization Factor (ng/intensity)"),
          Value = c(spike_in_sequence, spike_sample, control_sample,
                    spike_intensity, control_intensity, synthetic_intensity,
                    norm_factor)
        )
        
        addWorksheet(wb_norm, "Normalization Info")
        writeData(wb_norm, "Normalization Info", norm_info, startRow = 1, startCol = 1)
        addStyle(wb_norm, "Normalization Info", headerStyle, rows = 1, cols = 1:ncol(norm_info))
        
        # Add all normalized peptide intensities (limit rows to prevent Excel from crashing)
        addWorksheet(wb_norm, "All Normalized Peptides")
        writeData(wb_norm, "All Normalized Peptides", head(peptide_intensities, 10000), startRow = 1, startCol = 1)
        addStyle(wb_norm, "All Normalized Peptides", headerStyle, rows = 1, cols = 1:ncol(peptide_intensities))
        
        # Top abundance peptides overall
        top_peptides <- peptide_intensities %>%
          group_by(Peptide) %>%
          summarize(
            avg_intensity = mean(!!sym(intensity_column), na.rm = TRUE),
            avg_amount_ng = mean(normalized_amount, na.rm = TRUE),
            sample_count = n_distinct(SampleID),
            samples = paste(sort(unique(SampleID)), collapse = ", ")
          ) %>%
          arrange(desc(avg_amount_ng))
        
        addWorksheet(wb_norm, "Top Peptides")
        writeData(wb_norm, "Top Peptides", top_peptides, startRow = 1, startCol = 1)
        addStyle(wb_norm, "Top Peptides", headerStyle, rows = 1, cols = 1:ncol(top_peptides))
        
        # Save workbook
        saveWorkbook(wb_norm, paste0(excel_dir, "/normalized_peptide_intensities.xlsx"), overwrite = TRUE)
        
        # 11.7 Analyze sharing patterns with intensity taken into account
        # For robustly shared peptides, look at their intensity distribution
        intensity_by_sharing <- peptide_intensities %>%
          inner_join(peptide_report %>% select(Peptide, sample_count), by = "Peptide") %>%
          group_by(sample_count) %>%
          summarize(
            avg_intensity = mean(!!sym(intensity_column), na.rm = TRUE),
            avg_amount_ng = mean(normalized_amount, na.rm = TRUE),
            peptide_count = n_distinct(Peptide)
          )
        
        # Visualize intensity by sharing level
        p_intensity_sharing <- ggplot(intensity_by_sharing, aes(x = sample_count, y = avg_amount_ng,
                                                                text = paste0("Samples: ", sample_count,
                                                                              "<br>Avg peptide amount: ", round(avg_amount_ng, 4), " ng",
                                                                              "<br>Peptide count: ", peptide_count))) +
          geom_point(aes(size = peptide_count), color = "blue", alpha = 0.7) +
          geom_smooth(method = "loess", se = TRUE) +
          labs(
            title = "Average Peptide Amount by Sharing Level",
            x = "Number of Samples",
            y = "Average Peptide Amount (ng)",
            size = "Peptide Count"
          ) +
          theme_minimal()
        
        # Convert to interactive plotly plot
        p_intensity_sharing_interactive <- ggplotly(p_intensity_sharing, tooltip = "text")
        
        # Save as interactive HTML
        htmlwidgets::saveWidget(p_intensity_sharing_interactive,
                                paste0(viz_dir, "/interactive_intensity_by_sharing.html"),
                                selfcontained = TRUE)
        
        # 11.8 Analyze robustly shared peptides with intensity data
        # Create a special report for high-intensity shared peptides
        if (exists("robustly_shared_peptides")) {
          for (i in 1:length(robustly_shared_peptides)) {
            threshold <- names(robustly_shared_peptides)[i]
            robust_set <- robustly_shared_peptides[[threshold]]
            
            # Get intensity data for these peptides
            robust_intensity <- peptide_intensities %>%
              filter(Peptide %in% robust_set$Peptide) %>%
              group_by(Peptide) %>%
              summarize(
                avg_intensity = mean(!!sym(intensity_column), na.rm = TRUE),
                avg_amount_ng = mean(normalized_amount, na.rm = TRUE),
                sample_count = n_distinct(SampleID)
              ) %>%
              arrange(desc(avg_amount_ng))
            
            # Combine with original robust data
            robust_with_intensity <- robust_set %>%
              left_join(robust_intensity, by = "Peptide")
            
            # Clean up threshold name for Excel sheet name
            sheet_name <- paste0("Intensity_", i)
            
            # Add worksheet to the intensity workbook
            addWorksheet(wb_norm, sheet_name)
            writeData(wb_norm, sheet_name, robust_with_intensity, startRow = 1, startCol = 1)
            addStyle(wb_norm, sheet_name, headerStyle, rows = 1, cols = 1:ncol(robust_with_intensity))
          }
          
          # Create visualization showing the top shared peptides by intensity
          # Get the most abundant peptides from the highest sharing threshold
          high_sharing_threshold <- names(robustly_shared_peptides)[length(robustly_shared_peptides)]
          high_robust_set <- robustly_shared_peptides[[high_sharing_threshold]]
          
          if (nrow(high_robust_set) > 0) {
            top_shared_by_amount <- peptide_intensities %>%
              filter(Peptide %in% high_robust_set$Peptide) %>%
              group_by(Peptide) %>%
              summarize(
                avg_amount_ng = mean(normalized_amount, na.rm = TRUE),
                sample_count = n_distinct(SampleID)
              ) %>%
              arrange(desc(avg_amount_ng)) %>%
              head(20) # Top 20 peptides
            
            # Create a bar plot
            p_top_shared <- ggplot(top_shared_by_amount, aes(x = reorder(Peptide, avg_amount_ng), y = avg_amount_ng,
                                                             text = paste0("Peptide: ", Peptide,
                                                                           "<br>Avg amount: ", round(avg_amount_ng, 4), " ng",
                                                                           "<br>Found in ", sample_count, " samples"))) +
              geom_bar(stat = "identity", fill = "purple") +
              coord_flip() +
              labs(
                title = paste0("Top Abundant Peptides from ", high_sharing_threshold),
                x = "Peptide",
                y = "Average Amount (ng)"
              ) +
              theme_minimal()
            
            # Convert to interactive plotly plot
            p_top_shared_interactive <- ggplotly(p_top_shared, tooltip = "text")
            
            # Save as interactive HTML
            htmlwidgets::saveWidget(p_top_shared_interactive,
                                    paste0(viz_dir, "/interactive_top_shared_by_amount.html"),
                                    selfcontained = TRUE)
          }
        }
        
        # Save the updated workbook
        saveWorkbook(wb_norm, paste0(excel_dir, "/normalized_peptide_intensities.xlsx"), overwrite = TRUE)
      }
    } else {
      
      cat("WARNING: Could not find or match spike-in peptide. Falling back to TIC normalization.\n")
      
      # Fall back to TIC normalization (same code as in the earlier condition)
      total_intensity_by_sample <- peptides_data %>%
        group_by(SampleID) %>%
        summarize(TotalIntensity = sum(!!sym(intensity_column), na.rm = TRUE)) %>%
        arrange(SampleID)
      
      # Calculate TIC normalization factors
      median_tic <- median(total_intensity_by_sample$TotalIntensity)
      tic_norm_factors <- total_intensity_by_sample %>%
        mutate(norm_factor = median_tic / TotalIntensity)
      
      # Apply normalization to peptides
      peptide_intensities <- peptides_data %>%
        select(SampleID, Peptide, !!sym(intensity_column)) %>%
        left_join(tic_norm_factors %>% select(SampleID, norm_factor), by = "SampleID") %>%
        mutate(normalized_intensity = !!sym(intensity_column) * norm_factor,
               normalized_amount = normalized_intensity) # For consistency with other code paths
      
      cat("Applied TIC normalization as fallback. Results saved to Excel report.\n")
      
      # Create visualizations for TIC normalization
      # Plot original vs normalized total intensities
      tic_viz_data <- peptide_intensities %>%
        group_by(SampleID) %>%
        summarize(
          original_total = sum(!!sym(intensity_column), na.rm = TRUE),
          normalized_total = sum(normalized_intensity, na.rm = TRUE)
        ) %>%
        pivot_longer(cols = c(original_total, normalized_total),
                     names_to = "type",
                     values_to = "intensity")
      
      p_tic_norm <- ggplot(tic_viz_data, aes(x = SampleID, y = intensity, fill = type)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("original_total" = "darkred", "normalized_total" = "darkblue"),
                          labels = c("original_total" = "Original", "normalized_total" = "Normalized")) +
        labs(
          title = "TIC Normalization Effect (Fallback Method)",
          x = "Sample",
          y = "Total Intensity",
          fill = "Type"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Convert to interactive plotly plot  
      p_tic_norm_interactive <- ggplotly(p_tic_norm)
      
      # Save as interactive HTML
      htmlwidgets::saveWidget(p_tic_norm_interactive,
                              paste0(viz_dir, "/interactive_tic_normalization.html"),
                              selfcontained = TRUE)
      
      # Calculate sample summaries for reporting
      sample_amounts <- peptide_intensities %>%
        group_by(SampleID) %>%
        summarize(
          total_intensity = sum(!!sym(intensity_column), na.rm = TRUE),
          total_amount_ng = sum(normalized_amount, na.rm = TRUE),
          peptide_count = n_distinct(Peptide),
          avg_amount_per_peptide = total_amount_ng / peptide_count
        )
      
      # Create visualization of normalized amounts
      p_norm_amount <- ggplot(sample_amounts, aes(x = SampleID, y = total_amount_ng,
                                                  text = paste0("Sample: ", SampleID,
                                                                "<br>Total peptide amount: ", round(total_amount_ng, 2),
                                                                "<br>Peptide count: ", peptide_count,
                                                                "<br>Avg peptide amount: ", round(avg_amount_per_peptide, 4)))) +
        geom_bar(stat = "identity", fill = "darkgreen") +
        labs(
          title = "Normalized Peptide Amounts by Sample (TIC Method)",
          x = "Sample",
          y = "Total Peptide Amount"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Convert to interactive plotly plot
      p_norm_amount_interactive <- ggplotly(p_norm_amount, tooltip = "text")
      
      # Save as interactive HTML
      htmlwidgets::saveWidget(p_norm_amount_interactive,
                              paste0(viz_dir, "/interactive_normalized_amounts.html"),
                              selfcontained = TRUE)
      
      # Create Excel reports for normalized data
      wb_norm <- createWorkbook()
      
      # Add sample summary
      addWorksheet(wb_norm, "Sample Summary")
      writeData(wb_norm, "Sample Summary", sample_amounts, startRow = 1, startCol = 1)
      addStyle(wb_norm, "Sample Summary", headerStyle, rows = 1, cols = 1:ncol(sample_amounts))
      
      # Add normalization information
      norm_info <- data.frame(
        Parameter = c("Normalization Method", "Median TIC", "Spike Peptide"),
        Value = c("TIC Normalization (Fallback)", median_tic, "Not found")
      )
      
      addWorksheet(wb_norm, "Normalization Info")
      writeData(wb_norm, "Normalization Info", norm_info, startRow = 1, startCol = 1)
      addStyle(wb_norm, "Normalization Info", headerStyle, rows = 1, cols = 1:ncol(norm_info))
      
      # Add all normalized peptide intensities
      addWorksheet(wb_norm, "All Normalized Peptides")
      writeData(wb_norm, "All Normalized Peptides", head(peptide_intensities, 10000), startRow = 1, startCol = 1)
      addStyle(wb_norm, "All Normalized Peptides", headerStyle, rows = 1, cols = 1:ncol(peptide_intensities))
      
      # Top abundance peptides overall
      top_peptides <- peptide_intensities %>%
        group_by(Peptide) %>%
        summarize(
          avg_intensity = mean(!!sym(intensity_column), na.rm = TRUE),
          avg_amount_ng = mean(normalized_amount, na.rm = TRUE),
          sample_count = n_distinct(SampleID),
          samples = paste(sort(unique(SampleID)), collapse = ", ")
        ) %>%
        arrange(desc(avg_amount_ng))
      
      addWorksheet(wb_norm, "Top Peptides")
      writeData(wb_norm, "Top Peptides", top_peptides, startRow = 1, startCol = 1)
      addStyle(wb_norm, "Top Peptides", headerStyle, rows = 1, cols = 1:ncol(top_peptides))
      
      # Save workbook
      saveWorkbook(wb_norm, paste0(excel_dir, "/normalized_intensities_tic.xlsx"), overwrite = TRUE)
      
      # 11.7 Analyze sharing patterns with intensity taken into account
      # For robustly shared peptides, look at their intensity distribution
      intensity_by_sharing <- peptide_intensities %>%
        inner_join(peptide_report %>% select(Peptide, sample_count), by = "Peptide") %>%
        group_by(sample_count) %>%
        summarize(
          avg_intensity = mean(!!sym(intensity_column), na.rm = TRUE),
          avg_amount_ng = mean(normalized_amount, na.rm = TRUE),
          peptide_count = n_distinct(Peptide)
        )
      
      # Visualize intensity by sharing level
      p_intensity_sharing <- ggplot(intensity_by_sharing, aes(x = sample_count, y = avg_amount_ng,
                                                              text = paste0("Samples: ", sample_count,
                                                                            "<br>Avg peptide amount: ", round(avg_amount_ng, 4), " ng",
                                                                            "<br>Peptide count: ", peptide_count))) +
        geom_point(aes(size = peptide_count), color = "blue", alpha = 0.7) +
        geom_smooth(method = "loess", se = TRUE) +
        labs(
          title = "Average Peptide Amount by Sharing Level",
          x = "Number of Samples",
          y = "Average Peptide Amount (ng)",
          size = "Peptide Count"
        ) +
        theme_minimal()
      
      # Convert to interactive plotly plot
      p_intensity_sharing_interactive <- ggplotly(p_intensity_sharing, tooltip = "text")
      
      # Save as interactive HTML
      htmlwidgets::saveWidget(p_intensity_sharing_interactive,
                              paste0(viz_dir, "/interactive_intensity_by_sharing.html"),
                              selfcontained = TRUE)
      
      # 11.8 Analyze robustly shared peptides with intensity data
      # Create a special report for high-intensity shared peptides
      for (i in 1:length(robustly_shared_peptides)) {
        threshold <- names(robustly_shared_peptides)[i]
        robust_set <- robustly_shared_peptides[[threshold]]
        
        # Get intensity data for these peptides
        robust_intensity <- peptide_intensities %>%
          filter(Peptide %in% robust_set$Peptide) %>%
          group_by(Peptide) %>%
          summarize(
            avg_intensity = mean(!!sym(intensity_column), na.rm = TRUE),
            avg_amount_ng = mean(normalized_amount, na.rm = TRUE),
            sample_count = n_distinct(SampleID)
          ) %>%
          arrange(desc(avg_amount_ng))
        
        # Combine with original robust data
        robust_with_intensity <- robust_set %>%
          left_join(robust_intensity, by = "Peptide")
        
        # Clean up threshold name for Excel sheet name
        sheet_name <- paste0("Intensity_", i)
        
        # Add worksheet to the intensity workbook
        addWorksheet(wb_norm, sheet_name)
        writeData(wb_norm, sheet_name, robust_with_intensity, startRow = 1, startCol = 1)
        addStyle(wb_norm, sheet_name, headerStyle, rows = 1, cols = 1:ncol(robust_with_intensity))
      }
      
      # Save the updated workbook
      saveWorkbook(wb_norm, paste0(excel_dir, "/normalized_peptide_intensities.xlsx"), overwrite = TRUE)
    }
  }
}


################################################################################
# 12. Create an updated comprehensive HTML dashboard with all visualizations
################################################################################

# This should replace the original dashboard creation section (section 7)
# It includes all the new visualizations we've created

# List of all visualization files to include
viz_files <- list(
  "sample_count" = list(
    file = "interactive_sample_count_distribution.html",
    title = "Peptide Distribution by Sample Count",
    exists = TRUE  # Original visualization
  ),
  "unique_shared" = list(
    file = "interactive_unique_vs_shared.html",
    title = "Unique vs Shared Peptides",
    exists = TRUE  # Original visualization
  ),
  "length_dist" = list(
    file = "interactive_length_distribution.html",
    title = "Peptide Length Distribution",
    exists = TRUE  # Original visualization
  ),
  "avg_sharing" = list(
    file = "interactive_avg_sharing_by_length.html",
    title = "Average Sharing by Peptide Length",
    exists = TRUE  # Original visualization
  ),
  "gene_heatmap" = list(
    file = "interactive_gene_heatmap.html",
    title = "Peptides by Gene Across Samples",
    exists = TRUE  # Original visualization
  ),
  "intersections" = list(
    file = "interactive_sample_intersections.html",
    title = "Sample Intersections",
    exists = num_samples >= 2 && num_samples <= 4  # May not exist
  ),
  "robust_peptides" = list(
    file = "interactive_robust_peptides.html",
    title = "Robustly Shared Peptides",
    exists = TRUE  # New visualization
  ),
  "sample_coverage" = list(
    file = "interactive_sample_core_coverage.html",
    title = "Sample Coverage of Core Peptidome",
    exists = TRUE  # New visualization
  ),
  "sample_similarity" = list(
    file = "interactive_sample_similarity.html",
    title = "Sample Similarity (Jaccard Index)",
    exists = TRUE  # New visualization
  ),
  "sample_clustering" = list(
    file = "interactive_sample_clustering.html",
    title = "Sample Clustering",
    exists = TRUE  # New visualization
  ),
  "normalized_amounts" = list(
    file = "interactive_normalized_amounts.html",
    title = "Normalized Peptide Amounts",
    exists = FALSE  # May not exist if no intensity data
  ),
  "intensity_sharing" = list(
    file = "interactive_intensity_by_sharing.html",
    title = "Peptide Amount by Sharing Level",
    exists = FALSE  # May not exist if no intensity data
  ),
  "tic_normalization" = list(
    file = "interactive_tic_normalization.html",
    title = "TIC Normalization Effect",
    exists = FALSE  # May not exist
  )
)

# Check which files actually exist in the directory
for (i in 1:length(viz_files)) {
  viz_name <- names(viz_files)[i]
  file_path <- file.path(viz_dir, viz_files[[viz_name]]$file)
  
  if (viz_files[[viz_name]]$exists) {
    # For files that should exist, double-check they actually do
    viz_files[[viz_name]]$exists <- file.exists(file_path)
  } else {
    # For files that might exist, check if they do
    viz_files[[viz_name]]$exists <- file.exists(file_path)
  }
}

# Create HTML footer
html_footer <- '
  <div class="footer">
    <p>Analysis generated on %s | Input file: %s</p>
    <p>Peptide Sharing Analysis Tool v1.0</p>
  </div>
</body>
</html>
'

html_header <- '
<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>Peptide Sharing Analysis Dashboard</title>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <style>
    body {
      font-family: Arial, sans-serif;
      margin: 0;
      padding: 20px;
      background-color: #f5f5f5;
    }
    .dashboard-container {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
      gap: 20px;
      margin-bottom: 30px;
    }
    .dashboard-item {
      background-color: white;
      border-radius: 5px;
      box-shadow: 0 2px 5px rgba(0,0,0,0.1);
      padding: 15px;
      height: 450px;
    }
    .dashboard-wide {
      grid-column: 1 / -1;
      height: 550px;
    }
    iframe {
      width: 100%%;
      height: 100%%;
      border: none;
    }
    h1 {
      color: #333;
      text-align: center;
      margin-bottom: 30px;
    }
    h2 {
      color: #444;
      margin-top: 30px;
      margin-bottom: 20px;
      border-bottom: 1px solid #ddd;
      padding-bottom: 10px;
    }
    h3 {
      color: #555;
      margin-top: 0;
      margin-bottom: 15px;
      border-bottom: 1px solid #eee;
      padding-bottom: 8px;
    }
    .summary {
      background-color: white;
      border-radius: 5px;
      box-shadow: 0 2px 5px rgba(0,0,0,0.1);
      padding: 20px;
      margin-bottom: 20px;
    }
    .footer {
      text-align: center;
      margin-top: 30px;
      color: #777;
      font-size: 0.9em;
    }
    .tab-container {
      display: flex;
      flex-wrap: wrap;
      gap: 0px;
      margin-bottom: 20px;
    }
    .tab {
      padding: 10px 20px;
      background-color: #e0e0e0;
      border-radius: 5px 5px 0 0;
      cursor: pointer;
      border: 1px solid #ccc;
      border-bottom: none;
    }
    .tab.active {
      background-color: white;
      border-bottom: 1px solid white;
      margin-bottom: -1px;
      font-weight: bold;
    }
    .tab-content {
      display: none;
      padding: 20px;
      background-color: white;
      border: 1px solid #ccc;
      border-radius: 0 5px 5px 5px;
    }
    .tab-content.active {
      display: block;
    }
  </style>
  <script>
    function openTab(evt, tabName) {
      // Hide all tab content
      var tabcontent = document.getElementsByClassName("tab-content");
      for (var i = 0; i < tabcontent.length; i++) {
        tabcontent[i].className = tabcontent[i].className.replace(" active", "");
      }

      // Remove active class from all tabs
      var tabs = document.getElementsByClassName("tab");
      for (var i = 0; i < tabs.length; i++) {
        tabs[i].className = tabs[i].className.replace(" active", "");
      }

      // Show the current tab and add active class
      document.getElementById(tabName).className += " active";
      evt.currentTarget.className += " active";
    }
  </script>
</head>
<body>
  <h1>Peptide Sharing Analysis Dashboard</h1>

  <div class="summary">
    <h2>Analysis Summary</h2>
    <p><strong>Input file:</strong> %s</p>
    <p><strong>Analysis date:</strong> %s</p>
    <p>This dashboard presents an interactive analysis of peptide sharing across samples. The analysis is based on %d unique peptides across %d samples. Out of these peptides, %d (%d%%) are unique to a single sample, while %d (%d%%) are shared across multiple samples.</p>
    <p>The most prevalent peptide length is %dmer, making up %d%% of all peptides. Peptides of length %d amino acids show the highest rate of sharing across samples.</p>
  </div>
  
  <div class="tab-container">
    <button class="tab active" onclick="openTab(event, \'BasicAnalysis\')">Basic Peptide Analysis</button>
    <button class="tab" onclick="openTab(event, \'RobustSharing\')">Robust Sharing Analysis</button>
    <button class="tab" onclick="openTab(event, \'SampleSimilarity\')">Sample Similarity</button>
    <button class="tab" onclick="openTab(event, \'IntensityAnalysis\')">Intensity Analysis</button>
  </div>
'

# Function to create an iframe element for each plot
create_plot_iframe <- function(plot_info, title = NULL) {
  if (is.null(title)) {
    title <- plot_info$title
  }
  
  iframe_template <- '
  <div class="dashboard-item">
    <h3>%s</h3>
    <iframe src="%s"></iframe>
  </div>
  '
  
  # Format the iframe div
  sprintf(iframe_template, title, plot_info$file)
}

# Function to create a wide iframe (full width)
create_wide_plot_iframe <- function(plot_info, title = NULL) {
  if (is.null(title)) {
    title <- plot_info$title
  }
  
  iframe_template <- '
  <div class="dashboard-item dashboard-wide">
    <h3>%s</h3>
    <iframe src="%s"></iframe>
  </div>
  '
  
  # Format the iframe div
  sprintf(iframe_template, title, plot_info$file)
}

# Generate summary stats for the dashboard
total_peptides <- nrow(peptide_report)
unique_peptide_count <- sum(sample_count_dist$n[sample_count_dist$sample_count == 1])
shared_peptide_count <- sum(sample_count_dist$n[sample_count_dist$sample_count > 1])
unique_peptide_pct <- round(unique_peptide_count / total_peptides * 100)
shared_peptide_pct <- round(shared_peptide_count / total_peptides * 100)
most_common_length <- length_dist$`Peptide Length`[which.max(length_dist$n)]
most_common_length_pct <- round(max(length_dist$percentage))
best_sharing_length <- avg_sharing_by_length$`Peptide Length`[which.max(avg_sharing_by_length$avg_samples)]

# Format the header with summary stats
dashboard_header <- sprintf(html_header,
                            basename(combined_file),
                            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                            total_peptides,
                            num_samples,
                            unique_peptide_count,
                            unique_peptide_pct,
                            shared_peptide_count,
                            shared_peptide_pct,
                            most_common_length,
                            most_common_length_pct,
                            best_sharing_length)

# Format the footer with date and input file
dashboard_footer <- sprintf(html_footer, 
                            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                            basename(combined_file))

# Start building the dashboard HTML
dashboard_html <- dashboard_header

# Create the Basic Analysis tab content
dashboard_html <- paste0(dashboard_html, '<div id="BasicAnalysis" class="tab-content active">\n')
dashboard_html <- paste0(dashboard_html, '<div class="dashboard-container">\n')

# Add basic analysis plots
basic_plots <- c("sample_count", "unique_shared", "length_dist", "avg_sharing")
for (plot_name in basic_plots) {
  if (viz_files[[plot_name]]$exists) {
    dashboard_html <- paste0(dashboard_html, create_plot_iframe(viz_files[[plot_name]]))
  }
}

dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close container

# Add gene heatmap as wide plot
if (viz_files[["gene_heatmap"]]$exists) {
  dashboard_html <- paste0(dashboard_html, '<div class="dashboard-container">\n')
  dashboard_html <- paste0(dashboard_html, create_wide_plot_iframe(viz_files[["gene_heatmap"]]))
  dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close container
}

# Add sample intersections if they exist
if (viz_files[["intersections"]]$exists) {
  dashboard_html <- paste0(dashboard_html, '<div class="dashboard-container">\n')
  dashboard_html <- paste0(dashboard_html, create_plot_iframe(viz_files[["intersections"]]))
  dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close container
}

dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close tab content

# Create the Robust Sharing Analysis tab content
dashboard_html <- paste0(dashboard_html, '<div id="RobustSharing" class="tab-content">\n')

# Add robust sharing visualizations
dashboard_html <- paste0(dashboard_html, '<div class="dashboard-container">\n')
if (viz_files[["robust_peptides"]]$exists) {
  dashboard_html <- paste0(dashboard_html, create_plot_iframe(viz_files[["robust_peptides"]]))
}
dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close container

# Add sample core coverage as wide plot
if (viz_files[["sample_coverage"]]$exists) {
  dashboard_html <- paste0(dashboard_html, '<div class="dashboard-container">\n')
  dashboard_html <- paste0(dashboard_html, create_wide_plot_iframe(viz_files[["sample_coverage"]]))
  dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close container
}

dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close tab content

# Create the Sample Similarity tab content
dashboard_html <- paste0(dashboard_html, '<div id="SampleSimilarity" class="tab-content">\n')

# Add sample similarity visualizations
if (viz_files[["sample_similarity"]]$exists) {
  dashboard_html <- paste0(dashboard_html, '<div class="dashboard-container">\n')
  dashboard_html <- paste0(dashboard_html, create_wide_plot_iframe(viz_files[["sample_similarity"]]))
  dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close container
}

if (viz_files[["sample_clustering"]]$exists) {
  dashboard_html <- paste0(dashboard_html, '<div class="dashboard-container">\n')
  dashboard_html <- paste0(dashboard_html, create_wide_plot_iframe(viz_files[["sample_clustering"]]))
  dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close container
}

dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close tab content

# Create the Intensity Analysis tab content
dashboard_html <- paste0(dashboard_html, '<div id="IntensityAnalysis" class="tab-content">\n')

# Add intensity visualizations if they exist
intensity_plots <- c("normalized_amounts", "intensity_sharing", "tic_normalization")
has_intensity_plots <- FALSE

for (plot_name in intensity_plots) {
  if (viz_files[[plot_name]]$exists) {
    has_intensity_plots <- TRUE
    break
  }
}

if (has_intensity_plots) {
  dashboard_html <- paste0(dashboard_html, '<div class="dashboard-container">\n')
  
  for (plot_name in intensity_plots) {
    if (viz_files[[plot_name]]$exists) {
      dashboard_html <- paste0(dashboard_html, create_plot_iframe(viz_files[[plot_name]]))
    }
  }
  
  dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close container
} else {
  dashboard_html <- paste0(dashboard_html, 
                           '<div class="summary"><p>No intensity data was found in the input file. To enable intensity analysis, ensure your data includes intensity values.</p></div>\n')
}

dashboard_html <- paste0(dashboard_html, '</div>\n')  # Close tab content

# Add the footer
dashboard_html <- paste0(dashboard_html, dashboard_footer)

# Write the complete dashboard HTML to a file
cat("Writing dashboard HTML...\n")
writeLines(dashboard_html, paste0(viz_dir, "/peptide_sharing_dashboard.html"))
cat("Dashboard created successfully!\n")

################################################################################
# Print completion message with all new files
################################################################################

cat("\nEnhanced visualization and Excel report generation complete!\n")
cat("Input file:", combined_file, "\n")
cat("Output directory:", main_output_dir, "\n")
cat("Interactive visualizations saved to:", viz_dir, "\n")
cat("Excel reports saved to:", excel_dir, "\n")
cat("The following files were created:\n")

cat("\n1. Original visualizations:\n")
cat("  - interactive_sample_count_distribution.html\n")
cat("  - interactive_unique_vs_shared.html\n")
cat("  - interactive_length_distribution.html\n")
cat("  - interactive_avg_sharing_by_length.html\n")
cat("  - interactive_gene_heatmap.html\n")
if (num_samples >= 2 && num_samples <= 4) {
  cat("  - interactive_sample_intersections.html\n")
}

cat("\n2. New robust sharing analysis:\n")
cat("  - interactive_robust_peptides.html\n")
cat("  - interactive_sample_core_coverage.html\n")
cat("  - robust_shared_peptides.xlsx\n")
cat("  - sample_core_coverage.xlsx\n")

cat("\n3. New sample similarity analysis:\n")
cat("  - interactive_sample_similarity.html\n")
cat("  - interactive_sample_clustering.html\n")
cat("  - sample_sample_sharing.xlsx\n")

has_intensity <- FALSE
for (plot_name in intensity_plots) {
  if (viz_files[[plot_name]]$exists) {
    has_intensity <- TRUE
    break
  }
}

if (has_intensity) {
  cat("\n4. Intensity analysis:\n")
  for (plot_name in intensity_plots) {
    if (viz_files[[plot_name]]$exists) {
      cat("  - ", viz_files[[plot_name]]$file, "\n", sep = "")
    }
  }
  cat("  - normalized_peptide_intensities.xlsx\n")
}

cat("\n5. Comprehensive dashboard:\n")
cat("  - peptide_sharing_dashboard.html\n")

cat("\nComplete!\n")
