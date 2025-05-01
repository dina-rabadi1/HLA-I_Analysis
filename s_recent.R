# R script for interactive peptide sample distribution visualization

# Setting directory
setwd("~/Documents/Github/HLA-I_Analysis/")

# Path to the combined peptides file
combined_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/combined_peptides.tsv"

# Create the main output directory first
main_output_dir <- file.path(dirname(combined_file), "shared_peptide_visualizations")
dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)

# Create subdirectories within main output directory
viz_dir <- file.path(main_output_dir, "visualizations")
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

excel_dir <- file.path(main_output_dir, "excel_reports")
dir.create(excel_dir, recursive = TRUE, showWarnings = FALSE)

# Load required libraries
required_packages <- c("tidyverse", "plotly", "htmlwidgets", "openxlsx", 
                       "UpSetR", "VennDiagram", "RColorBrewer")

# Install and load necessary packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
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
# 7. Create a comprehensive HTML dashboard with all visualizations
################################################################################

# Create a comprehensive HTML dashboard that combines all plots
# First, create a function to wrap plotly objects for the dashboard
wrap_in_div <- function(plotly_obj, title) {
  div_template <- '
  <div class="dashboard-item">
    <h3>%s</h3>
    <div class="plot-container">
      %s
    </div>
  </div>
  '
  
  # Save the plot to a temporary file
  temp_file <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(plotly_obj, temp_file, selfcontained = TRUE)
  
  # Read the saved plot HTML
  plot_html <- readLines(temp_file)
  
  # Extract just the plot div (not the full HTML document)
  plot_start <- grep('<div class="plotly-html-widget', plot_html, fixed = TRUE)
  plot_end <- grep('</script>', plot_html, fixed = TRUE)
  plot_end <- plot_end[plot_end > plot_start][1]
  
  plot_content <- paste(plot_html[plot_start:plot_end], collapse = "\n")
  
  # Delete the temporary file
  unlink(temp_file)
  
  # Format the div
  sprintf(div_template, title, plot_content)
}

# Create HTML header with CSS styling
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
    .plot-container {
      height: 380px;
    }
    h1 {
      color: #333;
      text-align: center;
      margin-bottom: 30px;
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
  </style>
</head>
<body>
  <h1>Peptide Sharing Analysis Dashboard</h1>

  <div class="summary">
    <h2>Analysis Summary</h2>
    <p>This dashboard presents an interactive analysis of peptide sharing across samples. The analysis is based on %d unique peptides across %d samples. Out of these peptides, %d (%d%%) are unique to a single sample, while %d (%d%%) are shared across multiple samples.</p>
    <p>The most prevalent peptide length is %dmer, making up %d%% of all peptides. Peptides of length %d amino acids show the highest rate of sharing across samples.</p>
  </div>
'

# Create HTML footer
html_footer <- '
  <div class="footer">
    <p>Peptide Sharing Analysis Dashboard | Generated: %s</p>
  </div>
</body>
</html>
'

# Generate summary stats for the dashboard
total_peptides <- nrow(peptide_report)
unique_peptide_count <- sum(sample_count_dist$n[sample_count_dist$sample_count == 1])
shared_peptide_count <- sum(sample_count_dist$n[sample_count_dist$sample_count > 1])
unique_peptide_pct <- round(unique_peptide_count / total_peptides * 100)
shared_peptide_pct <- round(shared_peptide_count / total_peptides * 100)
most_common_length <- length_dist# R script for interactive peptide sample distribution visualization

# Setting directory
setwd("~/Documents/Github/HLA-I_Analysis/")

# Path to the combined peptides file
combined_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/combined_peptides.tsv"

# Create the main output directory first
main_output_dir <- file.path(dirname(combined_file), "shared_peptide_visualizations")
dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)

# Create subdirectories within main output directory
viz_dir <- file.path(main_output_dir, "visualizations")
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

excel_dir <- file.path(main_output_dir, "excel_reports")
dir.create(excel_dir, recursive = TRUE, showWarnings = FALSE)

# Load required libraries
required_packages <- c("tidyverse", "plotly", "htmlwidgets", "openxlsx", 
                       "UpSetR", "VennDiagram", "RColorBrewer")

# Install and load necessary packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
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
    most_common_length_pct <- round(max(length_dist$percentage))
    best_sharing_length <- avg_sharing_by_length# R script for interactive peptide sample distribution visualization
    
    # Setting directory
    setwd("~/Documents/Github/HLA-I_Analysis/")
    
    # Path to the combined peptides file
    combined_file <- "/Users/dinarabadi/Documents/Github/HLA-I_Analysis/combined_peptides.tsv"
    
    # Create the main output directory first
    main_output_dir <- file.path(dirname(combined_file), "shared_peptide_visualizations")
    dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create subdirectories within main output directory
    viz_dir <- file.path(main_output_dir, "visualizations")
    dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
    
    excel_dir <- file.path(main_output_dir, "excel_reports")
    dir.create(excel_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Load required libraries
    required_packages <- c("tidyverse", "plotly", "htmlwidgets", "openxlsx", 
                           "UpSetR", "VennDiagram", "RColorBrewer")
    
    # Install and load necessary packages
    for (pkg in required_packages) {
      if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
      }
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
    # 7. Create a comprehensive HTML dashboard with all visualizations
    ################################################################################
    
    # Create a comprehensive HTML dashboard that combines all plots
    # First, create a function to wrap plotly objects for the dashboard
    wrap_in_div <- function(plotly_obj, title) {
      div_template <- '
  <div class="dashboard-item">
    <h3>%s</h3>
    <div class="plot-container">
      %s
    </div>
  </div>
  '
      
      # Save the plot to a temporary file
      temp_file <- tempfile(fileext = ".html")
      htmlwidgets::saveWidget(plotly_obj, temp_file, selfcontained = TRUE)
      
      # Read the saved plot HTML
      plot_html <- readLines(temp_file)
      
      # Extract just the plot div (not the full HTML document)
      plot_start <- grep('<div class="plotly-html-widget', plot_html, fixed = TRUE)
      plot_end <- grep('</script>', plot_html, fixed = TRUE)
      plot_end <- plot_end[plot_end > plot_start][1]
      
      plot_content <- paste(plot_html[plot_start:plot_end], collapse = "\n")
      
      # Delete the temporary file
      unlink(temp_file)
      
      # Format the div
      sprintf(div_template, title, plot_content)
    }
    
    # Create HTML header with CSS styling
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
    .plot-container {
      height: 380px;
    }
    h1 {
      color: #333;
      text-align: center;
      margin-bottom: 30px;
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
  </style>
</head>
<body>
  <h1>Peptide Sharing Analysis Dashboard</h1>

  <div class="summary">
    <h2>Analysis Summary</h2>
    <p>This dashboard presents an interactive analysis of peptide sharing across samples. The analysis is based on %d unique peptides across %d samples. Out of these peptides, %d (%d%%) are unique to a single sample, while %d (%d%%) are shared across multiple samples.</p>
    <p>The most prevalent peptide length is %dmer, making up %d%% of all peptides. Peptides of length %d amino acids show the highest rate of sharing across samples.</p>
  </div>
'
    
    # Create HTML footer
    html_footer <- '
  <div class="footer">
    <p>Peptide Sharing Analysis Dashboard | Generated: %s</p>
  </div>
</body>
</html>
'
    
    # Generate summary stats for the dashboard
    total_peptides <- nrow(peptide_report)
    unique_peptide_count <- sum(sample_count_dist$n[sample_count_dist$sample_count == 1])
    shared_peptide_count <- sum(sample_count_dist$n[sample_count_dist$sample_count > 1])
    unique_peptide_pct <- round(unique_peptide_count / total_peptides * 100)
    shared_peptide_pct <- round(shared_peptide_count / total_peptides * 100)
    most_common_length <- length_dist$`Peptide Length`[which.max(length_dist$n)]
    most_common_length_pct <- round(max(length_dist$percentage))
    best_sharing_length <- avg_sharing_by_length$`Peptide Length`[which.max(avg_sharing_by_length$avg_samples)]
    
        ################################################################################
        # Print completion message
        ################################################################################
        
        cat("\nVisualization and Excel report generation complete!\n")
        cat("Interactive visualizations saved to:", viz_dir, "\n")
        cat("Excel reports saved to:", excel_dir, "\n")
        cat("The following files were created:\n")
        cat("- Interactive visualizations:\n")
        cat("  - interactive_sample_count_distribution.html\n")
        cat("  - interactive_unique_vs_shared.html\n")
        cat("  - interactive_length_distribution.html\n")
        cat("  - interactive_avg_sharing_by_length.html\n")
        cat("  - interactive_gene_heatmap.html\n")
        if (num_samples >= 2 && num_samples <= 4) {
          cat("  - interactive_sample_intersections.html\n")
        }
        cat("  - peptide_sharing_dashboard.html (comprehensive dashboard)\n")
        
        cat("\n- Excel reports:\n")
        cat("  - peptide_sample_count_distribution.xlsx\n")
        cat("  - unique_vs_shared_peptides.xlsx\n")
        cat("  - peptide_length_distribution.xlsx\n")
        cat("  - peptide_sharing_by_length.xlsx\n")
        cat("  - gene_peptide_distribution.xlsx\n")
        if (num_samples >= 2 && num_samples <= 4) {
          cat("  - sample_intersections.xlsx\n")
        }