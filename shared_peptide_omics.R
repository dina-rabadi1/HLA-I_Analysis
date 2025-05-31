# Load libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

# Set working directory
setwd("/Users/dinarabadi/Documents/Github/HLA-I_Analysis")

# Create output directory structure
if (!dir.exists("results/shared50_omic")) {
  dir.create("results/shared50_omic/")
}

# Load peptide detection data (Sheet3 from peptide file)
peptides <- read_excel("shared50_filtered_out_148N_manual.xlsx", sheet = "Sheet3")

# Load transcriptome/proteome data (sheet = "use")
omics <- read_excel("rawdata/dina_requena_2024.xlsx", sheet = "use")

# Clean and separate peptide gene mappings (handle multiple mappings)
peptides_clean <- peptides %>%
  rename_with(~str_replace_all(., "\\s+", "_")) %>%
  separate_rows(genes, sep = ";") %>%
  mutate(genes = str_trim(genes))  # clean whitespace

# Normalize gene names to uppercase for consistent joins
peptides_clean <- peptides_clean %>%
  mutate(genes = toupper(str_trim(genes)))

omics <- omics %>%
  mutate(Gene_Symbol = toupper(Gene_Symbol))

# Join omics data by gene symbol
peptide_long <- peptides_clean %>%
  left_join(omics, by = c("genes" = "Gene_Symbol")) %>%
  mutate(immuno_detected = TRUE)

# # Optional: create a filtered version based on expression
# # E.g., keep only peptides from genes with RNA/protein upregulation
# peptide_filtered <- peptide_long %>%
#   filter(!is.na(l2fc_transcriptome) & !is.na(lfc_proteome))
# # another filtering example
# peptide_filtered <- peptide_long %>%
#   filter(!is.na(l2fc_transcriptome) & !is.na(lfc_proteome)) %>%
#   filter(l2fc_transcriptome > 1 | lfc_proteome > 1)


# Export peptide-centric long table
write_xlsx(peptide_long, "results/shared50_omic/peptide_long_format.xlsx")
write_xlsx(peptide_filtered, "results/shared50_omic/peptide_filtered_by_expression.xlsx")






# scatter visualization
# Load libraries
library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr)

# Make sure we're working from raw Excel again
peptide_long <- readxl::read_xlsx("results/shared50_omic/peptide_long_format.xlsx")

nrow(peptide_long)  # total rows in original file

peptide_long %>%
  mutate(
    l2fc_transcriptome_num = as.numeric(l2fc_transcriptome),
    lfc_proteome_num = as.numeric(lfc_proteome)
  ) %>%
  filter(is.na(l2fc_transcriptome_num) | is.na(lfc_proteome_num)) %>%
  dplyr::select(Peptide, l2fc_transcriptome, lfc_proteome) %>%
  head(20)

# Explicit coercion of all relevant columns
peptide_plot <- peptide_long %>%
  filter(!is.na(`l2fc_transcriptome`), !is.na(`lfc_proteome`)) %>%
  mutate(
    Peptide = as.character(Peptide),
    sample_count = as.numeric(as.character(sample_count)),
    l2fc_transcriptome = as.numeric(l2fc_transcriptome),
    lfc_proteome = as.numeric(lfc_proteome)
  )

nrow(peptide_plot)  # how many passed the !is.na filter and coercion


# Break any lingering class inheritance
peptide_plot_clean <- as.data.frame(peptide_plot)

# Now extract top peptides
top_peptides <- peptide_plot_clean %>%
  dplyr::arrange(desc(sample_count)) %>%
  dplyr::distinct(Peptide, .keep_all = TRUE) %>%
  dplyr::slice(1:10)

nrow(peptide_plot_clean)  


# Plot
p <- ggplot(peptide_plot, aes(x = l2fc_transcriptome, y = lfc_proteome)) +
  geom_point(aes(color = sample_count), alpha = 0.7, size = 3) +
  scale_color_viridis_c(option = "plasma") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text_repel(data = top_peptides,
                  aes(label = Peptide),
                  size = 3.5,
                  box.padding = 0.3,
                  max.overlaps = 10) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Peptides by Transcriptomic and Proteomic Fold Change",
    x = "Log2 Fold Change (Transcriptome)",
    y = "Log Fold Change (Proteome)",
    color = "Sample Count"
  )

# Save the plot as PNG
ggsave(
  filename = "results/shared50_omic/peptide_scatterplot.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

summary(peptide_plot$sample_count)  # check if NAs remain
