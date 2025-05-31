#' HLA-I Analysis Core Utilities
#' peptide_core_utils.R
#' Contains shared functions used across all peptide analysis scripts
#' @author Your Name
#' @version 1.0

# Load required packages function
load_required_packages <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  cat("All required packages loaded successfully\n")
}

# Create directory structure with standardized subdirectories
create_output_directories <- function(base_dir, analysis_name = "peptide_analysis") {
  # Main output directory
  main_dir <- file.path(base_dir, analysis_name)
  dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create standard subdirectories
  viz_dir <- file.path(main_dir, "visualizations")
  dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
  
  excel_dir <- file.path(main_dir, "excel_reports")
  dir.create(excel_dir, recursive = TRUE, showWarnings = FALSE)
  
  data_dir <- file.path(main_dir, "processed_data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Return directory paths as a list
  return(list(
    main_dir = main_dir,
    viz_dir = viz_dir,
    excel_dir = excel_dir,
    data_dir = data_dir
  ))
}

# Function to read peptide TSV files with sample ID extraction
read_peptide_files <- function(file_list, sample_id_pattern = NULL, 
                               sample_id_replacements = NULL) {
  if (length(file_list) == 0) {
    stop("No files provided to read")
  }
  
  all_data <- list()
  
  for (file in file_list) {
    filename <- basename(file)
    
    # Extract sample ID - use pattern if provided, otherwise use filename
    if (!is.null(sample_id_pattern)) {
      sample_id <- gsub(sample_id_pattern, "\\1", filename)
      
      # Apply replacements if provided
      if (!is.null(sample_id_replacements) && sample_id %in% names(sample_id_replacements)) {
        sample_id <- sample_id_replacements[[sample_id]]
      }
    } else {
      # Default to filename without extension
      sample_id <- tools::file_path_sans_ext(filename)
    }
    
    cat("Reading file:", filename, "- Sample ID:", sample_id, "\n")
    
    # Read the file
    data <- tryCatch({
      df <- read.delim(file, stringsAsFactors = FALSE)
      # Add a column for sample ID and filename
      df$Sample_ID <- sample_id
      df$Filename <- filename
      df
    }, error = function(e) {
      cat("Error reading file:", file, "- Error:", e$message, "\n")
      NULL
    })
    
    if (!is.null(data)) {
      all_data[[length(all_data) + 1]] <- data
    }
  }
  
  if (length(all_data) == 0) {
    stop("No files could be read successfully")
  }
  
  # Combine all data frames
  combined_data <- bind_rows(all_data)
  
  return(combined_data)
}

# Create Excel output with consistent formatting
create_excel_report <- function(sheet_list, output_file, header_style = NULL) {
  # Create a workbook
  wb <- createWorkbook()
  
  # Create default header style if not provided
  if (is.null(header_style)) {
    header_style <- createStyle(textDecoration = "bold", fgFill = "#D9D9D9")
  }
  
  # Add each sheet with formatting
  for (sheet_name in names(sheet_list)) {
    # Add a worksheet
    addWorksheet(wb, sheet_name)
    
    # Write data
    writeData(wb, sheet_name, sheet_list[[sheet_name]], headerStyle = header_style)
    
    # Auto-adjust column widths
    setColWidths(wb, sheet_name, cols = 1:ncol(sheet_list[[sheet_name]]), widths = "auto")
    
    # Freeze the header row
    freezePane(wb, sheet_name, firstRow = TRUE)
  }
  
  # Save the workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  cat("Excel report saved to:", output_file, "\n")
}

# Standardized color palettes and plotting defaults
get_color_palette <- function(palette_type = "sequential", n = 10) {
  if (palette_type == "sequential") {
    return(colorRampPalette(c("white", "steelblue", "darkblue"))(n))
  } else if (palette_type == "diverging") {
    return(colorRampPalette(c("blue", "white", "red"))(n))
  } else if (palette_type == "categorical") {
    return(brewer.pal(min(n, 8), "Set1"))
  } else if (palette_type == "binary") {
    return(c("white", "darkblue"))
  } else {
    # Default
    return(colorRampPalette(c("white", "steelblue", "darkblue"))(n))
  }
}

# Save both PDF and PNG versions of a ggplot
save_plot <- function(plot, filename, width = 10, height = 8, 
                      png_res = 100, png_width = 800, png_height = 600) {
  # Save PDF
  pdf_file <- paste0(filename, ".pdf")
  pdf(pdf_file, width = width, height = height)
  print(plot)
  dev.off()
  
  # Save PNG
  png_file <- paste0(filename, ".png")
  png(png_file, width = png_width, height = png_height, res = png_res)
  print(plot)
  dev.off()
  
  cat("Plot saved as:", pdf_file, "and", png_file, "\n")
}

# Function to add standard metadata to a dataframe
add_metadata <- function(df, analysis_type, version = "1.0", timestamp = Sys.time()) {
  attr(df, "analysis_type") <- analysis_type
  attr(df, "version") <- version
  attr(df, "timestamp") <- timestamp
  return(df)
}