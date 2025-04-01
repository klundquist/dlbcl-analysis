# Setup and data acquisition script for lymphoma analysis project

# Check for packages and suggest using package_reinstall.R if missing
required_packages <- c(
  "TCGAbiolinks",
  "DESeq2", 
  "tidyverse",
  "SummarizedExperiment"
)

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(paste0(
    "The following required packages are missing: ", 
    paste(missing_packages, collapse = ", "), 
    "\nPlease run the package installation script first: source('scripts/package_reinstall.R')"
  ))
}

# Load required libraries with error handling
load_package <- function(pkg_name) {
  tryCatch({
    library(pkg_name, character.only = TRUE)
    message(paste0("Loaded package: ", pkg_name))
    return(TRUE)
  }, error = function(e) {
    message(paste0("Failed to load package: ", pkg_name, " - ", e$message))
    return(FALSE)
  })
}

# Load core packages
core_packages <- c("TCGAbiolinks", "tidyverse", "SummarizedExperiment")
core_loaded <- sapply(core_packages, load_package)

if (!all(core_loaded)) {
  stop("Failed to load core packages. Please check the installation.")
}

# Try to load DESeq2 (needed for later steps)
deseq_loaded <- load_package("DESeq2")
if (!deseq_loaded) {
  warning("DESeq2 could not be loaded. Some analysis steps may fail.")
}

# Create project directory structure
dirs <- c("data/raw", 
          "data/processed", 
          "results/figures", 
          "results/tables",
          "scripts")

for (dir in dirs) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# Download and save the TCGA-DLBC dataset
download_and_save_data <- function() {
  message("Downloading TCGA-DLBC RNA-seq dataset...")
  
  # Query TCGA data
  query <- GDCquery(
    project = "TCGA-DLBC",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  # Download data
  GDCdownload(query, method = "api", files.per.chunk = 10, directory = "data/raw")
  
  # Prepare data
  message("Preparing expression data...")
  data <- GDCprepare(query, directory = "data/raw")
  
  # Save the data
  message("Saving raw data...")
  saveRDS(data, "data/raw/TCGA_DLBC_data.rds")
  
  return(data)
}

# Process raw counts and save
process_raw_counts <- function(data) {
  message("Processing raw counts...")
  
  # Extract count matrix
  counts <- assay(data, "unstranded")
  
  # Save processed counts
  write.csv(counts, "data/processed/raw_counts_matrix.csv")
  
  return(counts)
}

# Create and save sample metadata
create_sample_metadata <- function(data) {
  message("Creating sample metadata...")
  
  # Extract essential clinical data
  sample_metadata <- colData(data) %>%
    as.data.frame() %>%
    dplyr::select(
      barcode,
      sample_type,
      gender,
      vital_status,
      age_at_diagnosis,
      ann_arbor_clinical_stage,
      primary_diagnosis
    ) %>%
    # Convert to simple data frame with character columns
    mutate(across(everything(), as.character)) %>%
    # Replace NA with "Unknown"
    mutate(across(everything(), ~replace_na(., "Unknown")))
  
  # Create stage factor based on Ann Arbor staging
  sample_metadata$stage <- factor(
    case_when(
      str_detect(sample_metadata$ann_arbor_clinical_stage, "^Stage I") ~ "Early",
      str_detect(sample_metadata$ann_arbor_clinical_stage, "^Stage II") ~ "Early",
      str_detect(sample_metadata$ann_arbor_clinical_stage, "^Stage III") ~ "Advanced",
      str_detect(sample_metadata$ann_arbor_clinical_stage, "^Stage IV") ~ "Advanced",
      TRUE ~ "Unknown"
    ),
    levels = c("Early", "Advanced", "Unknown")
  )
  
  # Save metadata
  message("Saving sample metadata...")
  write.csv(sample_metadata, 
            "data/processed/sample_metadata.csv", 
            row.names = FALSE)
  
  return(sample_metadata)
}

# Main processing function with fallback options
main_processing <- function(counts_matrix, sample_metadata) {
  message("Running main processing...")
  
  # Check if DESeq2 is available for normalization
  if (requireNamespace("DESeq2", quietly = TRUE)) {
    message("Using DESeq2 for normalization...")
    
    # Create DESeq2 object
    tryCatch({
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts_matrix,
        colData = sample_metadata,
        design = ~stage
      )
      
      # Filter low count genes
      keep <- rowSums(counts(dds)) >= 10
      dds <- dds[keep,]
      
      # Normalize counts
      dds <- DESeq2::estimateSizeFactors(dds)
      normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
      
      # Save processed data
      saveRDS(dds, "data/processed/deseq2_object.rds")
      
      result_object <- list(
        dds = dds,
        normalized_counts = normalized_counts
      )
    }, error = function(e) {
      message("Error in DESeq2 processing: ", e$message)
      message("Falling back to simple normalization...")
      
      # Fallback to simple normalization if DESeq2 fails
      result_object <- fallback_normalization(counts_matrix, sample_metadata)
    })
  } else {
    message("DESeq2 not available, using fallback normalization...")
    result_object <- fallback_normalization(counts_matrix, sample_metadata)
  }
  
  # Get normalized counts from whichever method succeeded
  normalized_counts <- result_object$normalized_counts
  
  # Save normalized counts
  write.csv(normalized_counts, "results/tables/normalized_counts.csv")
  
  # Save summary statistics
  summary_stats <- data.frame(
    total_genes = nrow(counts_matrix),
    filtered_genes = nrow(normalized_counts),
    samples = ncol(normalized_counts),
    tumor_samples = sum(sample_metadata$sample_type == "Primary Tumor", na.rm = TRUE),
    normal_samples = sum(sample_metadata$sample_type == "Solid Tissue Normal", na.rm = TRUE)
  )
  write.csv(summary_stats, "results/tables/summary_statistics.csv")
  
  message("Processing complete. Data saved in data/processed/ and results/")
  return(result_object)
}

# Fallback normalization when DESeq2 is not available
fallback_normalization <- function(counts_matrix, sample_metadata) {
  message("Performing simple CPM normalization...")
  
  # Filter low count genes
  keep <- rowSums(counts_matrix) >= 10
  filtered_counts <- counts_matrix[keep,]
  
  # Simple CPM normalization
  lib_sizes <- colSums(filtered_counts)
  norm_factors <- lib_sizes / mean(lib_sizes)
  normalized_counts <- t(t(filtered_counts) / norm_factors) * 1e6
  
  message("Simple normalization complete")
  
  # Return results
  return(list(
    normalized_counts = normalized_counts
  ))
}

# Main execution function
run_pipeline <- function() {
  message("Starting analysis pipeline...")
  
  # Download and save data
  data <- download_and_save_data()
  
  # Process counts
  counts <- process_raw_counts(data)
  
  # Create metadata
  metadata <- create_sample_metadata(data)
  
  # Run main processing
  processed_data <- main_processing(counts, metadata)
  
  message("Pipeline complete!")
  return(processed_data)
}