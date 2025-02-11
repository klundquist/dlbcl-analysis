# Setup and data acquisition script for lymphoma analysis project

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c(
  "TCGAbiolinks",
  "DESeq2",
  "edgeR",
  "pheatmap",
  "EnhancedVolcano",
  "SummarizedExperiment"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load required libraries
library(TCGAbiolinks)
library(DESeq2)
library(edgeR)
library(tidyverse)
library(SummarizedExperiment)

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

# Main processing function
main_processing <- function(counts_matrix, sample_metadata) {
  message("Running main processing...")
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = sample_metadata,
    design = ~stage
  )
  
  # Filter low count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Normalize counts
  dds <- estimateSizeFactors(dds)
  normalized_counts <- counts(dds, normalized=TRUE)
  
  # Save processed data
  saveRDS(dds, "data/processed/deseq2_object.rds")
  write.csv(normalized_counts, "results/tables/normalized_counts.csv")
  
  # Save summary statistics
  summary_stats <- data.frame(
    total_genes = nrow(counts_matrix),
    filtered_genes = nrow(dds),
    samples = ncol(dds),
    tumor_samples = sum(sample_metadata$condition == "Tumor"),
    normal_samples = sum(sample_metadata$condition == "Normal")
  )
  write.csv(summary_stats, "results/tables/summary_statistics.csv")
  
  message("Processing complete. Data saved in data/processed/ and results/")
  return(list(
    dds = dds,
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