# Helper script for downloading and working with pre-processed data
# Instead of storing large data files in the repo, this script downloads them from an external source

# Load basic libraries
message("Loading basic libraries...")
basic_pkgs <- c("utils") # Just need base R utils for download

for (pkg in basic_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
    message(paste0("Loaded: ", pkg))
  } else {
    warning(paste0("Package ", pkg, " not available. Some functionality may be limited."))
  }
}

# Create directories if they don't exist
create_directories <- function() {
  dirs <- c("data/processed", "results/tables")
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      message(paste0("Created directory: ", dir))
    }
  }
}

# Function to download pre-processed data
download_preprocessed_data <- function(force_download = FALSE) {
  message("Checking for pre-processed data...")
  
  # Define file paths
  counts_file <- "data/processed/raw_counts_matrix.csv"
  metadata_file <- "data/processed/sample_metadata.csv"
  normalized_file <- "results/tables/normalized_counts.csv"
  
  # Check if files already exist
  files_exist <- file.exists(counts_file) && 
                file.exists(metadata_file) && 
                file.exists(normalized_file)
  
  # Skip download if files exist and force_download is FALSE
  if (files_exist && !force_download) {
    message("Pre-processed data files already exist. Use force_download=TRUE to re-download.")
    return(TRUE)
  }
  
  # Create directories if needed
  create_directories()
  
  # In a real scenario, these would be URLs to a data repository
  # For now, we'll use placeholders and provide instructions
  # Example: "https://example.com/dlbcl-data/raw_counts_matrix.csv"
  
  message("\nNOTE: In a real implementation, this script would download data from URLs like:")
  message("  - https://example.com/dlbcl-data/raw_counts_matrix.csv")
  message("  - https://example.com/dlbcl-data/sample_metadata.csv")
  message("  - https://example.com/dlbcl-data/normalized_counts.csv")
  message("\nHowever, since we don't have an actual data repository set up, please:")
  message("1. Run the data processing pipeline OR")
  message("2. Generate a small sample dataset for testing\n")
  
  # Option to generate a small sample dataset for testing
  generate_sample <- readline(prompt = "Would you like to generate a small sample dataset for testing? (y/n): ")
  if (tolower(generate_sample) == "y") {
    message("Generating small sample dataset for testing...")
    generate_sample_data()
    return(TRUE)
  } else {
    message("Download aborted. Please run the data processing pipeline instead.")
    return(FALSE)
  }
}

# Function to generate a small sample dataset for testing
generate_sample_data <- function() {
  # Create a small gene expression matrix (100 genes x 20 samples)
  set.seed(42) # For reproducibility
  
  # Generate gene names
  gene_names <- paste0("GENE_", 1:100)
  
  # Generate sample names
  sample_names <- paste0("SAMPLE_", 1:20)
  
  # Generate raw counts (integers)
  raw_counts <- matrix(rpois(100 * 20, lambda = 100), nrow = 100, ncol = 20)
  rownames(raw_counts) <- gene_names
  colnames(raw_counts) <- sample_names
  
  # Generate sample metadata
  sample_types <- sample(c("Primary Tumor", "Solid Tissue Normal"), 20, replace = TRUE, prob = c(0.8, 0.2))
  stages <- sample(c("Early", "Advanced", "Unknown"), 20, replace = TRUE, prob = c(0.4, 0.4, 0.2))
  
  sample_metadata <- data.frame(
    barcode = sample_names,
    sample_type = sample_types,
    gender = sample(c("male", "female"), 20, replace = TRUE),
    vital_status = sample(c("Alive", "Dead"), 20, replace = TRUE, prob = c(0.7, 0.3)),
    age_at_diagnosis = round(rnorm(20, mean = 65, sd = 10)),
    ann_arbor_clinical_stage = sample(c("Stage I", "Stage II", "Stage III", "Stage IV"), 20, replace = TRUE),
    primary_diagnosis = rep("Diffuse Large B-Cell Lymphoma, NOS", 20),
    stage = stages
  )
  
  # Generate normalized counts (floating point)
  normalized_counts <- matrix(rnorm(100 * 20, mean = 8, sd = 2), nrow = 100, ncol = 20)
  rownames(normalized_counts) <- gene_names
  colnames(normalized_counts) <- sample_names
  normalized_counts <- 2^normalized_counts # Log2 back-transform to get positive values
  
  # Save the files
  message("Saving sample data files...")
  write.csv(raw_counts, "data/processed/raw_counts_matrix.csv")
  write.csv(sample_metadata, "data/processed/sample_metadata.csv", row.names = FALSE)
  write.csv(normalized_counts, "results/tables/normalized_counts.csv")
  
  message("Sample dataset generated successfully!")
  return(TRUE)
}

# Function to load pre-processed data
load_preprocessed_data <- function() {
  message("Loading pre-processed data...")
  
  # Check if files exist
  counts_file <- "data/processed/raw_counts_matrix.csv"
  metadata_file <- "data/processed/sample_metadata.csv"
  normalized_file <- "results/tables/normalized_counts.csv"
  
  # Check if files exist
  files_exist <- file.exists(counts_file) && 
                file.exists(metadata_file) && 
                file.exists(normalized_file)
  
  # If files don't exist, try to download them
  if (!files_exist) {
    message("Data files not found. Attempting to download...")
    if (!download_preprocessed_data()) {
      stop("Cannot load data. Please run the pipeline to generate the data first.")
    }
  }
  
  data <- list()
  
  # Load raw counts
  message("Loading raw counts matrix...")
  data$raw_counts <- read.csv(counts_file, row.names = 1)
  
  # Load sample metadata
  message("Loading sample metadata...")
  data$sample_metadata <- read.csv(metadata_file)
  
  # Load normalized counts
  message("Loading normalized counts...")
  data$normalized_counts <- read.csv(normalized_file, row.names = 1)
  
  # Return the loaded data
  return(data)
}

# Function to summarize dataset
summarize_data <- function(data) {
  message("\nData Summary:")
  
  if (!is.null(data$raw_counts)) {
    message(paste0("Raw counts matrix: ", nrow(data$raw_counts), " genes x ", 
                  ncol(data$raw_counts), " samples"))
  }
  
  if (!is.null(data$normalized_counts)) {
    message(paste0("Normalized counts matrix: ", nrow(data$normalized_counts), " genes x ", 
                  ncol(data$normalized_counts), " samples"))
  }
  
  if (!is.null(data$sample_metadata)) {
    message(paste0("Sample metadata: ", nrow(data$sample_metadata), " samples"))
    
    # Sample type breakdown
    if ("sample_type" %in% colnames(data$sample_metadata)) {
      sample_types <- table(data$sample_metadata$sample_type)
      message("\nSample types:")
      for (i in 1:length(sample_types)) {
        message(paste0("  ", names(sample_types)[i], ": ", sample_types[i]))
      }
    }
    
    # Stage breakdown if available
    if ("stage" %in% colnames(data$sample_metadata)) {
      stages <- table(data$sample_metadata$stage)
      message("\nDisease stages:")
      for (i in 1:length(stages)) {
        message(paste0("  ", names(stages)[i], ": ", stages[i]))
      }
    }
  }
  
  message("\nUse data$raw_counts, data$normalized_counts, and data$sample_metadata to access the data")
}

# Main function to load and work with pre-processed data
work_with_data <- function(force_download = FALSE) {
  # Attempt to download data if needed
  if (force_download) {
    download_preprocessed_data(force_download = TRUE)
  }
  
  # Load the data
  data <- load_preprocessed_data()
  
  # Summarize the data
  summarize_data(data)
  
  # Return the data for further analysis
  return(data)
}

# Message for users
message("To work with the data:")
message("1. Run 'data <- work_with_data()' to load existing pre-processed data")
message("2. Or run 'data <- work_with_data(force_download = TRUE)' to re-download or generate sample data")