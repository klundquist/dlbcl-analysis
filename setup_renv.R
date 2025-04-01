# Setup script for renv
# This creates an R environment using renv

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
library(renv)

# Initialize renv
renv::init()

# Install CRAN packages
cran_packages <- c(
  "tidyverse",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "data.table",
  "remotes"
)

# Install Bioconductor packages
bioc_packages <- c(
  "DESeq2",
  "TCGAbiolinks",
  "edgeR",
  "EnhancedVolcano", 
  "clusterProfiler",
  "org.Hs.eg.db",
  "SummarizedExperiment"
)

# Function to install packages with error handling
install_package <- function(pkg, bioc = FALSE) {
  tryCatch({
    message(paste0("Installing ", pkg, "..."))
    
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    
    return(TRUE)
  }, error = function(e) {
    message(paste0("Error installing ", pkg, ": ", e$message))
    return(FALSE)
  })
}

# Install packages
message("Installing CRAN packages...")
cran_results <- sapply(cran_packages, install_package)

message("Installing Bioconductor packages...")
bioc_results <- sapply(bioc_packages, function(pkg) install_package(pkg, bioc = TRUE))

# Snapshot the environment
renv::snapshot()

message("\nEnvironment setup complete. To activate this environment in a new R session, run:")
message("renv::restore()")