# Simple requirements file for installing needed R packages
# Run this directly in R to install all required packages:
# source("requirements.R")

# Function to install a package if it's not already installed
install_if_missing <- function(package_name, bioc = FALSE) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    message(paste("Installing", package_name))
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(package_name, update = FALSE)
    } else {
      install.packages(package_name)
    }
  } else {
    message(paste(package_name, "is already installed"))
  }
}

# CRAN packages
install_if_missing("tidyverse")
install_if_missing("ggplot2")
install_if_missing("pheatmap")
install_if_missing("RColorBrewer")
install_if_missing("data.table")

# Bioconductor packages
install_if_missing("BiocManager")
install_if_missing("DESeq2", bioc = TRUE)
install_if_missing("TCGAbiolinks", bioc = TRUE)
install_if_missing("edgeR", bioc = TRUE)
install_if_missing("EnhancedVolcano", bioc = TRUE)
install_if_missing("clusterProfiler", bioc = TRUE)
install_if_missing("org.Hs.eg.db", bioc = TRUE)
install_if_missing("SummarizedExperiment", bioc = TRUE)

message("All required packages have been installed or were already present.")