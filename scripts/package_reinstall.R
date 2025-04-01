# Package installation script for DLBCL RNA-seq analysis project

# Set repository options
options(repos = c(CRAN = "https://cloud.r-project.org"))

# First, install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)
}

# Function to install a package with error handling
install_package <- function(pkg_name, bioc = FALSE) {
  tryCatch({
    message(paste0("Attempting to install ", pkg_name, "..."))
    
    if (bioc) {
      BiocManager::install(pkg_name, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg_name, quiet = TRUE)
    }
    
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      message(paste0(pkg_name, " installed successfully"))
      return(TRUE)
    } else {
      message(paste0("Failed to install ", pkg_name))
      return(FALSE)
    }
  }, error = function(e) {
    message(paste0("Error installing ", pkg_name, ": ", e$message))
    return(FALSE)
  })
}

# CRAN packages
cran_packages <- c(
  "tidyverse",
  "ggplot2",
  "pheatmap",
  "RColorBrewer"
)

# Bioconductor packages
bioc_packages <- c(
  "DESeq2",
  "TCGAbiolinks",
  "edgeR",
  "EnhancedVolcano",
  "clusterProfiler",
  "org.Hs.eg.db",
  "SummarizedExperiment"
)

# Install CRAN packages
message("Installing CRAN packages...")
cran_results <- sapply(cran_packages, install_package)

# Install Bioconductor packages
message("Installing Bioconductor packages...")
bioc_results <- sapply(bioc_packages, function(pkg) install_package(pkg, bioc = TRUE))

# Check installation results
all_packages <- c(cran_packages, bioc_packages)
all_results <- c(cran_results, bioc_results)

# Print summary
message("\nPackage Installation Summary:")
for (i in seq_along(all_packages)) {
  status <- if (all_results[i]) "INSTALLED" else "FAILED"
  message(sprintf("  %s: %s", all_packages[i], status))
}

# Check if core packages are available
core_packages <- c("DESeq2", "TCGAbiolinks", "tidyverse")
if (all(sapply(core_packages, requireNamespace, quietly = TRUE))) {
  message("\nCore packages are installed. You can proceed with the analysis.")
} else {
  missing_core <- core_packages[!sapply(core_packages, requireNamespace, quietly = TRUE)]
  message("\nWARNING: Some core packages are missing: ", paste(missing_core, collapse = ", "))
  message("You may need to modify the analysis scripts or install these packages manually.")
}