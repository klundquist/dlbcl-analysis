# Check if the environment has all required packages
# Run this script to verify your R environment is properly set up

# Setup
cat("Checking R environment for DLBCL Analysis Project...\n\n")
cat("R version:", R.version.string, "\n\n")

# Function to check if package is installed and can be loaded
check_package <- function(package_name) {
  is_installed <- requireNamespace(package_name, quietly = TRUE)
  
  result <- list(
    name = package_name,
    installed = is_installed,
    loaded = FALSE,
    version = NA
  )
  
  if (is_installed) {
    tryCatch({
      # Try to load the package
      ns <- loadNamespace(package_name)
      result$loaded <- TRUE
      
      # Get version info
      if ("DESCRIPTION" %in% names(ns)) {
        result$version <- ns$DESCRIPTION$Version
      } else {
        desc <- packageDescription(package_name)
        result$version <- desc$Version
      }
    }, error = function(e) {
      result$error <- as.character(e)
    })
  }
  
  return(result)
}

# List of required packages
cran_packages <- c(
  "tidyverse",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "data.table"
)

bioc_packages <- c(
  "BiocManager",
  "DESeq2",
  "TCGAbiolinks",
  "edgeR",
  "EnhancedVolcano",
  "clusterProfiler",
  "org.Hs.eg.db",
  "SummarizedExperiment"
)

# Check all packages
cat("Checking CRAN packages:\n")
cran_results <- lapply(cran_packages, check_package)

cat("Checking Bioconductor packages:\n")
bioc_results <- lapply(bioc_packages, check_package)

all_results <- c(cran_results, bioc_results)

# Print results table
cat("\nPackage Status:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
cat(sprintf("%-20s %-10s %-10s %-15s\n", "Package", "Installed", "Loaded", "Version"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (result in all_results) {
  installed_str <- ifelse(result$installed, "Yes", "No")
  loaded_str <- ifelse(result$loaded, "Yes", "No")
  version_str <- ifelse(is.na(result$version), "-", result$version)
  
  cat(sprintf("%-20s %-10s %-10s %-15s\n", 
              result$name, installed_str, loaded_str, version_str))
}

cat(paste(rep("-", 60), collapse = ""), "\n\n")

# Check if there are any missing packages
missing_packages <- sapply(all_results, function(x) !x$installed)
missing_count <- sum(missing_packages)

if (missing_count > 0) {
  cat("WARNING: ", missing_count, " packages are missing!\n")
  cat("The following packages need to be installed:\n")
  for (i in which(missing_packages)) {
    cat("  -", all_results[[i]]$name, "\n")
  }
  cat("\nRun one of the environment setup scripts to install the missing packages.\n")
} else {
  cat("SUCCESS: All required packages are installed.\n")
  
  # Check for load failures
  load_fails <- sapply(all_results, function(x) x$installed && !x$loaded)
  load_fail_count <- sum(load_fails)
  
  if (load_fail_count > 0) {
    cat("WARNING: ", load_fail_count, " packages could not be loaded!\n")
    cat("The following packages have load issues:\n")
    for (i in which(load_fails)) {
      cat("  -", all_results[[i]]$name, "\n")
      if (!is.null(all_results[[i]]$error)) {
        cat("    Error:", all_results[[i]]$error, "\n")
      }
    }
  } else {
    cat("All packages can be loaded successfully.\n")
    cat("\nYour environment is ready for the DLBCL analysis!\n")
  }
}