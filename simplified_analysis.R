# Simplified DLBCL Analysis Script
# This script works with minimal dependencies for users who have trouble with package installation

# Check for the minimal required packages
minimal_packages <- c("stats", "graphics", "utils", "datasets", "base")

# Function to create synthetic data
generate_sample_data <- function() {
  # Create a small gene expression matrix (100 genes x 20 samples)
  set.seed(42) # For reproducibility
  
  # Generate gene names
  gene_names <- paste0("GENE_", 1:100)
  
  # Generate sample names
  sample_names <- paste0("SAMPLE_", 1:20)
  
  # Generate raw counts (integers)
  raw_counts <- matrix(stats::rpois(100 * 20, lambda = 100), nrow = 100, ncol = 20)
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
    age_at_diagnosis = round(stats::rnorm(20, mean = 65, sd = 10)),
    ann_arbor_clinical_stage = sample(c("Stage I", "Stage II", "Stage III", "Stage IV"), 20, replace = TRUE),
    primary_diagnosis = rep("Diffuse Large B-Cell Lymphoma, NOS", 20),
    stage = stages
  )
  
  # Create directory structure if needed
  if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
  if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)
  
  # Save data files
  message("Saving synthetic data...")
  write.csv(raw_counts, "data/processed/raw_counts_matrix.csv")
  write.csv(sample_metadata, "data/processed/sample_metadata.csv", row.names = FALSE)
  
  return(list(counts = raw_counts, metadata = sample_metadata))
}

# Basic differential expression analysis with minimal dependencies
basic_differential_expression <- function(counts, metadata) {
  message("Performing basic differential expression analysis...")
  
  # Get sample groups
  tumor_samples <- metadata$barcode[metadata$sample_type == "Primary Tumor"]
  normal_samples <- metadata$barcode[metadata$sample_type == "Solid Tissue Normal"]
  
  # Calculate means for each group
  tumor_means <- rowMeans(counts[, tumor_samples, drop = FALSE])
  normal_means <- rowMeans(counts[, normal_samples, drop = FALSE])
  
  # Calculate log2 fold change
  log2_fc <- log2((tumor_means + 1) / (normal_means + 1))
  
  # Calculate simple t-test p-values
  pvalues <- numeric(nrow(counts))
  for (i in 1:nrow(counts)) {
    if (length(normal_samples) > 1 && length(tumor_samples) > 1) {
      test_result <- try(stats::t.test(counts[i, tumor_samples], counts[i, normal_samples]), silent = TRUE)
      if (!inherits(test_result, "try-error")) {
        pvalues[i] <- test_result$p.value
      } else {
        pvalues[i] <- 1
      }
    } else {
      pvalues[i] <- NA
    }
  }
  
  # Multiple testing correction
  adjusted_pvalues <- stats::p.adjust(pvalues, method = "BH")
  
  # Create results table
  results <- data.frame(
    Gene = rownames(counts),
    log2FoldChange = log2_fc,
    pvalue = pvalues,
    padj = adjusted_pvalues
  )
  
  # Sort by p-value
  results <- results[order(results$pvalue), ]
  
  # Save results
  write.csv(results, "results/tables/differential_expression_results.csv", row.names = FALSE)
  
  # Return results
  return(results)
}

# Basic visualization function
create_basic_plots <- function(counts, metadata, results) {
  message("Creating basic visualizations...")
  
  # Check if ggplot2 is available for better plots
  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
  
  # 1. Histogram of p-values
  pdf("results/figures/pvalue_histogram.pdf")
  hist(results$pvalue, breaks = 20, main = "Histogram of P-values", 
       xlab = "P-value", col = "skyblue")
  dev.off()
  
  # 2. Volcano plot (with base R)
  pdf("results/figures/volcano_plot.pdf")
  plot(results$log2FoldChange, -log10(results$pvalue),
       main = "Volcano Plot",
       xlab = "Log2 Fold Change",
       ylab = "-Log10 P-value",
       pch = 16,
       col = ifelse(results$padj < 0.05 & abs(results$log2FoldChange) > 1, 
                    "red", "black"))
  abline(h = -log10(0.05), col = "blue", lty = 2)
  abline(v = c(-1, 1), col = "blue", lty = 2)
  dev.off()
  
  # 3. Top Genes Box Plot
  top_genes <- head(results[order(results$padj), ], 5)$Gene
  
  pdf("results/figures/top_genes_boxplot.pdf")
  par(mfrow = c(2, 3))
  for (gene in top_genes) {
    gene_idx <- which(rownames(counts) == gene)
    if (length(gene_idx) > 0) {
      boxplot(split(as.numeric(counts[gene_idx, ]), metadata$sample_type),
              main = gene,
              xlab = "Sample Type",
              ylab = "Expression",
              col = c("lightblue", "salmon"))
    }
  }
  dev.off()
  par(mfrow = c(1, 1))
  
  message("Basic plots saved to results/figures/")
}

# Main function to run the simplified analysis
run_simplified_analysis <- function() {
  message("Starting simplified DLBCL analysis...")
  
  # Generate sample data
  data <- generate_sample_data()
  
  # Run basic differential expression analysis
  results <- basic_differential_expression(data$counts, data$metadata)
  
  # Create visualizations
  create_basic_plots(data$counts, data$metadata, results)
  
  message("Analysis complete! Results saved in the 'results/' directory.")
  
  return(list(
    data = data,
    results = results
  ))
}

# Run the simplified analysis
cat("This script will run a simplified version of the DLBCL analysis\n")
cat("using minimal R dependencies and synthetic data.\n\n")
cat("Run 'analysis <- run_simplified_analysis()' to start the analysis.\n")