# Analysis and visualization script

# Reinstall and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install/reinstall required packages
required_packages <- c(
  "DESeq2", 
  "pheatmap", 
  "EnhancedVolcano",
  "clusterProfiler", 
  "org.Hs.eg.db", 
  "ggplot2",
  "tidyverse"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load libraries
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)

# Load the processed data
dds <- readRDS("data/processed/deseq2_object.rds")

# Create survival status factor
dds$survival_status <- factor(dds$vital_status, levels=c("Alive", "Dead"))

# Run differential expression analysis
run_deseq_analysis <- function(dds) {
  message("Running differential expression analysis...")
  
  # Print diagnostic information
  message("\nChecking survival status assignments:")
  print(table(dds$survival_status, useNA = "always"))
  
  # Run DESeq
  message("\nRunning DESeq...")
  dds$survival_status <- relevel(dds$survival_status, ref = "Alive")
  design(dds) <- ~survival_status
  dds <- DESeq(dds)
  
  # Get results
  message("\nGetting results...")
  res <- results(dds, name="survival_status_Dead_vs_Alive")
  
  # Add gene symbols
  message("\nAdding gene symbols...")
  clean_ids <- gsub("\\..*", "", row.names(res))
  res$symbol <- mapIds(org.Hs.eg.db,
                       keys=clean_ids,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
  
  # Save results
  write.csv(as.data.frame(res), "results/tables/differential_expression_results.csv")
  
  # Print summary
  message("\nSummary of differential expression results:")
  summary(res)
  
  return(res)
}

# Create basic barplot using base R
create_deg_barplot <- function(res) {
  message("Creating DEG barplot...")
  
  # Get significant DEGs
  sig_genes <- res[!is.na(res$padj) & res$padj < 0.05, ]
  up_genes <- sum(sig_genes$log2FoldChange > 0)
  down_genes <- sum(sig_genes$log2FoldChange < 0)
  
  # Create plot
  png("results/figures/deg_barplot.png", width=800, height=600)
  barplot(c(up_genes, down_genes), 
          names.arg=c("Up-regulated", "Down-regulated"),
          col=c("red", "blue"),
          main="Differentially Expressed Genes\n(Dead vs Alive)",
          ylab="Number of genes")
  dev.off()
  
  message(sprintf("Up-regulated genes: %d", up_genes))
  message(sprintf("Down-regulated genes: %d", down_genes))
}

# Create alternative heatmap using base R
create_heatmap <- function(dds, res) {
  message("Creating heatmap...")
  
  # Get normalized counts
  norm_counts <- counts(dds, normalized=TRUE)
  
  # Get top DEGs (lowest p-value)
  top_genes <- head(order(res$padj), 30)
  
  # Create matrix for heatmap
  mat <- norm_counts[top_genes, ]
  gene_names <- res$symbol[top_genes]
  
  # Scale rows
  mat_scaled <- t(scale(t(mat)))
  
  # Create plot
  png("results/figures/deg_heatmap.png", width=1200, height=800)
  par(mar=c(8, 10, 4, 2))  # Adjust margins for labels
  
  # Create heatmap
  heatmap(mat_scaled,
          labRow=gene_names,
          labCol=dds$survival_status,
          main="Top 30 Differentially Expressed Genes",
          cexRow=0.8,  # Adjust text size
          cexCol=0.8)
  
  dev.off()
}

# Run gene enrichment analysis
run_enrichment_analysis <- function(res) {
  message("Running gene enrichment analysis...")
  
  # Get significant up-regulated genes
  sig_up <- res[which(res$padj < 0.05 & res$log2FoldChange > 0), ]
  
  # Convert ENSEMBL IDs to ENTREZ IDs
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys=gsub("\\..*", "", rownames(sig_up)),
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
  
  # Run GO enrichment analysis
  ego <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)
  
  # Save results
  write.csv(as.data.frame(ego), "results/tables/go_enrichment_results.csv")
  message("Enrichment analysis results saved to results/tables/go_enrichment_results.csv")
  
  return(ego)
}

# Extract top DEGs
extract_top_degs <- function(res) {
  message("\nExtracting top differentially expressed genes...")
  
  # Get significant genes
  sig_genes <- res[!is.na(res$padj) & res$padj < 0.05, ]
  
  # Sort by absolute log2FoldChange
  sig_genes <- sig_genes[order(abs(sig_genes$log2FoldChange), decreasing=TRUE), ]
  
  # Create data frame with relevant information
  top_degs <- data.frame(
    Gene_Symbol = sig_genes$symbol,
    Log2_Fold_Change = round(sig_genes$log2FoldChange, 3),
    P_Value = format(sig_genes$pvalue, scientific=TRUE, digits=3),
    Adj_P_Value = format(sig_genes$padj, scientific=TRUE, digits=3),
    Expression_Change = ifelse(sig_genes$log2FoldChange > 0, "Up-regulated", "Down-regulated")
  )
  
  # Save to file
  write.csv(top_degs, "results/tables/top_degs.csv", row.names=FALSE)
  
  # Print top 10 genes
  message("\nTop 10 most differentially expressed genes:")
  print(head(top_degs, 10))
  
  return(top_degs)
}

# Update main analysis function
run_analysis <- function() {
  message("Starting analysis pipeline...")
  
  # Run differential expression analysis
  res <- run_deseq_analysis(dds)
  
  # Extract top DEGs
  top_degs <- extract_top_degs(res)
  
  # Create visualizations
  create_deg_barplot(res)
  create_heatmap(dds, res)
  
  # Run enrichment analysis
  enrich <- run_enrichment_analysis(res)
  
  message("Analysis complete! Check results/figures/ and results/tables/ for output")
  
  return(list(
    deseq_results = res,
    top_degs = top_degs,
    enrichment_results = enrich
  ))
}