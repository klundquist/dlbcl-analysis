# Enhanced visualization script
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)
library(tidyverse)
library(DESeq2)

# First, let's load our saved results
load_saved_results <- function() {
  message("Loading saved results...")
  
  # Load DESeq results from CSV
  deseq_results <- read.csv("results/tables/differential_expression_results.csv", row.names=1)
  
  # Load normalized counts
  normalized_counts <- read.csv("results/tables/normalized_counts.csv", row.names=1)
  
  # Load sample metadata
  sample_metadata <- read.csv("data/processed/sample_metadata.csv", row.names=1)
  
  # Load GO enrichment results
  enrichment_results <- read.csv("results/tables/go_enrichment_results.csv", row.names=1)
  
  # Print summary of loaded data
  message("Loaded data summary:")
  message("Number of genes: ", nrow(deseq_results))
  message("Number of samples: ", ncol(normalized_counts))
  message("Number of enriched terms: ", nrow(enrichment_results))
  
  return(list(
    deseq = deseq_results,
    counts = normalized_counts,
    metadata = sample_metadata,
    enrichment = enrichment_results
  ))
}

# 1. Create enhanced MA plot
create_ma_plot <- function(res_df) {
  message("Creating MA plot...")
  
  # Add significance column
  res_df$significant <- ifelse(res_df$padj < 0.05, "Significant", "Not Significant")
  
  # Create MA plot
  ma_plot <- ggplot(res_df, aes(x = log10(baseMean), y = log2FoldChange, color = significant)) +
    geom_point(size = 1, alpha = 0.6) +
    scale_color_manual(values = c("grey70", "red3")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
    annotate("text", x = max(log10(res_df$baseMean)), y = 2, 
             hjust = 1, label = "Higher in Deceased Patients", color = "darkred") +
    annotate("text", x = max(log10(res_df$baseMean)), y = -2, 
             hjust = 1, label = "Higher in Surviving Patients", color = "darkblue") +
    labs(x = "Log10 Mean Expression (average across all samples)", 
         y = "Log2 Fold Change (Deceased vs. Surviving)",
         title = "MA Plot: Gene Expression Changes in DLBCL",
         subtitle = "Red points: significantly different genes (adj. p-value < 0.05)\nPositive fold change: higher in deceased patients") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave("results/figures/ma_plot.png", ma_plot, width = 10, height = 8)
  return(ma_plot)
}

# 2. Create enhanced heatmap
create_enhanced_heatmap <- function(norm_counts, res_df, sample_metadata) {
  message("Creating enhanced heatmap...")
  
  # Get top 30 DEGs by p-value
  top_genes <- head(order(res_df$padj), 30)
  
  # Create matrix for heatmap
  mat <- as.matrix(norm_counts[top_genes, ])
  rownames(mat) <- res_df$symbol[top_genes]
  
  # Scale rows
  mat_scaled <- t(scale(t(mat)))
  
  # Create annotation dataframe
  anno_col <- data.frame(
    Patient_Outcome = factor(ifelse(sample_metadata$vital_status == "Alive", 
                                    "Surviving", "Deceased")),
    row.names = colnames(mat)
  )
  
  # Create annotation colors
  anno_colors <- list(
    Patient_Outcome = c(Surviving = "#2ECC71", Deceased = "#E74C3C")
  )
  
  # Create heatmap
  pheatmap(mat_scaled,
           annotation_col = anno_col,
           annotation_colors = anno_colors,
           show_rownames = TRUE,
           show_colnames = FALSE,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "complete",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           main = "Top 30 Differentially Expressed Genes",
           filename = "results/figures/enhanced_heatmap.png",
           width = 12,
           height = 10)
}

# 3. Create enhanced volcano plot
create_enhanced_volcano <- function(res_df) {
  message("Creating enhanced volcano plot...")
  
  volcano_plot <- EnhancedVolcano(res_df,
                                  lab = res_df$symbol,
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  title = 'Deceased vs. Surviving Patients',
                                  subtitle = 'Differential Expression Analysis',
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  pointSize = 3.0,
                                  labSize = 4.0,
                                  col = c('#808080', '#ADD8E6', '#FFB6C1', '#FF0000'),
                                  colAlpha = 0.5,
                                  legendPosition = 'right',
                                  drawConnectors = TRUE,
                                  widthConnectors = 0.2)
  
  ggsave("results/figures/enhanced_volcano.png", volcano_plot, width = 12, height = 10)
  return(volcano_plot)
}

# 4. Create expression boxplots
create_top_gene_boxplots <- function(norm_counts, res_df, sample_metadata) {
  message("Creating boxplots for top genes...")
  
  # Get top 6 genes by p-value
  top_genes <- head(order(res_df$padj), 6)
  gene_symbols <- res_df$symbol[top_genes]
  
  # Prepare data for plotting
  plot_data <- data.frame()
  for (i in seq_along(top_genes)) {
    gene_data <- data.frame(
      Gene = gene_symbols[i],
      Expression = as.numeric(norm_counts[top_genes[i], ]),
      Status = ifelse(sample_metadata$vital_status == "Alive", 
                      "Surviving", "Deceased")
    )
    plot_data <- rbind(plot_data, gene_data)
  }
  
  # Create boxplot
  boxplot <- ggplot(plot_data, aes(x = Status, y = log2(Expression + 1), fill = Status)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    facet_wrap(~Gene, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = c("#E74C3C", "#2ECC71")) +
    labs(title = "Expression of Top Differentially Expressed Genes",
         x = "Patient Outcome",
         y = "Log2 Normalized Counts") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("results/figures/top_genes_boxplot.png", boxplot, width = 12, height = 8)
  return(boxplot)
}

# 5. Create enrichment visualization
create_enrichment_plot <- function(enrichment_df) {
  message("Creating enrichment plot...")
  
  # Select top 15 terms by p-value
  top_terms <- head(enrichment_df[order(enrichment_df$pvalue), ], 15)
  
  # Create dot plot
  enrich_plot <- ggplot(top_terms, 
                        aes(x = reorder(Description, -log10(pvalue)), 
                            y = -log10(pvalue), 
                            size = Count, 
                            color = p.adjust)) +
    geom_point() +
    coord_flip() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "GO Term",
         y = "-log10(p-value)",
         title = "Top 15 Enriched GO Terms",
         size = "Gene Count",
         color = "Adjusted\np-value") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave("results/figures/enrichment_dotplot.png", enrich_plot, width = 12, height = 8)
  return(enrich_plot)
}

# Main visualization function
create_all_visualizations <- function() {
  message("Creating all visualizations...")
  
  # Load results
  results <- load_saved_results()
  
  # Create each visualization
  ma_plot <- create_ma_plot(results$deseq)
  create_enhanced_heatmap(results$counts, results$deseq, results$metadata)
  volcano_plot <- create_enhanced_volcano(results$deseq)
  boxplots <- create_top_gene_boxplots(results$counts, results$deseq, results$metadata)
  enrich_plot <- create_enrichment_plot(results$enrichment)
  
  message("All visualizations completed! Check results/figures/ directory.")
}

# Run visualizations
# create_all_visualizations()