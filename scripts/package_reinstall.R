# Remove problematic packages
remove.packages(c("ggplot2", "gtable", "pheatmap", "EnhancedVolcano"))

# Reinstall Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("DESeq2", "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"))

# Reinstall CRAN packages
install.packages(c("ggplot2", "pheatmap", "tidyverse", "RColorBrewer"))