FROM rocker/r-ver:4.3.2

LABEL maintainer="Karl Lundquist <klundquist@gmail.com>"
LABEL description="Docker image for the DLBCL RNA-seq Analysis Project"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN install2.r --error --deps TRUE \
    tidyverse \
    pheatmap \
    RColorBrewer

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('DESeq2', 'TCGAbiolinks', 'SummarizedExperiment', 'edgeR', 'EnhancedVolcano'), update = FALSE, ask = FALSE)"

# Create project directory
WORKDIR /project

# Create necessary directories
RUN mkdir -p data/raw data/processed results/figures results/tables

# Set entrypoint
ENTRYPOINT ["R"]