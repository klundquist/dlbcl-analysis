# DLBCL RNA-seq Analysis Project

## Overview
This project analyzes RNA-seq data from Diffuse Large B-Cell Lymphoma (DLBCL) patients to identify gene expression patterns associated with survival outcomes. Using data from the TCGA-DLBC project, the analysis includes differential expression analysis, pathway enrichment, and visualization of results.

## Project Structure
```
dlbcl-analysis/
├── data/
│   ├── raw/              # Raw TCGA-DLBC RNA-seq data
│   └── processed/        # Processed data files
├── results/
│   ├── figures/          # Generated plots and visualizations
│   └── tables/          # Analysis results in tabular format
├── scripts/
│   ├── 01_get_process_data.R    # Data acquisition and processing
│   ├── 02_analysis_visualization.R  # Main analysis script
│   └── 03_enhanced_visualization.R  # Enhanced visualization script
└── docs/
    └── index.html        # Project write-up and documentation
```

## Requirements
- R version 4.4.1 or higher
- Required R packages:
  - Bioconductor packages:
    - DESeq2
    - TCGAbiolinks
    - clusterProfiler
    - org.Hs.eg.db
    - EnhancedVolcano
  - CRAN packages:
    - tidyverse
    - pheatmap
    - RColorBrewer

## Installation
1. Clone this repository:
```bash
git clone https://github.com/klundquist/dlbcl-analysis.git
cd dlbcl-analysis
```

2. Install required R packages:
```R
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required packages
BiocManager::install(c("DESeq2", "TCGAbiolinks", "clusterProfiler", 
                      "org.Hs.eg.db", "EnhancedVolcano"))

# Install CRAN packages
install.packages(c("tidyverse", "pheatmap", "RColorBrewer"))
```

## Usage
1. Data Processing:
```R
source("scripts/01_get_process_data.R")
run_pipeline()
```

2. Analysis and Visualization:
```R
source("scripts/02_analysis_visualization.R")
results <- run_analysis()
```

3. Enhanced Visualizations:
```R
source("scripts/03_enhanced_visualization.R")
create_all_visualizations()
```

## Results
The analysis identifies:
- 54 significantly differentially expressed genes
- Key pathways associated with survival outcomes
- Potential prognostic markers

View the complete analysis and results in `docs/index.html`

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)

## Contact
Karl Lundquist - klundquist@gmail.com
