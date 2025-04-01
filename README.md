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

2. Set up the environment (choose one option):

### Option A: Using Docker (recommended)
This is the most reliable method that works across all platforms:

```bash
# Build and start the Docker container (this will take some time on first run)
./run_docker.sh
```

The Docker container comes with all required packages pre-installed. Once it starts, you can run the analysis directly:

```r
source("scripts/01_get_process_data.R")
run_pipeline()
```

If you have docker-compose installed:

```bash
# Alternative: Using docker-compose
docker-compose up -d
docker exec -it dlbcl-analysis R
```

### Option B: Using Conda
```bash
# Create and activate the conda environment
conda env create -f environment.yml
conda activate dlbcl-analysis

# Start R within the environment
R
```

### Option C: Using renv
```bash
# Start R and run the setup script
R
> source("setup_renv.R")
```

### Option D: Direct installation
If you prefer to install packages directly into your R environment:
```R
# Use our simple requirements script
source("requirements.R")
```

3. Verify your environment:
```R
# After setting up your environment, verify it's working correctly
source("check_environment.R")
```
This will display a table showing all required packages and their status.

## Usage

### Option 1: Ultra-simplified analysis
If you're having trouble with any of the environment setups or package installation, try our ultra-simplified analysis that uses only base R:

```R
source("simplified_analysis.R")
analysis <- run_simplified_analysis()
```

This will:
- Generate a small synthetic dataset
- Perform basic differential expression analysis with t-tests
- Create simple plots using base R graphics
- Save all results to the results directory

### Option 2: Using sample data
If your environment is set up but you want to avoid downloading large datasets:

```R
source("scripts/use_existing_data.R")
data <- work_with_data(force_download = TRUE)
```

This will generate a small synthetic dataset that you can use to test the workflow.

### Option 3: Running the full pipeline
If you want to run the full pipeline from scratch:

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

### Troubleshooting
If you encounter issues with package installation:

1. Try using the pre-processed data approach (Option 1)
2. Ensure you have the latest version of R (4.4.1 or higher recommended)
3. Run `source("scripts/package_reinstall.R")` which tries to install binary packages when possible
4. If you're using macOS, you might need additional development tools:
   ```bash
   xcode-select --install
   ```

## Data Management
This project works with large genomic datasets that are not stored in the repository. The data files are:
- Raw counts matrix (~8.5MB)
- Sample metadata (~8KB)
- Normalized counts (~26MB)

These files are excluded from version control using `.gitignore` to keep the repository size manageable. When running the pipeline:
- The download script will acquire data from TCGA
- Processed files will be saved locally
- Use the `use_existing_data.R` script to work with previously processed data

For testing purposes, you can generate a small synthetic dataset that mimics the structure of the real data.

## Results
The analysis identifies:
- 54 significantly differentially expressed genes
- Key pathways associated with survival outcomes
- Potential prognostic markers

View the complete analysis and results in `project-writeup.html`

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)

## Contact
Karl Lundquist - klundquist@gmail.com
