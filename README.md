# DESeq2 Workflow

A genomics workflow for differential gene expression analysis using DESeq2, specifically designed for comparing GTEX and TCGA data across different tissue types.

## Overview

This repository contains scripts for processing and analyzing RNA-seq data to identify differentially expressed genes between normal tissue samples (GTEX) and cancer samples (TCGA). The workflow consists of two main components:

1. **Data Preprocessing (Shell Script)**: Prepares gene expression matrices (GEM) for DESeq2 analysis
2. **Differential Expression Analysis (R Script)**: Performs the statistical analysis using DESeq2

## Requirements

### System Requirements
- Bash shell environment (Linux/macOS)
- R version 3.6.0 or higher

### R Packages
- DESeq2
- Other dependencies will be installed automatically when loading DESeq2

## Installation

Clone this repository to your local machine:

```bash
git clone https://github.com/gabby-grant/deseq-workflow
cd deseq-workflow
```

No additional installation is needed. The R script includes a command to set up a custom library path for R packages.

## Usage
### Step 1: Download Data
This example will use cervix from Wang et al. 
```bash
wget -O cervix_GTEX.gz https://figshare.com/ndownloader/files/9150199

wget -O cervix_TCGA.gz 
https://figshare.com/ndownloader/files/9150202
https://figshare.com/ndownloader/files/9150202
```
**Citation:**
Wang, Q., Armenia, J., Zhang, C. et al. Unifying cancer and normal RNA sequencing data from different sources. Sci Data 5, 180061 (2018). https://doi.org/10.1038/sdata.2018.61
### Step 2: Merge Data
```shell
module load anaconda3

git clone https://github.com/feltus/mergegem.git

python ./mergegem/mergegem.py cervix_GTEX cervix_TCGA cervix-gtex-tcga.txt
```
### Step 3: Data Preprocessing

The `process_gtex_tcga.sh` script converts raw expression data into the proper format for DESeq2 analysis.

```bash
./process_gtex_tcga.sh <tissue_type>
```

Example:
```bash
./process_gtex_tcga.sh cervix
```

If no tissue type is provided, the script will prompt you to enter one.

#### Input Files Required:
- `<tissue>-gtex-tcga.txt`: The raw gene expression matrix
- `<tissue>-gtex-tcga.labels.txt`: Sample labels file

#### Output Files:
- `<tissue>-gtex-tcga-integer.txt`: Expression values converted to integers
- `<tissue>-gtex-tcga-integer-unique.txt`: Duplicate genes removed
- `<tissue>-gtex-tcga-clean.txt`: Cleaned data with dashes replaced by underscores
- `<tissue>-gtex-tcga-clean.csv`: CSV version of the cleaned data
- `<tissue>-gtex-tcga-dash.labels.txt`: Labels with dashes replaced by underscores
- `<tissue>-gtex-tcga.comparison.csv`: Comparison file for DESeq2 analysis

### Step 4: Differential Expression Analysis

The `deseq_analysis.R` script performs differential expression analysis using DESeq2.

In R:
```R
source("deseq_analysis.R")
results <- main("your_counts_file.csv", "your_annotation_file.csv", "output_file.csv")
```

Example:
```R
results <- main("cervix-gtex-tcga-clean.csv", "cervix-gtex-tcga.comparison.csv", "cervix_results.csv")
```

#### Functions in the R Script:

- `subgem()`: Extracts a sub-matrix of counts for a specific group
- `subanot()`: Extracts a subset of the sample annotation matrix
- `run_deseq()`: Runs DESeq2 analysis with tissue type as parameter
- `main()`: Main function to execute the workflow

## Output

The analysis produces:
- A CSV file containing differentially expressed genes with adjusted p-values < 0.05
- Each gene row includes:
  - log2 fold change
  - p-value and adjusted p-value
  - Normalized expression values for all samples

## Data Format

### Expression Matrix Format:
- First column: Gene symbols (Hugo_Symbol)
- Subsequent columns: Expression values per sample

### Sample Annotation Format:
- Column 1: Sample identifiers
- Column 2: Group labels (e.g., cervix_GTEX, cervix_TCGA)
- Column 3: Comparison labels

## Troubleshooting

Common issues:
- **File not found errors**: Ensure that input files are in the correct location and named properly
- **R package installation errors**: Check if your R environment has internet access and permissions
- **Execution permission errors**: Make the shell script executable with `chmod +x process_gtex_tcga.sh`

## Contributing

Contributions to improve this workflow are welcome. Please feel free to submit a pull request or open an issue.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this workflow in your research, please cite:

```
Grant, G. (2025). DESeq Workflow: A pipeline for differential gene expression analysis.
GitHub repository: https://github.com/gabby-grant/deseq-workflow
```

## Contact

For questions or support, please open an issue in the GitHub repository.
