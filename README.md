# A Distinguished Roadmap of Fibroblast Senescence in Predicting Immunotherapy Response and Prognosis across Human Cancers


## Overview
This R script, developed by Dongjie Chen from the Pancreatic Disease Center, Ruijin Hospital, Shanghai, China, is used to construct a gene signature for pancreatic disease research based on single-cell RNA-sequencing (scRNA-seq) data. The core steps include calculating scores related to cellular senescence, identifying differential expression between cell types, and building a final gene signature based on these results.

### Author: Dongjie Chen  
### Date: November 3, 2024  
### Email: chen_dj@sjtu.edu.cn  

---

## Requirements

- **R version 4.2.2 or higher**
- **Libraries**:
  - `pbapply`: For parallel processing.
  - `future`: For parallel computing with multiple sessions.
  - `GSVA`: To compute Gene Set Variation Analysis.
  - `compositions`: For geometric mean calculations.
  - `dplyr`: For data manipulation.
  - `Seurat`: For scRNA-seq analysis.
  
You can install the required libraries using the following commands:
```r
install.packages(c("pbapply", "future", "dplyr", "Seurat"))
BiocManager::install(c("GSVA", "compositions"))
```

---

## Script Breakdown

### Step 1: Calculate FSx
- **Objective**: This step calculates a score for each single-cell RNA dataset based on cellular senescence gene sets (loaded from `cellular_senescence_gene.rds`).
- **Data**: It requires two RDS files:
  - `cellular_senescence_gene.rds`: Gene sets related to cellular senescence.
  - `40_pancancer_scRNA_seq_dat.rds`: List of scRNA-seq data from various pancreatic cancer samples.
  
The script performs Gene Set Variation Analysis (GSVA) to calculate the S_CAF score for each dataset, and correlates the expression of genes with this score using Spearman’s rank correlation. The result is stored as `FSx_list`, which includes coefficients and p-values for each gene.

### Step 2: Calculate FSy
- **Objective**: This step identifies differentially expressed genes between CAF (Cancer-Associated Fibroblasts) and control cells in each dataset.
- **Process**:
  - Creates a new column `group` to separate CAF and control cells.
  - Performs differential expression analysis using `FindMarkers` from Seurat, with thresholds for log fold-change and p-values.
  
The differentially expressed genes with an adjusted p-value less than 1e-05 are stored in `FSy_list`.

### Step 3: LM.SIG Construction
- **Objective**: This step constructs the final gene signature by combining the results from FSx and FSy.
  - **FSx filtering**: Only genes with a positive correlation coefficient (`coef > 0`) and a highly significant p-value (adjusted p-value < 1e-05) are selected.
  - **FSy filtering**: Only genes with a log fold-change greater than 0.25 are kept, and ribosomal proteins (denoted by gene names starting with "RP") are excluded.
  
The genes that appear in both FSx and FSy are merged to form the final gene list, with coefficients for each dataset. The geometric mean of the coefficients is calculated and used to filter out genes with low overall expression, resulting in the final signature `FSS`.

### Final Output
- The gene signature (`FSS`) is saved as an RDS file (`FSS.rds`) for later use.

---

## File Description
1. **`cellular_senescence_gene.rds`**: RDS file containing the gene set for cellular senescence.
2. **`40_pancancer_scRNA_seq_dat.rds`**: RDS file containing the single-cell RNA sequencing data from multiple pancreatic cancer samples.
3. **`FSS.rds`**: The final output containing the constructed gene signature, which is a list of genes that are significantly associated with the cellular senescence phenotype in the context of CAF.

---

## How to Run the Script
1. **Prepare your data**:
   - Ensure you have the required input RDS files (`cellular_senescence_gene.rds` and `40_pancancer_scRNA_seq_dat.rds`).
   
2. **Execute the script**:
   - Load the script into your R environment.
   - Run the script line by line or as a whole.

3. **Retrieve results**:
   - After the script finishes running, the gene signature will be saved in `FSS.rds`. You can load this file using the following command:
     ```r
     FSS <- readRDS("FSS.rds")
     ```

---

## Notes
- Ensure that the environment is set up to handle large-scale computations, especially for parallel processing with `future` and `pbapply`. The script is optimized for use on multi-core systems to speed up processing.
- The script assumes that all datasets in `scRNA_list` are compatible with Seurat objects and that the cell types are labeled in the `Celltype` metadata.

---

## Contact Information
For any questions or issues related to the script, please contact the author:
- **Email**: chen_dj@sjtu.edu.cn
- **Institution**: Pancreatic Disease Center, Ruijin Hospital, Shanghai, China


