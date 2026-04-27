# scRNA-seq Analysis (Seurat Workflow)
This repository contains a basic single-cell RNA-seq analysis workflow using Seurat.

## Requirements
R

R packages:
- Seurat
- dplyr

## Workflow
The script performs:
- Quality control (QC)
- Filtering of low-quality cells
- Normalization using SCTransform
- Dimensionality reduction (PCA)
- Clustering (Louvain)
- UMAP visualization
- Marker gene identification
- (Requires manual cell type annotation)

## Input Data
The script expects a 10x Genomics formatted matrix directory containing:

- matrix.mtx.gz
- barcodes.tsv.gz
- features.tsv.gz

## How to Run
1. Open R or RStudio
2. Update the `matrix_dir` path in the script:

   `matrix_dir <- "path/to/filtered_feature_bc_matrix"`

3. Run the script line-by-line adjusting parameters as needed
