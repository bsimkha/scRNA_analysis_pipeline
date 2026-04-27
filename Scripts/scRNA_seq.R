# -------------------------------
# Load required libraries
# -------------------------------
library(Seurat)
library(dplyr)

# -------------------------------
# Load 10x matrix data
# -------------------------------
matrix_dir <- "path/to/matrix/data"
counts <- Read10X(data.dir = matrix_dir)

# Quick sanity check
dim(counts)   # should be ~ genes x cells

# -------------------------------
# Create Seurat object
# -------------------------------
scObject <- CreateSeuratObject(counts = counts, project = "pbmc_1k")
head(scObject@meta.data)

# -------------------------------
# Add QC metrics
# -------------------------------
# % mitochondrial genes (common QC metric)
scObject[["percent.mt"]] <- PercentageFeatureSet(scObject, pattern = "^MT-")

# Visualize distributions
VlnPlot(scObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Check relationships between metrics
plot1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# -------------------------------
# Filter low-quality cells
# -------------------------------
# Basic thresholds — adjust depending on dataset
scObject <- subset(scObject, subset = nFeature_RNA > 200 &
                     nFeature_RNA < 2500 &
                     percent.mt < 10) #adjust ad needed

# Re-check after filtering
VlnPlot(scObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# -------------------------------
# Normalize + variance stabilization
# -------------------------------
# SCTransform handles normalization + scaling internally
scObject <- SCTransform(scObject)

# -------------------------------
# Dimensionality reduction
# -------------------------------
scObject <- RunPCA(scObject, assay = "SCT")

# Inspect PCs
ElbowPlot(scObject)
VizDimLoadings(scObject, dims = 1:2)

# -------------------------------
# Clustering
# -------------------------------
set.seed(123)  # helps with reproducibility

scObject <- FindNeighbors(scObject, dims = 1:6, assay = "SCT") #adjust dimensions as needed using elbow plot generated previously
scObject <- FindClusters(scObject, resolution = 0.5)

# -------------------------------
# UMAP visualization
# -------------------------------
scObject <- RunUMAP(scObject, dims = 1:6) #adjust dimensions as needed using elbow plot
DimPlot(scObject, reduction = "umap", label = TRUE)

# -------------------------------
# Identify marker genes
# -------------------------------
markers <- FindAllMarkers(
  scObject,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(markers) #sanity check

# Pull top markers per cluster
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

View(top_markers)

# -------------------------------
# Manual annotation
# -------------------------------
# Use marker genes + references (PanglaoDB, literature, proteinatlas.org etc.)
annotation_df <- data.frame(
  seurat_clusters = c("0","1","2","3","4"),
  cell_type = c("T cells","B cells","Monocytes","NK cells","Neutrophil")
) #THIS STEP REQUIRES MANUAL CURATION!!!

scObject$cell_type <- annotation_df$cell_type[
  match(as.character(scObject$seurat_clusters),
        annotation_df$seurat_clusters)
]

# -------------------------------
# Visualize annotated clusters
# -------------------------------
Idents(scObject) <- scObject$cell_type

DimPlot(scObject, reduction = "umap", label = TRUE, repel = TRUE)

# -------------------------------
# Heatmap of marker genes
# -------------------------------
heatmap_genes <- top_markers$gene
DoHeatmap(scObject, features = heatmap_genes)

# -------------------------------
# Save object for later use
# -------------------------------
saveRDS(scObject, file = "pbmc_1k_annotated.rds")