#Load libraries
library(Seurat)
library(dplyr)

#Load the directory containing matrix files generated from CellRanger
matrix_dir <- "~/Desktop/scRNA/1_Raw_data/filtered_feature_bc_matrix"
counts <- Read10X(data.dir = matrix_dir)

dim(counts) #sanity check to make sure the dimensions (Features*Cells) seem right

#Wrap counts into Seurat scObjectect
scObject <- CreateSeuratObject(counts = counts, project = "pbmc_1k")
head(scObject@meta.data) #Inspect metadata for sanity check

#Calculate mitochondrial percentage
scObject[["percent.mt"]] <- PercentageFeatureSet(scObject, pattern = "^MT-")
head(scObject@meta.data) #Should have a new percent.mt column

#Generate Violin plot
VlnPlot(scObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#These plots represent the relative distribution of features (genes per cell), total counts and mtDNA%
#Things to note
#1. Feature (genes per cell) should have unimodal distribution. Use the distribution to identify cutoff threshold
#2. Count (total) can sometimes have a slight bimodal distribution
#3. The mitochondrial contamination is generally expected to be < 25%

#Generate scatterplots to check for correlation between the three
plot1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mt") #No correlation expected
plot2 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #High correlation expected
scatterplot <- plot1 + plot2
scatterplot

#Use the violinplot to find the outliers and trim the data
#Baseline gold standard is 200 < x < 2500
scObject <- subset(scObject, nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<10)

#Plot the figures again to make sure the trimming worked
VlnPlot(scObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mt") #No correlation expected
plot2 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #High correlation expected
scatterplot <- plot1 + plot2
scatterplot

#Normalize the data using SCTransform
Normalized_scObject <- SCTransform(scObject)
head(Normalized_scObject)
head(scObject)

#Identify the highly variable genes across cell types
VariableFeaturePlot(Normalized_scObject) #plot without the variable genes labeled
top10 <- head(VariableFeatures(Normalized_scObject), 10) #Identify the top 10 variable genes
LabelPoints(plot = VariableFeaturePlot(Normalized_scObject), points = top10, repel = TRUE) #Plot with the top genes labeled

#Scale data with mean around 0
Normalized_scObject <- ScaleData(Normalized_scObject)

#Run the PCA
Normalized_scObject <- RunPCA(Normalized_scObject)

#Use elbow plot to identify the PC that don't add any further information
ElbowPlot(Normalized_scObject)

VizDimLoadings(Normalized_scObject, dims = 1:2)

#Find clusters of "like"cell types
Normalized_scObject <- FindNeighbors(Normalized_scObject, dims = 1:6) #the dimensions come from the eblow plot we plotted earlier
Normalized_scObject <- FindClusters(Normalized_scObject, resolution = 0.5)

#plot the UMAP to visaualize the data
Normalized_scObject <- RunUMAP(Normalized_scObject, dims = 1:6)
DimPlot(Normalized_scObject, reduction = "umap", label = TRUE)


#Find marker genes for cell type annotation
markers <- FindAllMarkers(Normalized_scObject, only.pos = TRUE)
head(markers)

#Find top markers
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

top_markers #check this df to manually curate a cell-type annotation file
#for example: Go to proteinatlas.org and search for the genes in the file
#Choose single cell and identify cell type. Do this for all genes to identify consensus for each cluster

#create annotation file
annotation_df <- data.frame(
  seurat_clusters = c("0","1","2","3","4"),
  cell_type = c("T cells","B cells","Monocytes","NK cells","Neutrophil")
)

Normalized_scObject$cell_type <- annotation_df$cell_type[
  match(as.character(Normalized_scObject$seurat_clusters),
        annotation_df$seurat_clusters)]

#Heat map for cell-type specific expression of top genes
DoHeatmap(Normalized_scObject, features = head(markers$gene, 20))

#Plot the UMAP with cell names labelled
Idents(Normalized_scObject) <- Normalized_scObject$cell_type
DimPlot(Normalized_scObject, reduction = "umap", repel = TRUE)
