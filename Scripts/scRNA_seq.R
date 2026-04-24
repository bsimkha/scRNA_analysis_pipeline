#Load the Seurat library
library(Seurat)

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
scObject <- subset(scObject, nFeature_RNA>1000 & nFeature_RNA<3000 & percent.mt<10)

#Plot the figures again to make sure the trimming worked
VlnPlot(scObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mt") #No correlation expected
plot2 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #High correlation expected
scatterplot <- plot1 + plot2
scatterplot

#Normalize the data
Normalized_scObject <- NormalizeData(scObject)
head(Normalized_scObject)
head(scObject)
