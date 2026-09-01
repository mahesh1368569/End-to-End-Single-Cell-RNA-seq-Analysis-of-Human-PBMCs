#Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
#Read the data (replace with your file path)
data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19")

# Create the Seurat object with minimal filtering
seurat_ob <- CreateSeuratObject(
  counts = data,               # Input count matrix
  project = "pbmc3k",          # Project name
  min.cells = 3,               # Include genes detected in at least 3 cells
  min.features = 200           # Include cells with at least 200 detected genes
)

# Calculate percentage of mitochondrial genes
# The pattern "^MT-" captures all genes starting with "MT-" (mitochondrial genes)
seurat_ob[["percent.mt"]] <- PercentageFeatureSet(seurat_ob, pattern = "^MT-")


# Create violin plots for QC metrics
VlnPlot(seurat_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter cells based on QC metrics
data_filter <- subset(seurat_ob,
                      subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)
#Visualize metrics after filtering
VlnPlot(data_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalize the data
data_norm <- NormalizeData(
  data_filter,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Find variable features
data_hv <- FindVariableFeatures(
  data_norm,
  selection.method = "vst",    # Variance stabilizing transformation
  nfeatures = 2000             # Number of features to select
)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data_hv), 10)
#Plot variable features
plot1 <- VariableFeaturePlot(data_hv)
plot1


#Add labels for top 10 variable features
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2 

# Scale the data
all.genes <- rownames(data_norm)
data_scaled <- ScaleData(data_norm, features = all.genes)

# Run PCA
data_pca <- RunPCA(data_scaled, features = VariableFeatures(object = data_hv))
#Visualize PCA results
print(data_pca[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(data_pca, reduction = "pca")

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# Visualize top features in each PC
DimHeatmap(data_pca, dims = 1:15, cells = 500, balanced = TRUE)

#Run JackStraw procedure
data_p <- JackStraw(data_pca, num.replicate = 100)
data_score <- ScoreJackStraw(data_p, dims = 1:20)
#Visualize JackStraw results
JackStrawPlot(data_score, dims = 1:20)

ElbowPlot(pbmc)

#Run UMAP
data_umap <- RunUMAP(data_pca, dims = 1:10)
#Visualize UMAP

DimPlot(data_umap, reduction = "umap")

saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

# Run t-SNE if desired
data_tsne <- RunTSNE(data_pca, dims = 1:10)
DimPlot(data_tsne, reduction = "tsne")

#Find nearest neighbors
data_nn <- FindNeighbors(data_umap, dims = 1:10)
#Find clusters
data_cluster <- FindClusters(data_nn, resolution = 0.5)

#Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# Visualize clusters on UMAP
DimPlot(data_cluster, reduction = "umap", label = TRUE)

#Find markers for all clusters
data.markers <- FindAllMarkers(
  data_cluster,
  only.pos = TRUE,          # Only find positive markers
  min.pct = 0.25,           # Minimum percentage of cells expressing the gene
  logfc.threshold = 0.25    # Minimum log fold-change
)
#View top markers per cluster
top_markers <- data.markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 5, order_by = avg_log2FC)
print(top_markers)

#Find markers for a specific cluster (e.g., cluster 2)
cluster2.markers <- FindMarkers(
  data_cluster,
  ident.1 = 2,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
head(cluster2.markers, n = 5)

#find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)


#===========================================================================================
# Manual Cell type assignment

#Define cell type labels based on canonical markers
new.cluster.ids <- c(
  "Naive CD4 T",  # Cluster 0
  "CD14+ Mono",   # Cluster 1
  "Memory CD4 T", # Cluster 2
  "B",            # Cluster 3
  "CD8 T",        # Cluster 4
  "FCGR3A+ Mono", # Cluster 5
  "NK",           # Cluster 6
  "DC",           # Cluster 7
  "Platelet"      # Cluster 8
)
names(new.cluster.ids) <- levels(data_cluster)
#Rename clusters
data_cluster_named <- RenameIdents(data_cluster, new.cluster.ids)


#Visualize named clusters
DimPlot(data_cluster_named, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#load libraries
library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

hdf5_obj <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19")

pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)

str(pbmc.seurat)

#QC and Filtering -----------
#explore QC
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 mitoPercent < 10)

#It is a good practice to filter out cells with non-sufficient genes identified and genes with non-sufficient expression across cells.
#pre-process standard workflow ---------------
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)

#running steps above to get clusters
View(pbmc.seurat.filtered@meta.data)

DimPlot(pbmc.seurat.filtered, reduction = 'umap')

#get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))

install.packages("rlang")
library(rlang)
