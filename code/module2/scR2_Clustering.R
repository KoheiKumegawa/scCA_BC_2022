#----------------------------------------------------------------------------
# scR2_Clustering.R
#----------------------------------------------------------------------------
#Kohei Kumegawa, Japanese Foundation for Cancer Research

### CAUTION!! THIS PROCESS NEEDS LARGE MEMORY ###
library(Seurat)
library(dplyr)

seu <- readRDS("rds/seu.rds")

#Clustering
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu)
seu <- FindNeighbors(seu, dims = 1:50)
seu <- FindClusters(seu, resolution = 0.2)
seu <- RunUMAP(seu, dims = 1:50, return.model = T)
  
saveRDS(seu, "rds/seu2.rds")
