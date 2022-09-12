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

#Rename clusters
seu$Clusters <- factor(paste0("C", as.numeric(as.character(seu$seurat_clusters)) + 1), 
                       levels = paste0("C", c(1:length(levels(seu$seurat_clusters)))))
saveRDS(seu, "rds/seu2.rds")

#marker visualization
marker <- read.table("ref/markers.tsv")
avg <- AverageExpression(seu, return.seurat = T, group.by = "Clusters")
mtx <- avg@assays$RNA@scale.data[marker,]

col_fun1 <- colorRamp2(c(-2,-1,0,1,2), viridis(5, option = "A"))
fh = function(x) hclust(dist(x), method="ward.D2")
ht1 <- Heatmap(mtx, name = "Scaled expression", cluster_columns = fh, cluster_rows = fh, 
               col = col_fun1, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))
p1 <- draw(ht1)
pdf("figures/GEx_markers_Heatmap.pdf", width = 7.5, height = 8)
p1
dev.off()

#manual cell type classification
TME <- paste0("C", c(1,15,2,31,32,9,26,6,27))

tme.seu <- subset(seu, subset = Clusters %in% TME)
epi.seu <- subset(seu, subset = Clusters %ni% TME)

saveRDS(tme.seu, "rds/tme_seu.rds")
saveRDS(epi.seu, "rds/epi_seu.rds")
