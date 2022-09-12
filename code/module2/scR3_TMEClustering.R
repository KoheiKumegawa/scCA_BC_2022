#----------------------------------------------------------------------------
# scR3_TMEClustering.R
#----------------------------------------------------------------------------
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridisLite)

tme_seu <- readRDS("rds/tme_seu.rds")

#Clustering
tme_seu <- NormalizeData(tme_seu, normalization.method = "LogNormalize", scale.factor = 10000)
tme_seu <- FindVariableFeatures(tme_seu, selection.method = "vst", nfeatures = 2000)
tme_seu <- ScaleData(tme_seu, features = rownames(tme_seu))
tme_seu <- RunPCA(tme_seu)
tme_seu <- FindNeighbors(tme_seu, dims = 1:50)
tme_seu <- FindClusters(tme_seu, resolution = 0.2)
tme_seu <- RunUMAP(tme_seu, dims = 1:50, return.model = T)
tme_seu$Clusters2 <- factor(paste0("TME", as.numeric(as.character(tme_seu$seurat_clusters)) + 1),
                            levels = paste0("TME", c(1:length(levels(tme_seu$seurat_clusters)))))

#plotting
clust_cols <- c(ArchR::ArchRPalettes$calm, ArchR::ArchRPalettes$kelly)[c(1:17)] %>% `names<-`(., paste0("TME",c(1:17)))
sampl_cols <- c(ArchR::ArchRPalettes$calm, ArchR::ArchRPalettes$kelly)[c(1:38)] %>% `names<-`(., unique(tme_seu$orig.ident))
p1 <- DimPlot(tme_seu, reduction = "umap", label = TRUE, group.by = "Clusters2", cols = clust_cols, raster = T) + 
      ggtitle("Clusters") + ArchR::theme_ArchR()
p2 <- DimPlot(tme_seu, reduction = "umap", label = F, group.by = "orig.ident", cols = sampl_cols, raster = T) + 
      ggtitle("Samples") + ArchR::theme_ArchR() + theme(legend.text=element_text(size = 4))
pdf("output/Plots/GEx_TME_UMAP.pdf", width = 5, height = 6.5)
p1
p2
dev.off()

#marker expression
p3 <- lapply(c("CD3D", "CD19", "ITGAM", "PECAM1", "FAP"), function(x) FeaturePlot(tme_seu, features = x, order = T) + 
               scale_colour_gradientn(colours = viridis(256, option = "A")) + ArchR::theme_ArchR())
pdf("output/Plots/GEx_TME_markerExp.pdf", width = 5, height = 6)
p3
dev.off()

saveRDS(tme_seu, "rds/tme_seu.rds")
