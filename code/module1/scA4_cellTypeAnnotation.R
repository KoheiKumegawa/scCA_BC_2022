#----------------------------------------------------------------------------
# scA4_cellTypeAnnotation.R
#----------------------------------------------------------------------------
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(ArchR)
library(ComplexHeatmap)
library(circlize)
addArchRThreads(threads = 8)

arc <- readRDS("rds/arc.rds")

#----- clusters
p1 <- plotEmbedding(arc, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", plotAs = "points", size = 1) 
p2 <- plotEmbedding(arc, colorBy = "cellColData", name = "Sample", embedding = "UMAP", plotAs = "points", size = 1)
plotPDF(p1, p2, name = "CA_UMAP.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)
cM <- table(arc$Sample, arc$Clusters)[, paste0("C", seq_along(unique(arc$Clusters)))]
write.csv(cM, "output/Tables/CA_SampleClusterMatrix.csv")

#----- gene score
markers <- read.csv("ref/markers.csv", header = F) %>% as.character()
arc <- addImputeWeights(arc)
GeneScoreClusters <- getMarkerFeatures(
  ArchRProj = arc, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

#heatmap
idx <- which(mcols(GeneScoreClusters)$name %in% markers)
mtx <- assays(GeneScoreClusters[idx,])$Mean %>% `rownames<-`(., mcols(GeneScoreClusters)$name[idx])
mtx <- t(scale(t(mtx)))

fh <- function(x) hclust(dist(x), method="ward.D2")
col_fun1 <- colorRamp2(c(-2,-1,0,1,2), paletteContinuous(set = "horizonExtra", n = 5))
ht1 <- Heatmap(mtx, name = "Mean GS z-score", cluster_rows = fh, cluster_columns = F, show_row_dend = F, col = col_fun1)
p3 <- draw(ht1)
plotPDF(p3, name = "CA_markersGS_Heatmap.pdf", ArchRProj = arc, addDOC = FALSE, width = 7, height = 10)

#genome track
p4 <- plotBrowserTrack(arc, 
                       groupBy = "Clusters", 
                       geneSymbol = c("EPCAM","VIM","KRT18","FLT4","THY1","CD3D","PAX5","IGLL5","ITGAX"), 
                       upstream = 10000, # 10kb
                       downstream = 10000, # 10kb
                       plotSummary = c("bulkTrack","geneTrack")
                       )
plotPDF(p4, name = "CA_markersGS_browserTrack.pdf", ArchRProj = arc, addDOC = FALSE, width = 7, height = 10)

#UMAP overlay
p5 <- plotEmbedding(arc, 
                    colorBy = "GeneScoreMatrix", 
                    name = markers, 
                    embedding = "UMAP", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = 1)
plotPDF(p5, name = "CA_markerGS_UMAPoverlay.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#----- motif score
VarMotifs <- getVarDeviations(arc, name = "homerMatrix", plot = F)
MotifScoreClusters <- getMarkerFeatures(arc, 
                                        useMatrix = "homerMatrix", 
                                        groupBy = "Clusters",
                                        bias = c("TSSEnrichment", "log10(nFrags)"),
                                        testMethod = "wilcoxon", 
                                        useSeqnames = "z")
idy <- which(mcols(MotifScoreClusters)$name %in% VarMotifs$name[c(1:50)])
mtx2 <- assays(MotifScoreClusters[idy,])$Mean %>% `rownames<-`(., mcols(MotifScoreClusters)$name[idy])

col_fun2 <- colorRamp2(c(-5,-2.5,0,2,5), paletteContinuous(set = "solarExtra", n = 5))
ht2 <- Heatmap(mtx2, name = "Mean motif score", cluster_rows = fh, cluster_columns = F, show_row_dend = F, col = col_fun2)
p6 <- draw(ht2)
plotPDF(p6, name = "CA_topVarMotifs_Heatmap.pdf", ArchRProj = arc, addDOC = FALSE, width = 8, height = 8)

#----- cell type assignment
#C1-13, epithelial
#C14, endothelial
#C15-16, fibroblast
#C17-21, T
#C22, B
#C23, Plasma
#C24-25, Myeloid
cellType <- list("Ep" = paste0("C", c(1:13)), 
                 "En" = paste0("C", 14),  
                 "Fb" = paste0("C", c(15:16)),
                 "Tc" = paste0("C", c(17:21)),
                 "Bc" = paste0("C", 22),
                 "Pc" = paste0("C", 23),
                 "My" = paste0("C", 24:25))
idz <- lapply(cellType, function(x) which(arc$Clusters %in% x))
arc$cellType <- rep(NA, nrow(arc@cellColData))
for(i in names(idz)){arc$cellType[idz[[i]]] <- i}
p7 <- plotEmbedding(arc, colorBy = "cellColData", name = "cellType", embedding = "UMAP", plotAs = "points", size = 1)
plotPDF(p7, name = "CA_UMAP_cellTypeAssign.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#----- cell type specific genes
GeneScoreCellType <- getMarkerFeatures(
  ArchRProj = arc, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "cellType",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

CellTypeFDR <- lapply(names(cellType), function(x){
  a <- assays(GeneScoreCellType)$FDR[,x]
  a <- -log10(a)
  names(a) <- mcols(GeneScoreCellType)$name
  out <- sort(a, decreasing = T)
}) %>% `names<-`(., names(cellType))
lapply(CellTypeFDR, head)

#----- peak profile
df <- reshape2::melt(table(arc@peakSet$peakType))
df$Var1 <- factor(df$Var1, levels = c("Exonic", "Distal", "Intronic", "Promoter"))
p8 <- ggplot(df, aes(x="", y= value, fill= Var1)) + 
       geom_bar(stat="identity", width=1, color="white") + 
       coord_polar("y", start=0) +
       theme_void() + scale_fill_brewer(palette="Dark2")
plotPDF(p8, name = "CA_peakGenomicPosition.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

tcga_bcPeaks <- import.bed("ref/TCGA-BRCA-PeaksHg19.bed")
o <- findOverlaps(arc@peakSet, tcga_bcPeaks)
jfcr_overlap_no <- unique(queryHits(o)) %>% length()
tcga_overlap_no <- unique(subjectHits(o)) %>% length()
mtx <- matrix(c(length(arc@peakSet), jfcr_overlap_no, length(tcga_bcPeaks), tcga_overlap_no), ncol = 2, nrow = 2)
mtx <- rbind(mtx, (mtx[1,] - mtx[2,])) %>% `colnames<-`(.,c("scATAC", "TCGA")) %>% `rownames<-`(., c("Total", "Overlaps", "Nonoverlaps"))
p9 <- ggplot(reshape2::melt(mtx[c(3,2),]), aes(x = Var2, y = value, fill = Var1)) + geom_bar(stat = "identity") + 
      theme_ArchR() + scale_fill_manual(values = c("lightgray", "darkorange"))
pdf("output/Plots/CA_peakTCGAoverlaps.pdf", width = 3, height = 5)
p9
dev.off()

p10 <- plotGroups(arc, groupBy = "Sample", name = "FRIP", plotAs = "violin", alpha = 0.5)
p11 <- plotGroups(arc, groupBy = "Clusters", name = "FRIP", plotAs = "violin", alpha = 0.5)
plotPDF(p10, p11, name = "CA_FRIP.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)

#----- identify DARs
PeakCellType <- getMarkerFeatures(arc, 
                                  useMatrix = "PeakMatrix", 
                                  groupBy = "cellType",
                                  bias = c("TSSEnrichment", "log10(nFrags)"),
                                  testMethod = "wilcoxon")


DAR <- getMarkers(PeakCellType, cutOff = "FDR < 0.01 & Log2FC >= 1")
DARmtx <- plotMarkerHeatmap(seMarker = PeakCellType, cutOff = "FDR < 0.01 & Log2FC >= 1", returnMatrix = T)
col_fun3 = colorRamp2(c(-2,-1,0,1,2), viridis::viridis(5))
ht3 <- Heatmap(DARmtx, name = "ATAC Z-score", cluster_rows = fh, cluster_columns = F, 
               show_row_dend = F, show_row_names = F, col = col_fun3, use_raster = T)
p12 <- draw(ht3)
plotPDF(p12, name = "CA_DAR_Heatmap.pdf", width = 6, height = 8, ArchRProj = arc, addDOC = FALSE)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = PeakCellType,
  ArchRProj = arc,
  peakAnnotation = "homer",
  cutOff = "FDR < 0.01 & Log2FC >= 1"
)
p13 <- plotEnrichHeatmap(enrichMotifs, n = 12, clusterCols = F, transpose = TRUE)
plotPDF(p27, name = "DAR-Motif-Heatmap.pdf", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)

DAR_gr <- getMarkers(PeakCellType, cutOff = "FDR < 0.01 & Log2FC >= 1", returnGR = T)
lapply(names(DAR_gr), function(x){
  gr <- DAR_gr[[x]]
  mcols(gr) <- NULL
  rtracklayer::export.bed(gr, paste0("output/output_bed/celltype_DAR/", x, "_DAR.bed"))
  return(NULL)
})

#----- epithelial cell id
epiCellNames <- arc$cellNames[arc$cellType == "Ep"]
saveRDS(epiCellNames, "rds/epiCellNames.rds")

saveRDS(arc, "rds/arc2.rds")
