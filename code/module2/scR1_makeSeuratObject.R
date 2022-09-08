#----------------------------------------------------------------------------
# scR1_makeSeuratObject.R
#----------------------------------------------------------------------------
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(Seurat)

#read scRNA-seq data
data_ls <- data.frame(name = list.files("data/GSE161529_RAW/"),
                      sample = stringr::str_split(list.files("data/GSE161529_RAW/"), pattern = "_", simplify = T)[,2],
                      detail = stringr::str_split(list.files("data/GSE161529_RAW/"), pattern = "_", simplify = T)[,4])
data <- as.list(data_ls[which(data_ls$sample %in% c("TN", "HER2", "ER")), 1])
d <- lapply(data, function(x){Read10X(paste0("data/GSE161529_RAW/", x))})
names(d) <- unlist(data)

#create seurat object
seu_ls <- parallel::mclapply(names(d), function(x){
  a <- CreateSeuratObject(counts = d[[x]])
  a$orig.ident <- x
  a <- RenameCells(a, add.cell.id = a$orig.ident)
  return(a)
  }, mc.cores = 8)
pre_seu <- merge(x = seu_ls[[1]], y = seu_ls[c(2:length(seu_ls))])
pre_seu[["percent.mt"]] <- PercentageFeatureSet(pre_seu, pattern = "^MT-")

p1 <- VlnPlot(pre_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
p2 <- FeatureScatter(pre_seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(pre_seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("output/Plots/GEx_QCplots.pdf")
p1
p2
p3
dev.off()

seu <- subset(pre_seu, subset = nFeature_RNA > 500 & percent.mt < 20)
saveRDS(seu, "rds/seu.rds")
