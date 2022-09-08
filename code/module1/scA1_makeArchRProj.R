#----------------------------------------------------------------------------
# scA1_makeArchRProj.R
#----------------------------------------------------------------------------
#Kohei Kumegawa, Japanese Foundation for Cancer Research
#system("mkdir data code rds ref output output/Plots output/Tables")
#Move tsv file into data directory
library(ArchR)
addArchRThreads(threads = 8)
addArchRGenome("hg19")

#assign sample tsv files
sampleName <- list.files(path = "data/", pattern = "tsv.gz")
names(sampleName) <- gsub(pattern = ".fragments.tsv.gz", replacement = "", sampleName)
outFile <- as.character(sampleName)

#make ArrowFiles
ArrowFiles <- character(length(sampleName))
ArrowFiles <- createArrowFiles(inputFiles = paste0("data//", outFile), 
                               sampleNames = names(sampleName),
                               minTSS = 4, 
                               minFrags = 1000, 
                               addTileMat = TRUE, addGeneScoreMat = TRUE, 
                               force = TRUE)

#infer doublets
doubScores <- addDoubletScores(ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)

#make pre-filtered ArchRProject
pre_arc <- ArchRProject(ArrowFiles, 
                        outputDirectory = "output", 
                        copyArrows = F)
pre_arc <- filterDoublets(pre_arc)

#save rds
saveRDS(pre_arc, "rds/pre_arc.rds")
