#----------------------------------------------------------------------------
# scC1_makeArchRProj.R
#----------------------------------------------------------------------------
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(ArchR)
addArchRThreads(threads = 8)
addArchRGenome("hg19")

#assign sample tsv files
sampleName <- list.files(path = "../data/", pattern = "tsv.gz$")
names(sampleName) <- gsub(pattern = ".fragments.tsv.gz", replacement = "", sampleName)
outFile <- as.character(sampleName)

#make ArrowFiles
ArrowFiles <- character(length(sampleName))
ArrowFiles <- createArrowFiles(inputFiles = paste0("../data//", outFile), 
                               sampleNames = names(sampleName),
                               minTSS = 4, 
                               minFrags = 1000, 
                               addTileMat = TRUE, addGeneScoreMat = TRUE, 
                               force = TRUE)
pre_arc <- ArchRProject(ArrowFiles, 
                        outputDirectory = "output", 
                        copyArrows = F)

#epithelial cell ids
epiCellNames <- readRDS("../rds/epiCellNames.rds")
epi_arc <- pre_arc[epiCellNames]

#save rds
saveRDS(epi_arc, "rds/arc.rds")
