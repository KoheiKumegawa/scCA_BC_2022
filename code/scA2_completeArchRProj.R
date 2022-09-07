#----------------------------------------------------------------------------
# scA2_completeArchRProj.R
#----------------------------------------------------------------------------
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(ArchR)
library(dplyr)
addArchRThreads(threads = 8)

#filter cells
pre_arc <- readRDS("rds/pre_arc.rds")
arc <- pre_arc[which(pre_arc$TSSEnrichment >= 8 & pre_arc$nFrags >= 3000)]

#clustering
arc <- addIterativeLSI(arc, useMatrix = "TileMatrix", name = "IterativeLSI") %>%
  addClusters(., reducedDims = "IterativeLSI", maxClusters = 30, force = T) %>%
  addUMAP(., reducedDims = "IterativeLSI", force = T)

#call peaks
pathToMacs2 <- findMacs2()
arc <- addGroupCoverages(arc, groupBy = "Clusters") %>%
  addReproduciblePeakSet(., groupBy = "Clusters", 
                         pathToMacs2 = pathToMacs2, 
                         method = "q", cutOff = 0.05) %>% 
  addPeakMatrix(.)

#calculate motif deviation
arc <- addMotifAnnotations(arc, motifSet = "homer", annoName = "homer")
arc <- addMotifAnnotations(arc, motifSet = "cisbp", annoName = "cisbp")
arc <- addBgdPeaks(arc) %>% addDeviationsMatrix(., peakAnnotation = "homer") %>% addDeviationsMatrix(., peakAnnotation = "cisbp")

#save rds
saveRDS(arc, "rds/arc.rds")
