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
arc <- addIterativeLSI(arc, useMatrix = "TileMatrix", name = "IterativeLSI", force = T) %>%
       addClusters(., reducedDims = "IterativeLSI", force = T) %>%
       addUMAP(., reducedDims = "IterativeLSI", force = T)

#call peaks
pathToMacs2 <- findMacs2()
arc <- addGroupCoverages(arc, groupBy = "Clusters") %>%
       addReproduciblePeakSet(., groupBy = "Clusters", 
                                 pathToMacs2 = pathToMacs2, 
                                 shift = -40, extsize = 80,
                                 method = "q", cutOff = 0.05) %>% 
       addPeakMatrix(.)

#calculate motif(homer) deviation
arc <- addMotifAnnotations(arc, motifSet = "homer", annoName = "homer")
arc <- addBgdPeaks(arc) %>% addDeviationsMatrix(., peakAnnotation = "homer")

#save rds
saveRDS(arc, "rds/arc.rds")
