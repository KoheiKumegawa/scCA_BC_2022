#----------------------------------------------------------------------------
# scAR1_Projection_TMEtoTotal.R
#----------------------------------------------------------------------------
library(ArchR)

arc <- readRDS("rds/arc2.rds")
tme_seu <- readRDS("rds/tme_seu.rds")

#Cross platform linkgage
arc <- addGeneIntegrationMatrix(
  ArchRProj = arc, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = tme_seu,
  addToArrow = FALSE,
  groupRNA = "cellType",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

saveRDS(arc, "rds/arc_projection1.rds")
