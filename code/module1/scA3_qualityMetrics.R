#----------------------------------------------------------------------------
# scA3_qualityMetrics.R
#----------------------------------------------------------------------------
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(ArchR)
library(dplyr)
addArchRThreads(threads = 8)

pre_arc <- readRDS("rds/pre_arc.rds")
arc <- readRDS("rds/arc.rds")

df <- data.frame(nFrags = log10(pre_arc$nFrags), TSSe = pre_arc$TSSEnrichment)
p1 <- ggplot(df, aes(x = nFrags, y = TSSe)) + geom_hex(bins = 100) + scale_fill_viridis_c(trans = "log") + 
      theme_ArchR() + geom_hline(yintercept = 8, lty = "dashed") + geom_vline(xintercept = log10(3000), lty = "dashed") +
      labs(x = "Log10 Unique Fragments", y = "TSS Enrichment")
p2 <- plotTSSEnrichment(arc)
p3 <- plotFragmentSizes(arc)
p4 <- plotGroups(arc, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p5 <- plotGroups(arc, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)",plotAs = "violin",alpha = 0.4, addBoxPlot = TRUE)

plotPDF(p1,p2,p3,p4,p5, name = "CA_QualityPlots.pdf", ArchRProj = arc, addDOC = FALSE, width = 5, height = 5)
