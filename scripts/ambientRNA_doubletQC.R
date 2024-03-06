library(SingleCellExperiment)
library(tidyverse)
library(singleCellTK)
#BiocManager::install("singleCellTK")



## Read preprocessed h5ad file as SCE 

inputdir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/cellranger/"

sce <- importCellRangerV3(
  cellRangerDirs = c(inputdir,inputdir),
  sampleDirs = c("drprnull_42D","w1118_42D"),
  sampleNames = c("drprnull_42D","w1118_42D"),
  dataType = "filtered")


#drprnull_42d <- importAnnData(sampleDirs = "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/preprocess",
#                     sampleNames = "drprnull_42D")

  cellRangerDirs = c(inputdir,inputdir),
  sampleDirs = c("drprnull_42D","w1118_42D"),
  sampleNames = c("drprnull_42D","w1118_42D"),
  dataType = "filtered")

samples <- colData(sce)$sample

#sce <- getUMAP(inSCE = sce, useAssay = "counts", logNorm = TRUE, sample = samples)

set.seed(12345)
sce <- runCellQC(sce, 
                 algorithms = c("QCMetrics", "scrublet", "doubletFinder", "decontX", "soupX"), 
                 sample = samples)

exportSCEtoAnnData(
  sce,
  useAssay = "counts",
  outputDir = inputdir,
  prefix = "combined",
  overwrite = TRUE,
  compression = "None",
  compressionOpts = NULL,
  forceDense = FALSE
)

sce <- sampleSummaryStats(sce, sample = samples, simple = FALSE)
getSampleSummaryStatsTable(sce, statsName = "qc_table")

decontxResults <- plotDecontXResults(
  inSCE = sce, sample = colData(sce)$sample, 
  reducedDimName = "decontX_drprnull_42D_UMAP", combinePlot = "all",
  titleSize = 8,
  axisLabelSize = 8,
  axisSize = 10,
  legendSize = 5,
  legendTitleSize = 7,
  relWidths = c(0.5, 1, 1),
  sampleRelWidths = c(0.5, 1, 1),
  labelSamples = TRUE,
  labelClusters = FALSE
)

soupxResults <- plotSoupXResults(
  inSCE = sce, sample = colData(sce)$sample, 
  reducedDimName = "SoupX_UMAP_drprnull_42D", combinePlot = "all",
  titleSize = 8,
  axisLabelSize = 8,
  axisSize = 10,
  legendSize = 5,
  legendTitleSize = 7,
  labelClusters = FALSE
)

soupxResults
doubletFinderResults <- plotDoubletFinderResults(
  inSCE = sce,
  sample = colData(sce)$sample, 
  reducedDimName = "doubletFinder_UMAP",
  combinePlot = "all",
  titleSize = 13,
  axisLabelSize = 13,
  axisSize = 13,
  legendSize = 13,
  legendTitleSize = 13
)

doubletFinderResults

qcresults <- as.data.frame(colData(sce))
ggplot(qcresults, aes(x=soupX_contamination,y=decontX_contamination))+
  geom_point()


ggplot(qcresults, aes(x=decontX_contamination))+
  geom_histogram()




