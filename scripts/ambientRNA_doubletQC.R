####################################################################################
# Run SCTK-QC on drprnull42d and w111842d cells separately 
# Estimate per-cell QC metrics, ambient RNA correction, and doublet detection using SCTK-QC
# Outputs AnnData .h5ad files for both samples separately
####################################################################################


library(SingleCellExperiment)
library(tidyverse)
library(singleCellTK)
library(zellkonverter)

library(SoupX)

library(reticulate)
adata <- import("anndata")
scrublet <- import("scrublet")

## Read preprocessed h5ad file as SCE 

inputdir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/cellranger/"
outputdir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/sctk_qc/"
dir.create(outputdir,showWarnings = FALSE)



drprnull_42d <- importCellRangerV3(
    cellRangerDirs = inputdir,
    sampleDirs = "drprnull_42D",
    sampleNames = "drprnull_42D",
    dataType = "filtered")


drprnull_42d_samples <- colData(drprnull_42d)$sample


set.seed(12345)
drprnull_42d <- runCellQC(drprnull_42d, 
                 algorithms = c("QCMetrics","doubletFinder", "decontX","cxds", "bcds", "cxds_bcds_hybrid"), 
                 sample = drprnull_42d_samples)


writeH5AD(
  drprnull_42d,
  paste0(outputdir,"drprnull_42d.h5ad"),
  X_name = "counts",
  skip_assays = FALSE,
  compression = "none"
)


##### QC W1118_42sd
w1118_42d <- importCellRangerV3(
  cellRangerDirs = inputdir,
  sampleDirs = "w1118_42D",
  sampleNames = "w1118_42D",
  dataType = "filtered")


w1118_42d_samples <- colData(w1118_42d)$sample


set.seed(12345)
w1118_42d <- runCellQC(w1118_42d, 
                          algorithms = c("QCMetrics","doubletFinder", "decontX","cxds", "bcds", "cxds_bcds_hybrid"), 
                          sample = w1118_42d_samples)


writeH5AD(
  w1118_42d,
  paste0(outputdir,"w1118_42d.h5ad"),
  X_name = "counts",
  skip_assays = FALSE,
  compression = "none"
)








