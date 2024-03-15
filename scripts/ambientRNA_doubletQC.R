####################################################################################
# Run QC on drprnull42d and w111842d cells separately 
# Estimate per-cell QC metrics, ambient RNA correction, and doublet detection using SCTK-QC and SoupX
# Outputs AnnData .h5ad files for both samples separately
# Code adpated from sanbomics, sc-bestpractices, and OSCA
####################################################################################

## Load required packages

library(SingleCellExperiment)
library(tidyverse)
library(singleCellTK)
library(zellkonverter)
library(SoupX)
library(reticulate)
library(Seurat)
library(viridis)
adata <- import("anndata")
#scrublet <- import("scrublet")
set.seed(12345)


## Set up I/O

inputdir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/cellranger/"
outputdir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/sctk_qc/"
dir.create(outputdir,showWarnings = FALSE)

mad_outlier <- function(sobj, metric, nmads){
  M <- sobj@meta.data[[metric]]
  median_M <- median(M, na.rm = TRUE)
  mad_M <- mad(M, na.rm = TRUE)
  outlier <- (M < (median_M - nmads * mad_M)) | (M > (median_M + nmads * mad_M))
  return(outlier)
}

read_raw_sobj <- function(sample_id){
  inputdir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/cellranger/"
  path <- paste0(inputdir,sample_id, "/outs/filtered_feature_bc_matrix/")
  sobj <- Read10X(data.dir = path)
  sobj <- CreateSeuratObject(counts = sobj, min.cells = 0, min.features = 200)
  sobj$sample_id <- sample_id
  
  #add QC metrics
  sobj$log1p_total_counts <- log1p(sobj@meta.data$nCount_RNA)
  sobj$log1p_n_genes_by_counts <- log1p(sobj@meta.data$nFeature_RNA)
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt:")

  return(sobj)
}

filter_by_counts <- function(sobj){
  #find outliers and subset
  bool_vector <- !mad_outlier(sobj, 'log1p_total_counts', 5) & !mad_outlier(sobj, 'log1p_n_genes_by_counts', 5) & !mad_outlier(sobj, 'percent.mt', 3) 
  sobj <- subset(sobj, cells = which(bool_vector))
  sobj <- subset(sobj, subset = percent.mt < 8)
  return(sobj)
}


samples <- c('w1118_42D', 'drprnull_42D')
data_list <- sapply(samples, read_raw_sobj)


ggplot(as.data.frame(data_list[1]$w1118_42D@meta.data), aes(x=nCount_RNA,y=nFeature_RNA, color=percent.mt))+
  geom_point()+
  scale_color_viridis_c(option = "magma")+
  theme_classic()

ggplot(as.data.frame(data_list[2]$drprnull_42D@meta.data), aes(x=nCount_RNA,y=nFeature_RNA, color=percent.mt))+
  geom_point()+
  scale_color_viridis_c(option = "magma")+
  theme_classic()

data_list <- sapply(data_list, filter_by_counts)

ggplot(as.data.frame(data_list[1]$w1118_42D@meta.data), aes(x=nCount_RNA,y=nFeature_RNA, color=percent.mt))+
  geom_point()+
  scale_color_viridis_c(option = "magma")+
  theme_classic()

ggplot(as.data.frame(data_list[2]$drprnull_42D@meta.data), aes(x=nCount_RNA,y=nFeature_RNA, color=percent.mt))+
  geom_point()+
  scale_color_viridis_c(option = "magma")+
  theme_classic()


FeatureScatter(data_list[1]$w1118_42D, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="percent.mt")
FeatureScatter(data_list[1]$w1118_42D, feature1 = "nCount_RNA", feature2 = "percent.mt")

get_soup_groups <- function(sobj){
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(object = sobj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
  sobj <- RunUMAP(sobj, dims = 1:10)
  
  #return(sobj@meta.data[['seurat_clusters']])
  return(sobj)
  
}

add_soup_groups <- function(sobj){
  sobj <- get_soup_groups(sobj)
  sobj$soup_group <- sobj@meta.data[['seurat_clusters']]
  return(sobj)
}

data_list <- sapply(data_list, add_soup_groups)


data_list[1]$w1118_42D[[]]

make_soup <- function(sobj){
  inputdir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/cellranger/"
  sample_id <- as.character(sobj$sample_id[1]) #e.g, Lung1
  path <- paste0(inputdir,sample_id, "/outs/raw_feature_bc_matrix/")
  raw <- Read10X(data.dir = path)
  
  sc = SoupChannel(raw,sobj@assays$RNA@counts)
  sc = setClusters(sc,sobj$soup_group)
  sc = autoEstCont(sc, doPlot=FALSE)
  out = adjustCounts(sc, roundToInt = TRUE)
  
  #optional keep original
  sobj[["original.counts"]] <- CreateAssayObject(counts = sobj@assays$RNA@counts)
  
  sobj@assays$RNA@counts <- out
  
  return(sobj)
  
}

data_list <- sapply(data_list, make_soup)

sum(data_list[1]$w1118_42D@assays$original.counts@counts)
sum(data_list[1]$w1118_42D@assays$RNA@counts)/sum(data_list[1]$w1118_42D@assays$original.counts@counts)*100

sum(data_list[2]$drprnull_42D@assays$original.counts@counts)
sum(data_list[2]$drprnull_42D@assays$RNA@counts)/sum(data_list[2]$drprnull_42D@assays$original.counts@counts)*100


DimPlot(data_list[2]$drprnull_42D)
data_list[2]$drprnull_42D@meta.data

## Load raw w1118_42d counts


w1118_42d <- importCellRangerV3(
  cellRangerDirs = inputdir,
  sampleDirs = "w1118_42D",
  sampleNames = "w1118_42D",
  dataType = "filtered")


w1118_42d_samples <- colData(w1118_42d)$sample

## Run SCTK-QC on w1118_42d cells

w1118_42d <- runCellQC(w1118_42d, 
                       algorithms = c("QCMetrics","doubletFinder", "decontX","cxds", "bcds", "cxds_bcds_hybrid"), 
                       sample = w1118_42d_samples)





drprnull_42d <- importCellRangerV3(
    cellRangerDirs = inputdir,
    sampleDirs = "drprnull_42D",
    sampleNames = "drprnull_42D",
    dataType = "filtered")


drprnull_42d_samples <- colData(drprnull_42d)$sample

## Normalize and log1p transform counts

drprnull_42d <- runNormalization(
  inSCE = drprnull_42d,
  useAssay = "counts",
  outAssayName = "logcounts",
  normalizationMethod = "LogNormalize",
  scale = FALSE,
  seuratScaleFactor = 10000,
  transformation = "log1p",
  pseudocountsBeforeNorm = NULL,
  pseudocountsBeforeTransform = NULL,
  trim = NULL,
  verbose = TRUE
  )


## Run SCTK-QC
drprnull_42d <- runCellQC(drprnull_42d, 
                 algorithms = c("QCMetrics","doubletFinder", "decontX","cxds", "bcds", "cxds_bcds_hybrid"), 
                 sample = drprnull_42d_samples)

drprnull_42d <- runCellQC(drprnull_42d, 
                          algorithms = c("soupX"), 
                          sample = drprnull_42d_samples)

## Run SoupX


drprnull_42d_soupxc <- Seurat::Read10X(file.path(inputdir, "drprnull_42D","outs","raw_feature_bc_matrix"))
drprnull_42d_soupxc <- autoEstCont(drprnull_42d_soupxc)
out = adjustCounts(sc)




writeH5AD(
  drprnull_42d,
  paste0(outputdir,"drprnull_42d.h5ad"),
  X_name = "counts",
  skip_assays = FALSE,
  compression = "none"
)






writeH5AD(
  w1118_42d,
  paste0(outputdir,"w1118_42d.h5ad"),
  X_name = "counts",
  skip_assays = FALSE,
  compression = "none"
)








