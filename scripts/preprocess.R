####################################################################################
# Run QC on drprnull42d and w111842d cells separately 
# Estimate per-cell QC metrics, ambient RNA correction, empty droplet and doublet detection using SCTK-QC and SoupX
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
library(robustbase)
library(scater)
library(scran)


#adata <- import("anndata")
#scrublet <- import("scrublet")
set.seed(12345)


## Set up I/O

basedir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/"
inputdir <- paste0(basedir,"cellranger/")
outputdir <- paste0(basedir,"preprocess/")
dir.create(outputdir,showWarnings = FALSE)

## Function definitions

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

flag_outliers <- function(sobj){
  ncount_outlier <- mad_outlier(sobj, 'log1p_total_counts', 5) 
  names(ncount_outlier) <- colnames(sobj)
  sobj <- AddMetaData(object = sobj,metadata = ncount_outlier,col.name = 'nCount_outlier')
  
  nfeature_outlier <- mad_outlier(sobj, 'log1p_n_genes_by_counts', 5) 
  names(nfeature_outlier) <- colnames(sobj)
  sobj <- AddMetaData(object = sobj,metadata = nfeature_outlier,col.name = 'nFeature_outlier')
  
  mito_outlier <- mad_outlier(sobj, 'percent.mt', 3) 
  names(mito_outlier) <- colnames(sobj)
  sobj <- AddMetaData(object = sobj,metadata = mito_outlier,col.name = 'mt_outlier')
}

plot_features <- function(sobjmetadata){
  p0 <- ggplot(as.data.frame(sobjmetadata), aes(x=nCount_RNA,y=nFeature_RNA, color=percent.mt))+
    geom_point()+
    scale_color_viridis_c(option = "magma")+
    theme_classic(base_size = 16)
  
  p1 <- ggplot(as.data.frame(sobjmetadata),aes(x=nCount_RNA))+
    geom_histogram()+
    theme_classic(base_size = 16)
  
  p2 <- ggplot(as.data.frame(sobjmetadata),aes(x=nFeature_RNA))+
    geom_histogram()+
    theme_classic(base_size = 16)
  
  p3 <- ggplot(as.data.frame(sobjmetadata),aes(x=log1p_total_counts))+
    geom_histogram()+
    theme_classic(base_size = 16)
  
  p4 <- ggplot(as.data.frame(sobjmetadata),aes(x=log1p_n_genes_by_counts))+
    geom_histogram()+
    theme_classic(base_size = 16)
  
  p5 <- ggplot(as.data.frame(sobjmetadata), aes(x=percent.mt))+
    geom_histogram(bins=20)+
    theme_classic(base_size = 16)
  
  #featureplots <- list(p0,p1,p2,p3,p4,p5)
  p0+p1+p2+p3+p4+p5 
  
  
}

plot_outliers <- function(sobjmetadata){
  ncount_outlier_plot <- ggplot(data = sobjmetadata, aes(x=sample_id,y=log1p_n_genes_by_counts,color=nFeature_outlier))+
    #geom_violin(alpha=0.5)+
    geom_jitter(alpha=0.5)+
    theme_classic(base_size = 16)
  nfeature_outlier_plot <- ggplot(data = sobjmetadata, aes(x=sample_id,y=log1p_total_counts,color=nCount_outlier))+
    #geom_violin(alpha=0.5)+
    geom_jitter(alpha=0.5)+
    theme_classic(base_size = 16)
  mt_outlier_plot <- ggplot(data = sobjmetadata, aes(x=sample_id,y=percent.mt,color=mt_outlier))+
    #geom_violin(alpha=0.5)+
    geom_jitter(alpha=0.5)+
    theme_classic(base_size = 16)
  return(ncount_outlier_plot+nfeature_outlier_plot+mt_outlier_plot)
}

filter_by_counts <- function(sobj){
  #find outliers and subset
  bool_vector <- !mad_outlier(sobj, 'log1p_total_counts', 5) & !mad_outlier(sobj, 'log1p_n_genes_by_counts', 5) & !mad_outlier(sobj, 'percent.mt', 3) 
  sobj <- subset(sobj, cells = which(bool_vector))
  sobj <- subset(sobj, subset = percent.mt < 1) ## Remove cells with percent.mt > 1% because mito counts in snRNAseq indicates incomplete cell lysis.
  return(sobj)
}


## Read raw counts and identify outliers
samples <- c('w1118_42D', 'drprnull_42D')
data_list <- sapply(samples, read_raw_sobj)
data_list <- sapply(data_list, flag_outliers)


## Plot features before filtering
plot_features(data_list[1]$w1118_42D@meta.data)
plot_features(data_list[2]$drprnull_42D@meta.data)


## Plot outliers before filtering
plot_outliers(data_list[1]$w1118_42D@meta.data)
plot_outliers(data_list[2]$drprnull_42D@meta.data)

## Filter low quality cells
data_list <- sapply(data_list, filter_by_counts)

## Plot features after filtering
plot_features(data_list[1]$w1118_42D@meta.data)
plot_features(data_list[2]$drprnull_42D@meta.data)

## Normalization
data_list <- sapply(data_list, SCTransform, vars.to.regress ="percent.mt")


ggplot(data=data_list[1]$w1118_42D@meta.data,aes(x=log10(nCount_SCT)))+
  geom_histogram()+
  theme_classic(base_size = 16)

ggplot(data=data_list[2]$drprnull_42D@meta.data,aes(x=nCount_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)

ggplot(data=data_list[1]$w1118_42D@meta.data,aes(x=nFeature_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)

ggplot(data=data_list[2]$drprnull_42D@meta.data,aes(x=nFeature_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)

ggplot(as.data.frame(data_list[1]$w1118_42D@meta.data), aes(x=nCount_SCT,y=nFeature_SCT, color=percent.mt))+
  geom_point()+
  ggtitle("W1118_42D")+
  scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)

ggplot(as.data.frame(data_list[2]$drprnull_42D@meta.data), aes(x=nCount_SCT,y=nFeature_SCT, color=percent.mt))+
  geom_point()+
  ggtitle("drprnull_42D")+
  scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)


## Clustering 

run_clustering <- function(sobj){
  sobj <- RunPCA(sobj, verbose = FALSE)
  sobj <- RunUMAP(sobj, dims = 1:30, verbose = FALSE)
  
  sobj <- FindNeighbors(sobj, dims = 1:30, verbose = FALSE)
  
  sobj <- FindClusters(sobj, verbose = FALSE, resolution=0.02)
  sobj <- FindClusters(sobj, verbose = FALSE, resolution=0.2)
  sobj <- FindClusters(sobj, verbose = FALSE, resolution=0.8)
  sobj <- FindClusters(sobj, verbose = FALSE, resolution=1)
  sobj <- FindClusters(sobj, verbose = FALSE, resolution=1.5)
  sobj <- FindClusters(sobj, verbose = FALSE, resolution=2.5)
  sobj <- FindClusters(sobj, verbose = FALSE, resolution=3.5)
  sobj <- FindClusters(sobj, verbose = FALSE, resolution=4.5)
  sobj <- FindClusters(sobj, verbose = FALSE, resolution=8)
  sobj <- FindClusters(sobj, verbose = FALSE, resolution=12)
  return(sobj)
}

data_list <- sapply(data_list, run_clustering)

cluster_res <- c("SCT_snn_res.0.02","SCT_snn_res.0.2","SCT_snn_res.0.8","SCT_snn_res.1","SCT_snn_res.1.5","SCT_snn_res.2.5",
                 "SCT_snn_res.3.5","SCT_snn_res.4.5","SCT_snn_res.8","SCT_snn_res.12")

markers <- c("elav","lncRNA:noe","VAChT","VGlut","Gad1","Vmat","SerT","Tdc2","ple", # neurons
             "repo","lncRNA:CR34335","alrm","wrapper","Indy","moody",#glia 
             "ninaC",	"trp",	"trpl", #photoreceptors
             "Hml", #hemocytes
             "ppl",#fatbody
             "drpr")


for(cr in cluster_res){
  p1 <- DimPlot(data_list[1]$w1118_42D,group.by = cr) &
    theme_classic(base_size = 16)
  
  ggsave(p1,filename=paste0(outputdir,"figures/","w1118_42d",cr,"_dimplot.pdf"),width=10,height=10)
  
  p2 <- DimPlot(data_list[2]$drprnull_42D,group.by = cr)&
    theme_classic(base_size = 16)
  
  ggsave(p2,filename=paste0(outputdir,"figures/","drprnull_42d",cr,"_dimplot.pdf"),width=10,height=10)
  
  p3 <- DotPlot(data_list[1]$w1118_42D, features = markers, group.by = cr)&
    scale_color_viridis_c(option = "magma",direction = -1)&
    theme_classic(base_size = 16) + RotatedAxis()
  
  ggsave(p3,filename=paste0(outputdir,"figures/","w1118_42d",cr,"_dotplot.pdf"), width=16,height=10)
  
  p4 <- DotPlot(data_list[2]$drprnull_42D, features = markers, group.by = cr)&
    scale_color_viridis_c(option = "magma",direction = -1)&
    theme_classic(base_size = 16)+ RotatedAxis()
  
  ggsave(p4,filename=paste0(outputdir,"figures/","drprnull_42d",cr,"_dotplot.pdf"), width=16,height=10)
  
}



## Ambient RNA correction using SoupX

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
  #sobj <- get_soup_groups(sobj)
  sobj$soup_group <- sobj@meta.data[['SCT_snn_res.0.8']] ## Use Clusters from SCT normalized assay
  return(sobj)
}

data_list <- sapply(data_list, add_soup_groups)

make_soup <- function(sobj){
  inputdir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/cellranger/"
  sample_id <- as.character(sobj$sample_id[1]) 
  path <- paste0(inputdir,sample_id, "/outs/raw_feature_bc_matrix/")
  raw <- Read10X(data.dir = path)
  
  sc = SoupChannel(raw,sobj@assays$RNA@counts)
  sc = setClusters(sc,sobj$soup_group)
  sc = autoEstCont(sc, doPlot=FALSE)
  out = adjustCounts(sc, roundToInt = TRUE)
  

  
  #sobj[["original.counts"]] <- CreateAssayObject(counts = sobj@assays$RNA@counts)
  #sobj@assays$RNA@counts <- out
  sobj[["soupXcounts"]] <- CreateAssayObject(counts = out)   #save SoupX corrected counts
  
  
  return(sobj)
  
}

data_list <- sapply(data_list, make_soup)

###  Check SoupX correction - % of RNA removed

sum(data_list[1]$w1118_42D@assays$RNA@counts)
sum(data_list[1]$w1118_42D@assays$soupXcounts@counts)/sum(data_list[1]$w1118_42D@assays$RNA@counts)*100

sum(data_list[2]$drprnull_42D@assays$RNA@counts)
sum(data_list[2]$drprnull_42D@assays$soupXcounts@counts)/sum(data_list[2]$drprnull_42D@assays$RNA@counts)*100

## Convert Seurat objects to SingleCellExperiment objects
data_list_sce <- sapply(data_list, as.SingleCellExperiment)



DefaultAssay(data_list[2]$drprnull_42D) <- "original.counts"
DefaultAssay(data_list[1]$w1118_42D) <- "original.counts"
FeaturePlot(data_list[2]$drprnull_42D, features=c("elav","repo","drpr"))
FeaturePlot(data_list[1]$w1118_42D, features=c("elav","repo","drpr"))

DefaultAssay(data_list[2]$drprnull_42D) <- "RNA"
DefaultAssay(data_list[1]$w1118_42D) <- "RNA"

FeaturePlot(data_list[1]$w1118_42D, features=c("repo","drpr"), blend = TRUE)
FeaturePlot(data_list[2]$drprnull_42D, features=c("repo","drpr"), blend = TRUE)



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








