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
library(SeuratDisk)


#adata <- import("anndata")
#scrublet <- import("scrublet")
set.seed(12345)


## Set up I/O

basedir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/"
inputdir <- paste0(basedir,"cellranger/")
outputdir <- paste0(basedir,"preprocess/")
dir.create(outputdir,showWarnings = FALSE)
dir.create(paste0(outputdir,"figures"),showWarnings = FALSE)

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

## Convert SCE object to Seurat for downstream analysis
w1118_42d <- setRowNames(w1118_42d,"feature_name") ## use gene symbols instead of flybase id
w1118_42d_seurat <- as.Seurat(w1118_42d,data = NULL)
w1118_42d_seurat[['decontXcounts']] <- CreateAssayObject(counts=assays(w1118_42d)$decontXcounts)
w1118_42d_seurat$nCount_RNA <- as.numeric(w1118_42d_seurat$nCount_originalexp)
w1118_42d_seurat$nFeature_RNA <-as.numeric( w1118_42d_seurat$nFeature_originalexp)

## Load raw drprnull_42d counts
drprnull_42d <- importCellRangerV3(
  cellRangerDirs = inputdir,
  sampleDirs = "drprnull_42D",
  sampleNames = "drprnull_42D",
  dataType = "filtered")

drprnull_42d_samples <- colData(drprnull_42d)$sample


## Run SCTK-QC on drprnull_42d cells
drprnull_42d <- runCellQC(drprnull_42d, 
                          algorithms = c("QCMetrics","doubletFinder", "decontX","cxds", "bcds", "cxds_bcds_hybrid"), 
                          sample = drprnull_42d_samples)

## Convert SCE object to Seurat for downstream analysis
drprnull_42d <- setRowNames(drprnull_42d,"feature_name") ## use gene symbols instead of flybase id
drprnull_42d_seurat <- as.Seurat(drprnull_42d,data = NULL)
drprnull_42d_seurat[['decontXcounts']] <- CreateAssayObject(counts=assays(drprnull_42d)$decontXcounts)
drprnull_42d_seurat$nCount_RNA <- as.numeric(drprnull_42d_seurat$nCount_originalexp)
drprnull_42d_seurat$nFeature_RNA <- as.numeric(drprnull_42d_seurat$nFeature_originalexp)

## Combine w1118_42d and drprnull_42d seurat objects into a list for easy analysis
data_list <- list("w1118_42d" = w1118_42d_seurat,"drprnull_42d" = drprnull_42d_seurat)

## Filter cells based on QC metrics 
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
  
  #add QC metrics
  sobj$log1p_total_counts <- log1p(sobj@meta.data$nCount_RNA)
  sobj$log1p_n_genes_by_counts <- log1p(sobj@meta.data$nFeature_RNA)
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt:")
  
  
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
  ncount_outlier_plot <- ggplot(data = sobjmetadata, aes(x=sample,y=log1p_n_genes_by_counts,color=nFeature_outlier))+
    #geom_violin(alpha=0.5)+
    geom_jitter(alpha=0.5)+
    theme_classic(base_size = 16)
  nfeature_outlier_plot <- ggplot(data = sobjmetadata, aes(x=sample,y=log1p_total_counts,color=nCount_outlier))+
    #geom_violin(alpha=0.5)+
    geom_jitter(alpha=0.5)+
    theme_classic(base_size = 16)
  mt_outlier_plot <- ggplot(data = sobjmetadata, aes(x=sample,y=percent.mt,color=mt_outlier))+
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
samples <- c('w1118_42d', 'drprnull_42d')
#data_list <- sapply(samples, read_raw_sobj)
data_list <- sapply(data_list, flag_outliers)


## Plot features before filtering
plot_features(data_list[1]$w1118_42d@meta.data)
plot_features(data_list[2]$drprnull_42d@meta.data)


## Plot outliers before filtering
plot_outliers(data_list[1]$w1118_42d@meta.data)
plot_outliers(data_list[2]$drprnull_42d@meta.data)

### Removing percent.mt > 1% cells in drpr sample removes more than half the samples. Finding a more relaxed threshold
drpr_percentmt_filter <- as.data.frame(data_list[2]$drprnull_42d@meta.data) %>%
  filter(percent.mt <5)
print(dim(drpr_percentmt_filter))
ggplot(drpr_percentmt_filter, aes(x=percent.mt))+
  geom_histogram()+
  theme_classic(base_size = 16)

## Filter low quality cells
data_list <- sapply(data_list, filter_by_counts)

## Plot features after filtering
plot_features(data_list[1]$w1118_42d@meta.data)
plot_features(data_list[2]$drprnull_42d@meta.data)

data_list[1]$w1118_42d ##9033 cells x 16507 genes 
data_list[2]$drprnull_42d ##6423 cells x 16507 genes 

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
  sobj <- get_soup_groups(sobj)
  sobj$soup_group <- sobj@meta.data[['seurat_clusters']]
  return(sobj)
}

data_list <- sapply(data_list, add_soup_groups)

make_soup <- function(sobj){
  inputdir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/cellranger/"
  sample_id <- as.character(sobj$sample[1]) 
  path <- paste0(inputdir,sample_id, "/outs/raw_feature_bc_matrix/")
  raw <- Read10X(data.dir = path)
  
  sc = SoupChannel(raw,sobj@assays$originalexp$counts)
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

sum(data_list[1]$w1118_42d@assays$originalexp$counts)
sum(data_list[1]$w1118_42d@assays$soupXcounts$counts)/sum(data_list[1]$w1118_42d@assays$originalexp@counts)*100

sum(data_list[2]$drprnull_42d@assays$originalexp$counts)
sum(data_list[2]$drprnull_42d@assays$soupXcounts@counts)/sum(data_list[2]$drprnull_42d@assays$originalexp@counts)*100


## Normalization
### Perform scTransform normalization on raw counts, soupX corrected counts, and decontX corrected counts
data_list <- sapply(data_list, SCTransform, vars.to.regress ="percent.mt",assay="originalexp", new.assay.name="SCT")
data_list <- sapply(data_list, SCTransform, vars.to.regress ="percent.mt",assay="soupXcounts", new.assay.name="soupX_SCT")
data_list <- sapply(data_list, SCTransform, vars.to.regress ="percent.mt",assay="decontXcounts", new.assay.name="decontX_SCT")

p1 <- ggplot(data=data_list[1]$w1118_42d@meta.data,aes(x=nCount_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)
p2 <- ggplot(data=data_list[1]$w1118_42d@meta.data,aes(x=nFeature_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)
p3 <- ggplot(as.data.frame(data_list[1]$w1118_42d@meta.data), aes(x=nCount_SCT,y=nFeature_SCT, color=percent.mt))+
  geom_point()+
  ggtitle("W1118_42D")+
  scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)

p1+p2+p3

p1 <- ggplot(data=data_list[1]$w1118_42d@meta.data,aes(x=nCount_soupX_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)
p2 <- ggplot(data=data_list[1]$w1118_42d@meta.data,aes(x=nFeature_soupX_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)
p3 <- ggplot(as.data.frame(data_list[1]$w1118_42d@meta.data), aes(x=nCount_soupX_SCT,y=nFeature_soupX_SCT, color=percent.mt))+
  geom_point()+
  ggtitle("W1118_42D")+
  scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)

p1+p2+p3


p1 <- ggplot(data=data_list[1]$w1118_42d@meta.data,aes(x=nCount_decontX_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)
p2 <- ggplot(data=data_list[1]$w1118_42d@meta.data,aes(x=nFeature_decontX_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)
p3 <- ggplot(as.data.frame(data_list[1]$w1118_42d@meta.data), aes(x=nCount_decontX_SCT,y=nFeature_decontX_SCT, color=percent.mt))+
  geom_point()+
  ggtitle("W1118_42D")+
  scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)

p1+p2+p3


p1 <- ggplot(data=data_list[2]$drprnull_42d@meta.data,aes(x=nCount_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)

p2 <- ggplot(data=data_list[2]$drprnull_42d@meta.data,aes(x=nFeature_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)

p3 <- ggplot(as.data.frame(data_list[2]$drprnull_42d@meta.data), aes(x=nCount_SCT,y=nFeature_SCT, color=percent.mt))+
  geom_point()+
  ggtitle("drprnull_42D")+
  scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)

p1+p2+p3

p1 <- ggplot(data=data_list[2]$drprnull_42d@meta.data,aes(x=nCount_soupX_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)

p2 <- ggplot(data=data_list[2]$drprnull_42d@meta.data,aes(x=nFeature_soupX_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)

p3 <- ggplot(as.data.frame(data_list[2]$drprnull_42d@meta.data), aes(x=nCount_soupX_SCT,y=nFeature_soupX_SCT, color=percent.mt))+
  geom_point()+
  ggtitle("drprnull_42D")+
  scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)

p1+p2+p3

p1 <- ggplot(data=data_list[2]$drprnull_42d@meta.data,aes(x=nCount_decontX_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)

p2 <- ggplot(data=data_list[2]$drprnull_42d@meta.data,aes(x=nFeature_decontX_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)

p3 <- ggplot(as.data.frame(data_list[2]$drprnull_42d@meta.data), aes(x=nCount_decontX_SCT,y=nFeature_decontX_SCT, color=percent.mt))+
  geom_point()+
  ggtitle("drprnull_42D")+
  scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)

p1+p2+p3


ggplot(as.data.frame(data_list[1]$w1118_42d@meta.data), aes(x=nCount_soupX_SCT,y=nFeature_soupX_SCT, color=scds_hybrid_call))+
  geom_point()+
  ggtitle("w1118_42d")+
  #scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)

ggplot(as.data.frame(data_list[2]$drprnull_42d@meta.data), aes(x=nCount_soupX_SCT,y=nFeature_soupX_SCT, color=scds_hybrid_call))+
  geom_point()+
  ggtitle("drprnull_42d")+
  #scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)


## Clustering 
### Run clustering with different resolutions
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

### Run clustering on scTransformed raw counts
DefaultAssay(data_list[1]$w1118_42d)="SCT"
DefaultAssay(data_list[2]$drprnull_42d)="SCT"
data_list <- sapply(data_list, run_clustering)

### Run clustering on scTransformed soupX corrected counts
DefaultAssay(data_list[1]$w1118_42d)="soupX_SCT"
DefaultAssay(data_list[2]$drprnull_42d)="soupX_SCT"
data_list <- sapply(data_list, run_clustering)

### Run clustering on scTransformed decontX corrected counts
DefaultAssay(data_list[1]$w1118_42d)="decontX_SCT"
DefaultAssay(data_list[2]$drprnull_42d)="decontX_SCT"
data_list <- sapply(data_list, run_clustering)

### Plot UMAP and marker gene dotplots for each correction/clustering resolution combination
norm_methods <- c("SCT","soupX_SCT","decontX_SCT")
cluster_res <- c("snn_res.0.02","snn_res.0.2","snn_res.0.8","snn_res.1","snn_res.1.5","snn_res.2.5",
                 "snn_res.3.5","snn_res.4.5","snn_res.8","snn_res.12")

markers <- c("elav","lncRNA:noe","VAChT","VGlut","Gad1","Vmat","SerT","Tdc2","ple", # neurons
             "repo","lncRNA:CR34335","alrm","wrapper","Indy","moody",#glia 
             "ninaC",	"trp",	"trpl", #photoreceptors
             "Hml", #hemocytes
             "ppl",#fatbody
             "drpr")

for(nm in norm_methods){
  for(cr in cluster_res){
    nmcr <- paste0(nm,"_",cr)
    p1 <- DimPlot(data_list[1]$w1118_42d,group.by = nmcr) &
      theme_classic(base_size = 16)
    
    ggsave(p1,filename=paste0(outputdir,"figures/","w1118_42d_",nmcr,"_dimplot.pdf"),width=10,height=10)
    
    p2 <- DimPlot(data_list[2]$drprnull_42d,group.by = nmcr)&
      theme_classic(base_size = 16)
    
    ggsave(p2,filename=paste0(outputdir,"figures/","drprnull_42d_",nmcr,"_dimplot.pdf"),width=10,height=10)
    
    p3 <- DotPlot(data_list[1]$w1118_42d, features = markers, group.by = nmcr)&
      scale_color_viridis_c(option = "magma",direction = -1)&
      theme_classic(base_size = 16) + RotatedAxis()
    
    ggsave(p3,filename=paste0(outputdir,"figures/","w1118_42d_",nmcr,"_dotplot.pdf"), width=16,height=10)
    
    p4 <- DotPlot(data_list[2]$drprnull_42d, features = markers, group.by = nmcr)&
      scale_color_viridis_c(option = "magma",direction = -1)&
      theme_classic(base_size = 16)+ RotatedAxis()
    
    ggsave(p4,filename=paste0(outputdir,"figures/","drprnull_42d_",nmcr,"_dotplot.pdf"), width=16,height=10)
    #break
  }
}


## Save seurat objects together

saveRDS(data_list, file=paste0(outputdir,"data_list.RDS"))



## Save seurat objects as h5ad files separately


SaveH5Seurat(data_list[1]$w1118_42d, filename = paste0(outputdir,"w1118_42d.h5Seurat"))
Convert(paste0(outputdir,"w1118_42d.h5Seurat"), dest = paste0(outputdir,"w1118_42d.h5ad"))

SaveH5Seurat(data_list[2]$drprnull_42d, filename = paste0(outputdir,"drprnull_42d.h5Seurat"))
Convert(paste0(outputdir,"drprnull_42d.h5Seurat"), dest = paste0(outputdir,"drprnull_42d.h5ad"))


## Save SCTK-QC intermediate SCE objects as h5ad files 

writeH5AD(
  drprnull_42d,
  paste0(outputdir,"drprnull_42d_SCTKSCE.h5ad"),
  X_name = "counts",
  skip_assays = FALSE,
  compression = "none"
)


writeH5AD(
  w1118_42d,
  paste0(outputdir,"w1118_42d_SCTKSCE.h5ad"),
  X_name = "counts",
  skip_assays = FALSE,
  compression = "none"
)




sessionInfo()



