
library(Seurat)
library(tidyverse)


input_dir <- file.path("/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/analysis/")
drprnull_42D_path <- file.path(paste0(input_dir,"drprnull_42D/","outs/filtered_feature_bc_matrix/"))
w1118_42D_path <- file.path(paste0(input_dir,"w1118_42D/","outs/filtered_feature_bc_matrix/"))

output_dir <- file.path(paste0(input_dir,"post_alignment_QC"))


drprnull_42D_raw <- Read10X(data.dir = drprnull_42D_path)
w1118_42D_raw <- Read10X(data.dir = w1118_42D_path)

drprnull_42D <- CreateSeuratObject(counts = drprnull_42D_raw, 
                                     project = "drprnull_42D", 
                                     min.cells = 3, 
                                     min.features = 200)
w1118_42D <- CreateSeuratObject(counts = w1118_42D_raw, 
                                   project = "w1118_42D", 
                                   min.cells = 3, 
                                   min.features = 200)


drprnull_42D <- PercentageFeatureSet(drprnull_42D, pattern = "-m", col.name = "percent.mt")
w1118_42D <- PercentageFeatureSet(w1118_42D, pattern = "-m", col.name = "percent.mt")

VlnPlot(drprnull_42D, group.by = "orig.ident", features = c("nFeature_RNA", 
                                                      "nCount_RNA", 
                                                      "percent.mt"), pt.size = 0.1, ncol = 3) + NoLegend()

VlnPlot(w1118_42D, group.by = "orig.ident", features = c("nFeature_RNA", 
                                                      "nCount_RNA", 
                                                      "percent.mt"), pt.size = 0.1, ncol = 3) + NoLegend()


drprnull_42D <- subset(drprnull_42D, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA < 25000 & percent.mt < 0.5)
w1118_42D <- subset(w1118_42D, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & nCount_RNA < 10000 & percent.mt < 0.5)

data_list <- c(drprnull_42D,w1118_42D)
data_list <- lapply(X = data_list, FUN = function(x) {
  x <- SCTransform(x,method = "glmGamPoi")
  
})
features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 2000)
data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = features)
data_list <- lapply(X = data_list, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(object.list = data_list, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
combined <- RunPCA(combined, verbose = FALSE)
ElbowPlot(combined)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)

combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

DimPlot(combined, reduction = "umap", group.by = "orig.ident")
DimPlot(combined, reduction = "umap", split.by = "orig.ident", label=TRUE)


DefaultAssay(combined) <- "integrated"
DefaultAssay(combined) <- "RNA"
markers <- FindAllMarkers(combined,
                                     assay = "RNA",
                                     min.pct = 0.25,
                                     only.pos = TRUE,
                                     logfc.threshold = 0.5) %>% 
  filter(p_val_adj <= 0.05) %>% arrange(cluster, desc(avg_log2FC))

top_markers <- markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top_noncg_markers <- markers %>% filter(avg_log2FC > 0) %>% group_by(cluster) %>% 
  filter(!grepl("CG", gene)) %>% top_n(n = 2, wt = avg_log2FC) %>% unique()

DotPlot(combined, features = unique(top_markers$gene), cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") +
  RotatedAxis()

celltype_markers <- c("VAChT",#cholinergic
                      "VGlut", #glutamatergic
                      "Gad1", #GABAergic
                      "Vmat",#monoaminergic
                      "SerT", #serotonergic
                      "Tdc2",#octopaminergic/tyraminergic
                      "ple"#dopaminergic
                      )

DotPlot(combined, features = celltype_markers, cols = c("blue", "red"), split.by = "orig.ident") +
  RotatedAxis()

FeaturePlot(combined, features=celltype_markers)








