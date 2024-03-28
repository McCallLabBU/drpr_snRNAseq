
p1 <- ggplot(data=data_list[1]$w1118_42D@meta.data,aes(x=nCount_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)
p2 <- ggplot(data=data_list[1]$w1118_42D@meta.data,aes(x=nFeature_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)
p3 <- ggplot(as.data.frame(data_list[1]$w1118_42D@meta.data), aes(x=nCount_SCT,y=nFeature_SCT, color=percent.mt))+
  geom_point()+
  ggtitle("W1118_42D")+
  scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)

p1+p2+p3

p1 <- ggplot(data=data_list[2]$drprnull_42D@meta.data,aes(x=nCount_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)

p2 <- ggplot(data=data_list[2]$drprnull_42D@meta.data,aes(x=nFeature_SCT))+
  geom_histogram()+
  theme_classic(base_size = 16)

p3 <- ggplot(as.data.frame(data_list[2]$drprnull_42D@meta.data), aes(x=nCount_SCT,y=nFeature_SCT, color=percent.mt))+
  geom_point()+
  ggtitle("drprnull_42D")+
  scale_color_viridis_c(option = "magma")+
  theme_classic(base_size = 16)

p1+p2+p3

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

