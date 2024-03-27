library(tidyverse)

basedir <- "/projectnb/mccall/sbandyadka/drpr42d_snrnaseq/reference/Fly head cell type markers-SCope/"
allfiles <- list.files(basedir)

find_markerfiles <-unlist(lapply(allfiles, grepl, pattern=".tsv"))
markerfiles <- allfiles[find_markerfiles]

celltype_markers <- c()
for(mf in markerfiles){
  #print(mf)
  removeext <- gsub(".tsv","",mf)
  removenum <- gsub("^[0-9]*.","",removeext)
  celltype <- str_replace_all(removenum, " ", "")
  #print(celltype)
  
  markers <- read.table(paste0(basedir,mf),header=TRUE)
  markers <- markers %>% filter(pval < 0.05) %>% arrange(desc(avg_logFC)) %>% head(15) %>% pull(gene)
  
  celltype_markers[[celltype]] <- markers
  
}

saveRDS(celltype_markers,file = paste0(basedir,"celltype_markers.RDS"))

