
Single-nucleus RNA-seq data analysis of drprnull 42D vs. w1118 42D 



```
> sessionInfo()
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: AlmaLinux 8.9 (Midnight Oncilla)

Matrix products: default
BLAS:   /share/pkg.7/r/4.2.1/install/lib64/R/lib/libRblas.so
LAPACK: /share/pkg.7/r/4.2.1/install/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scran_1.24.0                scater_1.24.0               scuttle_1.6.2              
 [4] robustbase_0.99-2           viridis_0.6.2               viridisLite_0.4.0          
 [7] SeuratObject_5.0.1          Seurat_4.1.1                reticulate_1.25            
[10] SoupX_1.6.2                 zellkonverter_1.6.5         singleCellTK_2.6.0         
[13] DelayedArray_0.22.0         Matrix_1.6-5                forcats_0.5.1              
[16] stringr_1.4.0               dplyr_1.0.9                 purrr_0.3.4                
[19] readr_2.1.2                 tidyr_1.2.0                 tibble_3.1.7               
[22] ggplot2_3.3.6               tidyverse_1.3.1             SingleCellExperiment_1.18.0
[25] SummarizedExperiment_1.26.1 Biobase_2.56.0              GenomicRanges_1.48.0       
[28] GenomeInfoDb_1.32.2         IRanges_2.30.0              S4Vectors_0.34.0           
[31] BiocGenerics_0.42.0         MatrixGenerics_1.8.0        matrixStats_0.62.0         

loaded via a namespace (and not attached):
  [1] utf8_1.2.2                R.utils_2.11.0            tidyselect_1.1.2          htmlwidgets_1.5.4        
  [5] grid_4.2.1                BiocParallel_1.30.3       Rtsne_0.16                DropletUtils_1.16.0      
  [9] ScaledMatrix_1.4.0        munsell_0.5.0             codetools_0.2-18          ica_1.0-2                
 [13] statmod_1.4.36            future_1.26.1             miniUI_0.1.1.1            withr_2.5.0              
 [17] spatstat.random_3.2-3     colorspace_2.1-0          progressr_0.10.1          filelock_1.0.2           
 [21] rstudioapi_0.13           ROCR_1.0-11               tensor_1.5                listenv_0.8.0            
 [25] labeling_0.4.2            GenomeInfoDbData_1.2.8    GSVAdata_1.32.0           polyclip_1.10-0          
 [29] farver_2.1.0              rhdf5_2.40.0              basilisk_1.8.1            parallelly_1.32.0        
 [33] vctrs_0.4.1               generics_0.1.3            fishpond_2.2.0            R6_2.5.1                 
 [37] ggbeeswarm_0.6.0          rsvd_1.0.5                locfit_1.5-9.5            bitops_1.0-7             
 [41] rhdf5filters_1.8.0        spatstat.utils_3.0-4      assertthat_0.2.1          promises_1.2.0.1         
 [45] scales_1.2.0              beeswarm_0.4.0            rgeos_0.5-9               gtable_0.3.0             
 [49] beachmat_2.12.0           globals_0.15.0            goftest_1.2-3             spam_2.8-0               
 [53] rlang_1.0.2               splines_4.2.1             lazyeval_0.2.2            spatstat.geom_3.2-9      
 [57] broom_0.8.0               reshape2_1.4.4            abind_1.4-7               modelr_0.1.8             
 [61] backports_1.4.1           httpuv_1.6.5              tools_4.2.1               ellipsis_0.3.2           
 [65] spatstat.core_2.4-4       RColorBrewer_1.1-3        ggridges_0.5.3            Rcpp_1.0.8.3             
 [69] plyr_1.8.7                sparseMatrixStats_1.8.0   zlibbioc_1.42.0           RCurl_1.98-1.7           
 [73] basilisk.utils_1.8.0      rpart_4.1.16              deldir_1.0-6              pbapply_1.5-0            
 [77] cowplot_1.1.1             zoo_1.8-10                haven_2.5.0               ggrepel_0.9.1            
 [81] cluster_2.1.3             fs_1.5.2                  svMisc_1.2.3              magrittr_2.0.3           
 [85] RSpectra_0.16-1           data.table_1.14.2         scattermore_1.2           lmtest_0.9-40            
 [89] reprex_2.0.1              RANN_2.6.1                fitdistrplus_1.1-8        hms_1.1.1                
 [93] patchwork_1.1.0.9000      mime_0.12                 xtable_1.8-6              readxl_1.4.0             
 [97] gridExtra_2.3             compiler_4.2.1            KernSmooth_2.23-20        crayon_1.5.1             
[101] R.oo_1.25.0               htmltools_0.5.2           mgcv_1.8-40               later_1.3.0              
[105] tzdb_0.3.0                lubridate_1.8.0           DBI_1.1.3                 dbplyr_2.2.0             
[109] MASS_7.3-57               cli_3.3.0                 R.methodsS3_1.8.2         metapod_1.4.0            
[113] parallel_4.2.1            dotCall64_1.0-1           igraph_1.3.2              pkgconfig_2.0.3          
[117] dir.expiry_1.4.0          sp_1.5-0                  plotly_4.10.0             spatstat.sparse_3.0-3    
[121] xml2_1.3.3                vipor_0.4.5               dqrng_0.3.0               XVector_0.36.0           
[125] rvest_1.0.2               digest_0.6.29             sctransform_0.4.1         RcppAnnoy_0.0.19         
[129] spatstat.data_3.0-4       cellranger_1.1.0          leiden_0.4.2              uwot_0.1.11              
[133] edgeR_3.38.1              DelayedMatrixStats_1.18.0 shiny_1.7.1               gtools_3.9.3             
[137] lifecycle_1.0.1           nlme_3.1-158              jsonlite_1.8.0            Rhdf5lib_1.18.2          
[141] BiocNeighbors_1.14.0      limma_3.52.2              fansi_1.0.3               pillar_1.7.0             
[145] lattice_0.20-45           DEoptimR_1.0-11           fastmap_1.1.0             httr_1.4.3               
[149] survival_3.3-1            glue_1.6.2                png_0.1-7                 bluster_1.6.0            
[153] stringi_1.7.6             HDF5Array_1.24.1          BiocSingular_1.12.0       irlba_2.3.5              
[157] future.apply_1.9.0  

```
