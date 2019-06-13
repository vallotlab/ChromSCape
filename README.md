# Analysis of single-cell Chromatine-ImmunoPrecipitation sequencing datasets
## What is scChIPseq_ShinyApp ?
This is a ready-to-launch user-friendly Shiny App for analysis of single-cell ChIP-seq datasets from count matrices to differential analysis & gene set enrichment analysis. ScChIPseq data can be produced using experimental protocol described in Grosselin et Al. (https://www.nature.com/articles/s41588-019-0424-9). The user should input one or many count matrices (in .txt or .tsv format). 

Note : Single-cell ATAC seq or scDNA-seq data under the same format should also work, as of same nature than scChIPseq data, but the application has not been tested for this type of data.

## Launch the App 

First download the repository in the location of your choice, either with `git clone https://github.com/vallotlab/scChIPseq.git scChIPseq` or by clicking on 'Clone or Download' -> 'Download ZIP' and unzip.

Go to the directory & modify the runApp.R script to the path to the directory & save file:
```
library(shiny)
runApp('path/to/scChIPseq_ShinyApp', launch.browser=TRUE)
```

Make sure to have all the libraries required (see ##Requirements) and start the App :

```
Rscript runApp.R
```
## Walkthrough of the App through screencast
![](www/img/scChIPseq_V1.gif)


## Sample datasets

The dataset from Grosselin et al. are PDX- triple negative breast cancer tumours resistant or not to chemotherapy (respectively HBCx_22 & HBCx_95).
Download the dataset of interest from GEO :https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117309. To run the app up to the differential analysis step, you need the count matrices. The peak calling and gene set enrichment parts require BAM files (available at https://figshare.com/s/fb04c2b17b234aa9d5eb). 


## Requirements

Before the first time you run the App, launch the installation script `Rscript ./installation_script.R` to install all the dependencies required.

```
  #Bioinfo
  library(scater)
  library(scran)
  library(IRanges)
  library(GenomicRanges)
  library(ConsensusClusterPlus)
  library(Rtsne)
  library(edgeR)
  
  #Data mining & utils
  library(tibble)
  library(dplyr)
  library(stringr)
  library(irlba)
  library(reshape2)
  library(DT)
  library(tidyr)
  library(splitstackshape)
  library(rlist)
  library(envDocument)
  library(rstudioapi)
  library(dplyr)
  
  #geco
  library(geco.utils)
  library(geco.visu)
  library(geco.unsupervised)
  library(geco.supervised)
  
  #Graphics
  library(RColorBrewer)
  library(colorRamps)
  library(colourpicker)
  library(kableExtra)
  library(knitr)
  library(viridis)
  library(ggplot2)
  library(gplots)
  library(png)
  library(grid)
  library(gridExtra)

  #Modules and functions
  source("Modules/geco.annotToCol2.R")
  source("Modules/geco.wilcox.R")
  
```

Install geco packages (under packages/) : 
```
install.packages("packages/geco.utils.tar.gz",repos = NULL,type = "source")
install.packages("packages/geco.visu.tar.gz",repos = NULL,type = "source")
install.packages("packages/geco.unsupervised.tar.gz",repos = NULL,type = "source")
install.packages("packages/geco.supervised.tar.gz",repos = NULL,type = "source")
```

Bash packages
```
  samtools 1.9 (Using htslib 1.9) (http://www.htslib.org/doc/samtools.html)
  bedtools v2.25.0 (https://github.com/arq5x/bedtools2/releases/tag/v2.25.0)
  macs2 2.1.2 (https://github.com/taoliu/MACS)
```

## Output

In the repo, the script should have created a directory **datasets** in which a new directory is created for each run with a different input name. Inside that directory are created a directory for each part of the analysis, containing RData and figures.
  
## Other

The Gene Set Enrichment Analysis is based on MSIG database (http://software.broadinstitute.org/gsea/msigdb).

# Authors
Please do not hesitate to post an issue or contact the authors :

Celine Vallot : celine.vallot@curie.fr

Pacome Prompsy : pacome.prompsy@curie.fr


# Session Info
```
R version 3.5.2 (2018-12-20)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=fr_FR.UTF-8      
 [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gridExtra_2.3               png_0.1-7                   gplots_3.0.1                viridis_0.5.1               viridisLite_0.3.0           knitr_1.22                 
 [7] kableExtra_1.1.0            colourpicker_1.0            colorRamps_2.3              RColorBrewer_1.1-2          shinyDirectoryInput_0.2.0   plotly_4.8.0               
[13] shinyjs_1.0                 shinydashboard_0.7.1        geco.supervised_1.0.0       geco.unsupervised_1.0.0     geco.visu_1.0.0             geco.utils_1.0.0           
[19] rlist_0.4.6.1               splitstackshape_1.4.6       tidyr_0.8.3                 DT_0.5                      ConsensusClusterPlus_1.46.0 Rtsne_0.15                 
[25] reshape2_1.4.3              edgeR_3.24.3                limma_3.38.3                irlba_2.3.3                 Matrix_1.2-15               stringr_1.4.0              
[31] dplyr_0.8.0.1               tibble_2.1.1                scran_1.10.2                scater_1.10.1               ggplot2_3.1.0               SingleCellExperiment_1.4.1 
[37] SummarizedExperiment_1.12.0 DelayedArray_0.8.0          BiocParallel_1.16.6         matrixStats_0.54.0          Biobase_2.42.0              GenomicRanges_1.34.0       
[43] GenomeInfoDb_1.18.2         IRanges_2.16.0              S4Vectors_0.20.1            BiocGenerics_0.28.0         shiny_1.2.0                

loaded via a namespace (and not attached):
 [1] ggbeeswarm_0.6.0         colorspace_1.4-1         dynamicTreeCut_1.63-1    rprojroot_1.3-2          XVector_0.22.0           BiocNeighbors_1.0.0      fs_1.2.7                
 [8] rstudioapi_0.10          remotes_2.0.2            xml2_1.2.0               pkgload_1.0.2            jsonlite_1.6             Cairo_1.5-10             cluster_2.0.7-1         
[15] HDF5Array_1.10.1         readr_1.3.1              BiocManager_1.30.4       compiler_3.5.2           httr_1.4.0               backports_1.1.3          assertthat_0.2.1        
[22] lazyeval_0.2.2           cli_1.1.0                later_0.8.0              htmltools_0.3.6          prettyunits_1.0.2        tools_3.5.2              igraph_1.2.4            
[29] gtable_0.3.0             glue_1.3.1               GenomeInfoDbData_1.2.0   Rcpp_1.0.1               gdata_2.18.0             crosstalk_1.0.0          DelayedMatrixStats_1.4.0
[36] xfun_0.6                 ps_1.3.0                 rvest_0.3.2              mime_0.6                 miniUI_0.1.1.1           gtools_3.8.1             devtools_2.0.1          
[43] statmod_1.4.30           zlibbioc_1.28.0          scales_1.0.0             hms_0.4.2                promises_1.0.1           rhdf5_2.26.2             yaml_2.2.0              
[50] memoise_1.1.0            stringi_1.4.3            desc_1.2.0               caTools_1.17.1.2         pkgbuild_1.0.3           rlang_0.3.3              pkgconfig_2.0.2         
[57] bitops_1.0-6             evaluate_0.13            lattice_0.20-38          purrr_0.3.2              Rhdf5lib_1.4.3           htmlwidgets_1.3          processx_3.3.0          
[64] tidyselect_0.2.5         plyr_1.8.4               magrittr_1.5             R6_2.4.0                 pillar_1.3.1             withr_2.1.2              RCurl_1.95-4.12         
[71] crayon_1.3.4             KernSmooth_2.23-15       rmarkdown_1.12           usethis_1.4.0            locfit_1.5-9.1           data.table_1.12.0        callr_3.2.0             
[78] webshot_0.5.1            digest_0.6.18            xtable_1.8-3             httpuv_1.5.0             munsell_0.5.0            beeswarm_0.2.3           vipor_0.4.5             
[85] sessioninfo_1.1.1
```
