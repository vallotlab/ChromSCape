# ChromSCape: Analysis of single-cell ChIP-seq data in a Shiny App

## Contents

- [What is ChromSCape ?](#What-is-ChromSCape-\?)
- [Launch the App](#Launch-the-App)
- [Docker version ](#Docker-version)
- [Test datasets](#Test-datasets)
- [Run Time](#Run-Time)
- [Walkthrough of the App with screencast](#Walkthrough-of-the-App-with-screencast)
- [Output](#Output)
- [Detailed walkthrough of the App](#Detailed-walkthrough-of-the-App)
- [Authors](#Authors)
- [Session Info](#Session-Info)

## What is ChromSCape ?

ChromSCape - Single-Cell Chromatin Landscape profiling - is a ready-to-launch user-friendly Shiny App for analysis of single-cell epigenomic datasets (scChIP-seq, scATAC-seq...) from count matrices to differential analysis & gene set enrichment analysis. ScChIPseq data can be produced using experimental protocol described in Grosselin et Al. (https://www.nature.com/articles/s41588-019-0424-9). The user should input one or many count matrices (in .txt or .tsv format). 

## Launch the App 

First download the repository in the location of your choice, either with `git clone https://github.com/vallotlab/ChromSCape.git ChromSCape` or by clicking on 'Clone or Download' -> 'Download ZIP' and unzip.

ChromSCape requires R version 3.5 or higher as well as bedtools v2.25.0 or higher (https://github.com/arq5x/bedtools2/releases/tag/v2.25.0). 

Make sure to have all the required libraries, to do so run in terminal :
```
Rscript ./installation_script.R
```

Go to the downloaded directory & modify the runApp.R script by indicating the path to ChromSCape directory & save file:
```
runApp('path/to/ChromSCape', launch.browser=TRUE)
```

and start the App from terminal:
```
Rscript runApp.R
```

## Docker version 

Docker is a software platform that allows you to build, test, and deploy applications quickly. In order for you to run the application without installing all R libraries on your computer, a docker container was created containing all libraries required to run ChromSCape. In order to use this docker environment, you need the admin rights of your computer. Please note that the total size of the docker environment is 3.0Gb.   

First, download & install Docker (https://hub.docker.com/?overlay=onboarding).  

Then, go to the branch "docker" of this repository, download or clone this repository on your local in the directory of your choice (e.g. /path/to/ChromSCape/ ).
Then create a data folder (e.g. /path/to/Data_ChromSCape/ ) that will be linked to the docker environment, and where you will be able to retrieve all of your data after closing the application.

### Change ownership of the data folders 
In order for the docker container to be able to write into the data folder you created, as well as the application folder for temporary files, the ownership of your directory must be changed for the container user 999 (shiny user) :

On Linux / Mac OS:
```
sudo chown -R 999:999 /path/to/Data_ChromSCape/
sudo chown -R 999:999 /path/to/ChromSCape/
sudo mkdir -p /path/to/Data_ChromSCape/bookmarks/shiny/
```



Then open the terminal and run the following command, replacing the path to application & data folder by your own :
```
sudo docker run --rm  -p 3838:3838  -v /path/to/ChromSCape/:/srv/shiny-server/ -v /path/to/Data_ChromSCape/:/var/lib/shiny-server/ -u shiny:shiny pacomito/chromscape:latest
```
Open a browser and go to http://localhost:3838/, you should see the application running.  

The user can't select a directory as of the data will be written directly in /path/to/Data_ChromSCape/.  
In the end of your analysis, if you want to retrieve your data run :

```
sudo cp /path/to/Data_ChromSCape/datasets/your_dataset /path/to/your_dataset_local
sudo chown username:username -R /path/to/your_dataset_local
```
Where your_dataset is the name you your dataset, /path/to/your_dataset_local is the path where you want to copy your dataset, username is your username.
  
Note that the peak calling is disabled for now in the docker application.

## Test datasets

The datasets correspond to mouse cells from 2 PDX models, luminal and triple negative breast cancer tumours resistant or not to cancer therapy (respectively HBCx_22 & HBCx_95, see Grosselin et al., 2019). 
Download count matrices from: https://figshare.com/projects/Single-Cell_ChIP-seq_of_Mouse_Stromal_Cells_in_PDX_tumour_models_of_resistance/66419 (theses count matrices have been processed using our latest data engineering pipeline, see https://github.com/vallotlab/scChIPseq_DataEngineering). The optional peak calling step requires  BAM files (also available on Figshare) to improve gene set enrichment analysis.

Alternatively, a ready-to-use pre-compiled analysis folder for HBCx22 & HBCx95 mouse H3K27me3 datasets is available at : https://figshare.com/articles/ChromSCape_scChIP_scATAC_compiled_datasets/11854371. A similar pre-compiled folder is available for the analysis of single-cell ATAC seq datasets from (Buenrostro et al., 2015, Corces et al., 2016, Schep et al., 2017). Download and uncompress the directory. Once in ChromSCape, select the directory containing the "dataset" folder and start exploring. 

## Run Time

On a Intel® Core™ i5-6500 CPU @ 3.20GHz × 4 with 31,3 Gio RAM, the installation took less than one hour. The running time of of scChIP_H3K27me3 test dataset was 25 minutes without peak calling and 35 minutes with peak calling.

## Walkthrough of the App with screencast

### 1 - Filtering, vizualiation & correlation
![](www/scChIPseq_App_V1.gif)

### 2 - Clustering, differential analysis & gene set enrichment
![](www/scChIPseq_App_V1_2.gif)

## Output

The app automatically creates a directory **datasets** in which a new directory is created for each analysis with a different input name. Inside that directory are created a directory for each part of the analysis, containing RData and figures.
  
## Other

The Gene Set Enrichment Analysis is based on MSIG database (http://software.broadinstitute.org/gsea/msigdb).


## Requirements

Before the first time you run the App, launch the **installation script** `Rscript ./installation_script.R` to install all the dependencies required.

```
  #Bioinfo
  library(scater)
  library(scran)
  library(IRanges)
  library(GenomicRanges)
  library(ConsensusClusterPlus)
  library(Rtsne)
  
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
  source("Modules/geco_functions.R")
  
```

Bash packages
```
  samtools 1.9 (Using htslib 1.9) (http://www.htslib.org/doc/samtools.html)
  bedtools v2.25.0 (https://github.com/arq5x/bedtools2/releases/tag/v2.25.0)
  macs2 2.1.2 (https://github.com/taoliu/MACS)
```

## Detailed walkthrough of the App
### 1. Upload your matrice(s)
### 2. Name & compile dataset
### 3. Fix filters to keep most covered regions & cells
### 4. Vizualize data in reduced dimension
### 5. Correlate cells and filter out cells with low correlation scores
### 6. Cluster cells
### 7. Peak call to refine signal (optional)
### 8. Find differentially bound regions in each cluster
### 9. Find enriched gene sets in differentially bound regions


# Authors
Please do not hesitate to post an issue or contact the authors :

Celine Vallot : celine.vallot@curie.fr

Pacome Prompsy : pacome.prompsy@curie.fr


# Session Info
```
R version 3.6.2 (2019-12-12)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.10

Matrix products: default
BLAS:   /usr/lib/libblas/libblas.so.3.6.1
LAPACK: /usr/lib/lapack/liblapack.so.3.6.1

locale:
 [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8     LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8    LC_PAPER=fr_FR.UTF-8      
 [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] DT_0.12                     gridExtra_2.3               png_0.1-7                   gplots_3.0.1.2              viridis_0.5.1               viridisLite_0.3.0           knitr_1.28                 
 [8] kableExtra_1.1.0            colourpicker_1.0            colorRamps_2.3              RColorBrewer_1.1-2          shinyDirectoryInput_0.2.0   plotly_4.9.2                shinyjs_1.1                
[15] shinydashboard_0.7.1        shiny_1.4.0                 rlist_0.4.6.1               splitstackshape_1.4.8       tidyr_1.0.2                 ConsensusClusterPlus_1.48.0 Rtsne_0.15                 
[22] reshape2_1.4.3              irlba_2.3.3                 Matrix_1.2-18               stringr_1.4.0               dplyr_0.8.4                 tibble_2.1.3                scran_1.12.1               
[29] scater_1.12.2               ggplot2_3.2.1               SingleCellExperiment_1.6.0  SummarizedExperiment_1.14.1 DelayedArray_0.10.0         BiocParallel_1.18.1         matrixStats_0.55.0         
[36] Biobase_2.44.0              GenomicRanges_1.36.1        GenomeInfoDb_1.20.0         IRanges_2.18.3              S4Vectors_0.22.1            BiocGenerics_0.30.0        

loaded via a namespace (and not attached):
 [1] ggbeeswarm_0.6.0         colorspace_1.4-1         dynamicTreeCut_1.63-1    XVector_0.24.0           BiocNeighbors_1.2.0      rstudioapi_0.11          RSpectra_0.16-0          xml2_1.2.2              
 [9] jsonlite_1.6.1           cluster_2.1.0            BiocManager_1.30.10      readr_1.3.1              compiler_3.6.2           httr_1.4.1               dqrng_0.2.1              assertthat_0.2.1        
[17] fastmap_1.0.1            lazyeval_0.2.2           limma_3.40.6             later_1.0.0              BiocSingular_1.0.0       htmltools_0.4.0          tools_3.6.2              rsvd_1.0.3              
[25] igraph_1.2.4.2           gtable_0.3.0             glue_1.3.1               GenomeInfoDbData_1.2.1   Rcpp_1.0.3               vctrs_0.2.3              gdata_2.18.0             DelayedMatrixStats_1.6.1
[33] xfun_0.12                rvest_0.3.5              mime_0.9                 miniUI_0.1.1.1           lifecycle_0.1.0          gtools_3.8.1             statmod_1.4.34           edgeR_3.26.8            
[41] zlibbioc_1.30.0          scales_1.1.0             hms_0.5.3                promises_1.1.0           stringi_1.4.6            caTools_1.18.0           rlang_0.4.4              pkgconfig_2.0.3         
[49] bitops_1.0-6             evaluate_0.14            lattice_0.20-38          purrr_0.3.3              htmlwidgets_1.5.1        tidyselect_1.0.0         plyr_1.8.5               magrittr_1.5            
[57] R6_2.4.1                 pillar_1.4.3             withr_2.1.2              RCurl_1.98-1.1           crayon_1.3.4             KernSmooth_2.23-16       rmarkdown_2.1            locfit_1.5-9.1          
[65] data.table_1.12.8        digest_0.6.25            webshot_0.5.2            xtable_1.8-4             httpuv_1.5.2             munsell_0.5.0            beeswarm_0.2.3           vipor_0.4.5
```
