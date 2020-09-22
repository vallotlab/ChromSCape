# ChromSCape: Analysis of single-cell epigenomic data in a Shiny App

## What is ChromSCape ?

**ChromSCape** - Single-Cell Chromatin Landscape profiling - is a ready-to-launch user-friendly **Shiny App** for analysis of single-cell epigenomic datasets (scChIP-seq, scATAC-seq...). It takes as input single-cell count matrices and let the user filter & cluster cells, run differential analysis & gene set enrichment analysis between epigenomic subpopulations, in an unsupervised manner.  
Various existing technologies allow to produce single-cell epigenomic datasets : scChIP-seq, scATAC-seq, scCUT&TAG, scChIL-seq, scChIC-seq ...

## Demo 

Checkout the application look & feel at : [Demo](https://vallotlab.shinyapps.io/ChromSCape/). 
On this demo application, you can follow analysis of Jurkat & Ramos scChIP H3K4me3 cells.

## Launch the App 

**ChromSCape** requires **R version 3.6.3** (does not work on R 4.0 yet).
To install **ChromSCape**, open **R** or **Rstudio** and copy the following commands : 

```
if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("vallotlab/ChromSCape")
```

Once the installation was sucessful, launch **ChromSCape** using the following command :

```
library(ChromSCape)
ChromSCape::launchApp()
```

It is recommended to use Chrome browser for optimal display of graphics & table.
If no browser opens, copy the url after 'Listening on ...' and paste in your browser.

## User guide

Take a look at the user guide before starting: 
[User guide](https://vallotlab.github.io/ChromSCape/ChromSCape_guide.html)

## Test datasets

ChromSCape takes as input one tab-separated count matrice (in .tsv or .txt) per sample. In order to upload multiple matrices, the matrices should be placed in the same folder of your computer. Before
you input your own matrices, it is recommended you try playing around and familiarize
with ChromSCape by downloading our example matrices and uploading them in ChromSCape :

Try out ChromSCape with various kind of dataset :
[Dropbox repository](https://www.dropbox.com/sh/vk7umx3ksgoez3x/AACEq9zn-rRbtwf_Al9uEUaQa?dl=0)

 

## Run Time

On a Intel® Core™ i5-6500 CPU @ 3.20GHz × 4 with 31,3 Gio RAM, the installation took less than one hour. The running time of of scChIP_H3K27me3 test dataset was 25 minutes without peak calling and 35 minutes with peak calling.

## Input

The matrix format should be tab-separated file, with Cells as column & Features 
as rows. The first line should be cell names, the first column should be feature 
names. Feature names can be either genomic coordinate in the format 'chr:start-end'
or 'chr_start_end' or gene symbols (e.g: A1BG, A1BG-AS1 for hg38 or Rab23, Bag2 
for mm10). 

## Output

The app automatically creates a directory **Chromscape_analysis** in which a new directory is created for each analysis with a different input name. Inside that directory are created a directory for each part of the analysis, containing RData and figures.
  
## Other

The Gene Set Enrichment Analysis is based on MSIG database (http://software.broadinstitute.org/gsea/msigdb).

## Advanced requirements for optional Peak Calling step

The peak calling step is important for Gene Set Enrichment Analysis particularly 
for features defined as genomic bins >= 20kbp or broad peaks. It will
aggregate signal of cells in each cluster ('in-silico cell sorting') and call peaks
separately for each cluster using MACS2 peak caller. Then the annotation of genes to
bins is refined and genes TSS not falling closer to 1000bp of any peaks are removed 
from annotation. This exclude any 'false' association of large genomic bins/regions to genes.  
This step requires **BAM files** of each sample (one BAM file must contains reads of all
 cells of a given sample) as input. 
The user should be on a Unix system (Mac, Linux) and have installed samtools & MACS2:

```
  samtools 1.9 (Using htslib 1.9) (http://www.htslib.org/doc/samtools.html)
  macs2 2.1.2 (https://github.com/taoliu/MACS)
```
The application will automatically check if these tools are available and will give
you a warning if they are not installed/available.

## Note for Windows users

Windows user are not able to run the peak calling step, as both samtools and macs2 are not yet available on windows.   

Also, if starting from a fresh installation of R3.6.3 on windows, you might encounter the following error   
```
WARNING: Rtools is required to build R packages, but is not currently installed.
Please download and install Rtools 3.5 from https://cran.r-project.org/bin/windows/Rtools/
```

Windows R version needs Rtools external software to install packages. Download Rtools 3.5 from https://cran.r-project.org/bin/windows/Rtools/history.html
and install it.  


# Authors
Please do not hesitate to post an issue or contact the authors :

Celine Vallot : celine.vallot@curie.fr

Pacome Prompsy : pacome.prompsy@curie.fr
