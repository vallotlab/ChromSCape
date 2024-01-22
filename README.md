# ChromSCape

## What is ChromSCape ?

**ChromSCape** - Chromatin landscape profiling for Single Cells - is a ready-to-launch user-friendly Shiny Application for the analysis of single-cell epigenomics datasets (scChIP-seq, scATAC-seq, scCUT&Tag, ...) from aligned data to differential analysis & gene set enrichment analysis. It is highly interactive, enables users to save their analysis and covers a wide range of analytical steps: QC, preprocessing, filtering, batch correction, dimensionality reduction, vizualisation, clustering, differential analysis and gene set analysis. 

[Prompsy, P., Kirchmeier, P., Marsolier, J. et al. Interactive analysis of single-cell epigenomic landscapes with ChromSCape.  
Nat Commun 11, 5702 (2020).](https://www.nature.com/articles/s41467-020-19542-x)

## Launching ChromSCape 

**ChromSCape** requires **R version 4.02**.
To install **ChromSCape**, open **R** or **Rstudio** and run the following commands : 

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("ChromSCape")
```

Once the installation was sucessful, launch **ChromSCape** using the following command :

```
library(ChromSCape)
launchApp()
```

It is recommended to use Chrome browser for optimal display of graphics & table.
If no browser opens, copy the url after 'Listening on ...' and paste in your browser.

## Demo 

Checkout the application look & feel at : [Demo](https://vallotlab.shinyapps.io/ChromSCape/). 
On this demo application, you can follow analysis of Jurkat & Ramos scChIP H3K4me3 cells.

## User guide

Take a look at the user guide before starting: 
[User guide](https://vallotlab.github.io/ChromSCape/articles/vignette.html)

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
The user should be on a Unix system (Mac, Linux) and have installed MACS2:

```
  macs2 2.1.2 (https://github.com/taoliu/MACS)
```
The application will automatically check if these tools are available and will give
you a warning if they are not installed/available.

## Note for Windows users

Windows user are not able to run the peak calling step, as both and macs2 are not yet available on windows.   
Therefore, if users want to run peak calling, they can use the docker version.

## Docker

A docker image with all dependencies is available at [DockerHub](https://hub.docker.com/repository/docker/pacomito/chromscape).
To run the docker image and launch ChromSCape, run :
```
sudo docker run --rm -p 4747:4747 -v ~/ChromSCape_analyses_docker:/root/output/ -t pacomito/chromscape:v0.0.9001
```
Explanation:

 * `sudo` run with admin rights, a password will be asked
 * `docker run -t pacomito/chromscape:v0.0.9001` download and run the image
 * `--rm` supress container when run ends
 * `-p 4747:4747` expose docker port 4747 to localhost:4747
 * `-v ~/ChromSCape_analyses_docker:/root/output/` output folder where ChromSCape_analyses folder will be created
 is linked to the container '/root/' folder. Change '~/ChromSCape_analyses_docker' to 
 your preferred output path
  
  
Optionally, if you want to input BAM, BED of Peak-Index-Barcode files, add another -v option from your local machine directory to the docker container:
```
-v ~/file_inputs:/root/input
```
After the downloading of the image and the loading of ChromSCape, navigate to : [http://127.0.0.1:4747](http://127.0.0.1:4747)

You can change the port number if it is already taken, e.g. port = 5858, by changing the -p option and adding '5858' as final argument :
```
sudo docker run --rm -p 4747:4747 -v ~/ChromSCape_analyses_docker:/root/output/ -t pacomito/chromscape:v0.0.9001 5858
```
  
  
MACS2 is installed on the Docker image so this is a way to run the 
peak calling on Windows.
  
# Authors
Please do not hesitate to post an issue or contact the authors :

Celine Vallot : celine.vallot@curie.fr

Pacome Prompsy : pacome.prompsy@curie.fr
