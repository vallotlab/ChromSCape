#0. Importing packages ####
library(scater)
library(scran)
library(tibble)
library(dplyr)
library(stringr)
library(irlba)
library(edgeR)
library(reshape2)
library(Rtsne)
library(ConsensusClusterPlus)
library(DT)
library(tidyr)
library(IRanges)
library(GenomicRanges)
library(splitstackshape)
library(rlist)

#geco
library(geco.utils)
library(geco.visu)
library(geco.unsupervised)
library(geco.supervised)

#Shiny
library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)
library(shinyDirectoryInput)

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
source("Modules/Filtering_and_Reduction.R")
source("Modules/geco.annotToCol2.R")
source("Modules/geco.wilcox.R")
source("Modules/subsetBam_fromAffectation.R")

# global variables
#options(scipen = 999)
jscode <- "shinyjs.closeWindow = function() { window.close(); }"

gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

choose_perplexity <- function(dataset){
  perplexity=30
  if(nrow(dataset) <= 200 ){
    perplexity = 20
  }
  if(nrow(dataset) <= 250 ){
    perplexity = 25
  }
  if (nrow(dataset) <= 150 ){
    perplexity= 15
  }
  if (nrow(dataset) <= 100 ){
    perplexity= 10
  } 
  if (nrow(dataset) <= 50 ) {
    perplexity= 5
  }
 perplexity
}

