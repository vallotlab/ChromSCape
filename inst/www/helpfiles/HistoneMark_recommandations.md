# Recommandations of features choice


## How to choose with which features to analyse your dataset ? 
 
 
Provided that you start with single-cell BED files, Fragment files or that you
have multiple features matrices (bin sizes, peaks, gene) available, **ChromSCape** 
lets you choos on which feature to analyse your dataset.  We do have recommandations
based on a benchmark done on the match analysis of multiple marks and either scRNA (Paired-Tag [1] ),
or surface proteins (CUT&Tag-pro [2] ) . However we still strongly recommand that
the user try out multiple features on the dataset.
  
As a general rule of thumb we recommend to first analyze the dataset with genomic bins,
which is the more 'unsupervised way', for the dimensionality reduction. Then for 
the supervised analysis, the user can add a more 'supervised' feature such as genes
or peaks to the analysis and perform differential analysis and gene set analysis
on these features.

For the various epigenetic features retrieved we recommend:

* H3K4me1, H3K4me3, ATAC  - bin sizes of 10,000bp - **20,000bp** - 50,000bp
* H3K27me3, H3K9me3, H3K27ac - bin sizes of 100,000bp-**200,000bp**

_____________________________

[1]: Zhu C, Zhang Y, Li YE, Lucero J, Behrens MM, Ren B. Joint profiling of histone modifications and transcriptome in single cells from mouse brain. Nat Methods. 2021;18(3):283-292. https://doi:10.1038/s41592-021-01060-3  

[2]: Zhang, B., Srivastava, A., Mimitou, E. et al. Characterizing cellular heterogeneity in chromatin state with scCUT&Tag-pro. Nat Biotechnol (2022). https://doi.org/10.1038/s41587-022-01250-0
