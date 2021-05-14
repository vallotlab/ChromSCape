# ChromSCape 1.1.3

## Major Changes

    * Support "multi-feature" analysis, e.g. parallel analysis of multiple 
    features (bins, peaks or gene) on the same object.
    
    * New "Coverage" tab & functions generate_coverage_tracks() and
    plot_coverage_BigWig() to generate cluster coverage tracks and interactively 
    visualise loci/genes of interest in the application.
    
    * New inter- and intra-correlation violin plots to vizualise cell correlation
    distribution between and within clusters.
    
    * New normalization method : TF-IDF combined with systematic removal of PC1
    strongly correlated with library size.
    
    * Simple 'Copy Number Alteration' approximation & visualization using 
    'calculate_CNA' function for genetically re-arranged samples, provided one 
    or more control samples.
    
    * New generate_analysis() & generate_report() functions to run a full-on 
    ChromSCape analysis and/or generate an HTML interactive report of an existing analysis.
    
    * Supports 'custom' differential analysis to find differential loci between
    a subset of samples and/or clusters.
    
    * New pathway overlay on UMAP to visualize cumulative pathways signal 
    directly on cells.
    
    * Now supports 'Fragment Files' input (e.g. from 10X cell ranger scATAC
    pipeline), using a wrapper around 'Signac' package FeatureMatrix() function.
    
    * New 'Contribution to PCA' plots showing most contributing features and 
    chromosome to PCA.
    
    * Restructuration of the ChromSCape directory & faster reading/saving of 
    S4 objects using package 'qs'.
    
    
## Minor Changes
 
    * RAM optimisation & faster pearson cell-to-cell correlations with 'coop'
    package, and use of 'Rcpp' for as_dist() RAM-efficient distance calculation.
    
    * Faster correlation filtering using multi-parallel processing.
    
    * plot_reduced_dim now supports gene input to color cells by gene signal.
    
    * All plots can now be saved in High Quality PDF files.
    
    * Changed 'geneTSS' to 'genebody' with promoter extension to better reflect
    the fact that mark spread in genebodies.
    
    * Possibility to rename samples in the application.
    
    * Downsampling of UMAPs & Heatmaps for fluider navigation.
    
    * Changed 'total cell percent based' feature selection to manual selection of
    top-covered features, as the previous was srongly dependent on the experiment size.
    
    * Faster sparse SVD calculation.
    
    * Faster differential analysis using pairWise Wilcoxon rank test from 'scran'
    package.




