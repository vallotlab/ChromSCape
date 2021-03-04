---
output: 
  html_document: 
    self_contained: no
---
## Differential Analysis

***

In this tab, investigate on the differences in levels of the mark of interest or
accessibility for scATAC-seq between clusters or samples. 
Groups can be compared in three kind of ways: 

- **One vs All**: This option will run comparison for each cluster, taking as 
reference all the cells in other clusters. This is interesting to understand
which features are specifically enriched / depleted in one cluster but not any
others.

- **Pairwise**: This option will compare each cluster two by two and consolidate
the multiple differential tests for each cluster using scran::combineMarkers 
function with pval.type = "any" : 'This approach does not explicitly favour [loci]
that are uniquely [differential] in a cluster. Rather, it focuses on combinations of 
[loci] that - together - drive separation of a cluster from the others.

- **Custom** Using this option, choose your own group and reference to compare,
either samples or clusters and do a single differential testing. This is useful
when you are only interested in the differences between two groups but not the 
others.

For each option, user can pick either Wilocoxon non-parametric testing or 
edgeR GLM parametric testing, relying on linear models.  

If running multiple comparisons using various options, all the comparisons will
be saved and accessible at any time.


