---
output: 
  html_document: 
    self_contained: no
---
## Inter & Intra correlation violin plots

***

In this box, the heterogeneity of cells in various groups (e.g. sample or 
cluster) are displayed with two plots:

- In the intra-correlation plot, the distribution of cell-to-cell Pearson
correlation scores **within** each cluster or sample is displayed. Heterogeneous 
groups will have lower correlation than homogeneous groups. 

- In the inter-correlation plot, the distribution of cell-to-cell Pearson
correlation scores **between** each cluster or sample is displayed. This allows to 
understand which groups are closer to each other than other groups.

For both plots, you can add a layer of information by adding the single-cell
as jitter points on the plots. This is interesting when combining samples and
clusters to see the repartition of conditions into the clusters.


