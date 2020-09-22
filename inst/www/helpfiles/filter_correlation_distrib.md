## Filter lowly correlated cells

***

 To improve the stability of the clustering 
 and to remove from the analysis isolated cells that do not belong to any subgroup, cells
 displaying a Pearson’s pairwise correlation score below a threshold t with at least **p % of cells**
 are filtered (p is set at 1 % by default). The **correlation threshold t** is calculated as a
 percentile of Pearson’s pairwise correlation scores for a randomized dataset (percentile is
 recommended to be set as the 99th percentile). The correlation threshold is represented by
 the red line on the distribution plot.