## Filtering parameters

***

Due to the current technical limitations of single-cell epigenomic technologies,
the signal is **sparse**, some cells are **lowly covered** and potential **doublets**
might be present in droplet-based technologies. To face these issues stringent 
filtering is necessary and the dataset can be filtered using 3 different parameters :  

* **Minimum coverage per cell** - represented
by the green line on the plot. Setting this filter will remove lowly covered cells. 
If you observe two peaks in the distribution, you can manually set this threshold
in between the two peaks as the first peak might be considered as empty droplets
with contaminant DNA.
* **Upper percentile of cells to remove (potential doublets)** - represented
by the red line on the plot. Cells with an outlier number of reads should be
discarded as there is a non-zero percentage of droplets containing 
2 cells or more (doublets).
* **Minimum percentage of cells to support a window**  - Regions not supported by this 
percentage of cells (with coverage greater than 1,000 reads) are filtered 
out. These regions might not be relevant to the analysis. The threshold at 1% is
recommended for genomic bins feature, but might be set to 0 if reads were counted
on predefined regions / features of intereset.
  
It is recommended to try out multiple combination of parameters to see what is 
the impact on clustering & differential analysis.
