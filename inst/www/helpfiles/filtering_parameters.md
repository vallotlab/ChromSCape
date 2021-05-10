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

* **Top covered features to keep**  - Number of features to keep for the PCA & 
downstream analysis. Usually not all of the features contain relevant information
and retrieving the most covered features is a simple way to filter for
the most relevant features. The violet line on the 'Feature coverage plot' indicates 
the minimum coverage in the selected features. 
  
It is recommended to try out multiple combination of parameters to see what is 
the impact on clustering & differential analysis.
