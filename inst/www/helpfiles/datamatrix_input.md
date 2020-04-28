## Count matrix

***

Browse here to select one or multiple single-cell count matrice(s) of epigenomic
features.  
These matrices are obtained after sequencing samples single-cell epigenomic 
experiment :  

* Targeting a specific protein mark: **scChIP-seq**, **CUT&Tag**, **scChIL-seq**,
**scChIC-seq**   
* Targeting open region of the chromatin: **snATAC-seq**, **scATAC-seq** 
  
The matrix format should be tab-separated file, with Cells as column & Features as rows. The first line should be cell names, the first column should be feature names. Feature names can be either genomic coordinate in the format 'chr:start-end' or 'chr_start_end' or gene symbols (e.g: A1BG, A1BG-AS1 for hg38 or Rab23, Bag2 for mm10). An example of such matrix is available
here [example matrix](https://github.com/vallotlab/ChromSCape/blob/package/inst/extdata/example_matrix.tsv.zip).
