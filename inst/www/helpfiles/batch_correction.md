## Batch correction

***

Sometimes datasets produced using different methods or in different laboratories
contain batch effect biasing interpretation of the differences between cells of 
various samples. In order to remove batch effect you can specify a number of batches
and the batch-sample correspondance to run a mutual nearest neigbour correction 
implemented by Haghverdi et al.^1^, in the package 'scran'. The advantage of this
approach is to not rely on predefined or equally sized cell populations accross
batches but it assumes that a subset of population of each batch to be similar.
We advice to always run the analysis without batch correction first and to rerun
with batch correction only if a clear batch effect is present.

For the differential analysis if the batch correction is 'ON', differential 
regions between two clusters are calculated using scran function 'pairwiseWilcox', 
using the batch ID as a blocking level.
  
1. Batch effects in single-cell RNA-sequencing data are corrected by matching 
mutual nearest neighbors. Haghverdi L, Lun ATL, Morgan MD, Marioni JC. Nat 
Biotechnol. 2018. 

