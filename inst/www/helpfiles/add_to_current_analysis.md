## Additional feature layer on the same cells

***
  
  
It is sometimes interesting to **appreciate the same data from multiple angles**. 
  
One might want to use *50kbp genomic bins* for a comprehensive 
dimensionality reduction and clustering, then run differential analysis 
using *genebodies*.
  
Or using different size of bins and peaks while being able to
visualize the different UMAPs on the same cells.
  
For genetically 
alterated biological samples one might also count reads at the cytoband level to
see if 'de-zooming' still allows to separate cells 
(see `?ChromSCape::calculate_CNA()`).
  
*** 

When ticking this checkbox, data will be added to the current analysis. The
data you wish to upload **must contain the same cells** (e.g. cell barcodes) than 
the current analysis. 
  
  
The additional data layers you upload will be kept in alternative slots
(see `?SingleCellExperiment::altExp()`). You can switch at any time between 
the different features you counted your data by clicking on '**Features**' in the
upper bar of the application.

