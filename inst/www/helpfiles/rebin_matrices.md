## Re-bin matrices on new features

***
  
  
It is sometimes interesting to **appreciate the same data from multiple angles**
and to re-count the features into new larger features. This is possible if you start
with very small genomic bins (e.g. 500-1,000bp, recommended or maximum 5,000bp) and 
you want to rebin on large bin sizes (20,000bp - 1,000,000bp) or onto peaks.  


To re count the bins into peak, upload a BED file containing your peaks in the 
slot 'Re-count on features (BED)'.  


**It is important to set the minimum overlap to half of the original bin size 
(e.g. 500bp for original bins of 1,000bp) so that no original bins is counted twice.**

*** 


