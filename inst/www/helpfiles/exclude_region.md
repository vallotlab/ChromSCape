## Exclude regions / features from analysis

***

If your samples contain known amplified / deleted regions, you might want to
exclude the amplified / deleted regions by uploading a BED file of those regions
in order to remove bias from the analysis.  
If your count matrices contains genomic coordinates (peaks or regions), upload
a tab-seperated '.bed' file under the format:  
```
chr1	880001	901000
chr1	913001	944000
chr1	959001	961000
...
```

Else if you want to remove specific genes / enhancers, provide a text file of 
one feature name per line:
```
CNR2
PNRC2
SRSF10
...
```


