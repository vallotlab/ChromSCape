## Input Data

***

ChromSCape allows user to input a variety of different format. Depending on the
output of the data-engineering/pre-processing pipeline used, the signal can
be already *summarized into features* :  

 * **Count matrix**: features x cells tab-separated matrix. Features as rownames and cells as column names 
 * **Peak-Index-Barcode** -3 files :  
    + '*_barcode.txt' : 1-column file of cell-barcode names
    + '*_peaks.bed': BED file of feature genomic location
    + '_index.txt': 3 column space-separated file containing row index, column index and value of non-zeroes entries in the sparse matrix  
    
Or the signal can be stored directly in *raw format* (**single-cell BAM**,**single-cell BED/BED.gz**) containing genomic location of deduplicated reads for one cell.

Anyhow the format, ChromSCape needs signal to be summarized into features. 
If inputing *raw signal* (**scBAM** or **scBED**), the application lets user
summarize signal of each cells into various features:

 * **Genomic features** (extended region around TSS of genes, enhancers)
 * **Peaks** called on bulk or single-cell signal (BED file must be provided by the user)
 * **Genomic bins** (windows of constant length, e.g. 100kbp, 50kbp, 5kbp...)  

**Important note 1:** For Peak-Index-Barcode, scBAM and scBED format, the user must precise the number of samples contained in the data.  

**Important note 2:** For Peak-Index-Barcode, scBAM and scBED format, the user must place all the files in one folder and select the folder containing the files, not directly the files.
