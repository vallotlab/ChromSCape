## Input Data

***

ChromSCape allows user to input a variety of different format. Depending on the
output of the data-engineering/pre-processing pipeline used, the signal can
be already *summarized into features* :  

 * **Count matrix**: features x cells tab-separated matrix. Features as rownames and cells as column names 
 * **Sparse Matrix (10X format)** - 3 files :  
    + '*barcodes.tsv' : 1-column file of cell-barcode names (or .gz)
    + '*features.tsv': Tab-separated file of feature genomic location (or .gz)
    + '*matrix.mtx': 3 column space-separated file containing row index,
    column index and value of non-zeroes entries in the sparse matrix  
    
Or the signal can be stored directly in *raw format* (**single-cell BAM**,**single-cell BED/BED.gz**) containing genomic location of deduplicated reads for one cell.

Anyhow the format, ChromSCape needs signal to be summarized into features. 
If inputing *raw signal* (**scBAM** or **scBED**), the application lets user
summarize signal of each cells into various features:

 * **Genomic features** (extended region around TSS of genes, enhancers)
 * **Peaks** called on bulk or single-cell signal (BED file must be provided by the user)
 * **Genomic bins** (windows of constant length, e.g. 100kbp, 50kbp, 5kbp...)  

**Important note:** For Sparse Matrix (10X format), scBAM and scBED format,
the user must have one folder per sample, and place all the folders in the same 
directory as follows :

Example of folder structure : <br>
      <u>scChIPseq</u><br>
         * <b>Sample_1</b><br>
                - cell_1.bed<br>
                - cell_2.bed<br>
                - cell_3.bed<br>
                - cell_4.bed<br>
         * <b>Sample_2</b><br>
                - cell_1.bed<br>
                - cell_2.bed<br>
                - cell_3.bed<br>
                - cell_4.bed<br>

In this case, select the <u>scChIPseq</u> directory and ChromSCape will recognize Sample_1 and 
Sample_2 folders.
