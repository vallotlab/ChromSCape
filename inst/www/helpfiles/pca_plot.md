## PCA - Principal Component Analysis

***

This linear dimensionality reduction method is applied to the normalized count matrix.
In this plot, each cell is a dot. The cells are placed according to their projection
in the Component_1 & Component_2 (you are able to select any components). Cells close 
together have similar epigenomic traits.  
Coloring by 'total_counts' is a good way to ensure there is no direct correlation between
the library size and the first components.   
Among other things, applying PCA to the large count matrix allows to:

* Overcome the [Curse of Dimensionality](https://deepai.org/machine-learning-glossary-and-terms/curse-of-dimensionality) for downstream clustering of cells.
* Eliminate redudancy between multiple features & reduce noise
* Vizualise high-dimensional dataset in lower 2D space

