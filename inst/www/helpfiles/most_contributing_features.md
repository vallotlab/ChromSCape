---
output: 
  html_document: 
    self_contained: no
---
## Contribution to PCA

***

In this box, the contribution of features to the PCA are displayed with two 
graphs:  

- In the **barplot**, you will get the top loci / closest 
genes associated to the top loci contributing to the Component of your choice.
The top contributing features are taken by sorting the values of the features
in the componenent of intereset (e.g. PC1, PC2...). If you observe a huge gap 
between the top 2-3 features and the rest, this means that very few features 
are driving most of the dimensionality reduction, and that other less contributing 
features might be hidden. This also gives you interesting features to look at 
in the *Peak Calling & Coverage* or *Differential Analysis* tabs.  

- In the **pie chart**, the chromosome repartition of the top 100 most contributing features
to the componenent of interest (e.g. PC1, PC2...) are displayed. This is used
to see if there is any strong imbalance in the "chromosome contribution", e.g.
if a chromosome is contributing more than it should be. Usually, if a 
single chromosome has > 30% representation, this might indicate that CNVs or 
genetic events are happening on this chromosome, and you might want to look into
more details at these events, or exclude them from the analysis (see
*Filter & Normalize* tab)


