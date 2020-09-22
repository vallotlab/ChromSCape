# This script is for testing reproducibility between
# the branch "package" & the branch "master" of ChromSCape

devtools::load_all("/media/pprompsy/LaCie/InstitutCurie/Documents/GitLab/ChromSCape/")
setwd("/media/pprompsy/LaCie/InstitutCurie/ChromSCape_inputs/scChIP_batch_effect_human_PDX/")

# Loading via import_scExp :
out = import_scExp(
    file_names = list.files("./",pattern = ".txt",full.names = T))

datamatrix = out$datamatrix
annot_raw = out$annot_raw

scExp = create_scExp(datamatrix,annot_raw)
gc()
table(scExp$sample_id)
scExp = filter_scExp(scExp)
table(scExp$sample_id)
gc()
dim(scExp)
scExp = normalize_scExp(scExp, "CPM")
summary(rowSums(normcounts(scExp)))
scExp = feature_annotation_scExp(scExp,ref="mm10")
tail(rowRanges(scExp))

#no batch cor
scExp = reduce_dims_scExp(scExp,dimension_reductions = c("PCA","UMAP"),n = 50, batch_correction = F)
scExp = colors_scExp(scExp)
plot_reduced_dim_scExp(scExp,reduced_dim = "UMAP")

# scExp_safe = scExp
#batch cor
scExp = reduce_dims_scExp(
    scExp,
    batch_correction = T,
    batch_list = list("batch_1" = c("HBCx_95_batch1","HBCx_95_CapaR_batch1"),
                      "batch_2" = c("HBCx_95_batch2"))
    )

reducedDim(scExp,"PCA")[1:3,1:3]
reducedDim(scExp,"UMAP")[1:10,1:2]

scExp = colors_scExp(scExp)
plot_reduced_dim_scExp(scExp,reduced_dim = "UMAP")
plot_reduced_dim_scExp(scExp,reduced_dim = "PCA")


pca = reducedDim(scExp,"PCA")
pca = as.matrix((pca * rowSums(abs(master$pca))) / rowSums(abs(as.matrix(reducedDim(scExp,"PCA")))))
colnames(pca) = paste0("PC",1:50)
expect_equal(cor(master$pca[,1], pca[,1]),1) # Is equal approximately
expect_equal(cor(master$pca[,17], pca[,17]),1) 
