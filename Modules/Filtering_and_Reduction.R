moduleFiltering_and_ReductionUI <- function(id, label = "Filtering_and_Reduction module") {
  ns <- NS(id)
}

moduleFiltering_and_Reduction <- function(input, output, session, raw_dataset_name, min_cov_cell, percentMin, quant_removal,
                                          datamatrix, annot_raw, data_folder,
                                          annotation_id, exclude_regions) {
  withProgress(message='Processing data set...', value = 0, {
    
    batch_string <- "uncorrected"
    
    ###############################################################
    # 1. Data loading
    ###############################################################
    
    incProgress(amount=0.3, detail=paste("Loading raw data"))
 
    umi <- SingleCellExperiment(assays = list(counts = as.matrix(datamatrix())), colData = annot_raw()) #, colData = annot_raw

    umi <- umi[rowSums(counts(umi)>0)>0, ] # remove windows that do not have any read in any cells
    umi <- scater::calculateQCMetrics(umi)
    thresh <- quantile(colSums(counts(umi)), probs=seq(0, 1, 0.01))
    
    ###############################################################
    # 2. Filtering & Window selection
    ###############################################################
    incProgress(amount=0.2, detail=paste("Filtering data"))
    
    #Cell selection based on number of total counts, between min_cov_cell and upper 5%
    sel1000 =  (colSums(counts(umi))>1000 & colSums(counts(umi))< thresh[quant_removal()+1]) 
    sel <- ( colSums(counts(umi))>min_cov_cell() & colSums(counts(umi)) < thresh[quant_removal()+1] )
    

    
    SelMatCov1000 <- counts(umi)[,sel1000]
    bina_counts <- SelMatCov1000  
    bina_counts[bina_counts<2] <-0
    bina_counts[bina_counts>1] <-1
    fixedWin <- names(which((rowSums(bina_counts) > (percentMin()*(dim(bina_counts)[2])) ))) # window selection
    
    SelMatCov <- counts(umi)[,sel]
   
    
    annot <- colData(umi)
    annot <- as.data.frame(annot[sel,])

      
    SelMatCov <- SelMatCov[fixedWin,]
    
    
    
    # Filtering based on exclude-regions from bed file, if provided
    if(!is.null(exclude_regions())){
      regions <- data.frame(loc=rownames(SelMatCov))
      regions <- separate(regions, loc, into=c("chr", "start", "stop"), sep="_", convert=TRUE)
      reg_gr <- makeGRangesFromDataFrame(regions, ignore.strand=TRUE, seqnames.field=c("chr"), start.field=c("start"), end.field="stop")
      excl_gr <- makeGRangesFromDataFrame(exclude_regions(), ignore.strand=TRUE, seqnames.field=c("chr"), start.field=c("start"), end.field="stop")
      ovrlps <- as.data.frame(findOverlaps(reg_gr, excl_gr))[, 1]
      SelMatCov <- SelMatCov[-unique(ovrlps), ]
    }
    

    mat <- NULL
    
    mat <- mean(colSums(SelMatCov))*t(t(SelMatCov)/colSums(SelMatCov))
    
    norm_mat <- mat
    mat <- mat-apply(mat, 1, mean)
    save(norm_mat, file=file.path(data_folder(), "datasets", raw_dataset_name(), "reduced_data", paste0(paste(raw_dataset_name(), min_cov_cell(), (percentMin()*100), quant_removal(), batch_string, sep="_"), "_normMat.RData"))) # used for supervised analysis
    
   
    
    ###############################################################
    # 3. Feature annotation
    ###############################################################
    
    incProgress(amount=0.1, detail=paste("Feature annotation"))
    ID <- rownames(norm_mat)
    chr <- unlist(lapply(strsplit(ID, split= "_"), FUN=function(x) unlist(x)[1]))
    start <- unlist(lapply(strsplit(ID, split= "_"), FUN=function(x) unlist(x)[2]))
    end <- unlist(lapply(strsplit(ID, split= "_"), FUN=function(x) unlist(x)[3]))
    feature <- data.frame(ID=rownames(norm_mat), chr=chr, start=start, end=end)
    write.table(feature[, 2:4], file=("feature.bed"), sep="\t", row.names=F, col.names=F, quote=F)
    system("bedtools sort -i feature.bed > featuresort.bed")
    system(paste0("bedtools closest -a featuresort.bed -b ", file.path("annotation", annotation_id(), "Gencode_TSS_pc_lincRNA_antisense.bed"), " -wa -wb -d> out.bed"))
    annotFeat <- read.table("out.bed", header=F)
    unlink(c("feature.bed", "featuresort.bed", "out.bed"))
    annotFeat <- annotFeat[, c(1:3, 7:8)]
    colnames(annotFeat) <- c("chr", "start", "end", "Gene", "distance")
    annotFeat$ID <- paste(annotFeat$chr, annotFeat$start, annotFeat$end, sep="_")
    annotFeat <- annotFeat %>% group_by(ID) %>% summarise_all(funs(paste(unique(.), collapse = ', '))) %>% as.data.frame()
    save(annotFeat, file=file.path(data_folder(), "datasets", raw_dataset_name(), "reduced_data", paste0(paste(raw_dataset_name(), min_cov_cell(), (percentMin()*100), quant_removal(), batch_string, sep="_"), "_annotFeat.RData"))) # used for supervised analysis
    
    # ###############################################################
    # # 4. Batch correction
    # ###############################################################
    # 
    # pca <- NULL
    # batches <- list()
    # 
    # if(doBatchCorr()){
    #   incProgress(amount=0.2, detail=paste("Performing batch correction and dimensionality reduction - batch Corr"))
    #   annot$batch_name <- "unspecified"
    #   for(i in 1:num_batches()){
    #     for(s_id in batch_sels()[[i]]){
    #       annot[annot$sample_id==s_id, "batch_name"] <- batch_names()[i]
    #     }
    #   }
    #   adj_annot <- data.frame()  # annot must be reordered based on batches
    #   b_names <- unique(annot$batch_name)
    #   for(i in 1:length(b_names)){
    #     b_name <- b_names[i]
    #     batches[[i]] <- mat[, annot$batch_name==b_name]
    #     adj_annot <- rbind(adj_annot, annot[annot$batch_name==b_name, ])
    #   }
    #   mnn.out <- do.call(fastMNN, c(batches, list(k=20, d=50, approximate=TRUE, auto.order=TRUE, cos.norm=FALSE)))
    #   pca <- mnn.out$corrected
    #   colnames(pca) <- paste0("PC", 1:dim(pca)[2])
    #   annot <- adj_annot
    # }
    #Extract the rotation to do the cross product
    # out.2 <- fastMNN(pcs[[1]], pcs[[2]], pc.input=TRUE)
    # all.equal(head(out,-1), out.2) # should be TRUE (no rotation)
    # 
    # # Obtaining corrected expression values for genes 1 and 10.
    # cor.exp <- tcrossprod(out$rotation[c(1,10),], out$corrected)
    # dim(cor.exp)
    
    ###############################################################
    # 5. PCA
    ###############################################################
    
    # batch_ids <- unique(annot$batch_id)
    # 
    # ## !!!!!!! Put Batch_ids = 1 for testing purpose !!!!!!
    # batch_ids = 1
    # ###################################################"###
    # 
    
    # if(!doBatchCorr()){
    #   print("No Batch Correction")
    #   incProgress(amount=0.2, detail=paste("Performing dimensionality reduction- NO batch Corr"))
    #   # if(length(batch_ids) > 1){
    #   #   print("Detecting  different samples")
    #   #   print(length(unique(annot$batch_id)))
    #   #   for(i in batch_ids){
    #   #     batches[[i]] <- mat[, annot$batch_id==i]
    #   #   }
    #   # }else{  # simulate two batches by splitting data
    # 
    #   rand <- runif(dim(mat)[2]) > 0.5
    #   batches[[1]] <- mat[, rand]
    #   batches[[2]] <- mat[, !rand]
    #   # }
    #   pca_batches <- do.call(multiBatchPCA, batches)
    #   pca <- do.call(rbind, pca_batches)
    #   colnames(pca) <- paste0('PC', 1:50)
    #   pca <- pca[rownames(annot),]
      
    
    #Original PCA
    print("Running pca ...")
    pca <- stats::prcomp(t(mat),center=F,scale.=F)
    pca = pca$x[,1:50]

    # }
    
    ###############################################################
    # 6. tSNE
    ###############################################################
    
    incProgress(amount=0.2, detail=paste("Performing tSNE"))
    
    #Reduce the perplexity if the number of samples is too low to avoid perplexity error
    tsne <- Rtsne(pca[,1:50], dims=2, pca=FALSE, theta=0.0, perplexity=choose_perplexity(pca), verbose=FALSE, max_iter=1000)
    
    ###############################################################
    # 7. Save data
    ###############################################################
    
    save(pca, mat,annot, tsne, file=file.path(data_folder(), "datasets", raw_dataset_name(), "reduced_data", paste0(paste(raw_dataset_name(), min_cov_cell(), (percentMin()*100), quant_removal(), batch_string, sep="_"), ".RData")))
  })
}