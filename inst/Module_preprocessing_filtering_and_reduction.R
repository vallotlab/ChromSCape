Module_preprocessing_filtering_and_reductionUI <- function(id, label = "Module_preprocessing_filtering_and_reduction")
{
    ns <- NS(id)
}

Module_preprocessing_filtering_and_reduction <- function(input, output, session, 
    raw_dataset_name, min_cov_cell, percentMin, quant_removal, datamatrix, annot_raw, 
    data_folder, annotation_id, exclude_regions, annotCol, doBatchCorr, batch_sels,
    run_tsne, subsample_n)
    {
    withProgress(message = "Processing data set...", value = 0, {
        
        print("Data loading...")
        ### 1. Data loading ###
        batch_string <- if (doBatchCorr()) 
            "batchCorrected" else "uncorrected"
        incProgress(amount = 0.1, detail = paste("Loading raw data..."))
        
        print(system.time({scExp = create_scExp(datamatrix(), annot_raw(), remove_zero_cells = T, remove_zero_features = T)}))
        
        ### 2. Filtering & Window selection ###
        
        incProgress(amount = 0.2, detail = paste("Filtering dataset..."))
        print("Filter...")
        print(system.time({ scExp = filter_scExp(scExp, min_cov_cell = min_cov_cell(), quant_removal = quant_removal(), 
            percentMin = percentMin())}))
        
        # Filtering based on exclude-regions from bed file, if provided
        if (!is.null(exclude_regions()))
        {
            scExp = exclude_features_scExp(scExp, exclude_regions(), by = "region")
        }
        
        ### 2.bis Optional subsampling ###
        if(length(subsample_n())>0){
            print("Doing subsampling")
            scExp = subsample_scExp(scExp, n_cells = subsample_n())
        }
        
        ### 3. Normalizing ###
        
        incProgress(amount = 0.1, detail = paste("Normalization..."))
        print("Normalize....")
        print(system.time({scExp = normalize_scExp(scExp, type = "CPM")}))
        
        ### 4. Feature annotation ###
        
        incProgress(amount = 0.3, detail = paste("Feature annotation..."))
        print("Feature annotation...")
        print(system.time({scExp = feature_annotation_scExp(scExp, ref = annotation_id())}))
        
        # Original PCA
        print("Running Dimensionality Reduction...")
        
        incProgress(amount = 0.3, detail = paste("Performing Dimensionality Reduction..."))
        print("Reduce dimensions annotation...")
        print(dim(scExp))
        print(class(SingleCellExperiment::normcounts(scExp)))
        if(run_tsne()) methods = c("PCA","TSNE","UMAP") else 
            methods = c("PCA","UMAP")
        print(system.time({scExp = reduce_dims_scExp(
            scExp, dimension_reductions = methods,
            batch_correction = doBatchCorr(), batch_list = batch_sels(), 
            verbose = F)}))
        
        ### 7. Add default colors ###
        print("Add colors ...")
        print(system.time({scExp = colors_scExp(scExp, annotCol())}))  # add colors 
        
        ### 8. Save data ###
        save(scExp, file = file.path(data_folder(), "ChromSCape_analyses", raw_dataset_name(), 
            "Filtering_Normalize_Reduce", paste0(paste(raw_dataset_name(), min_cov_cell(), percentMin(), 
                quant_removal(), batch_string, sep = "_"), ".RData")))
        
        gc()
        print("Filtering & Reduction done !")
    })
}
