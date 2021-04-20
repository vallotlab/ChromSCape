Module_preprocessing_filtering_and_reductionUI <- function(id, label = "Module_preprocessing_filtering_and_reduction")
{
    ns <- NS(id)
}

Module_preprocessing_filtering_and_reduction <- function(
    input, output, session, raw_dataset_name, min_cov_cell, n_top_features,
    quant_removal, datamatrix, annot_raw, data_folder, annotation_id,
    exclude_regions, doBatchCorr, batch_sels, run_tsne, subsample_n)
    {
    withProgress(message = "Processing data set...", value = 0, {
        
        print("Data loading...")
        ### 1. Data loading ###
        batch_string <- if (doBatchCorr()) 
            "batchCorrected" else "uncorrected"
        incProgress(amount = 0.1, detail = paste("Loading raw data..."))
        
        print(system.time({scExp = create_scExp(datamatrix(), annot_raw(),
                                                remove_zero_cells = TRUE,
                                                remove_zero_features = TRUE)}))
        gc()
        
        ### 2. Filtering & Window selection ###
        
        incProgress(amount = 0.2, detail = paste("Filtering dataset..."))
        print("Filter...")
        print(system.time({ scExp = filter_scExp(
            scExp,
            min_cov_cell = min_cov_cell(),
            quant_removal = quant_removal(),
            min_count_per_feature = 10)}))
        gc()
        
        ### 2.bis Finding top features ###
        
        incProgress(amount = 0.2, detail = paste("Filtering dataset..."))
        print("Finding top covered features...")
        print(system.time({ scExp = find_top_features(
            scExp,
            n =  n_top_features(),
            keep_others = FALSE)
        }))
        gc()
        
        # Filtering based on exclude-regions from bed file, if provided
        if (!is.null(exclude_regions()))
        {
            scExp = exclude_features_scExp(scExp,
                                           exclude_regions(), by = "region")
            gc()
        }
        
       ### 2.ter Optional subsampling ###
        if(length(subsample_n())>0){
            print("Doing subsampling")
            scExp = subsample_scExp(scExp, n_cells = subsample_n())
            gc()
            
        }
        
        ### 3. Normalizing ###
        
        incProgress(amount = 0.1, detail = paste("Normalization..."))
        print("Normalize....")
        print(system.time({scExp = normalize_scExp(scExp, type = "CPM")}))
        gc()
        
        ### 4. Feature annotation ###
        
        incProgress(amount = 0.3, detail = paste("Feature annotation..."))
        print("Feature annotation...")
        print(system.time(
            {scExp = feature_annotation_scExp(scExp, ref = annotation_id())}))
        gc()
        
        # Original PCA
        print("Running Dimensionality Reduction...")
        
        incProgress(amount = 0.3, detail = paste(
            "Performing Dimensionality Reduction..."))
        if(run_tsne()) methods = c("PCA","TSNE","UMAP") else 
            methods = c("PCA","UMAP")
        print(system.time({scExp = reduce_dims_scExp(
            scExp, dimension_reductions = methods,
            batch_correction = doBatchCorr(), batch_list = batch_sels(), 
            verbose = FALSE)}))
        gc()
        
        ### 7. Add default colors ###
        print("Add colors ...")
        if(doBatchCorr()){ annotCol. = c("sample_id","total_counts", "batch_name")} 
        else{annotCol. = c("sample_id","total_counts")}
        print(system.time({scExp = colors_scExp(scExp, annotCol.)}))
        
        ### 8. Running hierarchical clustering ###
        print(system.time(
            {scExp = correlation_and_hierarchical_clust_scExp(scExp)}))
        
        ### 8. Save data ###
        print("saving :")
        print(file.path(
            data_folder(), "ChromSCape_analyses", raw_dataset_name(), 
            "Filtering_Normalize_Reduce", paste0(
                paste(raw_dataset_name(), min_cov_cell(), n_top_features(), 
                      quant_removal(), batch_string, sep = "_"), ".qs")))
        qs::qsave(
            scExp, file = file.path(
                data_folder(), "ChromSCape_analyses", raw_dataset_name(), 
                "Filtering_Normalize_Reduce", paste0(
                    paste(raw_dataset_name(), min_cov_cell(), n_top_features(), 
                        quant_removal(), batch_string, sep = "_"), ".qs"))
        )
        
        gc()
        print("Filtering & Reduction done !")
    })
}
