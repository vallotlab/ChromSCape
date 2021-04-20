## Run a complete unsupervised analysis and create an HTML report
# ChromSCape_directory = "/media/pacome/LaCie/InstitutCurie/Documents/Data/Tests/Test_ChromSCape/Test_generate_analysis/ChromSCape_analyses/Mouse_ctrl_orga_50kbp"
# datamatrix = qs::qread(file.path(ChromSCape_directory, "datamatrix.qs"))
# scExp = qs::qread(file.path(ChromSCape_directory,"Filtering_Normalize_Reduce" ,"Mouse_ctrl_orga_50kbp_1300_0.5_95_uncorrected.qs"))
# scExp = qs::qread( "/media/pacome/LaCie/InstitutCurie/Documents/Data/Tests/Test_ChromSCape/Test_generate_analysis/ChromSCape_analyses/Mouse_ctrl_orga_50kbp_full/Filtering_Normalize_Reduce/Mouse_ctrl_orga_50kbp_full_1300_0.5_95_uncorrected.qs")
# scExp_cf = qs::qread(file.path(ChromSCape_directory,"Diff_Analysis_Gene_Sets" ,"Mouse_ctrl_orga_50kbp_1300_0.5_95_uncorrected_2_0.1_1_one_vs_rest.qs"))$scExp_cf

analysis_name = "Test"
# coverages = sapply(list.files(file.path(ChromSCape_directory,"coverage","Mouse_ctrl_orga_50kbp_1300_0.5_95_uncorrected_k2"),".bw", full.names = T), rtracklayer::import)
input_data_type = "scBED"
feature_count_on = "bins"
feature_count_parameter = 50000
min_reads_per_cell = 1300
min_percent_to_keep_feature = 0.5
max_quantile_read_per_cell = 99
maxK = 8
n_clust = 6
output_directory = "/media/pacome/LaCie/InstitutCurie/Documents/Data/Tests/Test_ChromSCape/Test_generate_analysis/"
input_data_folder = "/media/pacome/LaCie/InstitutCurie/ChromSCape_inputs/test_small_scBED/"
control_samples_CNA = c("C_a3_H3K27me3")
ref_genome = "mm10"
doBatchCorr = F
batch_sels = list("control" = c( "C_a3_H3K27me3","C_a4_5_H3K27me3"),
                  "juxta" = c("OJ_m6374_p2_H3K27me3","OT_m6374_p2_H3K27me3"))
subsample_n = NULL
exclude_regions = NULL
corr_threshold = 99
percent_correlation = 1
qval.th = 0.1
logFC.th = 1
enrichment_qval = 0.1
genes_to_plot = c("Krt8","Krt5","Tgfb1" ,"Foxq1", "Cdkn2b",
                  "Cdkn2a", "chr7:15000000-20000000")
# generate_analysis(
#     analysis_name = analysis_name,
#     input_data_folder = input_data_folder,
#     output_directory = output_directory,
#     input_data_type = input_data_type,
#     feature_count_on = feature_count_on,
#     feature_count_parameter = feature_count_parameter,
#     ref_genome = ref_genome,
#     control_samples_CNA = control_samples_CNA,
#     min_reads_per_cell = min_reads_per_cell,
#     max_quantile_read_per_cell = max_quantile_read_per_cell,
#     min_percent_to_keep_feature = min_percent_to_keep_feature,
#     maxK = maxK,
#     n_clust = n_clust,
#     doBatchCorr = doBatchCorr,
#     batch_sels = batch_sels)

generate_analysis <- function(input_data_folder,
                  analysis_name = "Analysis_1",
                  output_directory = "./",
                  input_data_type = c("DenseMatrix", "SparseMatrix", "scBED",
                                      "scBAM")[1],
                  feature_count_on = c("bins","geneTSS","peaks")[1],
                  feature_count_parameter = 50000,
                  ref_genome = c("hg38","mm10")[1],
                  run = c("filter", "CNA","cluster", "peak_call", "coverage", 
                          "DA", "GSA", "report"),
                  min_reads_per_cell = 1600,
                  max_quantile_read_per_cell = 95,
                  min_percent_to_keep_feature = 0.5,
                  subsample_n = NULL,
                  exclude_regions = NULL,
                  n_clust = NULL,
                  corr_threshold = 99,
                  percent_correlation = 1,
                  maxK = 10,
                  qval.th = 0.1,
                  logFC.th = 1,
                  enrichment_qval = 0.1,
                  doBatchCorr  = FALSE,
                  batch_sels  = NULL,
                  control_samples_CNA = NULL,
                  genes_to_plot = c("Krt8","Krt5","Tgfb1", "Foxq1", "Cdkn2b",
                                    "Cdkn2a", "chr7:15000000-20000000")
){
    ##### Run ChromSCape ####
    message("ChromSCape::generate_analysis - Running full analysis - ",
            analysis_name )
    
    #### Check parameters ####
    stopifnot(dir.exists(input_data_folder), is.character(analysis_name),
              is.character(output_directory),
              input_data_type %in% c("DenseMatrix", "SparseMatrix", "scBED", "scBAM"),
              ref_genome %in% c("hg38","mm10"), is.numeric(min_reads_per_cell),
              is.numeric(max_quantile_read_per_cell),
              is.numeric(min_percent_to_keep_feature),
              is.numeric(qval.th),
              is.numeric(logFC.th),
              is.numeric(enrichment_qval)
              )
    
    #### Create Analysis directory ####
    if(!dir.exists(output_directory)) dir.create(output_directory)
    ChromSCape_analyses = file.path(output_directory, "ChromSCape_analyses")
    if(!dir.exists(ChromSCape_analyses)) dir.create(ChromSCape_analyses)
    
    ChromSCape_directory = file.path(ChromSCape_analyses, analysis_name)
    if(!dir.exists(ChromSCape_directory)) dir.create(ChromSCape_directory)
    filt_dir <- file.path(ChromSCape_directory,"Filtering_Normalize_Reduce") 
    dir.create(filt_dir, showWarnings = FALSE)
    cor_dir <- file.path(ChromSCape_directory, "correlation_clustering")
    dir.create(cor_dir, showWarnings = FALSE)
    cor_plot_dir <- file.path(cor_dir,"Plots")
    dir.create(cor_plot_dir, showWarnings = FALSE)
    coverage_dir <- file.path(ChromSCape_directory,"coverage")
    dir.create(coverage_dir, showWarnings = FALSE)
    peak_dir <- file.path(ChromSCape_directory,"peaks")
    dir.create(peak_dir, showWarnings = FALSE)
    da_gsa_dir <- file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets")
    dir.create(da_gsa_dir, showWarnings = FALSE)
    
    write(ref_genome, file = file.path(ChromSCape_directory, "annotation.txt"))
    batch_string <- if (doBatchCorr) "batchCorrected" else "uncorrected"
    prefix = paste0(analysis_name, "_", min_reads_per_cell, "_",
           min_percent_to_keep_feature, "_",
           max_quantile_read_per_cell, "_", batch_string)
    
    time_analysis = system.time({
    #### Select & Import ####
    message("ChromSCape::generate_analysis - Importing datasets ...")
    out <- rawData_to_datamatrix_annot(input_data_folder, input_data_type,
                                       feature_count_on, feature_count_parameter,
                                       ref_genome
    )
    
    datamatrix = out$datamatrix
    annot_raw = out$annot_raw
    qs::qsave(datamatrix, file = file.path(ChromSCape_directory, "datamatrix.qs"))
    qs::qsave(annot_raw, file = file.path(ChromSCape_directory, "annot_raw.qs"))
    
    # Check that number of cells and features is enough or finish here
    n_cell = ncol(datamatrix)
    n_feature = nrow(datamatrix)
    
    if( (n_cell < 100) | (n_feature < 300) ){
        warning("Analyis ", analysis_name, " has stopped after raw data import,",
                " there were less than 100 cells or 300 features...")
        return()
    }
    #### Filter & Normalize ####
    if(!"filter" %in% run) {
        message("Finished generating analyis ", analysis_name, ". The analysis ",
            "is available at ", ChromSCape_directory, "...")
        return()
    } 
    message("ChromSCape::generate_analysis - Filtering, normalizing & ",
            "reducing dimensions ...")
    scExp = preprocessing_filtering_and_reduction(
        datamatrix = datamatrix,
        annot_raw = annot_raw,
        min_reads_per_cell = min_reads_per_cell,
        max_quantile_read_per_cell = max_quantile_read_per_cell,
        min_percent_to_keep_feature = min_percent_to_keep_feature,
        subsample_n = subsample_n,
        ref_genome = ref_genome,
        exclude_regions = exclude_regions,
        doBatchCorr  = doBatchCorr,
        batch_sels  = batch_sels
    )
    
    # Check that number of cells and features is enough or finish here
    n_cell = ncol(scExp)
    n_feature = nrow(scExp)
    
    if( (n_cell < 100) | (n_feature < 300) ){
        warning("Analyis ", analysis_name, " has stopped after QC filtering,",
                " there were less than 100 cells or 300 features...")
        return()
    }
    
    #### Calculate CNA
    if("CNA" %in% run) {
        message("ChromSCape::generate_analysis - running Copy Number Variant",
        " detection of cytobands on unfiltered features...")
        if(!is.null(control_samples_CNA)){
            scExp_full = create_scExp(datamatrix, annot_raw)
            scExp_full = scExp_full[,match(scExp$cell_id, scExp_full$cell_id)]
            scExp_full = calculate_CNA(scExp_full, control_samples = control_samples_CNA,
                                       ref_genome = ref_genome)
            
            SingleCellExperiment::reducedDim(scExp,"cytoBand") =
                SingleCellExperiment::reducedDim(scExp_full,"cytoBand")
            SingleCellExperiment::reducedDim(scExp,"logRatio_cytoBand") =
                SingleCellExperiment::reducedDim(scExp_full,"logRatio_cytoBand")
            SingleCellExperiment::reducedDim(scExp,"gainOrLoss_cytoBand") =
                SingleCellExperiment::reducedDim(scExp_full,"gainOrLoss_cytoBand")
            gc()
            rm(scExp_full)
            gc()
        }
    }
    qs::qsave(scExp, file = file.path(filt_dir, paste0(prefix, ".qs")))
    
    #### Correlation Filtering & Clustering ####
    if("cluster" %in% run) {
        message("ChromSCape::generate_analysis - Correlation Consensus Clustering ...")
        scExp_cf = filter_correlated_cell_scExp(scExp, random_iter = 50, corr_threshold = corr_threshold,
                                                percent_correlation = percent_correlation)
        
        # Check that number of cells and features is enough or finish here
        n_cell = ncol(scExp_cf)
        if( (n_cell < 100)){
            warning("Analyis ", analysis_name, " has too few cells after cell ",
                    "correlation filtering, skipping cell correlation filtering...")
            scExp_cf = scExp_cf@metadata$Unfiltered
        }
        
        scExp_cf = consensus_clustering_scExp(scExp_cf, reps = 100,
                                              maxK = maxK,
                                              clusterAlg = "hc",
                                              prefix = file.path(cor_plot_dir,prefix))
        
        ### Choose most robust cluster #####
        average_consensus_score = scExp_cf@metadata$icl$clusterConsensus %>% 
            as.data.frame %>% dplyr::group_by(k) %>% summarise(mean = mean(clusterConsensus))
        average_consensus_score$diff = 0 
        average_consensus_score$diff[2:(maxK-1)] = average_consensus_score$mean[2:(maxK-1)] - average_consensus_score$mean[1:(maxK-2)]
        nclust = average_consensus_score$k[which.max(abs(average_consensus_score$diff))-1]
        
        if(!is.null(n_clust)) nclust =n_clust
        
        message("ChromSCape::generate_analysis - Choosing k = ", nclust, " as the optimal cluster number...")
        scExp_cf = choose_cluster_scExp(scExp_cf, nclust = nclust, consensus = TRUE)
        
        data = list("scExp_cf" = scExp_cf)
        qs::qsave(data, file = file.path(cor_dir, paste0(prefix,".qs")))
    }
    
    #### Coverage #####
    if("coverage" %in% run) {
        coverages = NULL
        if(input_data_type %in% c("scBED")){
            message("ChromSCape::generate_analysis - Creating pseudo-bulk cluster ",
                    "coverage tracks...")
            
            sample_folders = list.dirs(input_data_folder, full.names = T, recursive = F)
            input_files_coverage = sapply(sample_folders, function(i) list.files(i, full.names = T, pattern = ".bed|.bed.gz"))
            names(input_files_coverage) = basename(sample_folders)
            
            coverage_dir_nclust = file.path(coverage_dir, paste0(prefix, "_k", nclust))
            if(!dir.exists(coverage_dir_nclust)) dir.create(coverage_dir_nclust)
            
            system.time({generate_coverage_tracks(scExp_cf,
                                                  input_files_coverage,
                                                  odir = coverage_dir_nclust,
                                                  ref_genome = ref_genome
            )
            })
            
            coverages = sapply(list.files(coverage_dir_nclust,".bw", full.names = T), rtracklayer::import)
        }
    }
    
    #### Peak Calling #####
    if("peak_call" %in% run) {
        if(input_data_type %in% c("scBED")){
            message("ChromSCape::generate_analysis - Calling peaks on pseudo-bulk ",
                    "clusters ...")
            
            sample_folders = list.dirs(input_data_folder, full.names = T, recursive = F)
            input_files = sapply(sample_folders, function(i) list.files(i, full.names = T, pattern = ".bed|.bed.gz"))
            names(input_files) = basename(sample_folders)
            
            peak_dir_nclust = file.path(peak_dir, paste0(prefix, "_k", nclust))
            if(!dir.exists(peak_dir_nclust)) dir.create(peak_dir_nclust)
            
            scExp_cf = subset_bam_call_peaks(scExp_cf, odir = peak_dir_nclust,
                                             input = input_files, format = "scBED",
                                             p.value = 0.05, ref =  ref_genome, 
                                             peak_distance_to_merge = 10000)
            refined_annotation = scExp_cf@metadata$refined_annotation
            qs::qsave(refined_annotation, file = file.path(peak_dir_nclust,
                                                           "refined_annotation.qs"))
            
        }
    }
    #### Differential Analysis #####
    if("DA" %in% run) {
        message("ChromSCape::generate_analysis - Running one vs rest differential  ",
                "analysis ...")
        
        if(doBatchCorr) block = TRUE else block = FALSE
        scExp_cf = differential_analysis_scExp(scExp = scExp_cf,
                                               method= "wilcox",
                                               de_type = "one_vs_rest",
                                               cdiff.th = logFC.th,
                                               qval.th = qval.th,
                                               block = block)
    }
    
    #### Gene Sets Analysis #####
    if("GSA" %in% run){
        message("ChromSCape::generate_analysis - Running Gene Set Analysis ...")
        
        if("refined_annotation" %in% names(scExp_cf@metadata)) use_peaks = TRUE else
            use_peaks = FALSE
        MSIG.classes <- c("c1_positional","c2_curated","c3_motif","c4_computational",
                          "c5_GO","c6_oncogenic","c7_immunologic","hallmark")
        scExp_cf = gene_set_enrichment_analysis_scExp(
            scExp_cf,
            enrichment_qval = 0.01, qval.th = qval.th,
            ref = ref_genome,
            cdiff.th = logFC.th,
            peak_distance = 1000,
            use_peaks = use_peaks,
            GeneSetClasses = MSIG.classes)
        
        data = list("scExp_cf" = scExp_cf)
        qs::qsave(data,
                  file = file.path(da_gsa_dir,
                                   paste0(prefix,"_",nclust,"_",qval.th,"_",
                                          logFC.th,"_","one_vs_rest",".qs")))
    }
    
    #### Creating HTML report ####
    if("report" %in% run){
        message("ChromSCape::generate_analysis - Running Gene Set Analysis ...")
        
        rmarkdown::render(
            input = file.path("/media/pacome/LaCie/InstitutCurie/Documents/GitLab/ChromSCape/inst/template.Rmd"),
            output_file = file.path(ChromSCape_directory, paste0(analysis_name,"_report.html")),
            params = list(
                analysis_name = analysis_name,
                datamatrix = datamatrix,
                scExp = scExp,
                scExp_cf = scExp_cf,
                genes_to_plot = genes_to_plot,
                ref_genome = ref_genome,
                coverages = coverages,
                control_samples_CNA = control_samples_CNA
            ))
    }
    })
    message("ChromSCape::generate_analysis - Done ! ...")
    message("ChromSCape::generate_analysis - finished complete analysis in ",
            round(time_analysis[3]/60,2), " minutes...")
    
    out = list("scExp" = scExp, "scExp_cf" = scExp_cf)
    return(out)
}

rawData_to_datamatrix_annot <- function(input_data_folder,
                                        input_data_type = c("DenseMatrix", "SparseMatrix", "scBED", "scBAM")[1],
                                        feature_count_on = c("bins","geneTSS","peaks")[1],
                                        feature_count_parameter = 50000,
                                        ref_genome = c("hg38","mm10")[1])
    {
    if(input_data_type == "DenseMatrix"){
        file_list = list.files(input_data_folder, pattern =  ".gz|.txt|.tsv")
        tmp_list = import_scExp(file_names = basename(as.character(file_list)),
                                path_to_matrix = file_list)
        datamatrix = tmp_list$datamatrix
        annot_raw = tmp_list$annot_raw
        out = list(datamatrix = datamatrix, annot_raw = annot_raw)
    } else {
        selected_sample_folders = list.dirs(input_data_folder, recursive = F)
        names(selected_sample_folders) = basename(selected_sample_folders)
        send_warning = FALSE
        if(input_data_type == "scBAM") if(length(list.files(selected_sample_folders[1],pattern = "*.bam$"))==0) send_warning = TRUE
        if(input_data_type == "scBED") if(length(list.files(selected_sample_folders[1],pattern = "*.bed$|.*.bed.gz"))==0) send_warning = TRUE
        if(input_data_type == "SparseMatrix") {
            combin = expand.grid(c(".*features", ".*barcodes", ".*matrix"), c(".mtx",".tsv",".txt",".bed",".*.gz"))[-c(1,2,6,9,11,12),]
            pattern = paste(combin$Var1,combin$Var2, sep="", collapse = "|")
            if(length(list.files(selected_sample_folders[1],
                                 pattern = pattern))!=3) send_warning = TRUE
        }
        if(send_warning) {
            stop("Error : Can't find any specified file types in the upload folder. 
                                    Select another upload folder or another data type.")
        }
        
        if(input_data_type == "SparseMatrix"){
            out =  ChromSCape:::raw_counts_to_feature_count_files(
                files_dir_list = selected_sample_folders,
                file_type = input_data_type,
                ref = ref_genome)
            
        } else if(input_data_type %in% c("scBAM","scBED") &
                  !is.null(selected_sample_folders)) {
            
            if(feature_count_on == "bins"){
                stopifnot(is.double(feature_count_parameter),
                          feature_count_parameter >= 500)
                
                out = ChromSCape:::raw_counts_to_feature_count_files(
                    files_dir_list = selected_sample_folders,
                    file_type = input_data_type,
                    bin_width = round(feature_count_parameter),
                    ref = ref_genome)
            } 
            if(feature_count_on == "peaks"){
                stopifnot(is.character(feature_count_parameter),
                          file.exists(feature_count_parameter))
                out = ChromSCape:::raw_counts_to_feature_count_files(
                    files_dir_list = selected_sample_folders,
                file_type = input_data_type,
                peak_file = as.character(feature_count_parameter),
                ref = ref_genome)
            }
            if(feature_count_on == "geneTSS"){
                stopifnot(is.double(feature_count_parameter),
                          feature_count_parameter >= 100)
                out = ChromSCape:::raw_counts_to_feature_count_files(
                    files_dir_list = selected_sample_folders,
                file_type = input_data_type,
                geneTSS = TRUE,
                aroundTSS = feature_count_parameter,
                ref = ref_genome)
            }
        } else {
            stop("No data folder or data files selected.")
        }
    }
    
    return(out)
}


preprocessing_filtering_and_reduction <- function(
    datamatrix, annot_raw,
    min_reads_per_cell = 1600,
    max_quantile_read_per_cell = 95,
    min_percent_to_keep_feature = 0.5,
    subsample_n = NULL,
    ref_genome,
    exclude_regions = NULL,
    doBatchCorr  = FALSE,
    batch_sels  = NULL)
{
    
    scExp = create_scExp(datamatrix, annot_raw,
                         remove_zero_cells = TRUE,
                         remove_zero_features = TRUE)
    gc()
    
    scExp = filter_scExp(
        scExp,
        min_cov_cell = min_reads_per_cell,
        quant_removal = max_quantile_read_per_cell,
        percentMin = min_percent_to_keep_feature)
    gc()
    
    # Filtering based on exclude-regions from bed file, if provided
    if (!is.null(exclude_regions))
    {
        exclude_regions = import(exclude_regions)
        scExp = exclude_features_scExp(scExp, exclude_regions, by = "region")
        gc()
    }
    
    ### 2.bis Optional subsampling ###
    if(length(subsample_n)>0){
        scExp = subsample_scExp(scExp, n_cells = subsample_n)
        gc()
        
    }
    
    ### 3. Normalizing ###
    scExp = normalize_scExp(scExp, type = "CPM")
    gc()
    
    ### 4. Feature annotation ###
    scExp = feature_annotation_scExp(scExp, ref = ref_genome)
    gc()
    
    # Original PCA
    print("Running Dimensionality Reduction...")
    
    scExp = reduce_dims_scExp(
        scExp, dimension_reductions = c("PCA","UMAP"),
        batch_correction = doBatchCorr, batch_list = batch_sels, 
        verbose = FALSE)
    gc()
    
    ### 7. Add default colors ###
    
    if(doBatchCorr){ annotCol. = c("sample_id","total_counts", "batch_name")} 
    else{annotCol. = c("sample_id","total_counts")}
    scExp = colors_scExp(scExp, annotCol.)
    
    ### 8. Running hierarchical clustering ###
    scExp = correlation_and_hierarchical_clust_scExp(scExp)
    gc()
    
    return(scExp)
}
