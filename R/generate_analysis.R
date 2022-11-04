## Run a complete unsupervised analysis and create an HTML report

#' Generate a complete ChromSCape analysis
#' @usage generate_analysis(input_data_folder,
#' analysis_name = "Analysis_1",
#' output_directory = "./",
#' input_data_type = c("scBED", "DenseMatrix", "SparseMatrix", "scBAM")[1],
#' rebin_sparse_matrix = FALSE,
#' feature_count_on = c("bins","genebody","peaks")[1],
#' feature_count_parameter = 50000,
#' ref_genome = c("hg38","mm10")[1],
#' run = c("filter", "CNA","cluster", "consensus","peak_call", "coverage", 
#'        "DA", "GSA", "report")[c(1,3,6,7,8,9)],
#' min_reads_per_cell = 1000,
#' max_quantile_read_per_cell = 99,
#' n_top_features = 40000,
#' norm_type = "CPM",
#' subsample_n = NULL,
#' exclude_regions = NULL,
#' n_clust = NULL,
#' corr_threshold = 99,
#' percent_correlation = 1,
#' maxK = 10,
#' qval.th = 0.1,
#' logFC.th = 1,
#' enrichment_qval = 0.1,
#' doBatchCorr  = FALSE,
#' batch_sels  = NULL,
#' control_samples_CNA = NULL,
#' genes_to_plot = c("Krt8","Krt5","Tgfb1", "Foxq1", "Cdkn2b",
#'                  "Cdkn2a", "chr7:15000000-20000000")
#' )
#' 
#' @param input_data_folder Directory containing the input data.
#' @param analysis_name Name given to the analysis.
#' @param output_directory Directory where to create the analysis and the 
#' HTML report.
#' @param input_data_type The type of input data. 
#' @param feature_count_on For raw data type, on which features to count the 
#' cells.
#' @param feature_count_parameter Additional parameter corresponding to the 
#' 'feature_count_on' parameter. E.g. for 'bins' must be a numeric, e.g. 50000, 
#' for 'peaks' must be a character containing path towards a BED peak file.
#' @param rebin_sparse_matrix A boolean specifying if the SparseMatrix should
#' be rebinned on features (see feature_count_on and feature_count_parameter).
#' @param ref_genome The genome of reference.
#' @param run What steps to run. By default runs everything. Some steps are
#' required in order to run downstream steps.
#' @param min_reads_per_cell Minimum number of reads per cell.
#' @param max_quantile_read_per_cell Upper quantile above which to consider 
#' cells doublets.
#' @param n_top_features Number of features to keep in the analysis. 
#' @param norm_type Normalization type.
#' @param subsample_n Number of cells per condition to downsample to, for 
#' performance principally.
#' @param exclude_regions Path towards a BED file containing CNA to exclude from
#' the analysis (optional).
#' @param n_clust Number of clusters to force choice of clusters.
#' @param corr_threshold Quantile of correlation above which two cells are 
#' considered as correlated. 
#' @param percent_correlation Percentage of the total cells that a cell must be
#' correlated with in order to be kept in the analysis.
#' @param maxK Upper cluster number to rest for ConsensusClusterPlus.
#' @param qval.th Adjusted p-value below which to consider features 
#' differential.
#' @param logFC.th Log2-fold-change above/below which to consider a feature 
#' depleted/enriched. 
#' @param enrichment_qval Adjusted p-value below which to consider a gene set as
#' significantly enriched in differential features. 
#' @param doBatchCorr Logical indicating if batch correction using fastMNN
#'  should be run.
#' @param batch_sels If doBatchCorr is TRUE, a named list containing the samples
#' in each batch.
#' @param control_samples_CNA If running CopyNumber Analysis, a character vector
#' of the sample names that are 'normal'. 
#' @param genes_to_plot A character vector containing genes of interest of which
#' to plot the coverage.
#'
#' @return Creates a ChromSCape-readable directory and saved objects, as well as
#' a multi-tabbed HTML report resuming the analysis.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' generate_analysis("/path/to/data/", "Analysis_1")
#' }
generate_analysis <- function(input_data_folder,
                  analysis_name = "Analysis_1",
                  output_directory = "./",
                  input_data_type = c("scBED", "DenseMatrix", "SparseMatrix",
                                      "scBAM")[1],
                  feature_count_on = c("bins","genebody","peaks")[1],
                  feature_count_parameter = 50000,
                  rebin_sparse_matrix = FALSE,
                  ref_genome = c("hg38","mm10")[1],
                  run = c("filter", "CNA","cluster", "consensus", "coverage", 
                          "DA", "GSA", "report")[c(1,3,5,6,7,8)],
                  min_reads_per_cell = 1000,
                  max_quantile_read_per_cell = 99,
                  n_top_features = 40000,
                  norm_type = "CPM",
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
              is.numeric(n_top_features),
              is.numeric(qval.th),
              is.numeric(logFC.th),
              is.numeric(enrichment_qval),
              is.character(norm_type)
              )
    
    #### Create Analysis directory ####
    ChromSCape_directory <- create_project_folder(output_directory,
                                                  analysis_name,
                                                  ref_genome)
    
    batch_string <- if (doBatchCorr) "batchCorrected" else "uncorrected"
    prefix = paste0(analysis_name, "_", min_reads_per_cell, "_",
           n_top_features, "_",
           max_quantile_read_per_cell, "_", batch_string)
    
    time_analysis = system.time({
    #### Select & Import ####
    message("ChromSCape::generate_analysis - Importing datasets ...")
    if(input_data_type == "DenseMatrix"){
        out <- import_scExp(file_paths = list.files(input_data_folder, full.names = TRUE))
    } else {
        out <- rawData_to_datamatrix_annot(input_data_folder, input_data_type,
                                           feature_count_on, feature_count_parameter,
                                           ref_genome
        )
    }
        
    
    datamatrix = out$datamatrix
    if(input_data_type == "SparseMatrix" & rebin_sparse_matrix ){
      regions = head(rownames(datamatrix), 1000)
      start = as.numeric(gsub(".*_","", gsub("_\\d*$","", regions)))
      end = as.numeric(gsub(".*_","", regions))
      overlap = ceiling(mean((end - start) / 2))
      message("Rebinning the sparse matrix into ", feature_count_on, " with parameter ", feature_count_parameter, " and overlap of ", overlap)
      message("Saving raw matrix for later usage for coverage or additional feature engineering...")
      qs::qsave(datamatrix, file.path(ChromSCape_directory, "raw_mat.qs"))

      if(feature_count_on == "bins") {
        datamatrix = rebin_matrix(datamatrix,
                                  bin_width = feature_count_parameter,
                                  minoverlap = overlap, ref = ref_genome)
      } else{
        custom_feature = rtracklayer::import(feature_count_parameter)
        datamatrix = rebin_matrix(datamatrix,
                                  custom_annotation = custom_feature,
                                  minoverlap = overlap, ref = ref_genome)
      }

    }
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
        n_top_features = n_top_features,
        norm_type = norm_type,
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
    rm(annot_raw, datamatrix)
    gc()
    qs::qsave(scExp, file = file.path(ChromSCape_directory,
                                      "Filtering_Normalize_Reduce",
                                      paste0(prefix, ".qs")))
    
    #### Correlation Filtering & Clustering ####
    if("cluster" %in% run) {
        message("ChromSCape::generate_analysis - Correlation Consensus Clustering ...")
        system.time({
        scExp_cf = filter_correlated_cell_scExp(scExp, random_iter = 50, corr_threshold = corr_threshold,
                                                percent_correlation = percent_correlation)
        })
        
        # Check that number of cells and features is enough or finish here
        n_cell = ncol(scExp_cf)
        if( (n_cell < 100)){
            warning("Analyis ", analysis_name, " has too few cells after cell ",
                    "correlation filtering, skipping cell correlation filtering...")
            scExp_cf = scExp_cf@metadata$Unfiltered
        }
        rm(scExp)
        gc()
        
        if("consenus" %in% run){
        scExp_cf = consensus_clustering_scExp(scExp_cf, reps = 100,
                                              maxK = maxK,
                                              clusterAlg = "hc",
                                              prefix = file.path(
                                                  ChromSCape_directory,
                                                  "correlation_clustering",
                                                  "Plots", prefix))
        
        ### Choose most robust cluster #####
        average_consensus_score = scExp_cf@metadata$icl$clusterConsensus %>% 
            as.data.frame %>% dplyr::group_by(k) %>% summarise(mean = mean(clusterConsensus))
        average_consensus_score$diff = 0 
        average_consensus_score$diff[2:(maxK-1)] = average_consensus_score$mean[2:(maxK-1)] - average_consensus_score$mean[1:(maxK-2)]
        nclust = average_consensus_score$k[which.max(abs(average_consensus_score$diff))-1]
        
        if(!is.null(n_clust)) nclust = n_clust
        }
        if(is.null(n_clust) && !'consensus' %in% run) nclust = 5 else nclust = n_clust
            
        message("ChromSCape::generate_analysis - Choosing k = ", nclust, " as the optimal cluster number...")
        scExp_cf = choose_cluster_scExp(scExp_cf, nclust = nclust, consensus = 'consensus' %in% run)
        
        data = list("scExp_cf" = scExp_cf)
        qs::qsave(data, file = file.path(ChromSCape_directory,
                                         "correlation_clustering",
                                         paste0(prefix,".qs")))
        rm(data)
        gc()
    }
    
    #### Coverage #####
    if("coverage" %in% run) {
        coverages = NULL
        format = input_data_type
        if(input_data_type %in% c("scBED", "SparseMatrix")){
            message("ChromSCape::generate_analysis - Creating pseudo-bulk cluster ",
                    "coverage tracks...")
            
            sample_folders = list.dirs(input_data_folder, full.names = TRUE,
                                       recursive = FALSE)
            if(input_data_type == "scBED"){
              input = sapply(sample_folders, function(i)
                list.files(i, full.names = TRUE, pattern = ".bed|.bed.gz"))
            names(input) = basename(sample_folders)
            } else {
              format = "raw_mat"
              input = qs::qread(file.path(ChromSCape_directory, "raw_mat.qs"))
            }
            coverage_dir_nclust = file.path(ChromSCape_directory,
                                            "coverage", paste0(prefix, "_k", nclust))
            if(!dir.exists(coverage_dir_nclust)) dir.create(coverage_dir_nclust)
            
            system.time({
              generate_coverage_tracks(
                scExp_cf,
                input = input,
                odir = coverage_dir_nclust,
                format = format,
                ref_genome = ref_genome,
                bin_width = 150,
                n_smoothBin = 5,
                read_size = 101,
                quantile_for_peak_calling = 0.85,
                by = "cell_cluster"
              )
              
            })
        }
    }
   
    #### Differential Analysis #####
    if("DA" %in% run) {
        message("ChromSCape::generate_analysis - Running one vs rest differential  ",
                "analysis ...")
        
        if(doBatchCorr) block = TRUE else block = FALSE
        scExp_cf = differential_analysis_scExp(scExp = scExp_cf,
                                               method= "wilcox",
                                               de_type = "one_vs_rest_fast",
                                               block = block)
        gc()
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
            enrichment_qval = 0.01,
            qval.th = qval.th,
            ref = ref_genome,
            logFC.th = logFC.th,
            peak_distance = 1000,
            use_peaks = use_peaks,
            GeneSetClasses = MSIG.classes)
        
        data = list("scExp_cf" = scExp_cf)
        qs::qsave(data,
                  file = file.path(ChromSCape_directory,
                                   "Diff_Analysis_Gene_Sets",
                                   paste0(prefix,"_",nclust,"_",qval.th,"_",
                                          logFC.th,"_","one_vs_rest",".qs")))
        rm(data)
        rm(scExp_cf)
        gc()
    }
    
    #### Creating HTML report ####
    if("report" %in% run){

        message("ChromSCape::generate_analysis - Running Gene Set Analysis ...")
        generate_report(ChromSCape_directory = ChromSCape_directory,
                        prefix = prefix,
                        run = run,
                        genes_to_plot = genes_to_plot,
                        control_samples_CNA = control_samples_CNA
                        )
    }
    gc()
    })
    
    message("ChromSCape::generate_analysis - Done ! ...")
    message("ChromSCape::generate_analysis - finished complete analysis in ",
            round(time_analysis[3]/60,2), " minutes...")
    scExp = qs::qread(file = file.path(ChromSCape_directory,
                                       "Filtering_Normalize_Reduce",
                                       paste0(prefix, ".qs")))
    gc()
    scExp_cf = qs::qread(file.path(ChromSCape_directory,
                                   "Diff_Analysis_Gene_Sets",
                                   paste0(prefix,"_",nclust,"_",qval.th,"_",
                                          logFC.th,"_","one_vs_rest",".qs")))$scExp_cf
    gc()
    out = list("scExp" = scExp, "scExp_cf" = scExp_cf)
    return(out)
}

rawData_to_datamatrix_annot <- function(input_data_folder,
                                        input_data_type = c("DenseMatrix", "SparseMatrix", "scBED", "scBAM")[1],
                                        feature_count_on = c("bins","genebody","peaks")[1],
                                        feature_count_parameter = 50000,
                                        ref_genome = c("hg38","mm10")[1])
    {
    if(input_data_type == "DenseMatrix"){
        file_list = list.files(input_data_folder, pattern =  ".gz|.txt|.tsv")
        tmp_list = import_scExp(file_paths = basename(as.character(file_list)),
                                temp_path = file_list)
        datamatrix = tmp_list$datamatrix
        annot_raw = tmp_list$annot_raw
        out = list(datamatrix = datamatrix, annot_raw = annot_raw)
    } else {
        selected_sample_folders = list.dirs(input_data_folder, recursive = FALSE)
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
            out =  read_sparse_matrix(
                files_dir_list = selected_sample_folders,
                ref = ref_genome)
            
        } else if(input_data_type %in% c("scBAM","scBED") &
                  !is.null(selected_sample_folders)) {
            
            if(feature_count_on == "bins"){
                stopifnot(is.double(feature_count_parameter),
                          feature_count_parameter >= 500)
                
                out = raw_counts_to_sparse_matrix(
                    files_dir_list = selected_sample_folders,
                    file_type = input_data_type,
                    bin_width = round(feature_count_parameter),
                    ref = ref_genome)
            } 
            if(feature_count_on == "peaks"){
                stopifnot(is.character(feature_count_parameter),
                          file.exists(feature_count_parameter))
                out = raw_counts_to_sparse_matrix(
                    files_dir_list = selected_sample_folders,
                file_type = input_data_type,
                peak_file = as.character(feature_count_parameter),
                ref = ref_genome)
            }
            if(feature_count_on == "genebody"){
                stopifnot(is.double(feature_count_parameter),
                          feature_count_parameter >= 100)
                out = raw_counts_to_sparse_matrix(
                    files_dir_list = selected_sample_folders,
                file_type = input_data_type,
                genebody = TRUE,
                extendPromoter = feature_count_parameter,
                ref = ref_genome)
            }
        } else {
            stop("No data folder or data files selected.")
        }
    }
    
    return(out)
}

#' Create ChromSCape project folder 
#' 
#' @description 
#' Creates a project folder that will be recognizable by ChromSCape Shiny
#' application. 
#' 
#' @param output_directory Path towards the directory to create the
#'  'ChromSCape_Analyses' folder and the analysis subfolder. If this path 
#'  already contains the 'ChromSCape_Analyses' folder, will only create the 
#'  analysis subfolder.
#' @param analysis_name Name of the analysis. Must only contain alphanumerical
#' characters or '_'. 
#' @param ref_genome  Reference genome, either 'hg38' or 'mm10'.
#'
#' @return Creates the project folder and returns the root of the project.
#' 
#' @export
#'
#' @examples 
#' dir = tempdir()  
#' create_project_folder(output_directory = dir,
#'  analysis_name = "Analysis_1")
#' list.dirs(file.path(dir))
create_project_folder <- function(output_directory,
                                  analysis_name = "Analysis_1",
                                  ref_genome = c("hg38","mm10")[1]){
    stopifnot(is.character(output_directory), dir.exists(output_directory),
              is.character(analysis_name), is.character(ref_genome),
              ref_genome %in% c("hg38","mm10"))

    if(!grepl("^[A-Za-z0-9_]+$", analysis_name))
        stop("ChromSCape::create_project_folder - The analysis name must ",
        "contain only alphanumerical characters or '_'.")
    
    # Dont create "ChromSCape_analyses" directory if it is created already
    if(!grepl("ChromSCape_analyses/$|ChromSCape_analyses$", output_directory)) {
        ChromSCape_analyses = file.path(output_directory, "ChromSCape_analyses")
        if(!dir.exists(ChromSCape_analyses)) dir.create(ChromSCape_analyses) 
    } else{
        ChromSCape_analyses = file.path(output_directory)
    }
    
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
    return(ChromSCape_directory)
}

#' Preprocess and filter matrix annotation data project folder to SCE
#'
#' @param datamatrix A sparse count matrix of features x cells.
#' @param annot_raw A data.frame with barcode, cell_id, sample_id, batch_id, total_counts 
#' @param min_reads_per_cell Minimum read per cell to keep the cell
#' @param max_quantile_read_per_cell Upper count quantile threshold above which cells are removed
#' @param n_top_features Number of features to keep
#' @param norm_type Normalization type c("CPM", "TFIDF", "RPKM", "TPM", "feature_size_only")
#' @param n_dims An integer specifying the number of dimensions to keep for PCA
#' @param subsample_n Number of cells to subsample.
#' @param ref_genome Reference genome ("hg38" or "mm10").
#' @param exclude_regions GenomicRanges with regions to remove from the object.
#' @param remove_PC A vector of string indicating which principal components to 
#' remove before downstream analysis as probably correlated to 
#' library size. Should be under the form : 'Component_1', 'Component_2', ...
#' Recommended when using 'TFIDF' normalization method. (NULL)
#' @param doBatchCorr Run batch correction ? TRUE or FALSE
#' @param batch_sels If doBatchCorr is TRUE, List of characters.
#'  Names are batch names, characters are sample names. 
#'
#' @return A SingleCellExperiment object containing feature spaces.
#' @export
#'
#' @examples
#' raw <- create_scDataset_raw()
#' scExp = preprocessing_filtering_and_reduction(raw$mat, raw$annot)
preprocessing_filtering_and_reduction <- function(
    datamatrix, annot_raw,
    min_reads_per_cell = 1600,
    max_quantile_read_per_cell = 95,
    n_top_features = 40000,
    norm_type = "CPM",
    n_dims = 10,
    remove_PC = NULL,
    subsample_n = NULL,
    ref_genome = "hg38",
    exclude_regions = NULL,
    doBatchCorr  = FALSE,
    batch_sels  = NULL)
{
    
    scExp = create_scExp(datamatrix, annot_raw,
                         remove_zero_cells = TRUE,
                         remove_zero_features = TRUE)
    gc()
    
    # Filtering based on exclude-regions from bed file, if provided
    if (!is.null(exclude_regions))
    {
        scExp = exclude_features_scExp(scExp, exclude_regions, by = "region")
        gc()
    }
    
    scExp = filter_scExp(
        scExp,
        min_cov_cell = min_reads_per_cell,
        quant_removal = max_quantile_read_per_cell)
    gc()
    
    scExp = find_top_features(
        scExp,
        n = n_top_features,
        keep_others = FALSE)
    
    ### 2.bis Optional subsampling ###
    if(length(subsample_n)>0){
        scExp = subsample_scExp(scExp,n_cell_per_sample = subsample_n)
        gc()
        
    }
    
    ### 3. Normalizing ###
    scExp = normalize_scExp(scExp, type = norm_type)
    gc()
    
    ### 4. Feature annotation ###
    scExp = feature_annotation_scExp(scExp, ref = ref_genome)
    gc()
    
    # Original PCA
    print("Running Dimensionality Reduction...")
    
    scExp = reduce_dims_scExp(
        scExp, n = n_dims,
        dimension_reductions = c("PCA","UMAP"),
        batch_correction = doBatchCorr, batch_list = batch_sels,  
        remove_PC = remove_PC,
        verbose = FALSE)
    gc()
    
    ### 7. Add default colors ###
    
    if(doBatchCorr){
      annotCol. = c("sample_id","total_counts", "batch_name")
    }  else{
      annotCol. = c("sample_id","total_counts")
    }
    scExp = colors_scExp(scExp, annotCol.)
    
    ### 8. Running hierarchical clustering ###
    scExp = correlation_and_hierarchical_clust_scExp(scExp)
    gc()
    
    return(scExp)
}

#' From a ChromSCape analysis directory, generate an HTML report.
#'
#' @param ChromSCape_directory Path towards the ChromSCape directory of which 
#' you want to create the report. The report will be created at the root of this
#' directory.
#' 
#' @param prefix Name of the analysis with the filtering parameters
#'  (e.g. Analysis_3000_100000_99_uncorrected). You will find the prefix in the
#'  Filtering_Normalize_Reduce subfolder.
#' @param run Which steps to report ("filter", "CNA","cluster", "consensus", "peak_call",
#' "coverage", "DA", "GSA", "report"). Only indicate steps that were done in the 
#' analysis. By default do not report CNA, consensus and peak calling.
#' @param genes_to_plot For the UMAP, which genes do you want to see in the 
#' report.
#' @param control_samples_CNA If running the Copy Number Alteration (CNA) part,
#' which samples are the controls
#'
#' @return Generate an HTML report at the root of the analysis directory.
#'
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' generate_analysis("/path/to/data/", "Analysis_1")
#' }
generate_report <- function(ChromSCape_directory,
                            prefix = NULL,
                            run = c("filter", "CNA","cluster", "consensus", "peak_call",
                                    "coverage", "DA", "GSA", "report")[c(1,3,6,7,8,9)],
                            genes_to_plot =  c(
                                "Krt8","Krt5","Tgfb1", "Foxq1", "Cdkn2b",
                                "Cdkn2a", "chr7:15000000-20000000"),
                            control_samples_CNA = NULL
                            ){
    
    ChromSCape_directory = ChromSCape_directory
    filt_dir <- file.path(ChromSCape_directory,"Filtering_Normalize_Reduce") 
    coverage_dir <- file.path(ChromSCape_directory,"coverage")
    cor_dir <- file.path(ChromSCape_directory, "correlation_clustering")
    da_gsa_dir <- file.path(ChromSCape_directory, "Diff_Analysis_Gene_Sets")
    ref_genome = read.table(file = file.path(ChromSCape_directory, "annotation.txt"))
    
    analysis_name = basename(ChromSCape_directory)
    datamatrix = qs::qread(file.path(ChromSCape_directory, "datamatrix.qs"))
    
    if("filter" %in% run)
        scExp_files = list.files(file.path(filt_dir), full.names = TRUE) else 
            scExp_files = ""
    if("DA" %in% run)
        scExp_cf_files = list.files(file.path(da_gsa_dir),
                            full.names = TRUE) else 
                            if("cluster" %in% run)
                            scExp_cf_files = list.files(file.path(cor_dir),
                                full.names = TRUE) else scExp_cf_files = ""
    if("coverage" %in% run) 
        coverage_dirs = list.dirs(file.path(coverage_dir),
                recursive = FALSE, full.names = TRUE) else coverage_dirs = ""
    
    if(!is.null(prefix)){
        scExp = qs::qread(scExp_files[grep(prefix,scExp_files)][1])
        scExp_cf = qs::qread(
            scExp_cf_files[grep(prefix,scExp_cf_files)][1])$scExp_cf
        coverage_dir_nclust = coverage_dirs[grep(prefix,coverage_dirs)][1]
        coverages = sapply(list.files(coverage_dir_nclust,".bw",
                                      full.names = TRUE), rtracklayer::import)
    } else{
        scExp = qs::qread(scExp_files[1])
        scExp_cf = qs::qread(scExp_cf_files[1])$scExp_cf
        coverage_dir_nclust = coverage_dirs[1]
        coverages = sapply(list.files(coverage_dir_nclust,".bw",
                                      full.names = TRUE), rtracklayer::import)
    }
    if(shiny::is.reactive(scExp)) scExp = shiny::isolate(scExp())
    if(shiny::is.reactive(scExp_cf)) scExp = shiny::isolate(scExp_cf())
    
    rmarkdown::render(
        input = file.path(system.file(package="ChromSCape","template.Rmd")),
        output_file = file.path(ChromSCape_directory, paste0(analysis_name, "_report.html")),
        params = list(
            analysis_name = analysis_name,
            datamatrix = datamatrix,
            run = run,
            scExp = scExp,
            scExp_cf = scExp_cf,
            genes_to_plot = genes_to_plot,
            ref_genome = ref_genome,
            coverages = coverages,
            control_samples_CNA = control_samples_CNA
        ))
}
