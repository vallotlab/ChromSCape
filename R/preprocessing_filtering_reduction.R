## Authors : Pacôme Prompsy Title : Wrappers & functions to preprocess & reduce
## dimensionality of single-cell Matrix Description : Funtions to load, filter,
## normalize & reduce single-cell Epigenetic Matrices prior to analysis

#' Heuristic discovery of samples based on cell labels
#'
#' Identify a fixed number of common string (samples) in a set of varying
#' strings (cells). E.g. in the set
#' "Sample1_cell1","Sample1_cell2","Sample2_cell1","Sample2_cell2" and with
#' nb_samples=2, the function returns "Sample1","Sample1","Sample2","Sample2".
#'
#'
#' @param barcodes Vector of cell barcode names (e.g. Sample1_cell1,
#'   Sample1_cell2...)
#' @param nb_samples Number of samples to find
#'
#' @return character vector of sample names the same length as cell labels
#' @export
#'
#' @importFrom stringdist stringdistmatrix
#' @importFrom qualV LCS
#' @examples
#'
#' barcodes = c(paste0("HBCx22_BC_",seq_len(100)),
#' paste0("mouse_sample_XX",208:397))
#' samples = detect_samples(barcodes, nb_samples=2)
#' 
detect_samples <- function(barcodes, nb_samples = 1) {
    t = system.time({
        mat = stringdist::stringdistmatrix(
            barcodes, barcodes, useNames = "strings",
            method = "lv", weight = c(1, 1, 1))
        
        dist = stats::as.dist(1 - stats::cor(mat))
        hc = stats::hclust(dist, method = 'ward.D')
        hc$labels = rep("", ncol(mat))
    })
    message("ChromSCape::detect_samples - found samples in ", t[3], " secs.\n")
    sample_groups = stats::cutree(hc, k = nb_samples)
    samp_name = c()
    for (i in seq_len(nb_samples)) {
        x = barcodes[which(sample_groups == i)]
        samp = as.character(sample(x, min(length(x), 10), replace = TRUE))
        lcs = strsplit(samp[1], "")[[1]]
        for (j in 2:min(length(x), 10)) {
            lcs = unlist(qualV::LCS(lcs, strsplit(samp[j], "")[[1]])[4])
        }
        samp_name = c(samp_name, paste(lcs, collapse = ""))
    }
    samples_names = gsub("[[:punct:]]$", "", samp_name)
    samples_names = gsub("^[[:punct:]]", "", samples_names)
    if (nb_samples > 1) {
        mat <- create_sample_name_mat(nb_samples,samples_names)
        longest_common_between_samples = 
            mat[which.min(unlist(lapply(mat, nchar)))][1]
        if (length(longest_common_between_samples) > 2) {
            lc <- grep(pattern = longest_common_between_samples, samples_names)
            if (length(lc) == nb_samples)
                samples_names = gsub(
                    longest_common_between_samples, "", samples_names)
        }
    }
    samples_names = samples_names[sample_groups]
    return(samples_names)
}

#' Create a sample name matrix 
#'
#' @param nb_samples Number of samples
#' @param samples_names Character vector of sample names
#'
#' @return A matrix
#'
create_sample_name_mat <- function(nb_samples, samples_names){
    mat = matrix("", nrow = nb_samples, ncol = nb_samples)
    for (i in seq_len(nb_samples)) {
        for (j in seq_len(nb_samples)) {
            x = as.character(samples_names[i])
            y = as.character(samples_names[j])
            x = strsplit(x, "")[[1]]
            y = strsplit(y, "")[[1]]
            lcs = qualV::LCS(x, y)
            ls = split(lcs[6]$va, cumsum(c(TRUE, diff(lcs[6]$va) != 1)))
            ls = ls[[which.max(lengths(ls))]]
            x = as.character(samples_names[i])
            x = strsplit(x, "")[[1]]
            lcs = paste0(x[ls], collapse = "")
            if (length(lcs) == 0) lcs = ""
            mat[i, j] = lcs
        }
    }
    
    mat = as.character(mat)
    return(mat)
}


#' Read in one or multiple sparse matrices (10X format)
#' 
#' @description 
#' Given one or multiple directories, look in each directory for a combination 
#' of the following files :
#' - A 'features' file containing unique feature genomic locations 
#' -in tab separated format ( *_features.bed / .txt / .tsv / .gz), 
#' e.g. chr, start and end
#' - A 'barcodes' file containing unique barcode names 
#' ( _barcode.txt / .tsv / .gz)
#' - A 'matrix' A file containing indexes of non zero entries 
#' (_matrix.mtx / .gz)
#'
#' @param files_dir_list A named character vector containing the full path 
#' towards folders. Each folder should contain only the Feature file, the 
#' Barcode file and the Matrix file (see description). 
#' @param verbose Print ?
#' @param ref Reference genome (used to filter non-canonical chromosomes).
#'
#' @return Returns a list containing a datamatrix and cell annotation
#' @export
#'
#' @importFrom Matrix readMM
#' @importFrom GenomicRanges seqnames start end
#' @importFrom rtracklayer import.bed
#' @examples
#' \dontrun{
#' sample_dirs = c("/path/to/folder1/", "/path/to/folder2/")
#' names(sample_dirs) = c("sample_1", "sample_2")
#' out <- read_sparse_matrix(sample_dirs, ref = "hg38")
#' head(out$datamatrix)
#' head(out$annot_raw)
#' }
#' 
read_sparse_matrix <- function(files_dir_list,
                               ref = c("hg38","mm10")[1],
                               verbose = TRUE
                              ){
    
    feature_file <- barcode_file <- matrix_file <- NULL
    
    pattern = ".*features.tsv|.*features.txt|.*features.bed|.*features.*.gz"
    feature_file = list.files(path = files_dir_list[1], full.names = TRUE,
                              pattern = pattern)
    pattern = ".*barcodes.tsv|.*barcodes.txt|.*barcodes.*.gz"
    barcode_file = list.files(path = files_dir_list[1], full.names = TRUE,
                              pattern = pattern)
    pattern = ".*matrix.mtx|.*matrix.*.gz"
    matrix_file = list.files(path = files_dir_list[1], full.names = TRUE,
                             pattern = pattern)
    
    if (length(c(feature_file, matrix_file, barcode_file)) != 3)
        stop(paste0(
            "ChromSCape::read_sparse_matrix - For ",
            "SparseMatrix Count type, the folder must contain exactly two files ",
            "matching respectively *index.txt, peaks.bed, barcodes.txt"))
    
    feature_indexes = list()
    samples_ids = barcodes = c()
    for(i in seq_along(files_dir_list)){
        t1 = system.time({
            sample_id = names(files_dir_list)[i]
            files = list.files(files_dir_list[i], full.names = TRUE)
            feature_file = files[grep("feature",files)][1]
            barcode_file = files[grep("barcode",files)][1]
            matrix_file = files[grep("matrix",files)][1]
            
            mat = Matrix::readMM(matrix_file)
            barcode = as.character(read.table(barcode_file)[,1])
            name_cell <- paste0(sample_id, "_", barcode)
            barcodes <- c(barcodes,barcode)
            
            which = tryCatch({rtracklayer::import.bed(feature_file)},
                             error = function(e) NULL) 
            if(is.null(which)){
                separator <- separator_count_mat(feature_file)
                format_test = read.table(feature_file, header = TRUE,
                                         sep = separator, nrows = 5)
                separated_chr_start_end = c(
                    grep("chr", colnames(format_test)[seq_len(3)]),
                    grep("start|begin", colnames(format_test)[seq_len(3)]),
                    grep("end|stop", colnames(format_test)[seq_len(3)]))
                if (length(separated_chr_start_end) > 0 &&
                    all.equal(separated_chr_start_end, c(1, 2, 3)) == TRUE) {
                    val = c("character",
                            "integer",
                            "integer",
                            rep("NULL", ncol(format_test) -3))
                    which = read.table(
                        feature_file,
                        header = FALSE,
                        sep = separator,
                        colClasses = val
                    )
                    which = as(which, "GRanges")
                } else{
                    which = read.table(feature_file,
                                       sep = separator,
                                       header = FALSE,
                                       )
                    if(ncol(which) > 3 ){
                        which = as(
                        setNames(which[,1:3],c("chr","start","end")), "GRanges")
                    } else {
                        which = as(which[,1, drop=FALSE] %>% tidyr::separate(
                        col = 1, into = c("chr","start","end"),sep =":|-"),
                        "GRanges")
                    }
                }
                
            }
            rownames(mat) = paste0(GenomicRanges::seqnames(which),"_",
                                   GenomicRanges::start(which),"_",
                                   GenomicRanges::end(which))
            colnames(mat) = name_cell
            feature_indexes[[sample_id]] = mat
            samples_ids = c(samples_ids, rep(sample_id, length(name_cell)))
        })
        if (verbose)
            message("ChromSCape::read_sparse_matrix - ",
                    "Count matrix ", sample_id ,
                    " created from Index-Peak-Barcodes files in ",
                    round(t1[3],3), " sec.")
    }

    tryCatch({mat <- do.call("cbind", feature_indexes)}, error = function(e){
        paste0("ChromSCape::read_sparse_matrix - All feature file do not have ",
               "the same number of features", e)})
    
    eval(parse(text = paste0("data(", ref, ".chromosomes)")))
    chr <- eval(parse(text = paste0("", ref, ".chromosomes")))

    regions_to_remove = which(!as.character(which@seqnames) %in% chr$chr)
    if(length(regions_to_remove) > 0) mat = mat[-regions_to_remove,]
    annot_raw = data.frame("barcode" = barcodes,
                           "cell_id" = colnames(mat),
                           "sample_id" = samples_ids,
                           "batch_id" = factor(rep(1, ncol(mat)))
    )
    out = list("datamatrix" = mat, "annot_raw" = annot_raw)
    return(out)
}

#' Create a sparse count matrix from various format of input data.
#'
#' This function takes three different type of single-cell input: - Single cell
#' BAM files (sorted) - Single cell BED files (gzipped) - A combination of an
#' index file, a peak file and cell barcode file (The index file is composed of
#' three column: index i, index j and value x for the non zeroes entries in the
#' sparse matrix.)
#'
#' This functions re-counts signal on either fixed genomic bins, a set of
#' user-defined peaks or around the TSS of genes.
#'
#' @param ref reference genome to use (hg38)
#' @param verbose Verbose (TRUE)
#' @param use_Signac Use Signac wrapper function 'FeatureMatrix' if the Signac
#' package is installed (TRUE).
#' @param files_dir_list A named character vector of directories containing
#'  the files. The names correspond to sample names.
#' @param peak_file A file containing genomic location of peaks (NULL)
#' @param n_bins The number of bins to tile the genome (NULL)
#' @param bin_width The size of bins to tile the genome (NULL)
#' @param genebody Count on genes (body + promoter) ? (NULL)
#' @param extendPromoter If counting on genes, number of base pairs to extend up or
#'  downstream of TSS (2500).
#' @param file_type Input file(s) type(s) ('scBED','scBAM','FragmentFile')
#' @param progress Progress object for Shiny
#' @param BPPARAM BPPARAM object for multiprocessing. See
#'  \link[BiocParallel]{bpparam} for more informations. Will take the default
#'  BPPARAM set in your R session.
#'  
#' @return A sparse matrix of features x cells
#' @source 
#' @references Stuart el al.,  Multimodal single-cell chromatin analysis with 
#' Signac bioRxiv \url{https://doi.org/10.1101/2020.11.09.373613}
#' @importFrom IRanges IRanges
#' @importFrom parallel detectCores mclapply
#' @importFrom GenomicRanges GRanges tileGenome width seqnames GRangesList
#' sort.GenomicRanges
#' 
raw_counts_to_sparse_matrix <- function(
    files_dir_list, file_type = c("scBED", "scBAM", "FragmentFile"), use_Signac = TRUE,
    peak_file = NULL, n_bins = NULL, bin_width = NULL, genebody = NULL,
    extendPromoter = 2500, verbose = TRUE, ref = c("hg38","mm10")[1],
    progress = NULL, BPPARAM = BiocParallel::bpparam()) {
    warning_raw_counts_to_sparse_matrix(
        files_dir_list, file_type, peak_file, n_bins, bin_width, genebody,
        extendPromoter, verbose, ref)
    
    if (!is.null(progress)) progress$set(detail = "Defining features...",
                                         value = 0.1)
    which <- define_feature(ref, peak_file, bin_width, genebody,
                            extendPromoter)
    
    if(use_Signac && requireNamespace("Signac", quietly=TRUE) &&
       file_type == "FragmentFile"){
        out <- wrapper_Signac_FeatureMatrix(files_dir_list = files_dir_list,
                                            which = which,
                                            ref = ref,
                                            verbose = verbose,
                                            progress = progress)
    } else {
        out <- import_count_input_files(
            files_dir_list, file_type, which, ref, verbose, progress,
            BPPARAM = BPPARAM)
        
        if (!is.null(progress)) progress$set(detail = "Combining matrices...",
                                             value = 0.05)
        mat = list()
        samples_ids = barcodes = c()
        for(i in seq_along(out$feature_indexes)){
            sample_id = names(out$feature_indexes)[i]
            feature_indexes <- out$feature_indexes[[i]]
            name_cells <- paste0(sample_id, "_", out$name_cells[[i]])
            barcodes <- c(barcodes,out$name_cells[[i]])
            samples_ids = c(samples_ids, rep(sample_id, length(name_cells)))
            
            mat[[sample_id]] = Matrix::sparseMatrix(
                i = as.numeric(feature_indexes$feature_index),
                j = feature_indexes$barcode_index,
                x = feature_indexes$counts,
                dims = c(length(which), length(name_cells)),
                dimnames = list(
                    rows = gsub(":|-", "_", as.character(which)),
                    cols = name_cells
                ))
        }
        mat <- do.call("cbind", mat)
        
        # If Gene information is present in metadata, add it to the rownames
        is_gene = ifelse("Gene" %in% colnames(GenomicRanges::elementMetadata(which)),
                         TRUE, FALSE)
        if(is_gene) rownames(mat) = paste0(which$Gene,":",rownames(mat))
        
        eval(parse(text = paste0("data(", ref, ".chromosomes)")))
        chr <- eval(parse(text = paste0("", ref, ".chromosomes")))
        mat = mat[which(as.character(which@seqnames) %in% chr$chr),]
        
        annot_raw = data.frame(barcode = barcodes,
                               cell_id = colnames(mat),
                               sample_id = samples_ids,
                               batch_id = factor(rep(1, ncol(mat)))
        )
        out = list("datamatrix" = mat, "annot_raw" = annot_raw)
    }
    return(out)
}


#' Warning for raw_counts_to_sparse_matrix
#'
#' @param ref reference genome to use (hg38)
#' @param verbose Verbose (TRUE)
#' @param files_dir_list A named character vector of directory containing
#'  the raw files
#' @param peak_file A file containing genomic location of peaks (NULL)
#' @param n_bins The number of bins to tile the genome (NULL)
#' @param bin_width The size of bins to tile the genome (NULL)
#' @param genebody Count on genes (body + promoter) ? (NULL)
#' @param extendPromoter If counting on genes, number of base pairs to extend up or
#'  downstream of TSS (2500).
#' @param file_type Input file(s) type(s) ('scBED','scBAM','SparseMatrix')
#'
#' @return Error or warnings if the input are not correct
warning_raw_counts_to_sparse_matrix <- function(
    files_dir_list, file_type = c("scBAM", "scBED", "SparseMatrix"),
    peak_file = NULL, n_bins = NULL, bin_width = NULL, genebody = NULL,
    extendPromoter = 2500, verbose = TRUE, ref = "hg38"){
        stopifnot(dir.exists(files_dir_list), is.numeric(extendPromoter),
                ref %in% c("mm10", "hg38"))
        
        if (!is.null(peak_file) && !file.exists(peak_file))
            stop("ChromSCape::raw_counts_to_sparse_matrix - 
            Can't find peak file.")
        
        if (!is.null(n_bins) && !is.numeric(n_bins))
            stop("ChromSCape::raw_counts_to_sparse_matrix - 
            n_bins must be a number.")
        
        if (!is.null(bin_width) && !is.numeric(bin_width))
            stop("ChromSCape::raw_counts_to_sparse_matrix - 
            bin_width must be a number.")
        
        if (!is.null(genebody) && !is.logical(genebody))
            stop(paste0(
                "ChromSCape::raw_counts_to_sparse_matrix - genebody ",
                "must be a TRUE or FALSE"))
    }

#' Define the features on which reads will be counted
#' 
#' @usage define_feature(ref = c("hg38","mm10")[1],
#'  peak_file = NULL,
#'  bin_width  = NULL,
#'  genebody = FALSE,
#'  extendPromoter = 2500)
#' 
#' @param ref Reference genome
#' @param peak_file A bed file if counting on peaks
#' @param bin_width A number of bins if divinding genome into fixed width bins
#' @param genebody A logical indicating if feature should be counted in 
#' genebodies and promoter.
#' @param extendPromoter Extension length before TSS (2500).
#'
#' @return A GRanges object
#'   
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges duplicated tileGenome width seqnames
#' @importFrom IRanges subsetByOverlaps
#' 
#' @export
#' @examples
#' gr_bins = define_feature("hg38", bin_width = 50000)
#' gr_genes = define_feature("hg38", genebody = TRUE, extendPromoter = 5000)
#' 
define_feature <- function(ref = c("hg38","mm10")[1], peak_file = NULL,
                           bin_width  = NULL, genebody = FALSE,
                           extendPromoter = 2500){
    eval(parse(text = paste0("data(", ref, ".chromosomes)")))
    chr <- eval(parse(text = paste0("", ref, ".chromosomes")))
    chr <- GenomicRanges::GRanges(chr)
    if (!is.null(peak_file)) {
        message('ChromSCape::define_feature - ',
                'Reading in peaks file...')
        features <- rtracklayer::import(peak_file, format="bed")
        which <- IRanges::subsetByOverlaps(features,chr)
        if(length(which(GenomicRanges::duplicated(which)))>0) 
            message("ChromSCape::define_feature - Removing ",
                    length(which(GenomicRanges::duplicated(which))), 
                    " duplicated regions from peak file.")
        which = which[!GenomicRanges::duplicated(which)]
    } else if (!is.null(bin_width)) {
        message('ChromSCape::define_feature - ',
                'Counting on genomic bins...')
            which <- unlist(GenomicRanges::tileGenome(
                setNames(
                    GenomicRanges::width(chr),
                    GenomicRanges::seqnames(chr)
                ),
                tilewidth = bin_width
            ))
    } else if (genebody == TRUE) {
        # Retrieve gene TSS from ref and create GRanges
        eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
        geneTSS_df <- eval(parse(text = paste0("", ref, ".GeneTSS")))

        geneTSS_df$start = ifelse(geneTSS_df$strand == "+",
                                  geneTSS_df$start - extendPromoter,
                                  geneTSS_df$start)
        
        geneTSS_df$end = ifelse(geneTSS_df$strand == "-",
                                  geneTSS_df$end + extendPromoter,
                                  geneTSS_df$end)
    
        geneTSS_df$strand = NULL
        geneTSS_df$distanceToTSS = 0
        which = GenomicRanges::GRanges(geneTSS_df)
        which = which[!GenomicRanges::duplicated(which)]
        names(which@ranges) = NULL
    }
    return(which)
}

#' Wrapper around 'FeatureMatrix' function from Signac Package
#' 
#' @param files_dir_list A named character vector of directories containing
#'  the files. The names correspond to sample names.
#' @param which A GenomicRanges containing the features to count on.
#' @param ref Reference genome to use (hg38).Chromosomes that are not present in 
#' the canonical chromosomes of the given reference genome will be excluded from
#'  the matrix.
#' @param process_n Number of regions to load into memory at a time, per thread.
#'  Processing more regions at once can be faster but uses more memory. (2000)
#' @param set_future_plan Set 'multisession' plan within the function (TRUE). 
#' If TRUE, the previous plan (e.g. future::plan()) will be set back on exit.
#' @param verbose Verbose (TRUE).
#' @param progress Progress object for Shiny.
#'
#' @return A sparse matrix of features x cells
#' 
#' @details Signac & future are not required packages for ChromSCape as they are
#' required only for the fragment matrix calculations. To use this function,
#' install Signac package first (future will be installed as a dependency).  
#' For the simplicity of the application & optimization, the function
#' by defaults sets future::plan("multisession") with workers = 
#' future::availableCores() - 1 in order to allow parallel processing
#' with Signac. On exit the plan is re-set to the previously set future plan.
#' Note that future multisession may have trouble running when VPN is on. To 
#' run in parallel, first deactivate your VPN if you encounter long runtimes.
#' 
#' @references Stuart el al.,  Multimodal single-cell chromatin analysis with 
#' Signac bioRxiv \url{https://doi.org/10.1101/2020.11.09.373613}
#' 
#' @export
#' @examples 
#' \dontrun{
#' gr_bins = define_feature("hg38", bin_width = 50000)
#' wrapper_Signac_FeatureMatrix("/path/to/dir_containing_fragment_files",
#'  gr_bins, ref = "hg38")
#' }
wrapper_Signac_FeatureMatrix <- function(files_dir_list, which, ref = "hg38",
                                         process_n = 2000,
                                         set_future_plan = TRUE, verbose = TRUE, 
                                         progress = NULL){
    
    fragment_files = lapply(files_dir_list, function(dir) 
        list.files(dir, full.names = TRUE, recursive = FALSE,
                   include.dirs = FALSE, pattern = ".tsv$|.tsv.gz$")
    )
    
    # Number of workers equal to number of workers in BPPARAM
    if(set_future_plan){
        cl <- future::availableCores() - 1
        oplan <- future::plan("multisession", workers = cl,  gc = TRUE)
    }
    
    message("ChromSCape::wrapper_Signac_FeatureMatrix - Running ",
            "Signac::FeatureMatrix with ", cl,
            " workers...")
    
    if(any(lapply(fragment_files, length) > 1)){
        message("ChromSCape::wrapper_Signac_FeatureMatrix - Multiple files ",
                "detected per folder, assigning sample name based on files",
                " names.")
        fragment_files = unlist(fragment_files)
        names(fragment_files) = gsub(".tsv|.gz","",basename(fragment_files))
    } else{
        names(fragment_files) = basename(files_dir_list)
    }
    
    mat_list = list()
    samples_ids = barcodes = c()
    for(sample_name in names(fragment_files)){
        message("ChromSCape::wrapper_Signac_FeatureMatrix - Running Signac for",
                " sample ", sample_name, "...")
        
        if(!is.null(progress)) progress$inc(
            amount = 0.6/(length(fragment_files)),
            detail = paste0("Counting sample - ", sample_name))
        
        fragment_file = fragment_files[[sample_name]]
        dir = dirname(fragment_file)
        # if needed index file using samtools tabix
        if(!file.exists(file.path(dir, paste0(basename(fragment_file),".tbi")))){
            message("ChromSCape::wrapper_Signac_FeatureMatrix - creating tabix ",
                    "index...")
            idx = Rsamtools::indexTabix(fragment_file, format="bed")
        }
        # create Signac Fragment 
        fragments = Signac::CreateFragmentObject(fragment_file)
       
        
        # create count matrix
        system.time({
            mat = Signac::FeatureMatrix(fragments = fragments,
                                        features = which,
                                        process_n = process_n,
                                        verbose = TRUE)
        })
        rownames(mat) = gsub("-","_", rownames(mat))
        sample_id = sample_name
        barcodes <- c(barcodes, colnames(mat))
        colnames(mat) = paste0(sample_id,"_", colnames(mat))
        samples_ids = c(samples_ids, rep(sample_id, ncol(mat)))
        mat_list[[sample_id]] = mat
        
        gc()
    }
    if(set_future_plan) future::plan(oplan)
    
    if (!is.null(progress)) progress$set(detail = "Combining matrices...",
                                         value = 0.05)
    mat <- do.call("cbind", mat_list)
    
    # If Gene information is present in metadata, add it to the rownames
    is_gene = ifelse(
        "Gene" %in% colnames(GenomicRanges::elementMetadata(which)),TRUE, FALSE)
    if(is_gene) rownames(mat) = paste0(which$Gene,":",rownames(mat))
    
    
    eval(parse(text = paste0("data(", ref, ".chromosomes)")))
    chr <- eval(parse(text = paste0("", ref, ".chromosomes")))

    mat = mat[which(as.character(which@seqnames) %in% chr$chr),]
    
    annot_raw = data.frame(barcode = barcodes,
                           cell_id = colnames(mat),
                           sample_id = samples_ids,
                           batch_id = factor(rep(1, ncol(mat)))
                           
    )
    
    out <- list("datamatrix" = mat,
                "annot_raw" = annot_raw)
    
    return(out)
}

#' Import and count input files depending on their format
#'
#' @param files_dir_list A named list of directories containing the input files.
#' @param which A GRanges object of features. 
#' @param ref Reference genome.
#' @param verbose Print ?
#' @param file_type Input file type.
#' @param progress A progress object for Shiny.
#' @param BPPARAM BPPARAM object for multiprocessing. See
#'  \link[BiocParallel]{bpparam} for more informations. Will take the default
#'  BPPARAM set in your R session.
#'  
#' @return A list with the 
#' feature indexes data.frame containing non-zeroes entries in the count matrix
#' and the cell names
#'
import_count_input_files <- function(files_dir_list, file_type,
                                    which, ref, verbose, progress,
                                    BPPARAM = BiocParallel::bpparam()){

    feature_indexes = list()
    name_cells = list()
    for(i in seq_along(files_dir_list)){
        dir = files_dir_list[[i]]
        sample_id = names(files_dir_list)[i]
        if(!is.null(progress)) progress$inc(
            detail=paste0("Computing ", sample_id, "...",
            "\nLook for progress in the console."),
            amount = 0.6/length(files_dir_list))
        message("ChromSCape::import_count_input_files - Computing ",
                sample_id, "...")
        if (file_type == "scBAM") {
            t1 = system.time({
                l = bams_to_matrix_indexes(dir, which, BPPARAM = BPPARAM)
            })
            if (verbose)
                message("ChromSCape::raw_counts_to_sparse_matrix - ",
                        "Count matrix ", sample_id ,"created from BAM files in "
                        , t1[3], " sec.\n")
        }
        else if (file_type == "scBED") {
            t1 = system.time({
                l = beds_to_matrix_indexes(dir, which, BPPARAM = BPPARAM)
            })
            if (verbose)
                message("ChromSCape::raw_counts_to_sparse_matrix - ",
                        "Count matrix ", sample_id ,"
                        created from BED files in ", t1[3], " sec.")}
        
        feature_indexes[[i]] = l[[1]]
        names(feature_indexes)[i] = sample_id
        name_cells[[i]] = l[[2]]
        names(name_cells)[i] = sample_id
    }
    
    out <- list("feature_indexes" = feature_indexes,
                "name_cells" = name_cells)

    return(out)
}

#' Count bam files on interval to create count indexes
#'
#' @param dir A directory containing single cell  BAM files and BAI files
#' @param which Genomic Range on which to count
#' @param BPPARAM BPPARAM object for multiprocessing. See
#'  \link[BiocParallel]{bpparam} for more informations. Will take the default
#'  BPPARAM set in your R session.
#'  
#' @return A list containing a "feature index" data.frame and a 
#' count vector for non 0 entries, both used to form the sparse matrix
#'
#' @importFrom Rsamtools BamFileList indexBam ScanBamParam countBam
#' @importFrom BiocParallel bplapply
bams_to_matrix_indexes = function(dir, which, 
                                  BPPARAM = BiocParallel::bpparam()) {
    single_cell_bams = list.files(dir,
                                full.names = TRUE,
                                pattern = paste0(".*.bam$"))
    bam_files = Rsamtools::BamFileList(single_cell_bams)
    name_cells = gsub("*.bam$", "", basename(as.character(single_cell_bams)))
    names(bam_files) = name_cells
    
    indexes = list.files(dir,
                        full.names = TRUE,
                        pattern = paste0(".*.bai$"))
    if (length(indexes) < length(bam_files)) {
        message("ChromSCape::bams_to_matrix_indexes - indexing BAM files...")
        BiocParallel::bplapply(bam_files, Rsamtools::indexBam)
        indexes = list.files(files_dir_list,
                            full.names = TRUE,
                            pattern = paste0(".*.bai$"))
    }
    names_index = gsub("*.bam.bai$", "", basename(as.character(indexes)))
    names(indexes) = names_index
    if (!all.equal(names_index, name_cells))
        stop("Different number of BAM and indexes files. Stopping.")
    param = Rsamtools::ScanBamParam(which = which)
    
    BiocParallel::bpprogressbar(BPPARAM) <- TRUE
    if(BiocParallel::bpworkers(BPPARAM) > 1) 
        BiocParallel::bptasks(BPPARAM) <- ceiling(length(single_cell_bams) /
                                    (10*BiocParallel::bpworkers(BPPARAM)))
    system.time({
        feature_list = BiocParallel::bplapply(BPPARAM = BPPARAM,
            names(bam_files), function(bam_name, bam_files, param) {
                bam_files[[bam_name]]$index = indexes[[bam_name]]
                tmp = Rsamtools::countBam(
                    file = bam_files[[bam_name]], param = param)
                tmp$feature_index = rownames(tmp)
                tmp$cell_id = bam_name
                sel = which(tmp$records > 0)
                if (length(sel) > 0)
                    tmp = tmp[sel, c("cell_id", "feature_index", "records")]
                tmp
            }, bam_files = bam_files, param = param)
    })
    
    feature_indexes = do.call(rbind, feature_list)
    feature_indexes$barcode_index = 
        as.numeric(as.factor(feature_indexes$cell_id))
    colnames(feature_indexes)[3] = "counts"
    return(list(feature_indexes, name_cells))
}

#' Read index-peaks-barcodes trio files on interval to create count indexes
#'
#'
#' @param feature_file A file containing the features genomic locations
#' @param matrix_file A file containing the indexes of non-zeroes values and
#'   their value (respectively i,j,x,see sparseMatrix)
#' @param barcode_file A file containing the barcode ids
#' @param binarize Binarize matrix ?
#'
#' @importFrom GenomicRanges GRanges
#'
#' @return A list containing a "feature index" data.frame, name_cells, and a
#'  region GenomicRange object used to form the sparse matrix
#'   
index_peaks_barcodes_to_matrix_indexes = function(
    feature_file, matrix_file, barcode_file, binarize = FALSE) {
    regions = GenomicRanges::GRanges(setNames(
        read.table(feature_file, sep = "\t", quote = "")[, seq_len(3)],
        c("chr", "start", "end")
    ))
    
    if(any(GenomicRanges::duplicated(regions))){
        stop("ChromSCape::index_peaks_barcodes_to_matrix_indexes - The ",
             "features are not unique, check your file : ", feature_file)
    }
    
    name_cells = read.table(barcode_file, sep = "", quote = "",
                            stringsAsFactors = FALSE)[, 1]
    if(length(unique(name_cells)) != length(name_cells)){
        stop("ChromSCape::index_peaks_barcodes_to_matrix_indexes - The ",
             "barcode IDs are not unique, check your file : ",barcode_file)
    }
    feature_indexes = setNames(
        Matrix::readMM(matrix_file),
        c("feature_index", "barcode_index", "counts")
    )
    
    feature_indexes$cell_id = name_cells[feature_indexes$barcode_index]
    feature_indexes = feature_indexes[, c("cell_id", "feature_index",
                                        "counts", "barcode_index")]
    return(list(feature_indexes, name_cells, regions))
}

#' Count bed files on interval to create count indexes
#'
#' @param dir A directory containing the single cell BED files
#' @param which Genomic Range on which to count
#' @param BPPARAM BPPARAM object for multiprocessing. See
#'  \link[BiocParallel]{bpparam} for more informations. Will take the default
#'  BPPARAM set in your R session.
#'
#' @importFrom BiocParallel bplapply
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges countOverlaps
#' 
#' @return A list containing a "feature index" data.frame and a 
#' names of cells as vector both used to form the sparse matrix
#' 
beds_to_matrix_indexes <- function(dir, which,
                                   BPPARAM = BiocParallel::bpparam()) {
    single_cell_beds = list.files(
        dir, full.names = TRUE, pattern = ".*.bed$|.*.bed.gz$")
    names_cells = gsub(".bed$|.bed.gz$", 
                    "", basename(as.character(single_cell_beds)))
    names(single_cell_beds) = names_cells

    BiocParallel::bpprogressbar(BPPARAM) <- TRUE
    # BiocParallel::bptasks(BPPARAM) <- ceiling(length(single_cell_beds) /
    #                                 (10*BiocParallel::bpworkers(BPPARAM)))

    system.time({
    feature_list = BiocParallel::bplapply( BPPARAM = BPPARAM,
            names(single_cell_beds),
            function(bed_name, single_cell_beds, which) {
                bed = rtracklayer::import(single_cell_beds[[bed_name]],
                                          format = "bed")
                suppressWarnings({
                    tmp = data.frame(
                        counts = GenomicRanges::countOverlaps(
                            which, bed, minoverlap = 1))
                })
                tmp$feature_index = as.numeric(rownames(tmp))
                tmp = tmp[which(tmp$counts > 0), ]
                if(nrow(tmp)>0) {
                    tmp$cell_id = bed_name 
                    tmp = tmp[, c("cell_id", "feature_index", "counts")]
                } else { tmp = NULL}
                tmp
            }, single_cell_beds = single_cell_beds,
            which = which)
    })
    
    gc()
    message("ChromSCape::import_count_input_files - Merging single-cells into ",
            "indexes...")
    feature_indexes = do.call(rbind, feature_list)
    feature_indexes$barcode_index = as.numeric(
        as.factor(feature_indexes$cell_id))
    out = list(feature_indexes, names_cells)
    gc()
    return(out)
}

#' Transforms a peaks x cells count matrix into a bins x cells count matrix.
#'
#' This functions is best used to re-count large number of small peaks (e.g. <=
#' 5000bp) into equal or larger bins. The genome is either cut in fixed bins
#' (e.g. 50,000bp) or into an user defined number of bins. Bins are calculated
#' based on the canconical chromosomes. Note that if peaks are larger than bins,
#' or if peaks are overlapping multiple bins, the signal is added to each bin.
#' Users can increase the minimum overlap to consider peaks overlapping bins (by
#' default 150bp, size of a nucleosome) to disminish the number of peaks
#' overlapping multiple region. Any peak smaller than the minimum overlapp
#' threshold will be dismissed. Therefore, library size might be slightly
#' different from peaks to bins if signal was duplicated into multiple bins or
#' ommitted due to peaks smaller than minimum overlap.
#'
#' @param mat A matrix of peaks x cells
#' @param bin_width width of bins to produce in base pairs (minimum 500) (50000)
#' @param ref reference genome to use (hg38)
#' @param n_bins number of bins (exclusive with bin_width)
#' @param minoverlap Minimum overlap between a peak and a bin to consider the
#'   peak as overlapping the bin (150).
#' @param verbose Verbose
#'
#' @return A sparse matrix of bins instead of peaks
#' @export
#'
#' @importFrom IRanges IRanges
#' @importFrom BiocParallel bpaggregate
#' @importFrom GenomicRanges GRanges tileGenome width seqnames GRangesList
#'   sort.GenomicRanges findOverlaps
#'
#' @examples
#' mat = create_scDataset_raw()$mat
#' binned_mat = peaks_to_bins(mat,bin_width = 10e6)
#' dim(binned_mat)
#' 
peaks_to_bins <- function(mat, bin_width = 50000, n_bins = NULL,
                        minoverlap = 150, verbose = TRUE,ref = "hg38"){
    stopifnot(!is.null(mat), ref %in% c("mm10", "hg38"))
    if (is.matrix(mat)) mat = as(mat, "dgCMatrix")
    if (is.null(n_bins) & is.null(bin_width)) 
        stop("One of bin_width or n_bins must be set")
    eval(parse(text = paste0("data(", ref, ".chromosomes)")))
    chr <- eval(parse(text = paste0("", ref, ".chromosomes")))
    
    chr <- GenomicRanges::GRanges(chr)
    if (!is.null(n_bins)) {
        bin_ranges <- unlist(GenomicRanges::tileGenome(
            setNames(GenomicRanges::width(chr), GenomicRanges::seqnames(chr)),
            ntile = n_bins))} else {
        bin_ranges <- unlist(GenomicRanges::tileGenome(
            setNames( GenomicRanges::width(chr), GenomicRanges::seqnames(chr)),
            tilewidth = bin_width )) }
    peaks = rownames(mat)
    peaks_chr = as.character(lapply(strsplit(peaks, ":|-|_"), function(x) x[1]))
    peaks_start = as.numeric(lapply(strsplit(peaks, ":|-|_"), function(x) x[2]))
    peaks_end = as.numeric(lapply(strsplit(peaks, ":|-|_"), function(x) x[3]))
    if (anyNA(peaks_start) | anyNA(peaks_end))
        stop("The rows of mat should be regions in format chr_start_end or",
            "chr:start-end, without non canonical chromosomes")
    if (verbose) {
        message("ChromSCape::peaks_to_bins - converting ", dim(mat)[1],
                " peaks into ", length(bin_ranges)[1], " bins of ",
                mean(bin_ranges@ranges@width), " bp in average.") }
    peaks = GenomicRanges::GRanges(
        seqnames = peaks_chr, ranges = IRanges::IRanges(peaks_start, peaks_end))
    hits <- GenomicRanges::findOverlaps(
        bin_ranges, peaks, minoverlap = minoverlap)
    bins_names = paste0(bin_ranges@seqnames, ":", 
                        GenomicRanges::start(bin_ranges), "-",
                        GenomicRanges::end(bin_ranges))
    bin_mat = NULL
    hits = as.matrix(hits)
    hits = cbind(hits, mat[hits[, "subjectHits"], ])
    print("Running aggregation of peaks to bins in parallel")
    hits = as.matrix(hits)
    t = system.time({ bin_mat = BiocParallel::bpaggregate(
        x = hits[, 3:dim(hits)[2]],
        by = list(bins = hits[, 1, drop = FALSE]), FUN = sum) })
    message("ChromSCape::peaks_to_bins - From peaks to bins in ",t[3]," sec.")
    bin_mat = as(as.matrix(bin_mat[, 2:ncol(bin_mat)]), "dgCMatrix")
    rownames(bin_mat) = bins_names[unique(hits[, 1])]
    gc()
    if (verbose) message("ChromSCape::peaks_to_bins - removed ",
                    length(bin_ranges) - nrow(bin_mat),
                    " empty bins from the binned matrix.")
    return(bin_mat)}

#' Create a simulated single cell datamatrix & cell annotation
#'
#' @param cells Number of cells (300)
#' @param features Number of features (600)
#' @param featureType Type of feature (window)
#' @param sparse Is matrix sparse ? (TRUE)
#' @param nsamp Number of samples (4)
#' @param ref Reference genome ('hg38')
#' @param batch_id Batch origin (factor((1,1,1,1))
#'
#' @return A list composed of * mat : a sparse matrix following an approximation
#'   of the negative binomial law (adapted to scChIPseq) * annot : a data.frame
#'   of cell annotation * batches : an integer vector with the batch number for
#'   each cell
#' @export
#'
#' @importFrom GenomicRanges GRanges tileGenome width seqnames GRangesList
#'   sort.GenomicRanges
#'
#' @examples
#' # Creating a basic sparse 600 genomic bins x 300 cells matrix and annotation
#' l = create_scDataset_raw()
#' head(l$mat)
#' head(l$annot)
#' head(l$batches)
#'
#' # Specifying number of cells, features and samples
#' l2 = create_scDataset_raw(cells = 500, features = 500, nsamp=2)
#'
#' # Specifying species
#' mouse_l = create_scDataset_raw(ref="mm10")
#' 
#' # Specifying batches
#' batch_l = create_scDataset_raw(nsamp=4, batch_id = factor(c(1,1,2,2)))
#' 
#' # Peaks of different size as features
#' peak_l = create_scDataset_raw(featureType="peak")
#' head(peak_l$mat)
#' 
#' # Genes as features
#' gene_l = create_scDataset_raw(featureType="gene")
#' head(gene_l$mat)
create_scDataset_raw <- function(
    cells = 300, features = 600, featureType = c("window", "peak", "gene"),
    sparse = TRUE, nsamp = 4, ref = "hg38", batch_id = factor(rep(1, nsamp))){
    stopifnot( featureType %in% c("window", "peak", "gene"),
            ref %in% c("mm10", "hg38"), nsamp >= 1, cells >= nsamp,
            features >= 1, length(batch_id) == nsamp)
    stopifnot(is.factor(batch_id))
    # Create cell names
    cell_counts <-
        as.numeric(lapply(split(
            seq_len(cells), sample(nsamp, cells, repl = TRUE)), length))
    cell_names <- sample <- batches <- list()
    for (i in seq_along(cell_counts))
    {
        cell_names[[i]] <- paste0("sample_", i, "_c", seq_len(cell_counts[i]))
        sample[[i]] <- rep(paste0("sample_", i), cell_counts[i])
        batches[[i]] <- rep(batch_id[i], cell_counts[i])
    }
    cell_names <- as.character(unlist(cell_names))
    sample <- as.character(unlist(sample))
    batches <- factor(as.numeric(unlist(batches)))
    
    feature_names <- generate_feature_names(featureType, ref, features)
    mat <- generate_count_matrix(cells, features, sparse,
                                cell_names, feature_names)
    if (length(levels(batches)) > 1)
    {
        mat <- mat %*% as(Matrix::diag(as.numeric(batches) * 
                                            as.numeric(batches)), "dgCMatrix")
    }
    colnames(mat) <- cell_names
    annot <- data.frame(
        barcode = paste0("BC", seq_len(cells)),
        cell_id = cell_names,
        sample_id = sample,
        batch_id = batches,
        total_counts = Matrix::colSums(mat)
    )
    return(list(
        mat = mat,
        annot = annot,
        batches = batches
    ))
}


#' Generate feature names
#'
#' @param featureType Type of feature
#' @param ref Reference genome
#' @param features Number of features to generate
#'
#' @return A character vector of feature names
#'
generate_feature_names <- function(featureType, ref,
                                features){
    
    eval(parse(text = paste0("data(", ref, ".chromosomes)")))
    chr <- eval(parse(text = paste0("", ref, ".chromosomes")))
    chr <- GenomicRanges::GRanges(chr)
    
    if (featureType[1] == "window"){
        chr_ranges <- unlist(GenomicRanges::tileGenome(
            setNames(GenomicRanges::width(chr), GenomicRanges::seqnames(chr)),
            ntile = features
        ))[seq_len(features)]
        
        feature_names <- paste(as.data.frame(chr_ranges)$seqnames,
                                as.data.frame(chr_ranges)$start,
                                as.data.frame(chr_ranges)$end, sep = "_")
    }
    if (featureType[1] == "peak"){
        size_peaks <- c(1000, 2500, 7999, 10000, 150000, 10 ^ 6)
        
        peaks <- as.numeric(lapply(split(seq_len(features), sample(
            length(size_peaks), features, replace = TRUE)), length))
        chr_ranges_list <- GenomicRanges::GRangesList()
        for (i in seq_along(peaks)){
            chr_ranges <- unlist(
                GenomicRanges::tileGenome(
                    setNames(GenomicRanges::width(chr),
                            GenomicRanges::seqnames(chr)),
                    tilewidth = size_peaks[i],
                    cut.last.tile.in.chrom = FALSE))
            chr_ranges_list[[i]] <-
                chr_ranges[sample(seq_along(chr_ranges),
                                size = peaks[i]),]
        }
        chr_ranges <-
            GenomicRanges::sort.GenomicRanges(
                unlist(chr_ranges_list))[seq_len(features)]
        
        feature_names <- paste(as.data.frame(chr_ranges)$seqnames,
                                as.data.frame(chr_ranges)$start,
                                as.data.frame(chr_ranges)$end,
                                sep = "_")
    }
    if (featureType[1] == "gene"){
        eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
        chr <- eval(parse(text = paste0("", ref, ".GeneTSS")))
        feature_names <- as.character(sample(
            chr$Gene, features, replace = FALSE))
    }
    return(feature_names)
}

#' Generate count matrix
#'
#' @param cells Number of cells
#' @param features Number of features
#' @param sparse Is matrix sparse ?
#' @param cell_names Cell names
#' @param feature_names Feature names
#'
#' @return A matrix or a sparse matrix
#'
generate_count_matrix <- function(cells, features, sparse,
                                cell_names, feature_names){
    vec <- stats::rpois(cells * features, 0.5)
    # Add count to values > 0, iteratively
    for (i in seq_len(10))
    {
        vec[vec >= i] <- vec[vec >= i] + i ^ 2 *
            stats::rpois(length(vec[vec >= i]),
                        0.5)
    }
    
    indices_vec <- which(vec > 0)
    j <- ceiling(indices_vec / features)
    i <- ceiling(indices_vec %% (features))
    i[i == 0] <- features
    if (sparse)
    {
        mat <- Matrix::sparseMatrix(i,
                                    j,
                                    x = vec[indices_vec],
                                    dimnames = list(feature_names,
                                                    cell_names))
    } else
    {
        mat <- matrix(
            vec,
            nrow = features,
            ncol = cells,
            dimnames = list(feature_names,
                            cell_names)
        )
    }
    return(mat)
}


#' Read single-cell matrix(ces) into scExp
#'
#' Combine one or multiple matrices together to create a sparse matrix and cell
#' annotation data.frame.
#'
#' @param file_paths A character vector of file names towards single cell
#' epigenomic matrices (features x cells) (must be .txt / .tsv)
#' @param temp_path In case matrices are stored in temporary folder,
#' a character vector of path towards temporary files. (NULL)
#' @param remove_pattern A string pattern to remove from the sample names. Can
#' be a regexp.
#' 
#' @return A list containing:  
#' * datamatrix: a sparseMatrix of features x cells  
#' * annot_raw: an annotation of cells as data.frame  
#'  
#' @export
#' @examples 
#' mat1 = mat2 = create_scDataset_raw()$mat
#' tmp1 = tempfile(fileext = ".tsv")
#' tmp2 = tempfile(fileext = ".tsv")
#' write.table(as.matrix(mat1),file=tmp1,sep = "\t",
#' row.names = TRUE,col.names = TRUE,quote = FALSE)
#' write.table(as.matrix(mat2),file=tmp2, sep = "\t",
#' row.names = TRUE,col.names = TRUE,quote = FALSE)
#' file_paths = c(tmp1,tmp2)
#' out = import_scExp(file_paths)
#' 
#' @importFrom scater readSparseCounts
#' @md
import_scExp <- function(file_paths,
                        remove_pattern = "",
                        temp_path = NULL) {
    stopifnot(is.character(file_paths))
    if (length(grep("(.tsv$)|(.txt$)|(.csv$)|(.gz$)", file_paths)) 
        < length(file_paths))
        stop(paste0("ChromSCape::import_scExp - Matrix files must be in",
                    " .txt, .csv or .tsv or .gz format."))
    if (is.null(temp_path)) temp_path = file_paths
    if (FALSE %in% as.logical(lapply(temp_path, file.exists))) 
        stop("ChromSCape::import_scExp - can't find one of the matrix files.")
    datamatrix = annot_raw = NULL
    for (i in seq_along(file_paths)) {
        sample_name <- gsub('(.tsv$)|(.txt$)|(.csv$)', "", gsub('(.gz$)', '',
                                               basename(file_paths[i])))
        if(remove_pattern != "") sample_name = gsub(remove_pattern, "",
                                                    sample_name)
        
        separator <- separator_count_mat(temp_path[i])
        format_test = read.table(temp_path[i], header = TRUE,
                                sep = separator, nrows = 5)
        separated_chr_start_end = c(
            grep("chr", colnames(format_test)[seq_len(3)]),
            grep("start|begin", colnames(format_test)[seq_len(3)]),
            grep("end|stop", colnames(format_test)[seq_len(3)]))
        if (length(separated_chr_start_end) > 0 &&
            all.equal(separated_chr_start_end, c(1, 2, 3))) {
            datamatrix_single = read_count_mat_with_separated_chr_start_end(
                temp_path[i], format_test, separator)
        } else{
            datamatrix_single <- scater::readSparseCounts(
                temp_path[i], sep = separator, chunk = 1000L)
        }
        gc()
        datamatrix_single <- check_correct_datamatrix(datamatrix_single,
                                                      sample_name)
        gc()
        total_cell <- length(datamatrix_single[1, ])
        annot_single <- data.frame(
            barcode = colnames(datamatrix_single),
            cell_id = paste0(sample_name, "_", colnames(datamatrix_single)),
            sample_id = rep(sample_name, total_cell), batch_id = 1)
        colnames(datamatrix_single) <- annot_single$cell_id
        datamatrix <- combine_datamatrix(datamatrix,
                                        datamatrix_single, file_paths, i)
        rm(datamatrix_single)
        if (is.null(annot_raw)) annot_raw <- annot_single
        else annot_raw <- rbind(annot_raw, annot_single)
        rm(annot_single)
        gc()
    }
    annot_raw$batch_id <- as.factor(annot_raw$batch_id)
    out = list("datamatrix" = datamatrix, "annot_raw" = annot_raw)
    return(out)
}



#' Determine Count matrix separator ("tab" or ",")
#'
#' @param path_to_matrix A path towards the count matrix to check
#'
#' @return A character separator
#'
separator_count_mat <- function(path_to_matrix){
    format_test = as.character(
        read.table(path_to_matrix, header = TRUE, sep = "\t", nrows = 5)[4, ])
    
    if (length(format_test) > 3)
        separator = "\t"
    else
        separator = ","
    
    return(separator)
}

#' Read a count matrix with three first columns (chr,start,end)
#'
#' @param path_to_matrix Path to the count matrix
#' @param format_test Sample of the read.table
#' @param separator Separator character
#'
#' @return A sparseMatrix with rownames in the form "chr1:1222-55555"
#'
read_count_mat_with_separated_chr_start_end <- function(
    path_to_matrix,
    format_test,
    separator
){
    val = c("NULL", "NULL", "NULL",
            rep("integer", length(colnames(
                format_test
            )) - 3))
    
    datamatrix_single = as(as.matrix(
        read.table(
            path_to_matrix,
            header = TRUE,
            sep = separator,
            colClasses = val
        )
    ), "dgCMatrix")
    gc()
    
    val = c("character",
            "integer",
            "integer",
            rep("NULL", ncol(datamatrix_single)))
    regions = read.table(
        path_to_matrix,
        header = TRUE,
        sep = separator,
        colClasses = val
    )
    regions = paste0(regions[, 1], ":", regions[, 2], "-", regions[, 3])
    rownames(datamatrix_single) = regions
    
    return(datamatrix_single)
}


#' Check if matrix rownames are well formated and correct if needed
#'
#' Throws warnings / error if matrix is in the wrong format
#' 
#' @param datamatrix_single A sparse matrix
#' @param sample_name Matrix sample name for warnings
#' 
#' 
#' @return A sparseMatrix in the right rownames format 
#'
check_correct_datamatrix <- function(datamatrix_single,sample_name=""){
    matchingRN <-
        grep("[[:alnum:]]+(:|_|-)[[:digit:]]+(-|_)[[:digit:]]+",
            rownames(datamatrix_single)) # check rowname format
    if (length(matchingRN) < length(rownames(datamatrix_single))) {
        warning(paste0(
            "ChromSCape::import_scExp - ", sample_name, " contains ",(length(
                rownames(datamatrix_single)) - length(matchingRN)),
            " rownames that do not conform to the required format.
                    Please check your data matrix and try again."))
        if (length(matchingRN) < 5) {
            stop(
                "ChromSCape::import_scExp - Maybe your rownames are 
                contained in the first column instead? In this case, remove the 
                header of this column so that they are interpreted as rownames."
            )
        }
        return()
    }
    numericC <-
        apply(datamatrix_single[seq_len(5), seq_len(5)], MARGIN = 2, is.numeric)
    if (sum(numericC) < 5) {
        stop(
            paste0(
                "ChromSCape::import_scExp - ", sample_name,
                " contains non-numeric columns at the following indices: ",
                which(numericC == FALSE), 
                ". Please check your data matrix and try again."))
    }
    if (all.equal(rownames(datamatrix_single)[seq_len(5)],
                c("1", "2", "3", "4", "5")) == TRUE) {
        warning(
            paste0("ChromSCape::import_scExp - ", sample_name, " contains ",
                "numeric rownames, taking first column as rownames."))
        names = datamatrix_single$X0
        datamatrix_single = datamatrix_single[, -1]
        rownames(datamatrix_single) = names
    }
    datamatrix_single <-
        datamatrix_single[!duplicated(rownames(datamatrix_single)), ]
    if (length(grep("chr", rownames(datamatrix_single)[seq_len(10)], 
                    perl = TRUE)) >= 9) {
        rownames(datamatrix_single) <-
            gsub(":|-", "_", rownames(datamatrix_single))
        rownames(datamatrix_single) <-
            gsub("-", "_", rownames(datamatrix_single))
    }
    return(datamatrix_single)
}

#' Combine two matrices and emit warning if no regions are in common
#'
#' @param datamatrix A sparse matrix or NULL if empty
#' @param datamatrix_single Another sparse matrix
#' @param file_names File name corresponding to the matrix for warnings
#' @param i file number
#'
#' @return A combined sparse matrix
#'
combine_datamatrix <- function(datamatrix, datamatrix_single,
                            file_names, i){
    if (is.null(datamatrix)) {datamatrix <- datamatrix_single
    } else {
        
        common_regions <- intersect(rownames(datamatrix),
                                    rownames(datamatrix_single))
        if (length(common_regions) > 0) {
            datamatrix <-
                Matrix::cbind2(datamatrix[common_regions, ],
                            datamatrix_single[common_regions, ])
        } else {
            stop(paste0(
                "ChromSCape::import_scExp - ", file_names[i],
                " contains no common regions with ",file_names[i - 1]))
        }
        if (length(common_regions) < nrow(datamatrix)) {
            warning(paste0(
                "ChromSCape::import_scExp - ", file_names[i],
                " contains less than ", ceiling(
                    100 * length(common_regions) / (nrow(datamatrix))),
                " common regions with ", file_names[i - 1]))
        }
    }
    return(datamatrix)
}

#' Wrapper to create the single cell experiment from count matrix and feature
#' dataframe
#'
#' Create the single cell experiment from (sparse) datamatrix and feature
#' dataframe containing feature names and location. Also optionally removes zero
#' count Features, zero count Cells, non canconical chromosomes, and chromosome
#' M. Calculates QC Metrics (scran).
#'
#' @param datamatrix A matrix or sparseMatrix of raw counts. Features x Cells
#'   (rows x columns).
#' @param annot A data.frame containing informations on cells. Should have the
#'   same number of rows as the number of columns in datamatrix.
#' @param remove_zero_cells remove cells with zero counts ? (TRUE)
#' @param remove_zero_features remove cells with zero counts ? (TRUE)
#' @param remove_non_canonical remove non canonical chromosomes ?(TRUE)
#' @param remove_chr_M remove chromosomes M ? (TRUE)
#' @param mainExpName Name of the mainExpName e.g. 'bins', 'peaks'... 
#' ("default")
#' @param verbose (TRUE)
#'
#' @return Returns a SingleCellExperiment object.
#' @export
#'
#'
#' @importFrom SingleCellExperiment SingleCellExperiment counts colData
#' @importFrom SummarizedExperiment rowRanges colData
#' @importFrom Matrix rowSums colSums
#'
#' @examples
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp
#' 
create_scExp <- function(
    datamatrix, annot, remove_zero_cells = TRUE, remove_zero_features = TRUE,
    remove_non_canonical = TRUE, remove_chr_M = TRUE, mainExpName = "main",
    verbose = TRUE)
{
    stopifnot(is.data.frame(annot), remove_zero_cells %in% c(TRUE, FALSE),
            remove_zero_features %in% c(TRUE, FALSE))
    if (ncol(datamatrix) != nrow(annot)) 
        stop(
            "ChromSCape::create_scExp - datamatrix and annot should contain",
            "the same number of cells")
    if (length(match(c("cell_id", "sample_id"), colnames(annot))) < 2) 
        stop("ChromSCape::create_scExp - annot should contain cell_id &
            sample_id as column names")
    if (is(datamatrix, "data.frame")) datamatrix <- as.matrix(datamatrix)
    message("ChromSCape::create_scExp - the matrix has ",
        dim(datamatrix)[2], " cells and ", dim(datamatrix)[1], " features.")
    scExp <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = datamatrix), colData = annot)
    
    # If rownames is formatted as gene:chr_start_end, put genes into rowRanges
    contains_genes = all(grepl(":chr",rownames(scExp)[1:5]))
    if(contains_genes){
        message("ChromSCape::create_scExp - Genes detected in rownames...")
        Genes = gsub(":.*","",rownames(scExp))
        rownames(scExp) = gsub(".*:","",rownames(scExp))
        SummarizedExperiment::rowRanges(scExp) = get_genomic_coordinates(scExp)
        SummarizedExperiment::rowRanges(scExp)$Gene = Genes
        SummarizedExperiment::rowRanges(scExp)$distanceToTSS = 0
        rownames(scExp) = SummarizedExperiment::rowRanges(scExp)$ID
    }
    if (has_genomic_coordinates(scExp)) {
        if (remove_non_canonical) scExp <- remove_non_canonical_fun(scExp,
                                                                    verbose)
        if (remove_chr_M) scExp <- remove_chr_M_fun(scExp, verbose)}
    dim_b <- dim(scExp)
    if (remove_zero_features)
        scExp <- scExp[(Matrix::rowSums(
            SingleCellExperiment::counts(scExp) > 0) > 0), ]
    if (remove_zero_cells) 
        scExp <- scExp[, (Matrix::colSums(
            SingleCellExperiment::counts(scExp) > 0) > 0)]
    if (dim(scExp)[2] != dim_b[2]){
        message("ChromSCape::create_scExp - ", dim_b[2] - dim(scExp)[2],
            " cells with 0 signals were removed.")
    }
    if (dim(scExp)[1] != dim_b[1]){
        message("ChromSCape::create_scExp - ",
            dim_b[1] - dim(scExp)[1], " features with 0 signals were removed.")
    }
    if (has_genomic_coordinates(scExp) && !contains_genes){
        rows <- rownames(scExp)
        SummarizedExperiment::rowRanges(scExp) <- get_genomic_coordinates(scExp)
        rownames(scExp) <- rows
    }
    SummarizedExperiment::colData(scExp)$total_counts = 
        colSums(SingleCellExperiment::counts(scExp))
    # if(!is.na(ncol(scExp) * nrow(scExp))){
    # SummarizedExperiment::colData(scExp)$detected = 
    #     apply(SingleCellExperiment::counts(scExp), 2,
    #           function(i) length(which(i>0)))
    # } else{
    #     SummarizedExperiment::colData(scExp)$detected = 
    #         SummarizedExperiment::colData(scExp)$total_counts
    # }
    SingleCellExperiment::mainExpName(scExp) <- mainExpName
    return(scExp)
}

#' Remove non canonical chromosomes from scExp
#'
#' @param scExp A SingleCellExperiment
#' @param verbose Print ?
#'
#' @return A SingleCellExperiment without non canonical chromosomes 
#' (random,unknown, contigs etc...)
#'
remove_non_canonical_fun <- function(scExp, verbose){
    # Removing non-canonical chromosomes
    splitID <- lapply(rownames(scExp), function(x)
    {
        unlist(strsplit(as.character(x),
                        split = "_|:|-",
                        fixed = FALSE))
    })
    normal_chr <- which(as.numeric(lapply(splitID, length)) <= 3)
    # weird chromosomes contain _/:/- in the name
    nrow_init <- nrow(scExp)
    scExp <- scExp[normal_chr, ]
    if (length(normal_chr) < nrow_init && verbose)
    {
        message("ChromSCape::create_scExp - ",nrow(scExp) - length(normal_chr),
                " non canonical regions were removed.\n")
    }
    return(scExp)
}

#' Remove chromosome M from scExprownames
#'
#' @param scExp A SingleCellExperiment
#' @param verbose Print ?
#'
#' @return A SingleCellExperiment without chromosome M (mitochondrial chr)
#'
remove_chr_M_fun <- function(scExp, verbose){
    # Remove chrM from mat if it is inside
    chrM_regions <- grep("chrM", rownames(scExp))
    if (length(chrM_regions) > 0)
    {
        scExp <- scExp[-chrM_regions, ]
        if (verbose)
        {
            message("ChromSCape::create_scExp - ", length(chrM_regions),
                " chromosome M regions were removed.\n")
        }
    }
    return(scExp)
}

#' Filter cells and features
#'
#' Function to filter out cells & features from SingleCellExperiment based on
#' total count per cell, number of cells 'ON' in features and top covered cells
#' that might be doublets.
#'
#' @param scExp A SingleCellExperiment object.
#' @param min_cov_cell Minimum counts for each cell. (1600)
#' @param quant_removal Centile of cell counts above which cells are removed.
#'   (95)
#' @param min_count_per_feature Minimum number of reads per feature (10).
#' @param verbose (TRUE)
#'
#' @return Returns a filtered SingleCellExperiment object.
#'
#' @export
#'
#' @examples
#'
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp. = filter_scExp(scExp)
#'
#' # No feature filtering (all features are valuable)
#' scExp. = filter_scExp(scExp,min_count_per_feature=30)
#'
#' # No cell filtering (all features are valuable)
#' scExp. = filter_scExp(scExp,min_cov_cell=0,quant_removal=100)
#'
#' @importFrom SingleCellExperiment SingleCellExperiment counts colData
#' @importFrom Matrix colSums rowSums
filter_scExp =  function (
    scExp, min_cov_cell = 1600, quant_removal = 95, min_count_per_feature = 10,
    verbose = TRUE){
    stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(min_cov_cell),
            is.numeric(quant_removal),
            is.numeric(min_count_per_feature), verbose %in% c(FALSE, TRUE))
    if (is.null(scExp)) warning(
        "ChromSCape::filter_scExp - Please specify a  SingleCellExperiment")
    
    cellCounts <- Matrix::colSums(SingleCellExperiment::counts(scExp))
    thresh <- stats::quantile(cellCounts, probs = seq(0, 1, 0.01))
    sel <- (cellCounts > min_cov_cell & cellCounts <= thresh[quant_removal + 1])
    if (verbose) message(
        "ChromSCape::filter_scExp - ", length(which(sel)), " cells pass the ",
        " threshold of ", min_cov_cell, " minimum reads and are lower than ",
        "the ", quant_removal, "th centile of library size ~= ",
        round(thresh[quant_removal + 1]), " reads.")
    
    scExp <- scExp[, sel]
    counts <- SingleCellExperiment::counts(scExp)
    sel_feature <- (Matrix::rowSums(counts) >= min_count_per_feature)
    
    if (verbose) message(
        "ChromSCape::filter_scExp - ", length(which(sel_feature))," features pass ",
        "the threshold of ", min_count_per_feature," count per feature.")

    scExp <- scExp[sel_feature,]
    empty_cells = (Matrix::colSums(SingleCellExperiment::counts(scExp)) < min_cov_cell)
    if(any(empty_cells)) scExp <- scExp[,!empty_cells]
    empty_features = (Matrix::rowSums(counts) < min_count_per_feature)
    if(any(empty_features)) scExp <- scExp[!empty_features,]
    
    SummarizedExperiment::colData(scExp)$total_counts = 
        colSums(SingleCellExperiment::counts(scExp))
    return(scExp)
}

#' Find most covered features
#'
#' @description  
#' Find the top most covered features that will be used for dimensionality 
#' reduction. Optionally remove non-top features. 
#'
#' @param scExp A SingleCellExperiment.
#' @param n Either an integer indicating the number of top covered regions to 
#' find or a character vector of the top percentile of features to keep (e.g.
#' 'q20' to keep top 20% features). 
#' @param keep_others Logical indicating if non-top regions are to be removed 
#' from the SCE or not (FALSE). 
#' @param prioritize_genes First filter by loci being close to genes ? E.g. for
#' differential analysis, it is more relevant to keep features close to genes
#' @param max_distanceToTSS If prioritize_genes is TRUE, the maximum distance to 
#' consider a feature close to a gene.
#' 
#' @param verbose Print ?
#' 
#' @return A SCE with top features
#' @export
#'
#' @examples
#' data(scExp)
#' scExp_top = find_top_features(scExp, n = 4000, keep_others = FALSE)
#' 
find_top_features <- function (scExp, n = 20000, keep_others = FALSE,
                               prioritize_genes = FALSE, 
                               max_distanceToTSS = 10000,
                               verbose = TRUE){
    scExp. = scExp
    if(prioritize_genes){
        # filter by loci close to genes
        scExp. = scExp.[which(
            SingleCellExperiment::rowData(scExp)$distanceToTSS <= max_distanceToTSS),]
        if(verbose) message("ChromScape::find_top_features - ", 
                            nrow(scExp.), " features kept as ",
                            " being closer than ", max_distanceToTSS,
                            "bp to genes TSS...")
    }
    n = min(n, nrow(scExp.))
    feature_counts = Matrix::rowSums(counts(scExp.))
    
    if(is.numeric(n)){
        feature_counts_top = sort(feature_counts,decreasing = TRUE)[seq_len(n)]   
    } else{
        n = as.numeric(gsub("q","",n))
        feature_counts_top = feature_counts[which(
            feature_counts>quantile(feature_counts,1-(n/100) ) )]
    }
    SummarizedExperiment::rowData(scExp)[,"top_feature"] = FALSE
    SummarizedExperiment::rowData(scExp)[names(feature_counts_top),
                                         "top_feature"] = TRUE
    if(keep_others == FALSE) scExp = scExp[names(feature_counts_top),]
    
    if(verbose) message("ChromScape::find_top_features - ", 
                        length(feature_counts_top), " features kept as ",
                        " most covered features...")
    return(scExp)
}

#' Does SingleCellExperiment has genomic coordinates in features ?
#'
#' @param scExp A SingleCellExperiment object
#'
#' @return TRUE or FALSE
#' @export
#' @examples 
#' 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' has_genomic_coordinates(scExp)
#' scExp_gene = create_scExp(create_scDataset_raw(featureType="gene")$mat,
#'   create_scDataset_raw(featureType="gene")$annot)
#' has_genomic_coordinates(scExp_gene)
#' 
has_genomic_coordinates <- function(scExp)
{
    stopifnot(is(scExp, "SingleCellExperiment"), !is.null(rownames(scExp)))
    ID <- rownames(scExp)[seq_len(min(10, length(rownames(scExp))))]
    chr <- unlist(lapply(
        strsplit(ID, split = "_|-|:"),
        FUN = function(x)
            unlist(x)[1]
    ))
    
    if (length(grep("chr|(^\\d+$|^X$|^Y$)",
                    chr[seq_len(min(10, length(chr)))], ignore.case = TRUE)) >=
        min(10, length(chr)))
    {
        return(TRUE)
    } else
    {
        return(FALSE)
    }
}

#' Get SingleCellExperiment's genomic coordinates
#'
#' @param scExp A SingleCellExperiment object.
#'
#' @return A GRanges object of genomic coordinates.
#' @importFrom GenomicRanges GRanges
#' @export
#' @examples 
#' 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' feature_GRanges = get_genomic_coordinates(scExp)
#' 
get_genomic_coordinates <- function(scExp)
{
    stopifnot(is(scExp, "SingleCellExperiment"))
    if (!has_genomic_coordinates(scExp))
    {
        stop("Feature names are not genomic coordinates")
    }
    
    ID <- rownames(scExp)
    chr <- unlist(lapply(
        strsplit(ID, split = "_|-|:"),
        FUN = function(x)
            unlist(x)[1]
    ))
    start <- unlist(lapply(
        strsplit(ID, split = "_|-|:"),
        FUN = function(x)
            unlist(x)[2]
    ))
    end <- unlist(lapply(
        strsplit(ID, split = "_|-|:"),
        FUN = function(x)
            unlist(x)[3]
    ))
    feature <-
        GRanges(data.frame(
            ID = ID,
            chr = chr,
            start = start,
            end = end
        ))
    return(feature)
}

#' Remove specific features (CNA, repeats)
#'
#' @param scExp A SingleCellExperiment object.
#' @param features_to_exclude A GenomicRanges object or data.frame containing
#' genomic regions or features to exclude or path towards a BED file containing
#' the features to exclude.
#' @param by Type of features. Either 'region' or 'feature_name'. If 'region',
#'  will look for genomic coordinates in columns 1-3 (chr,start,stop).
#' If 'feature_name', will look for a genes in first column. ('region')
#' @param verbose (TRUE)
#'
#' @return A SingleCellExperiment object without features to exclude.
#' @export
#'
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame findOverlaps
#'  intersect
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom rtracklayer import
#'
#' @examples 
#' 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' features_to_exclude = data.frame(chr=c("chr4","chr7","chr17"),
#' start=c(50000,8000000,2000000),
#' end=c(100000,16000000,2500000))
#' features_to_exclude = as(features_to_exclude,"GRanges")
#' scExp = exclude_features_scExp(scExp,features_to_exclude)
#' scExp
#' 
exclude_features_scExp <-
    function(scExp, features_to_exclude, by = "region", verbose = TRUE){
        stopifnot(is(scExp, "SingleCellExperiment"), is.character(by[1]))
        if(!is(features_to_exclude,"GenomicRanges") & 
           !is(features_to_exclude,"data.frame")) {
            if(!file.exists(features_to_exclude)) {
                stop(
                "ChromSCape::exclude_features_scExp - features_to_exclude must be
           either GenomicRanges or data.frame containg feature names.")
            } else {
                features_to_exclude = rtracklayer::import(features_to_exclude)
            }
        }
        if (!by[1] %in% c("region", "feature_name")) 
            stop("ChromSCape::exclude_features_scExp - by must be either
            'region' or 'feature_name'")
        if (by[1] == "region")
        {
            if (!has_genomic_coordinates(scExp))
                stop(paste0("ChromSCape::exclude_features_scExp -",
                    "Feature names are not genomic coordinates"))
            regions <- SummarizedExperiment::rowRanges(scExp)
            
            suppressWarnings({
                ovrlps <- as.data.frame(
                    GenomicRanges::findOverlaps(regions, features_to_exclude))[,1]
            })
            if (length(unique(ovrlps) > 0)) scExp <- scExp[-unique(ovrlps), ]
            if (verbose) 
                message("ChromSCape::exclude_features_scExp - Removed ",
                    length(unique(ovrlps)), " regions from the analysis.")}
        if (by[1] == "feature_name"){
            if (has_genomic_coordinates(scExp)) 
                warning(
                    "ChromSCape::exclude_features_scExp - Excluding by feature
                    name while object feature names are genomic coordinates !")
            features <- rownames(scExp)
            features_to_exclude <-
                as.character(features_to_exclude[, 1])
            suppressWarnings({
                ovrlps <-
                    GenomicRanges::intersect(features, features_to_exclude)
            })
            if (length(unique(ovrlps) > 0)){
                scExp <- scExp[-which(rownames(scExp) %in% ovrlps), ]
            }
            if (verbose)
                message("ChromSCape::exclude_features_scExp - Removed ",
                    length(unique(ovrlps)), " features from the analysis.")
        }
        return(scExp)
    }

#' Preprocess scExp - Transcripts per Million (TPM)
#'
#' @param scExp A SingleCellExperiment Object
#'
#' @return A SingleCellExperiment object.
#' @importFrom GenomicRanges width
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom SingleCellExperiment counts normcounts
#' @importFrom Matrix t colSums
#' @export
#' @examples 
#' 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp = preprocess_TPM(scExp)
#' head(SingleCellExperiment::normcounts(scExp))
#' 
preprocess_TPM <- function(scExp)
{
    size <- GenomicRanges::width(SummarizedExperiment::rowRanges(scExp))
    SummarizedExperiment::assay(scExp, "normcounts", withDimnames = FALSE) <-
        SingleCellExperiment::counts(scExp) / size
    SummarizedExperiment::assay(scExp, "normcounts", withDimnames = FALSE) <-
        10 ^ 6 * Matrix::t(Matrix::t(SingleCellExperiment::normcounts(scExp)) /
                    Matrix::colSums(SingleCellExperiment::normcounts(scExp)))
    return(scExp)
}

#' Preprocess scExp - Read per Kilobase Per Million (RPKM)
#'
#' @param scExp A SingleCellExperiment Object
#'
#' @return A SingleCellExperiment object.
#' @importFrom GenomicRanges width
#' @importFrom SummarizedExperiment rowRanges assay
#' @importFrom SingleCellExperiment counts normcounts
#' @importFrom Matrix t colSums
#' @export
#' @examples 
#' 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp = preprocess_RPKM(scExp)
#' head(SingleCellExperiment::normcounts(scExp))
preprocess_RPKM <- function(scExp)
{
    SummarizedExperiment::assay(scExp, "normcounts", withDimnames = FALSE) <-
        10 ^ 9 * Matrix::t(
            Matrix::t(SingleCellExperiment::counts(scExp)) /
                Matrix::colSums(SingleCellExperiment::counts(scExp))
        )
    size <-
        GenomicRanges::width(SummarizedExperiment::rowRanges(scExp))
    SummarizedExperiment::assay(scExp, "normcounts", withDimnames = FALSE) <-
        SingleCellExperiment::normcounts(scExp) / size
    
    return(scExp)
}

#' Preprocess scExp - Counts Per Million (CPM)
#'
#' @param scExp A SingleCellExperiment Object
#'
#' @return A SingleCellExperiment object.
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment counts normcounts
#' @importFrom Matrix t colSums
#' @export
#' @examples 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp = preprocess_CPM(scExp)
#' head(SingleCellExperiment::normcounts(scExp))
#' 
preprocess_CPM <- function(scExp)
{
    SummarizedExperiment::assay(scExp, "normcounts", withDimnames = FALSE) <- 
        10 ^ 6 * Matrix::t(Matrix::t(SingleCellExperiment::counts(scExp)) /
                            Matrix::colSums(SingleCellExperiment::counts(scExp))
        )
    return(scExp)
}


#' Preprocess scExp - TF-IDF
#'
#' @param scExp A SingleCellExperiment Object
#' @param scale A numeric to multiply the matrix in order to have human readeable
#' numbers. Has no impact on the downstream analysis
#' @param log Wether to use neperian log on the TF-IDF normalized data or not.
#'
#' @return A SingleCellExperiment object.
#' @importFrom GenomicRanges width
#' @importFrom SummarizedExperiment rowRanges assay
#' @importFrom SingleCellExperiment counts normcounts
#' @importFrom Matrix t colSums
#' @export
#' @examples 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp = preprocess_TFIDF(scExp)
#' head(SingleCellExperiment::normcounts(scExp))
#'
preprocess_TFIDF <- function(scExp, scale = 10000, log = TRUE)
{
        counts = SingleCellExperiment::counts(scExp)
        n_features <- Matrix::colSums(counts)
        tf <- Matrix::t(Matrix::t(counts) / n_features)
        idf <- 1 + ncol(counts) / Matrix::rowSums(counts)
        normcounts <- Matrix::Diagonal(length(idf), idf) %*% tf
        if(log) normcounts = log1p(normcounts * scale) else normcounts = 
            normcounts * scale
        SummarizedExperiment::assay(scExp, "normcounts",
                                    withDimnames = FALSE) <- normcounts
        return(scExp)
}

#' Preprocess scExp - size only
#'
#' @param scExp A SingleCellExperiment Object
#'
#' @return A SingleCellExperiment object.
#' @importFrom GenomicRanges width
#' @importFrom SummarizedExperiment rowRanges assay
#' @importFrom SingleCellExperiment counts normcounts
#' @importFrom Matrix t colSums
#' @export
#' @examples 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp = preprocess_feature_size_only(scExp)
#' head(SingleCellExperiment::normcounts(scExp))
#'
preprocess_feature_size_only <- function(scExp)
{
    size <- GenomicRanges::width(SummarizedExperiment::rowRanges(scExp))
    SummarizedExperiment::assay(scExp, "normcounts", withDimnames = FALSE) <-
        SingleCellExperiment::counts(scExp) / size
    return(scExp)
}

#' Normalize counts
#'
#' @param scExp A SingleCellExperiment object.
#' @param type Which normalization to apply. Either 'CPM', 'TFIDF','RPKM', 'TPM' or
#'  'feature_size_only'. Note that for all normalization by size
#'  (RPKM, TPM, feature_size_only), the features must have defined
#'  genomic coordinates.
#'
#' @return A SingleCellExperiment object containing normalized counts.
#'  (See ?normcounts())
#' @export
#'
#' @examples 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp = normalize_scExp(scExp)
#' head(SingleCellExperiment::normcounts(scExp))
#'
normalize_scExp <- function(scExp,
                            type = c("CPM", "TFIDF", "RPKM", "TPM",
                                    "feature_size_only"))
{
    stopifnot(
        type[1] %in% c("CPM", "TFIDF", "RPKM", "TPM", "feature_size_only"),
        is(scExp, "SingleCellExperiment")
    )
    if (!is(SingleCellExperiment::counts(scExp), "dgCMatrix")) {
        SingleCellExperiment::counts(scExp) <- as(
            SingleCellExperiment::counts(scExp), "dgCMatrix")
    }
    
    if (!has_genomic_coordinates(scExp))
    {
        warning(
            "ChromSCape::normalize_scExp - Switching to CPM normalization
            as features are not genomic coordinates."
        )
        type <- "CPM"
    }
    switch(
        type[1],
        RPKM = return(preprocess_RPKM(scExp)),
        TPM = return(preprocess_TPM(scExp)),
        feature_size_only = return(preprocess_feature_size_only(scExp)),
        CPM = return(preprocess_CPM(scExp)),
        TFIDF = return(preprocess_TFIDF(scExp))
    )
}


#' Add gene annotations to features
#'
#' @param scExp A SingleCellExperiment object.
#' @param ref Reference genome. Either 'hg38' or 'mm10'. ('hg38')
#' @param reference_annotation A data.frame containing gene (or else) annotation
#'   with genomic coordinates.
#'
#' @return A SingleCellExperiment object with annotated rowData.
#' @export
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame distanceToNearest
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom dplyr mutate select group_by summarise_all
#'
#' @examples
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp = feature_annotation_scExp(scExp)
#' head(SummarizedExperiment::rowRanges(scExp))
#'
#' # Mouse
#' scExp = create_scExp(create_scDataset_raw(ref="mm10")$mat,
#'   create_scDataset_raw(ref="mm10")$annot)
#' scExp = feature_annotation_scExp(scExp,ref="mm10")
#' head(SummarizedExperiment::rowRanges(scExp))
feature_annotation_scExp <- function(scExp, ref = "hg38",
                                    reference_annotation = NULL)
{
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(ref))
    if (is.null(SummarizedExperiment::rowRanges(scExp)))
        stop("ChromSCape::feature_annotation_scExp - The object doesn't have
            ranges of coordinates as rowData")
    if (is.null(reference_annotation) & !(ref %in% c("hg38", "mm10")))
        stop("ChromSCape::feature_annotation_scExp - If reference_annotation is
            null, ref must be either 'hg38' or 'mm10' to automatically load
            reference gene annotation.")
    if (is.null(reference_annotation) & (ref %in% c("hg38", "mm10")))
    {
        message("ChromSCape::feature_annotation_scExp - Selecting ",
                ref, " genes from Gencode.")
        eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
        reference_annotation <- eval(parse(text = paste0("", ref, ".GeneTSS")))
        start = reference_annotation$start
        reference_annotation$start = ifelse(reference_annotation$strand == "+",
                                            reference_annotation$start,
                                            reference_annotation$end)
        reference_annotation$end = ifelse(reference_annotation$strand == "+",
                                start + 1,
                                reference_annotation$end + 1)
        reference_annotation$strand = NULL
    }
    if (is.data.frame(reference_annotation)) 
        reference_annotation <-
            GenomicRanges::makeGRangesFromDataFrame(reference_annotation,
                                                    keep.extra.columns = TRUE)
    
    feature_ranges <- SummarizedExperiment::rowRanges(scExp)[, 1]
    hits <- GenomicRanges::distanceToNearest(feature_ranges, 
                                                reference_annotation,
                                                ignore.strand = TRUE,
                                                select = "all")
    q_hits <- S4Vectors::queryHits(hits)
    s_hits <- S4Vectors::subjectHits(hits)
    annotFeat <- data.frame(
        "chr" = as.character(GenomicRanges::seqnames(feature_ranges[q_hits])),
        "start" = as.character(GenomicRanges::start(feature_ranges[q_hits])),
        "end" = as.character(GenomicRanges::end(feature_ranges[q_hits])),
        "Gene" = as.character(
            reference_annotation@elementMetadata$Gene)[s_hits],
        "distanceToTSS" = hits@elementMetadata$distance)
    annotFeat <- annotFeat %>% dplyr::mutate(
        "ID" = paste(.data$chr, .data$start, .data$end, sep = "_")) %>%
        dplyr::select("ID", "chr", "start", "end", "Gene", "distanceToTSS")
    annotFeat <- annotFeat %>% dplyr::group_by(.data$ID, .data$chr,
                                                .data$start, .data$end) %>% 
        dplyr::summarise("Gene" = paste(.data$Gene, collapse = ", "),
                            "distanceToTSS" = max(.data$distanceToTSS)) %>%
        as.data.frame()
    annotFeat <- annotFeat[match(rownames(scExp), annotFeat$ID), ]
    rownames(annotFeat) <- annotFeat$ID
    if(all(c("ID","Gene","distanceToTSS") %in%
           colnames(SummarizedExperiment::rowData(scExp)))){
        SummarizedExperiment::rowData(scExp) <- annotFeat[,c("ID","Gene","distanceToTSS")]
    } else{
        SummarizedExperiment::rowData(scExp) <- cbind(
            SummarizedExperiment::rowData(scExp),
            annotFeat[,c("Gene","distanceToTSS")]
        )
    }
    
    return(scExp)
}

#' Choose perplexity depending on number of cells for Tsne
#'
#' @param dataset A matrix of features x cells (rows x columns)
#'
#' @return A number between 5 and 30 to use in Rtsne function
choose_perplexity <- function(dataset)
{
    stopifnot(!is.null(dataset), !is.null(dim(dataset)))
    perplexity <- 30
    if (nrow(dataset) <= 200)
    {
        perplexity <- 20
    }
    if (nrow(dataset) <= 250)
    {
        perplexity <- 25
    }
    if (nrow(dataset) <= 150)
    {
        perplexity <- 15
    }
    if (nrow(dataset) <= 100)
    {
        perplexity <- 10
    }
    if (nrow(dataset) <= 50)
    {
        perplexity <- 5
    }
    perplexity
}

#' Reduce dimensions (PCA, TSNE, UMAP)
#'
#' @param scExp A SingleCellExperiment object.
#' @param dimension_reductions A character vector of methods to apply.
#' (c('PCA','TSNE','UMAP'))
#' @param n Numbers of dimensions to keep for PCA. (50)
#' @param batch_correction Do batch correction ? (FALSE)
#' @param batch_list List of characters. Names are batch names, characters are
#'  sample names.
#' @param remove_PC1 Remove PC1 before UMAP & T-SNE, as probably correlated to 
#' library size ? Recommended when using 'TFIDF' normalization method. (FALSE)
#' @param verbose Print messages ?(TRUE)
#'
#' @return A SingleCellExperiment object containing feature spaces. See
#'  ?reduceDims().
#'
#' @export
#'

#' @importFrom Rtsne Rtsne
#' @importFrom umap umap umap.defaults
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment counts normcounts reducedDims
#' @importFrom batchelor fastMNN
#' @importFrom Matrix t
#' 
#' @examples 
#' 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp = reduce_dims_scExp(scExp,dimension_reductions=c("PCA","UMAP"))
#' scExp = normalize_scExp(scExp)
#' scExp = reduce_dims_scExp(scExp,dimension_reductions=c("PCA","UMAP"))
reduce_dims_scExp <-
    function(scExp, dimension_reductions = c("PCA", "UMAP"), n = 50,
            batch_correction = FALSE, batch_list = NULL,
            remove_PC1 = FALSE, verbose = TRUE)
    {
        stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(n),
            dimension_reductions[1] %in% c("PCA", "TSNE", "UMAP"),
            is.logical(remove_PC1))
        if (!"normcounts" %in% names(SummarizedExperiment::assays(scExp))){
            warning("ChromSCape::reduce_dims_scExp - The raw counts are not
                normalized, running dimensionality reduction on raw counts.")
            mat <- SingleCellExperiment::counts(scExp)
        } else{
            mat <- SingleCellExperiment::normcounts(scExp)
            }
        if (batch_correction && !is.list(batch_list)){
            stop("ChromSCape::reduce_dims_scExp - If doing batch correction,
                batch_list must be a list of the samples IDs of each batch.")}
        if (batch_correction){
            out <- reduce_dim_batch_correction(scExp, mat, batch_list, n)
            scExp <- out$scExp
            pca <- out$pca
        } else{
            scExp$batch_id <- factor(1)
            if ( is(mat, "dgCMatrix") | is(mat, "dgTMatrix")) {
                pca <- pca_irlba_for_sparseMatrix(Matrix::t(mat), n)
            } else{
                pca <- stats::prcomp(Matrix::t(mat),center = TRUE,
                                    scale. = FALSE)
                pca <- pca$x[, seq_len(n)]}}
        pca <- as.data.frame(as.matrix(pca))
        colnames(pca) <- paste0("Component_", seq_len(ncol(pca)))
        if(remove_PC1) {
            pca <- pca[,2:n]
            if(verbose) message("ChromSCape::reduce_dims_scExp - removing ",
                                "PC1... (probably correlated to library size)")
        }
        rownames(pca) <- colnames(scExp)
        if ("TSNE" %in% dimension_reductions){
            tsne <- Rtsne::Rtsne(pca, dims = 2, pca = FALSE, theta = 0,
                perplexity = choose_perplexity(pca), verbose = verbose,
                max_iter = 1000)
            tsne <- as.data.frame(tsne$Y)
            colnames(tsne) <- c("Component_1", "Component_2")
        }
        if ("UMAP" %in% dimension_reductions){
            config <- umap::umap.defaults
            config$metric <- "cosine"
            umap <- umap::umap(pca, config = config, method = "naive")
            umap <- as.data.frame(umap$layout)
            colnames(umap) <- c("Component_1", "Component_2")
        }
        listReducedDim <- list(PCA = pca)
        if ("TSNE" %in% dimension_reductions) listReducedDim$TSNE <- tsne
        if ("UMAP" %in% dimension_reductions) listReducedDim$UMAP <- umap
        SingleCellExperiment::reducedDims(scExp) <- listReducedDim
        return(scExp)
    }

#' Reduce dimension with batch corrections
#'
#' @param scExp SingleCellExperiment
#' @param batch_list List of batches
#' @param mat The normalized count matrix
#' @param n Number of PCs to keep
#'
#' @return A list containing the SingleCellExperiment with batch info and
#' the corrected pca
reduce_dim_batch_correction <- function(scExp, mat, batch_list, n){
    print("Running Batch Correction ...")
    num_batches <- length(batch_list)
    batch_names <- names(batch_list)
    batches <- list()
    scExp. <- scExp
    scExp.$batch_name <- "unspecified"
    for (i in seq_len(num_batches))
    {
        for (s_id in batch_list[[i]])
        {
            SummarizedExperiment::colData(
                scExp.)[scExp.$sample_id == s_id,
                        "batch_name"] <- batch_names[i]
        }
    }
    adj_annot <- data.frame()
    b_names <- unique(scExp.$batch_name)
    scExp.$batch_id <- as.factor(as.numeric(as.factor(scExp.$batch_name)))
    if (is(mat,"dgCMatrix") | is(mat, "dgTMatrix"))
    {
        pca <- pca_irlba_for_sparseMatrix(Matrix::t(mat), n)
    } else
    {
        pca <- stats::prcomp(Matrix::t(mat),
                            center = TRUE,
                            scale. = FALSE)
        pca <- pca$x[, seq_len(n)]
    }
    for (i in seq_along(b_names))
    {
        b_name <- b_names[i]
        batches[[i]] <-
            as.matrix(mat[,which(scExp.$batch_name == b_name)])
        adj_annot <- rbind(adj_annot, as.data.frame(
            SummarizedExperiment::colData(
                scExp.)[scExp.$batch_name ==b_name, ]))
    }
    mnn.out <- batchelor::fastMNN(batches, k = 25, d = n, ndist = 3,
                auto.merge = FALSE, cos.norm = FALSE)
    pca <- SingleCellExperiment::reducedDim(mnn.out,"corrected")
    SummarizedExperiment::colData(scExp.) <- as(adj_annot,
                                            "DataFrame")
    out = list("scExp" = scExp., "pca"= pca)
    return(out)
}

#' Run sparse PCA using irlba SVD
#' 
#' @description 
#' This function allows to run a PCA using IRLBA Singular Value Decomposition
#' in a fast & memory efficient way. The increamental Lanczos bidiagonalisation
#' algorithm allows to keep the matrix sparse as the "loci" centering is 
#' implicit. 
#' The function then multiplies by the approximate singular values (svd$d) in 
#' order to get more importance to the first PCs proportionnally to their 
#' singular values. This step is crucial for downstream approaches, e.g. UMAP or
#' T-SNE.
#' 
#' @param x A sparse normalized matrix (features x cells)
#' @param n_comp The number of principal components to keep
#' @param work Working subspace dimension, larger values can speed convergence
#'  at the cost of more memory use.
#'
#' @return The rotated data, e.g. the cells x PC column in case of sc data.
#'
#' @importFrom irlba irlba
#' @importFrom Matrix colMeans
#'
pca_irlba_for_sparseMatrix <- function(x, n_comp, work = 3 * n_comp)
{
    x.means <- Matrix::colMeans(x)
    svd.0 <- irlba::irlba(x, center = x.means, nv = n_comp, work = work)
    pca <- svd.0$u %*% diag(svd.0$d)
    return(pca)
}

sweep_sparse <- function(x, margin, stats, fun = "*") {
    f <- match.fun(fun)
    if (margin == 1) {
        idx <- x@i + 1
    } else {
        idx <- x@p + 1
    }
    x@x <- f(x@x, stats[idx])
    return(x)
}

#' Table of cells
#'
#' @param annot An annotation of cells. Can be obtain through 'colData(scExp)'. 
#' @param datamatrix A matrix of cells per regions  before filtering.
#' 
#' @export
#' @return A formatted kable in HTML.
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr bind_rows tibble left_join n summarise
#' @importFrom kableExtra kable kable_styling group_rows
#' 
#' @examples 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' \dontrun{num_cell_scExp(SingleCellExperiment::colData(scExp))}
num_cell_scExp <- function(annot, datamatrix)
{
    stopifnot(!is.null(annot), !is.null(datamatrix))
    annot$total_counts = colSums(datamatrix)
    table <- as.data.frame(
        annot %>% dplyr::group_by(sample_id) %>% 
            dplyr::summarise("#Cells" = n(),
                             "Median"= round(median(total_counts),0),
                             "Std"= round(stats::sd(total_counts),0)))
    
    rownames(table) <- NULL
    
    table[, 1] <- as.character(table[, 1])
    table[nrow(table) + 1, ] = c(
        "", sum(table[, 2]), round(median(annot$total_counts),0),
        round(stats::sd(annot$total_counts),0))
    
    table %>% kableExtra::kable(escape = FALSE, align = "c") %>%
        kableExtra::kable_styling(c("striped",
                                    "condensed"), full_width = TRUE) %>%
        kableExtra::group_rows("Total",
                            dim(table)[1], dim(table)[1])
}

#' Table of cells before / after QC
#'
#' @param scExp A SingleCellExperiment object.
#' @param annot A raw annotation data.frame of cells before filtering.
#' @param datamatrix A matrix of cells per regions  before filtering.
#'
#' @export
#' @return A formatted kable in HTML.
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr bind_rows tibble left_join n summarise
#' @importFrom kableExtra kable kable_styling group_rows
#' @importFrom kableExtra kable kable_styling group_rows
#'
#' @examples 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp_filtered = filter_scExp(scExp)
#' \dontrun{ num_cell_after_QC_filt_scExp(
#' scExp_filtered,SingleCellExperiment::colData(scExp))}
#' 
num_cell_after_QC_filt_scExp <- function(scExp, annot, datamatrix)
{
    stopifnot(is(scExp, "SingleCellExperiment"), !is.null(annot))
    annot$total_counts = colSums(datamatrix)
    table <- as.data.frame(
        annot %>% dplyr::group_by(sample_id) %>% 
        dplyr::summarise("#Cells" = n(),
                         "Median"= round(median(total_counts),0),
                         "Std"= round(stats::sd(total_counts),0)))
    table_filtered <-
        as.data.frame(SingleCellExperiment::colData(scExp))
    table_filtered <- as.data.frame(
        table_filtered %>% dplyr::group_by(sample_id) %>% 
            dplyr::summarise("#Cells (filt.)" = n(),
                             "Median (filt.)"= round(median(total_counts),0),
                             "Std (filt.)"= round(stats::sd(total_counts),0)))
    
    colnames(table)[1] <- c("Sample")
    rownames(table) <- NULL
    colnames(table_filtered) [1] <- c("Sample")
    rownames(table_filtered) <- NULL
    
    table_both <-
        dplyr::left_join(table, table_filtered, by = c("Sample"))
    table_both[, 1] <- as.character(table_both[, 1])
    table_both <- dplyr::bind_rows(
        table_both,
        dplyr::tibble(
            "Sample" = "",
            "#Cells" = sum(table_both[,2]),
            "Median"= round(median(annot$total_counts),0),
            "Std"= round(stats::sd(annot$total_counts),0),
            "#Cells (filt.)" = sum(table_both[, 4]),
            "Median (filt.)"= round(median(scExp$total_counts),0),
            "Std (filt.)"= round(stats::sd(scExp$total_counts),0)))
    table_both %>% kableExtra::kable(escape = FALSE, align = "c") %>%
        kableExtra::kable_styling(c("striped",
                                    "condensed"), full_width = TRUE) %>%
        kableExtra::group_rows("Total",
                            dim(table_both)[1], dim(table_both)[1])
}

#' Subsample scExp
#'
#' Randomly sample x cells from each sample in a SingleCellExperiment to return
#' a subsampled SingleCellExperiment with all samples having maximum n cells. If
#' n is higher than the number of cell in a sample, this sample will not be
#' subsampled.
#'
#' @param scExp A SingleCellExperiment
#' @param n_cells An integer number of cells to subsample for each sample (500)
#'
#' @return A subsampled SingleCellExperiment
#' @export
#'
#' @importFrom SingleCellExperiment colData
#' 
#' @examples 
#' raw <- create_scDataset_raw()
#' scExp = create_scExp(raw$mat, raw$annot)
#' scExp_sub = subsample_scExp(scExp,50)
#' \dontrun{num_cell_scExp(scExp_sub)}
#' 
subsample_scExp <- function(scExp, n_cells = 500) {
    stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(n_cells))
    
    annot = as.data.frame(SingleCellExperiment::colData(scExp))
    counts = SingleCellExperiment::counts(scExp)
    samples = as.character(unique(annot$sample_id))
    counts. = NULL
    annot. = NULL
    
    for (samp in samples) {
        cells = as.character(annot$cell_id[which(annot$sample_id == samp)])
        cells = sample(cells, min(n_cells, length(cells)), replace = FALSE)
        if (is.null(counts.))
            counts. = counts[, cells]
        else
            counts. = Matrix::cbind2(counts., counts[, cells])
        if (is.null(annot.))
            annot. = annot[cells,]
        else
            annot. = Matrix::rbind2(annot., annot[cells,])
    }
    message("ChromSCape::subsample_scExp - Subsampling each sample to ",n_cells)
    ord = order(colnames(counts.))
    ord2 = order(rownames(annot.))
    counts. = counts.[, ord]
    annot. = annot.[ord2,]
    scExp. = create_scExp(
        counts.,
        annot.,
        remove_zero_cells = FALSE,
        remove_zero_features = FALSE,
        remove_non_canonical = FALSE,
        remove_chr_M = FALSE
    )
    
    return(scExp.)
}
