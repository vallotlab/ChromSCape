## Authors : Pacôme Prompsy Title : Wrappers & functions to preprocess & reduce
## dimensionality of single-cell Matrix Description : Funtions to load, filter,
## normalize & reduce single-cell Epigenetic Matrices prior to analysis


#' Create a simulated single cell datamatrix & cell annotation
#'
#' @param cells
#' @param features
#' @param featureType
#' @param sparse
#' @param nsamp
#' @param ref
#' @param batch_id
#'
#' @return
#' @export
#'
#' @importFrom GenomicRanges GRanges tileGenome width seqnames GRangesList
#' sort.GenomicRanges
#'
create_scDataset_raw <- function(cells = 300,
                                 features = 600,
                                 featureType = c("window", "peak", "gene"),
                                 sparse = T,
                                 nsamp = 4,
                                 ref = "hg38",
                                 batch_id = rep(1, nsamp))
{
    stopifnot(
        featureType %in% c("window", "peak", "gene"),
        ref %in% c("mm10", "hg38"),
        nsamp >= 1,
        cells >= nsamp,
        features >= 1,
        length(batch_id) == nsamp
    )
    
    set.seed(47)
    
    # Create cell names
    cell_counts <-
        sapply(split(1:cells, sample(nsamp, cells, repl = TRUE)),
               length)
    cell_names <- sample <- batches <- list()
    for (i in 1:length(cell_counts))
    {
        cell_names[[i]] <- paste0("sample_", i, "_c", 1:cell_counts[i])
        sample[[i]] <- rep(paste0("sample_", i), cell_counts[i])
        batches[[i]] <- rep(batch_id[i], cell_counts[i])
    }
    cell_names <- as.character(unlist(cell_names))
    sample <- as.character(unlist(sample))
    batches <- as.numeric(unlist(batches))
    
    # Create feature names
    eval(parse(text = paste0(
        "chr <- ChromSCape::", ref, ".chromosomes"
    )))
    # load species chromsizes
    chr <- GenomicRanges::GRanges(chr)
    
    if (featureType[1] == "window")
    {
        chr_ranges <- unlist(GenomicRanges::tileGenome(
            setNames(
                GenomicRanges::width(chr),
                GenomicRanges::seqnames(chr)
            ),
            ntile = features
        ))[1:features]
        # ~constant window size
        features_names <- paste(
            as.data.frame(chr_ranges)$seqnames,
            as.data.frame(chr_ranges)$start,
            as.data.frame(chr_ranges)$end,
            sep = "_"
        )
    }
    if (featureType[1] == "peak")
    {
        size_peaks <- c(1000, 2500, 7999, 10000, 150000, 10 ^ 6)
        # Different size of peaks
        peaks <- sapply(split(1:features, sample(
            length(size_peaks),
            features, repl = TRUE
        )),
        length)
        chr_ranges_list <- GenomicRanges::GRangesList()
        for (i in 1:length(peaks))
        {
            chr_ranges <- unlist(
                GenomicRanges::tileGenome(
                    setNames(
                        GenomicRanges::width(chr),
                        GenomicRanges::seqnames(chr)
                    ),
                    tilewidth = size_peaks[i],
                    cut.last.tile.in.chrom = F
                )
            )
            chr_ranges_list[[i]] <-
                chr_ranges[sample(1:length(chr_ranges),
                                  size = peaks[i]), ]
        }
        chr_ranges <- GenomicRanges::sort.GenomicRanges(unlist(chr_ranges_list))[1:features]
        
        features_names <- paste(
            as.data.frame(chr_ranges)$seqnames,
            as.data.frame(chr_ranges)$start,
            as.data.frame(chr_ranges)$end,
            sep = "_"
        )
    }
    if (featureType[1] == "gene")
    {
        eval(parse(text = paste0(
            "chr <- ChromSCape::", ref, ".GeneTSS"
        )))
        # load species chromsizes
        features_names <-
            as.character(sample(chr$gene, features, replace = F))
    }
    
    
    vec <- stats::rpois(cells * features, 0.5)
    # Add count to values > 0, iteratively
    for (i in 1:10)
    {
        vec[vec >= i] <- vec[vec >= i] + i ^ 2 *
            stats::rpois(length(vec[vec >= i]),
                         0.5)
    }
    
    indices_vec <- which(vec > 0)
    j <- ceiling(indices_vec / features)
    i <- ceiling(indices_vec %% (features))
    i[i == 0] <- 600
    if (sparse)
    {
        mat <- Matrix::sparseMatrix(i,
                                    j,
                                    x = vec[indices_vec],
                                    dimnames = list(features_names,
                                                    cell_names))
    } else
    {
        mat <- matrix(
            vec,
            nrow = features,
            ncol = cells,
            dimnames = list(features_names,
                            cell_names)
        )
    }
    
    if (length(unique(batches)) > 1)
    {
        mat <- mat %*% as(Matrix::diag(batches * batches), "dgCMatrix")
    }
    colnames(mat) <- cell_names
    annot <- data.frame(
        barcode = paste0("BC", 1:cells),
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

#' Wrapper to create the single cell experiment from count matrix and feature
#'  dataframe
#'
#' Create the single cell experiment from (sparse) datamatrix and feature 
#' dataframe containing feature names and location.
#' Also optionally removes zero count Features, zero count Cells, non canconical
#'  chromosomes, and chromosome M. Calculates QC Metrics (scran).
#'
#' @param datamatrix A matrix or sparseMatrix of raw counts. Features x Cells
#'  (rows x columns).
#' @param annot A data.frame containing informations on cells. Should have the
#' same number of rows as the number of columns in datamatrix.
#' @param remove_zero_cells remove cells with zero counts ? [T]
#' @param remove_zero_features remove cells with zero counts ? [T]
#' @param remove_non_canonical remove non canonical chromosomes ?[T]
#' @param remove_chr_M remove chromosomes M ? [T]
#' @param verbose [TRUE]
#'
#' @return Returns a SingleCellExperiment object.
#' @export
#'
#' @examples
#'
#' @importFrom SingleCellExperiment SingleCellExperiment counts colData
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom Matrix rowSums colSums
#' @importFrom scater calculateQCMetrics
#'
#'
create_scExp <- function(datamatrix,
                         annot,
                         remove_zero_cells = TRUE,
                         remove_zero_features = TRUE,
                         remove_non_canonical = TRUE,
                         remove_chr_M = TRUE,
                         verbose = TRUE)
{
    stopifnot(
        is.data.frame(annot),
        remove_zero_cells %in% c(T, F),
        remove_zero_features %in%
            c(T, F)
    )
    
    if (ncol(datamatrix) != nrow(annot))
    {
        stop(
            "ChromSCape::create_scExp - datamatrix and annot should contain
             the same number of cells"
        )
    }
    if (length(match(c("cell_id", "sample_id"), colnames(annot))) < 2)
    {
        stop(
            "ChromSCape::create_scExp - annot should contain cell_id &
             sample_id as column names"
        )
    }
    if (class(datamatrix) == "data.frame")
    {
        datamatrix <- as.matrix(datamatrix)
    }
    
    scExp <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = datamatrix),colData = annot)
    
    if (has_genomic_coordinates(scExp))
    {
        if (remove_non_canonical)
        {
            # Removing non-canonical chromosomes
            splitID <- sapply(rownames(scExp), function(x)
            {
                strsplit(as.character(x),
                         split = "_|:|-",
                         fixed = F)
            })
            normal_chr <- which(sapply(splitID, length) <= 3)
            # weird chromosomes contain _/:/- in the name
            nrow_init <- nrow(scExp)
            scExp <- scExp[normal_chr,]
            if (length(normal_chr) < nrow_init && verbose)
            {
                cat(
                    "ChromSCape::create_scExp -",
                    nrow(scExp) - length(normal_chr),
                    "non canonical regions were removed.\n"
                )
            }
        }
        if (remove_chr_M)
        {
            # Remove chrM from mat if it is inside
            chrM_regions <- grep("chrM", rownames(scExp))
            if (length(chrM_regions) > 0)
            {
                scExp <- scExp[-chrM_regions,]
                if (verbose)
                {
                    cat(
                        "ChromSCape::create_scExp -",
                        length(chrM_regions),
                        "chromosome M regions were removed.\n"
                    )
                }
            }
        }
    }
    dim_b <- dim(scExp)
    if (remove_zero_features)
    {
        scExp <- scExp[(Matrix::rowSums(counts(scExp) > 0) > 0),]
    }  # remove windows that do not have any read in any cells
    if (remove_zero_cells)
    {
        scExp <- scExp[, (Matrix::colSums(counts(scExp) > 0) > 0)]
    }  # remove cells that do not have any read in any cells
    
    if (dim(scExp)[2] != dim_b[2])
    {
        cat(
            "ChromSCape::create_scExp -",
            dim_b[2] - dim(scExp)[2],
            "cells with 0 signals were removed.\n"
        )
    }
    if (dim(scExp)[1] != dim_b[1])
    {
        cat(
            "ChromSCape::create_scExp -",
            dim_b[1] - dim(scExp)[1],
            "features with 0 signals were removed.\n"
        )
    }
    
    
    if (has_genomic_coordinates(scExp))
    {
        rows <- rownames(scExp)
        SummarizedExperiment::rowRanges(scExp) <-
            get_genomic_coordinates(scExp)
        rownames(scExp) <- rows
    }
    
    scExp <- scater::calculateQCMetrics(scExp)
    
    return(scExp)
}


#' Filter cells and features
#'
#' Function to filter out cells & features from SingleCellExperiment based on 
#' total count per cell, number of cells 'ON' in features and top covered
#' cells that might be doublets.
#'
#' @param scExp A SingleCellExperiment object.
#' @param min_cov_cell Minimum counts for each cell. [1600]
#' @param quant_removal Centile of cell counts above which cells are removed.
#'  [95]
#' @param percentMin Minimum percent of cells 'ON' in feature. [1]
#' @param bin_min_count Minimum number of counts to define if cell is 'ON'. [2]
#' @param verbose [T]
#'
#' @return Returns a filtered SingleCellExperiment object.
#'
#' @export
#'
#' @examples
#'
#' @importFrom SingleCellExperiment SingleCellExperiment counts colData
#' @importFrom Matrix colSums rowSums
#' @importFrom scater calculateQCMetrics
filter_scExp <- function(scExp,
                         min_cov_cell = 1600,
                         quant_removal = 95,
                         percentMin = 1,
                         bin_min_count = 2,
                         verbose = T)
{
    stopifnot(
        is(scExp, "SingleCellExperiment"),
        is.numeric(min_cov_cell),
        is.numeric(quant_removal),
        is.numeric(percentMin),
        is.numeric(bin_min_count),
        verbose %in% c(F, T)
    )
    
    if (is.null(scExp))
    {
        warning("ChromSCape::filter_scExp -
                Please specify a SingleCellExperiment")
    }
    
    cellCounts <- Matrix::colSums(counts(scExp))
    
    thresh <- stats::quantile(cellCounts, probs = seq(0, 1, 0.01))
    
    sel1000 <-
        (cellCounts > 1000 & cellCounts <= thresh[quant_removal + 1])
    sel <-
        (cellCounts > min_cov_cell &
             cellCounts <= thresh[quant_removal + 1])
    
    if (verbose)
    {
        cat(
            "ChromSCape::filter_scExp -",
            length(which(sel)),
            "cells pass the threshold of",
            min_cov_cell,
            "minimum reads and are lower than the ",
            quant_removal,
            "th centile of library size ~=",
            round(thresh[quant_removal + 1]),
            "reads.\n"
        )
    }
    
    # Create the binary matrix of cells >= 1000 counts if count >= bin_min_count,
    # count = 1 else
    bina_counts <- counts(scExp)[, sel1000]
    sel_below_2 <- (bina_counts < bin_min_count)
    bina_counts[sel_below_2] <- 0
    sel_above_1 <- (bina_counts >= bin_min_count)
    bina_counts[sel_above_1] <- 1
    
    nCells_in_feature <- Matrix::rowSums(bina_counts)
    fixedFeature <- names(which(nCells_in_feature >
                                    ((percentMin / 100) * (
                                        ncol(bina_counts)
                                    ))))
    if (verbose)
    {
        cat(
            "ChromSCape::filter_scExp -",
            length(fixedFeature),
            "features pass the threshold of",
            percentMin,
            " % of total cells 'ON', representing a minimum of",
            round((percentMin / 100) * (ncol(
                bina_counts
            ))),
            "cells.\n"
        )
    }
    
    scExp <- scExp[, sel]
    scExp <- scExp[fixedFeature,]
    scExp <- scater::calculateQCMetrics(scExp)
    return(scExp)
}

#' Does SingleCellExperiment has genomic coordinates in features ?
#'
#' @param scExp
#'
#' @return TRUE or FALSE
#'
has_genomic_coordinates <- function(scExp)
{
    stopifnot(is(scExp, "SingleCellExperiment"),!is.null(rownames(scExp)))
    ID <- rownames(scExp)[1:min(10, length(rownames(scExp)))]
    chr <- unlist(lapply(
        strsplit(ID, split = "_|-|:"),
        FUN = function(x)
            unlist(x)[1]
    ))
    
    if (length(grep("chr|(^\\d+$|^X$|^Y$)",
                    chr[1:min(10, length(chr))], ignore.case = T)) >=
        min(10, length(chr)))
    {
        return(T)
    } else
    {
        return(F)
    }
}

#' Get SingleCellExperiment's genomic coordinates
#'
#' @param scExp A SingleCellExperiment object.
#'
#' @return A GRanges object of genomic coordinates.
#' @importFrom GenomicRanges GRanges

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
#' @param features_to_exclude A data.frame containing features to exclude.
#' @param by Type of features. Either 'region' or 'feature_name'. If 'region',
#'  will look for genomic coordinates in columns 1-3 (chr,start,stop). 
#' If 'feature_name', will look for a genes in first column. ['region']
#' @param verbose [T]
#'
#' @return A SingleCellExperiment object without features to exclude.
#' @export
#'
#' @examples
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame findOverlaps
#'  intersect
#' @importFrom SummarizedExperiment rowRanges
#'
exclude_features_scExp <-
    function(scExp,
             features_to_exclude,
             by = "region",
             verbose = T)
    {
        stopifnot(
            is(scExp, "SingleCellExperiment"),
            is.data.frame(features_to_exclude),
            is.character(by[1])
        )
        if (!by[1] %in% c("region", "feature_name"))
        {
            stop(
                "ChromSCape::exclude_features_scExp - by must be either 'region'
             or 'feature_name'"
            )
        }
        
        if (by[1] == "region")
        {
            if (!has_genomic_coordinates(scExp))
            {
                stop(
                    "ChromSCape::exclude_features_scExp -
                 Feature names are not genomic coordinates"
                )
            }
            regions <- SummarizedExperiment::rowRanges(scExp)
            colnames(features_to_exclude)[1:3] <-
                c("chr", "start", "stop")
            excl_gr <-
                GenomicRanges::makeGRangesFromDataFrame(
                    features_to_exclude,
                    ignore.strand = TRUE,
                    seqnames.field = c("chr"),
                    start.field = c("start"),
                    end.field = c("stop")
                )
            ovrlps <-
                as.data.frame(GenomicRanges::findOverlaps(regions, excl_gr))[,
                                                                             1]
            if (length(unique(ovrlps) > 0))
            {
                scExp <- scExp[-unique(ovrlps),]
            }
            if (verbose)
            {
                cat(
                    "ChromSCape::exclude_features_scExp - Removed",
                    length(unique(ovrlps)),
                    "regions from the analysis.\n"
                )
            }
        }
        if (by[1] == "feature_name")
        {
            if (has_genomic_coordinates(scExp))
            {
                warning(
                    "ChromSCape::exclude_features_scExp - Excluding by feature
                    name while object feature names are genomic coordinates !"
                )
            }
            features <- rownames(scExp)
            features_to_exclude <-
                as.character(features_to_exclude[, 1])
            ovrlps <-
                GenomicRanges::intersect(features, features_to_exclude)
            if (length(unique(ovrlps) > 0))
            {
                scExp <- scExp[-which(rownames(scExp) %in% ovrlps),]
            }
            if (verbose)
            {
                cat(
                    "ChromSCape::exclude_features_scExp - Removed",
                    length(unique(ovrlps)),
                    " features from the analysis.\n"
                )
            }
        }
        return(scExp)
    }

#' Preprocess scExp - Transcripts per Million (TPM)
#'
#' @param scExp
#'
#' @return A SingleCellExperiment object.
#' @importFrom GenomicRanges width
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom SingleCellExperiment counts normcounts
#' @importFrom Matrix t colSums
preprocess_TPM <- function(scExp)
{
    size <- GenomicRanges::width(SummarizedExperiment::rowRanges(scExp))
    SingleCellExperiment::normcounts(scExp) <-
        SingleCellExperiment::counts(scExp) / size
    SingleCellExperiment::normcounts(scExp) <- 10 ^ 6 *
        Matrix::t(Matrix::t(SingleCellExperiment::normcounts(scExp)) /
                      Matrix::colSums(normcounts(scExp)))
    return(scExp)
}

#' Preprocess scExp - Read per Kilobase Per Million (RPKM)
#'
#' @param scExp
#'
#' @return A SingleCellExperiment object.
#' @importFrom GenomicRanges width
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom SingleCellExperiment counts normcounts
#' @importFrom Matrix t colSums
preprocess_RPKM <- function(scExp)
{
    SingleCellExperiment::normcounts(scExp) <- 10 ^ 9 *
        Matrix::t(
            Matrix::t(SingleCellExperiment::counts(scExp)) /
                Matrix::colSums(SingleCellExperiment::counts(scExp))
        )
    size <-
        GenomicRanges::width(SummarizedExperiment::rowRanges(scExp))
    SingleCellExperiment::normcounts(scExp) <-
        SingleCellExperiment::normcounts(scExp) / size
    
    return(scExp)
}

#' Preprocess scExp - Counts Per Million (CPM)
#'
#' @param scExp
#'
#' @return A SingleCellExperiment object.
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom SingleCellExperiment counts normcounts
#' @importFrom Matrix t colSums
preprocess_CPM <- function(scExp)
{
    SingleCellExperiment::normcounts(scExp) <- 10 ^ 6 *
        Matrix::t(
            Matrix::t(SingleCellExperiment::counts(scExp)) /
                Matrix::colSums(SingleCellExperiment::counts(scExp))
        )
    return(scExp)
}

#' Preprocess scExp - size only
#'
#' @param scExp
#'
#' @return A SingleCellExperiment object.
#' @importFrom GenomicRanges width
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom SingleCellExperiment counts normcounts
#' @importFrom Matrix t colSums
preprocess_feature_size_only <- function(scExp)
{
    size <- GenomicRanges::width(SummarizedExperiment::rowRanges(scExp))
    SingleCellExperiment::normcounts(scExp) <-
        SingleCellExperiment::counts(scExp) / size
    return(scExp)
}

#' Normalize counts
#'
#' @param scExp A SingleCellExperiment object.
#' @param type Which normalization to apply. Either 'RPKM', 'CPM', 'TPM' or
#'  'feature_size_only'. Note that for all normalization by size
#'  (RPKM, TPM, feature_size_only), the features must have defined
#'  genomic coordinates.
#'
#' @return A SingleCellExperiment object containing normalized counts.
#'  (See ?normcounts())
#' @export
#'
#' @examples
normalize_scExp <- function(scExp,
                            type = c("RPKM", "CPM", "TPM",
                                     "feature_size_only"))
{
    stopifnot(
        type[1] %in% c("RPKM", "CPM", "TPM", "feature_size_only"),
        is(scExp, "SingleCellExperiment")
    )
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
        CPM = return(preprocess_CPM(scExp))
    )
}

#' Add gene annotations to features
#'
#' @param scExp A SingleCellExperiment object.
#' @param ref Reference genome. Either 'hg38' or 'mm10'. ['hg38']
#' @param reference_annotation A data.frame containing gene (or else) annotation with
#' genomic coordinates.
#'
#' @return A SingleCellExperiment object with annotated rowData.
#' @export
#'
#' @examples
#' @importFrom GenomicRanges makeGRangesFromDataFrame distanceToNearest
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom dplyr mutate select group_by summarise_all
feature_annotation_scExp <- function(scExp,
                                     ref = "hg38",
                                     reference_annotation = NULL)
{
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(ref))
    
    if (is.null(SummarizedExperiment::rowRanges(scExp)))
    {
        stop(
            "ChromSCape::feature_annotation_scExp - The object doesn't have
             ranges of coordinates as rowData"
        )
    }
    
    if (is.null(reference_annotation) &
        !(ref %in% c("hg38", "mm10")))
    {
        stop(
            "ChromSCape::feature_annotation_scExp - If reference_annotation is
        null, ref must be either 'hg38' or 'mm10' to automatically load
             reference gene annotation."
        )
    }
    
    if (is.null(reference_annotation) &
        (ref %in% c("hg38", "mm10")))
    {
        message(
            paste0(
                "ChromSCape::feature_annotation_scExp - Selecting ",
                ref,
                " genes from Gencode."
            )
        )
        eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
        reference_annotation <-
            eval(parse(text = paste0("", ref, ".GeneTSS")))
    }
    
    if (is.data.frame(reference_annotation))
    {
        reference_annotation <- GenomicRanges::makeGRangesFromDataFrame(reference_annotation,
                                                                        keep.extra.columns = T)
    }
    
    feature_ranges <- SummarizedExperiment::rowRanges(scExp)
    hits <-
        GenomicRanges::distanceToNearest(
            feature_ranges,
            reference_annotation,
            ignore.strand = T,
            select = "all"
        )
    
    annotFeat <- data.frame(
        chr = as.character(GenomicRanges::seqnames(feature_ranges[S4Vectors::queryHits(hits)])),
        start = as.character(GenomicRanges::start(feature_ranges
                                                  [S4Vectors::queryHits(hits)])),
        end = as.character(GenomicRanges::end(feature_ranges[S4Vectors::queryHits(hits)])),
        Gene = as.character(reference_annotation@elementMetadata$gene)[S4Vectors::subjectHits(hits)],
        distance = hits@elementMetadata$distance
    ) %>%
        dplyr::mutate(ID = paste(chr, start, end, sep = "_")) %>%
        dplyr::select(ID, chr, start, end, Gene, distance)
    
    system.time({
        annotFeat <- annotFeat %>% dplyr::group_by(ID, chr, start, end) %>%
            dplyr::summarise(Gene = paste(Gene,
                                          collapse = ", "),
                             distance = max(distance)) %>% as.data.frame()
    })
    
    annotFeat <- annotFeat[match(rownames(scExp), annotFeat$ID),]
    SummarizedExperiment::rowData(scExp) <- annotFeat
    return(scExp)
}

#' Choose perplexity depending on number of cells for Tsne
#'
#' @param dataset A matrix of features x cells (rows x columns)
#'
#' @return A number between 5 and 30 to use in Rtsne function
#'
choose_perplexity <- function(dataset)
{
    stopifnot(!is.null(dataset),!is.null(dim(dataset)))
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
#' [c('PCA','TSNE','UMAP')]
#' @param n Numbers of dimensions to keep for PCA. [50]
#' @param batch_correction Do batch correction ? [F]
#' @param batch_list List of characters. Names are batch names, characters are
#'  sample names.
#' @param verbose [T]
#'
#' @return A SingleCellExperiment object containing feature spaces. See
#'  ?reduceDims().
#'
#' @export
#'
#' @examples
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap umap.defaults
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment counts normcounts reducedDims
#' @importFrom Matrix t
reduce_dims_scExp <-
    function(scExp,
             dimension_reductions = c("PCA", "TSNE", "UMAP"),
             n = 50,
             batch_correction = F,
             batch_list = NULL,
             verbose = T)
    {
        stopifnot(
            is(scExp, "SingleCellExperiment"),
            is.numeric(n),
            dimension_reductions[1] %in%
                c("PCA", "TSNE", "UMAP")
        )
        
        if (!"normcounts" %in% names(SummarizedExperiment::assays(scExp)))
        {
            warning(
                "ChromSCape::reduce_dims_scExp - The raw counts are not
        normalized, running dimensionality reduction on raw counts."
            )
            mat <- SingleCellExperiment::counts(scExp)
        } else
        {
            mat <- SingleCellExperiment::normcounts(scExp)
        }
        
        if (batch_correction && !is.list(batch_list))
        {
            stop(
                "ChromSCape::reduce_dims_scExp - If doing batch correction,
        batch_list must be a list containing the samples IDs of each batch."
            )
        }
        
        batches <- list()
        
        if (batch_correction)
        {
            print("Running Batch Correction ...")
            num_batches <- length(batch_list)
            batch_names <- names(batch_list)
            
            scExp$batch_name <- "unspecified"
            
            for (i in 1:num_batches)
            {
                for (s_id in batch_list[[i]])
                {
                    SummarizedExperiment::colData(scExp)[scExp$sample_id == s_id,
                                                         "batch_name"] <- batch_names[i]
                }
            }
            
            adj_annot <- data.frame()
            b_names <- unique(scExp$batch_name)
            
            if (class(mat) %in% c("dgCMatrix", "dgTMatrix"))
            {
                pca <- pca_irlba_for_sparseMatrix(Matrix::t(mat), n)
            } else
            {
                pca <- stats::prcomp(Matrix::t(mat),
                                     center = T,
                                     scale. = F)
                pca <- pca$x[, 1:n]
            }
            for (i in 1:length(b_names))
            {
                b_name <- b_names[i]
                batches[[i]] <- pca[scExp$batch_name == b_name,]
                adj_annot <- rbind(adj_annot,
                                   SummarizedExperiment::colData(scExp)[scExp$batch_name ==
                                                                            b_name,])
            }
            
            mnn.out <-
                do.call(scran::fastMNN, c(
                    batches,
                    list(
                        k = 25,
                        d = 50,
                        ndist = 3,
                        pc.input = T,
                        auto.order = T,
                        cos.norm = F,
                        compute.variances = T
                    )
                ))
            pca <- mnn.out$corrected
            SummarizedExperiment::colData(scExp) <- adj_annot
        } else
        {
            scExp$batch_id <- "batch_1"
            if (class(mat) %in% c("dgCMatrix", "dgTMatrix"))
            {
                pca <- pca_irlba_for_sparseMatrix(Matrix::t(mat), n)
            } else
            {
                pca <- stats::prcomp(Matrix::t(mat),
                                     center = T,
                                     scale. = F)
                pca <- pca$x[, 1:n]
            }
        }
        pca <- as.data.frame(as.matrix(pca))
        colnames(pca) <- paste0("Component_", 1:n)
        rownames(pca) <- colnames(scExp)
        
        set.seed(47)
        # Reduce the perplexity if the number of samples is too low to avoid perplexity
        # error
        if ("TSNE" %in% dimension_reductions)
        {
            tsne <- Rtsne::Rtsne(
                pca,
                dims = 2,
                pca = FALSE,
                theta = 0,
                perplexity = choose_perplexity(pca),
                verbose = verbose,
                max_iter = 1000
            )
            tsne <- as.data.frame(tsne$Y)
            colnames(tsne) <- c("Component_1", "Component_2")
        }
        
        if ("UMAP" %in% dimension_reductions)
        {
            config <- umap::umap.defaults
            config$metric <- "cosine"
            umap <- umap::umap(pca, config = config, method = "naive")
            umap <- as.data.frame(umap$layout)
            colnames(umap) <- c("Component_1", "Component_2")
        }
        
        # save PCA & T-SNE in scExp object
        listReducedDim <- list(PCA = pca)
        if ("TSNE" %in% dimension_reductions)
        {
            listReducedDim$TSNE <- tsne
        }
        if ("UMAP" %in% dimension_reductions)
        {
            listReducedDim$UMAP <- umap
        }
        SingleCellExperiment::reducedDims(scExp) <- listReducedDim
        
        return(scExp)
    }

#' Run sparse PCA using irlba SVD
#'
#' @param x
#' @param n_comp
#'
#' @return
#'
#' @importFrom irlba irlba
#' @importFrom Matrix colMeans
#'
#' @examples
pca_irlba_for_sparseMatrix <- function(x, n_comp)
{
    system.time({
        x.means <- Matrix::colMeans(x)
        svd.0 <- irlba::irlba(x, center = x.means, nv = n_comp)
        x. <- sweep(x, 2, x.means, "-")
        pca <- x. %*% svd.0$v
    })
    return(pca)
}

#' Table of cells before / after QC
#'
#' @param scExp A SingleCellExperiment object.
#' @param annot A raw annotation data.frame of cells before filtering.
#'
#' @export
#' @return A formatted kable in HTML.
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr bind_rows tibble left_join
#' @importFrom kableExtra kable kable_styling group_rows
#'
num_cell_after_QC_filt_scExp <- function(scExp, annot)
{
    stopifnot(is(scExp, "SingleCellExperiment"),!is.null(annot))
    
    table <- as.data.frame(table(annot$sample_id))
    table_filtered <- as.data.frame(table(SingleCellExperiment::colData(scExp)$sample_id))
    
    colnames(table) <- c("Sample", "#Cells Before Filtering")
    rownames(table) <- NULL
    colnames(table_filtered) <-
        c("Sample", "#Cells After Filtering")
    rownames(table_filtered) <- NULL
    
    table_both <-
        dplyr::left_join(table, table_filtered, by = c("Sample"))
    table_both[, 1] <- as.character(table_both[, 1])
    table_both <-
        table_both %>% dplyr::bind_rows(
            .,
            dplyr::tibble(
                Sample = "",
                `#Cells Before Filtering` = sum(table_both[,
                                                           2]),
                `#Cells After Filtering` = sum(table_both[, 3])
            )
        )
    
    table_both %>% kableExtra::kable(escape = F, align = "c") %>%
        kableExtra::kable_styling(c("striped",
                                    "condensed"), full_width = T) %>%
        kableExtra::group_rows("Total cell count",
                               dim(table_both)[1], dim(table_both)[1])
}
