## Authors : Pacôme Prompsy Title : Wrappers & functions to preprocess & reduce
## dimensionality of single-cell Matrix Description : Funtions to load, filter,
## normalize & reduce single-cell Epigenetic Matrices prior to analysis

#' Heuristic discovery of samples based on cell labels 
#'
#' @param barcodes 
#' @param nb_samples 
#'
#' @return character vector of sample names the same length as cell labels
#' @export
#'
#' @importFrom stringdist stringdistmatrix
#' @importFrom qualV LCS
#' @examples
detect_samples <- function(barcodes, nb_samples=1){
  
  t = system.time({
    mat = stringdist::stringdistmatrix(barcodes,barcodes,useNames="strings", method = "lv",weight = c(1,1,1))
    dist = stats::as.dist(1-stats::cor(mat))
    hc = stats::hclust(dist,method = 'ward.D')
    hc$labels=rep("",ncol(mat))
  })
  cat("ChromSCape::detect_samples - found samples in",t[3], "secs.")
  sample_groups = stats::cutree(hc,k = nb_samples)
  samp_name = c()
  
  for(i in 1:nb_samples){
    x = barcodes[which(sample_groups==i)]
    set.seed(47)
    samp = as.character(sample(x,min(length(x),10),replace = T))
    
    lcs = strsplit(samp[1], "")[[1]]
    for(j in 2:min(length(x),10)){
      lcs = unlist(qualV::LCS(lcs, strsplit(samp[j], "")[[1]])[4])
    }
    samp_name = c(samp_name,paste(lcs,collapse = ""))   
  }
  samples_names = samp_name
  
  # remove trailing special char
  samples_names = gsub("[[:punct:]]$", "", samples_names)
  # remove starting special char
  samples_names = gsub("^[[:punct:]]", "", samples_names)
  
  # remove longest common string
  if(nb_samples > 1){
    mat = matrix("",nrow=nb_samples,ncol=nb_samples)
    for(i in 1:(nb_samples)){
      for(j in 1:(nb_samples)){
        x = as.character(samples_names[i])
        y = as.character(samples_names[j])
        
        x = strsplit(x, "")[[1]]
        y = strsplit(y, "")[[1]]
        
        lcs = qualV::LCS(x,y)
        ls = split(lcs[6]$va, cumsum(c(TRUE, diff(lcs[6]$va) != 1)))
        ls =ls[[which.max(lengths(ls))]]
        
        x = as.character(samples_names[i])
        x = strsplit(x, "")[[1]]
        lcs = paste0(x[ls],collapse = "")
        
        if(length(lcs)==0) lcs =""
        mat[i,j] = lcs
      }
    }
    
    mat = as.character(mat)
    longest_common_between_samples = mat[which.min(unlist(lapply(mat, nchar)))][1]
    if(length(longest_common_between_samples)>0){
      if(length(grep(pattern = longest_common_between_samples, samples_names)) 
         == nb_samples) samples_names = gsub(longest_common_between_samples,"",
                                             samples_names)
    }
  }
  samples_names = samples_names[sample_groups]
  
  return(samples_names)
}

#' split_bam_file_into_single_cell_bams
#'
#' @param bam_file 
#' @param sample_id 
#' @param outdir 
#' @param min_coverage 
#' @param barcode_tag 
#' @param verbose 
#'
#' @return
#'
#' @importFrom BiocParallel bpparam
#' @examples
split_bam_file_into_single_cell_bams <- function(bam_file,
                                                 sample_id,
                                                 outdir,
                                                 min_coverage=500,
                                                 barcode_tag="XB",
                                                 verbose = T)
{
    stopifnot(
        file.exists(bam_file),
        dir.exists(outdir),
        is.character(sample_id),
        is.numeric(min_coverage),
        is.character(barcode_tag)
        )
    if(length(barcode_tag)>2) stop("ChromSCape::split_bam_file_into_single_cell_bams 
                                   The barcode tag ID should be two characters (e.g. 'XB')")
    try({
        samtools = system2("which",args = "samtools",stdout = TRUE)
    })
    
    if(length(samtools)==0) stop("samtools is not found in path. Reinstall samtools or
                                 launch application from terminal.")
    num_cores = BiocParallel::bpparam()
    num_cores = as.numeric(num_cores@.xData$workers)
    
    system2(command = paste0("bash ",file.path(system.file(package = "ChromSCape"),
                                               "split_bam.sh"))
                             , args =paste0(bam_file,sample_id,
                                            outdir,min_coverage,
                                            barcode_tag, num_cores))

}



#' Transforms a peaks x cells count matrix into a bins x cells count matrix. 
#' 
#' This functions is best used to re-count large number of small peaks (e.g. <= 5000bp)
#' into equal or larger bins.
#' The genome is either cut in fixed bins (e.g. 50,000bp) or into an user defined 
#' number of bins. Bins are calculated based on the canconical chromosomes. Note
#' that if peaks are larger than bins, or if peaks are overlapping multiple bins, 
#' the signal is added to each bin. Users can increase the minimum overlap to 
#' consider peaks overlapping bins (by default 150bp, size of a nucleosome) to 
#' disminish the number of peaks overlapping multiple region. Any peak smaller 
#' than the minimum overlapp threshold will be dismissed. Therefore, library size
#' might be slightly different from peaks to bins if signal was duplicated into 
#' multiple bins or ommitted due to peaks smaller than minimum overlap.
#'
#' @param ref reference genome to use [hg38]
#' @param verbose 
#' @param files_dir 
#' @param peak_file 
#' @param n_bins 
#' @param bin_width 
#' @param geneTSS 
#' @param aroundTSS 
#' @param file_type 
#'
#' @return
#' @export
#'
#' @importFrom IRanges IRanges
#' @importFrom parallel detectCores mclapply
#' @importFrom GenomicRanges GRanges tileGenome width seqnames GRangesList
#' sort.GenomicRanges
raw_counts_to_feature_count_files <- function(files_dir,
                                              file_type = c("BAM","BED","Index_Peak_Barcode"),
                                              peak_file = NULL,
                                              n_bins = NULL,
                                              bin_width = NULL,
                                              geneTSS = NULL,
                                              aroundTSS = 2500,
                                              verbose = T,
                                              ref = "hg38"){
    stopifnot(
        dir.exists(files_dir),
        is.numeric(aroundTSS),
        ref %in% c("mm10", "hg38")
    )
    
    if(!is.null(peak_file) && !file.exists(peak_file)) 
        stop("ChromSCape::raw_counts_to_feature_count_files - Can't find peak file.")
    
    if(!is.null(n_bins) && !is.numeric(n_bins)) 
        stop("ChromSCape::raw_counts_to_feature_count_files - n_bins must be a number.")
    
    if(!is.null(bin_width) && !is.numeric(bin_width)) 
        stop("ChromSCape::raw_counts_to_feature_count_files - bin_width must be a number.")
    
    if(!is.null(geneTSS) && !is.logical(geneTSS)) 
        stop("ChromSCape::raw_counts_to_feature_count_files - geneTSS must be a TRUE or FALSE")
    
    if(file_type == "Index_Peak_Barcode"){
      peak_file_2 = list.files(path = files_dir,full.names = T, pattern = ".*peaks.bed$")
      index_file = list.files(path = files_dir,full.names = T, pattern = ".*index.txt$")
      barcode_file = list.files(path = files_dir,full.names = T, pattern = ".*barcodes.txt$")
      if(length(c(peak_file_2,index_file,barcode_file))!=3) stop("ChromSCape::raw_counts_to_feature_count_files - For Index Count type, the folder 
                 must contain exactly two files matching respectively [*index.txt$],[*peaks.bed$],
               [*barcodes.txt$]")
    }

    type = length(which(c(!is.null(peak_file),
                   !is.null(n_bins),
                   !is.null(bin_width),
                   !is.null(geneTSS))))
   
    if(type>1 | type==0 ) {
        stop("ChromSCape::raw_counts_to_feature_count_files - peak_file, n_bins,
    bin_width & geneTSS are required and mutually exclusive.")
    }
    
    # Create genomic ranges of bins
    eval(parse(text = paste0(
        "chr <- ChromSCape::", ref, ".chromosomes"
    )))
    
    if(!is.null(peak_file)){
        cat('ChromSCape::raw_counts_to_feature_count_files - Reading in peaks file...\n')
        features = read.table(peak_file,sep = "\t",
                              header = F,quote = "")
        features = features[which(features$V1 %in% chr$chr),]
        features$V1 = as.character(features$V1)
        which <- GenomicRanges::GRanges(seqnames = features$V1,ranges = IRanges::IRanges(features$V2, features$V3)) 
    } else if(!is.null(n_bins) | !is.null(bin_width)){
        cat('ChromSCape::raw_counts_to_feature_count_files - Counting on genomic bins...\n')
        # load species chromsizes
        chr <- GenomicRanges::GRanges(chr)
        if(!is.null(n_bins)){
            which <- unlist(GenomicRanges::tileGenome(
                setNames(
                    GenomicRanges::width(chr),
                    GenomicRanges::seqnames(chr)
                ),
                ntile = n_bins
            ))
        } else if(!is.null(bin_width)){
            which <- unlist(GenomicRanges::tileGenome(
                setNames(
                    GenomicRanges::width(chr),
                    GenomicRanges::seqnames(chr)
                ),
                tilewidth = bin_width
            ))
        }
    } else if(geneTSS == TRUE) {
        # Retrieve gene TSS from ref and create GRanges
        eval(parse(text = paste0(
            "geneTSS_df <- ChromSCape::", ref, ".GeneTSS"
        )))
        geneTSS_df$start = geneTSS_df$start - aroundTSS
        geneTSS_df$end = geneTSS_df$end + aroundTSS
        which = GenomicRanges::GRanges(geneTSS_df)
    }

    if(file_type == "BAM"){
        t1 = system.time({l = bams_to_matrix_indexes(files_dir, which)})
        feature_indexes = l[[1]]
        name_cells = l[[2]]
        if(verbose) cat("ChromSCape::raw_counts_to_feature_count_files - Count matrix created from BAM files in",t1[3],"sec.\n")
    }
    else if(file_type == "BED"){
        t1 = system.time({l = beds_to_matrix_indexes(files_dir, which) })
        feature_indexes = l[[1]]
        name_cells = l[[2]]
        if(verbose) cat("ChromSCape::raw_counts_to_feature_count_files - Count matrix created from BED files in",t1[3],"sec.\n")
    }
    else if(file_type == "Index_Peak_Barcode"){
        name_cells = read.table(barcode_file,sep="",quote="",stringsAsFactors = F)[,1]
        t1 = system.time({
            l = index_peaks_barcodes_to_matrix_indexes(peak_file = peak_file_2,
                                                       index_file = index_file,
                                                       name_cells = name_cells,
                                                       ref = ref)
        })
        feature_indexes=l[[1]]
        which = l[[2]]
        
        if(verbose) cat("ChromSCape::raw_counts_to_feature_count_files - Count matrix created from Index-Peak-Barcodes files in",t1[3],"sec.\n")
    }
    
    mat = Matrix::sparseMatrix(i = as.numeric(feature_indexes$feature_index),
                               j = feature_indexes$barcode_index,
                               x = feature_indexes$counts,
                               dims = c(length(which),length(name_cells)),
                               dimnames = list(rows=gsub(":|-","_",as.character(which)),
                                               cols=name_cells))
    
    
    return(mat)
}

#' Count bam files on interval to create count indexes
#'
#' @return feature_indexes
#'
#' @importFrom Rsamtools BamFileList indexBam ScanBamParam countBam
#' @importFrom BiocParallel bplapply
#' @examples
bams_to_matrix_indexes = function(files_dir,which){
    single_cell_bams = list.files(files_dir,full.names = T,pattern = paste0(".*.bam$"))
    bam_files = Rsamtools::BamFileList(single_cell_bams)
    name_cells = gsub("*.bam$","",basename(as.character(single_cell_bams)))
    names(bam_files) = name_cells

    indexes = list.files(files_dir,full.names = T,pattern = paste0(".*.bai$"))
    if(length(indexes) < length(bam_files)){
      cat("ChromSCape::bams_to_matrix_indexes - indexing BAM files...")
      BiocParallel::bplapply(bam_files,Rsamtools::indexBam)
      indexes = list.files(files_dir,full.names = T,pattern = paste0(".*.bai$"))
    }
    names_index = gsub("*.bam.bai$","",basename(as.character(indexes)))
    names(indexes) = names_index
    if(!all.equal(names_index,name_cells)) stop("Different number of BAM and indexes files. Stopping.")
    param = Rsamtools::ScanBamParam(which = which)
    system.time({
        feature_list = BiocParallel::bplapply(names(bam_files), function(bam_name){
            bam_files[[bam_name]]$index = indexes[[bam_name]]
            tmp = Rsamtools::countBam(file = bam_files[[bam_name]], param=param)
            tmp$feature_index = rownames(tmp)
            tmp$cell_id = bam_name
            sel=which(tmp$records>0)
            if(length(sel)>0) tmp = tmp[sel,c("cell_id","feature_index","records")]
            tmp
        }) 
    })
    
    feature_indexes = do.call(rbind,feature_list)
    feature_indexes$barcode_index = as.numeric(as.factor(feature_indexes$cell_id))
    colnames(feature_indexes)[3] = "counts"
    return(list(feature_indexes,name_cells))
}

#' Read index-peaks-barcodes trio files on interval to create count indexes
#'
#' @return feature_indexes
#' 
#' @importFrom GenomicRanges GRanges
#' @examples
index_peaks_barcodes_to_matrix_indexes = function(peak_file,
                                                  index_file,
                                                  name_cells,
                                                  binarize = F,
                                                  ref = "hg38"){
    
    regions = GenomicRanges::GRanges(setNames(read.table(peak_file,sep="\t",quote="")[,1:3],
                               c("chr","start","end")))
    eval(parse(text = paste0(
        "chr <- ChromSCape::", ref, ".chromosomes"
    )))
    
    regions = regions[which(as.character(regions@seqnames) %in% chr$chr),]
    feature_indexes = setNames(read.table(index_file,sep="",quote="")[,1:3],
                       c("feature_index","barcode_index","counts"))
    
    
    feature_indexes$cell_id = name_cells[feature_indexes$barcode_index]
    feature_indexes = feature_indexes[,c("cell_id","feature_index","counts","barcode_index")]
    return(list(feature_indexes,regions))
}

#' Count bed files on interval to create count indexes
#'
#' @return feature_indexes
#' @importFrom BiocParallel bplapply
#' @importFrom GenomicRanges GRanges countOverlaps
#' 
#' @examples
beds_to_matrix_indexes = function(files_dir,which){
    single_cell_beds = list.files(files_dir,full.names = T,pattern = ".*.bed$|.*.bed.gz$")
    names_cells = gsub(".bed$|.bed.gz$","",basename(as.character(single_cell_beds)))
    names(single_cell_beds) = names_cells
    
    system.time({
        feature_list = BiocParallel::bplapply(names(single_cell_beds), function(bed_name){
            
            gzipped = (length(grep(".gz$",single_cell_beds[[bed_name]])) > 0)
            if(gzipped) {
                bed = read.table(gzfile(single_cell_beds[[bed_name]]),
                           sep="\t",quote = "",stringsAsFactors = F)
            }else {bed = read.table(single_cell_beds[[bed_name]], sep="\t", 
                                  quote = "",stringsAsFactors = F)}
            bed = GenomicRanges::GRanges(setNames(bed[,1:3],c("chr","start","end")))       
            suppressWarnings({
              tmp = data.frame(counts=GenomicRanges::countOverlaps(which,bed, minoverlap = 1))
            })
            
            tmp$feature_index = as.numeric(rownames(tmp))
            tmp = tmp[which(tmp$counts>0),]
            tmp$cell_id = bed_name
            tmp = tmp[,c("cell_id","feature_index","counts")]
            tmp
        }) 
    })
    gc()
    feature_indexes = do.call(rbind,feature_list)
    feature_indexes$barcode_index = as.numeric(as.factor(feature_indexes$cell_id))
    return(list(feature_indexes,names_cells))
}

#' Transforms a peaks x cells count matrix into a bins x cells count matrix. 
#' 
#' This functions is best used to re-count large number of small peaks (e.g. <= 5000bp)
#' into equal or larger bins.
#' The genome is either cut in fixed bins (e.g. 50,000bp) or into an user defined 
#' number of bins. Bins are calculated based on the canconical chromosomes. Note
#' that if peaks are larger than bins, or if peaks are overlapping multiple bins, 
#' the signal is added to each bin. Users can increase the minimum overlap to 
#' consider peaks overlapping bins (by default 150bp, size of a nucleosome) to 
#' disminish the number of peaks overlapping multiple region. Any peak smaller 
#' than the minimum overlapp threshold will be dismissed. Therefore, library size
#' might be slightly different from peaks to bins if signal was duplicated into 
#' multiple bins or ommitted due to peaks smaller than minimum overlap.
#'
#' @param mat A matrix of peaks x cells 
#' @param bin_width width of bins to produce in base pairs (minimum 500) [50000]
#' @param ref reference genome to use [hg38]
#' @param n_bins number of bins (exclusive with bin_width)
#' @param minoverlap Minimum overlap between a peak and a bin to consider the 
#' peak as overlapping the bin [150].
#' @param verbose 
#'
#' @return
#' @export
#'
#' @importFrom IRanges IRanges
#' @importFrom BiocParallel bpaggregate
#' @importFrom GenomicRanges GRanges tileGenome width seqnames GRangesList
#' sort.GenomicRanges findOverlaps
peaks_to_bins <- function(mat,
                          bin_width = 50000,
                          n_bins = NULL,
                          minoverlap = 150,
                          verbose = T,
                          ref = "hg38")
{
    stopifnot(
        !is.null(mat),
        ref %in% c("mm10", "hg38")
    )
    # Render sparse
    if(is.matrix(mat)) mat = as(mat,"dgCMatrix")
    if(is.null(n_bins) & is.null(bin_width)) 
        stop("One of bin_width or n_bins must be set")
    
    # Create genomic ranges of bins
    eval(parse(text = paste0(
        "chr <- ChromSCape::", ref, ".chromosomes"
    )))
    # load species chromsizes
    chr <- GenomicRanges::GRanges(chr)
    
    if(!is.null(n_bins)){
        bin_ranges <- unlist(GenomicRanges::tileGenome(
            setNames(
                GenomicRanges::width(chr),
                GenomicRanges::seqnames(chr)
            ),
            ntile = n_bins
        ))
    } else {
        bin_ranges <- unlist(GenomicRanges::tileGenome(
            setNames(
                GenomicRanges::width(chr),
                GenomicRanges::seqnames(chr)
            ),
            tilewidth = bin_width
        ))
    }
    
    peaks = rownames(mat)
    peaks_chr = sapply(strsplit(peaks,":|-|_"),function(x) x[1])
    peaks_start = as.numeric(sapply(strsplit(peaks,":|-|_"),function(x) x[2]))
    peaks_end = as.numeric(sapply(strsplit(peaks,":|-|_"),function(x) x[3]))
    
    if(anyNA(peaks_start) | anyNA(peaks_end))
        stop("The rows of mat should be regions in format chr_start_end or 
             chr:start-end, without non canonical chromosomes")
    if(verbose) {
        cat(
            "ChromSCape::peaks_to_bins - converting",
            dim(mat)[1],
            "peaks into",length(bin_ranges)[1], "bins of",
            mean(bin_ranges@ranges@width), "bp in average.\n"
        )
    }
    peaks = GenomicRanges::GRanges(seqnames = peaks_chr,
                                   ranges = IRanges::IRanges(peaks_start, peaks_end))
    
    hits <- GenomicRanges::findOverlaps(bin_ranges,peaks,minoverlap = minoverlap)
    bins_names = paste0(bin_ranges@seqnames, ":",
                        GenomicRanges::start(bin_ranges), "-",
                        GenomicRanges::end(bin_ranges))
    
    bin_mat = NULL
    hits = as.matrix(hits)
    hits = cbind(hits, mat[hits[,"subjectHits"],])
    print("Running aggregation of peaks to bins in parallel")
    hits = as.matrix(hits)
    system.time({
    bin_mat = BiocParallel::bpaggregate(x = hits[,3:dim(hits)[2]], 
                                        by = list(bins = hits[,1,drop=F]), 
                                        FUN = sum)
    })
    
    bin_mat = as(as.matrix(bin_mat[,2:ncol(bin_mat)]),"dgCMatrix")
    rownames(bin_mat) = bins_names[unique(hits[,1])]
    gc()
    
    if(verbose){
        cat(
            "ChromSCape::peaks_to_bins - removed",
            length(bin_ranges) - nrow(bin_mat),
            "empty bins from the binned matrix.\n"
        )
    }
    if(length(empty_bins)>0) bin_mat=bin_mat[-empty_bins,]
    
}

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


#' Read single-cell matrix(ces) into scExp
#'
#' Combine one or multiple matrices together to create a sparse matrix and cell
#' annotation data.frame.
#' 
#' @param file_names A character vector of file names towards single cell 
#' epigenomic matrices (features x cells) (must be .txt / .tsv)
#' @param path_to_matrix In case matrices are stored in temporary folder, 
#' a character vector of path towards temporary files. [NULL]
#' @param ref Reference genome. Either 'hg38' or 'mm10'. (only matters if
#'  rows are genes)['hg38']
#'
#' @return a list containing spa
#'  * datamatrix: a sparseMatrix of features x cells
#'  * annot_raw: an annotation of cells as data.frame
#'  
#' @export
#'
#' @importFrom scater readSparseCounts
#' @examples
import_scExp <- function(file_names,
                         path_to_matrix = NULL,
                         ref = "hg38"
                         ){
    
    stopifnot(is.character(file_names))

    if(length(grep("(.tsv$)|(.txt$)|(.csv$)",file_names)) < length(file_names))
        stop(paste0("ChromSCape::import_scExp - Matrix files must be in .txt or .tsv format."))
    
    if(is.null(path_to_matrix)) path_to_matrix = file_names
    
    if( FALSE %in% sapply(path_to_matrix,file.exists)) 
        stop(paste0("ChromSCape::import_scExp - can't find one of the matrix files."))
    
    datamatrix = annot_raw = NULL
    
    for(i in 1:length(file_names)){

        sample_name <- gsub('.{4}$', '', basename(file_names[i]))
        
        format_test = as.character(read.table(path_to_matrix[i],
                                 header = F,
                                 sep = "\t", nrows = 5)[4,])
        
        if(length(format_test)>3) separator = "\t" else separator = ","
        
        format_test = read.table(path_to_matrix[i],
                                              header = T,
                                              sep = separator, nrows = 5)
        separated_chr_start_end = c(grep("chr",colnames(format_test)[1:3]),
          grep("start|begin", colnames(format_test)[1:3]),
          grep("end|stop",colnames(format_test)[1:3]))
        
        
        if(length(separated_chr_start_end)>0 && all.equal(separated_chr_start_end,c(1,2,3))){
          
          val = c("NULL","NULL","NULL",
                  rep("integer", length(colnames(format_test))-3))
          datamatrix_single = as(as.matrix(read.table(path_to_matrix[i],
                                                   header = T,
                                                   sep = separator,
                                                   colClasses = val)), "dgCMatrix")
          gc()
          
          val = c("character","integer","integer",rep("NULL", ncol(datamatrix_single)))
          regions = read.table(path_to_matrix[i],
                               header = T,
                               sep = separator,
                               colClasses = val)
          regions = paste0(regions[,1],":",regions[,2],"-",regions[,3])
          rownames(datamatrix_single) = regions
        
        } else{
          datamatrix_single <- scater::readSparseCounts(path_to_matrix[i],
                                                      sep=separator,
                                                      chunk = 1000L)
        }
        gc()
        #perform some checks on data format
        matchingRN <- grep("[[:alnum:]]+(:|_)[[:digit:]]+(-|_)[[:digit:]]+",
                           rownames(datamatrix_single)) # check rowname format
        
        if(length(matchingRN) < length(rownames(datamatrix_single))){
          eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
          reference_annotation <-
            eval(parse(text = paste0("", ref, ".GeneTSS")))
          matching_genes = intersect(rownames(datamatrix_single),
                                     reference_annotation$gene)
          if(length(matching_genes) < length(rownames(datamatrix_single))){
            warning(paste0("ChromSCape::import_scExp - ",sample_name, " contains ",
                           (length(rownames(datamatrix_single))-length(matchingRN)),
                           " rownames that do not conform to the required format.
                           Must be either genomic coordinates (chr:start-end) or 
                           conform gene name, check data(hg38.GeneTSS) or
                           data(mm10.GeneTSS).
                           Please check your data matrix and try again."))
            return()
          }
        }
        
        numericC <- apply(datamatrix_single[1:5,1:5],MARGIN=2,is.numeric) # check if matrix is numeric
        if(sum(numericC) < 5){
            stop(paste0("ChromSCape::import_scExp - ",sample_name,
                        " contains non-numeric columns at the following indices: ",
                        which(numericC==FALSE), ". Please check your data matrix and try again."))
        }
        
        if(all.equal(rownames(datamatrix_single)[1:5], c("1","2","3","4","5")) == TRUE){
            warning(paste0("ChromSCape::import_scExp - ", sample_name,
                        " contains numeric rownames, taking first column as rownames."))
            names = datamatrix_single$X0
            datamatrix_single = datamatrix_single[,-1]
            rownames(datamatrix_single) = names
        }
        
        datamatrix_single <- datamatrix_single[!duplicated(rownames(datamatrix_single)),] 
        
        if(length(grep("chr",rownames(datamatrix_single)[1:10],perl = T)) >= 9){
            rownames(datamatrix_single) <- gsub(":", "_", rownames(datamatrix_single))
            rownames(datamatrix_single) <- gsub("-", "_", rownames(datamatrix_single))
        }
        gc()
        total_cell <- length(datamatrix_single[1,])
        annot_single <- data.frame(barcode = colnames(datamatrix_single),
                                   cell_id = paste0(sample_name, "_", colnames(datamatrix_single)),
                                   sample_id = rep(sample_name, total_cell),
                                   batch_id = i)
        colnames(datamatrix_single) <- annot_single$cell_id
        if(is.null(datamatrix)){
            datamatrix <- datamatrix_single
        }else{
            common_regions <- intersect(rownames(datamatrix), rownames(datamatrix_single))

            if(length(common_regions)>0){
                datamatrix <- Matrix::cbind2(datamatrix[common_regions,], datamatrix_single[common_regions,])
            } else {
               stop(paste0("ChromSCape::import_scExp - ",file_names[i], " contains no common regions with ",file_names[i-1])) 
            }
            if(length(common_regions) < nrow(datamatrix)){
              warning(paste0("ChromSCape::import_scExp - ",file_names[i],
                             " contains less than ",
                             ceiling(100 * length(common_regions) / (nrow(datamatrix)))
                             ," common regions with ",file_names[i-1])) 
            }
            
        }
        rm(datamatrix_single)
        gc()
        if(is.null(annot_raw)){
            annot_raw <- annot_single
        } else {
            annot_raw <- rbind(annot_raw, annot_single)
        }
        rm(annot_single)
        gc()
    }
    out = list("datamatrix" = datamatrix, "annot_raw" = annot_raw)
    return(out)
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
#' @param ref Reference genome. Either 'hg38' or 'mm10'. (only matters if rows are genes)['hg38']
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
                         ref = "hg38",
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
    
    cat(
        "ChromSCape::create_scExp - the matrix has",
        dim(datamatrix)[2],
        "cells and ",dim(datamatrix)[1], "features.\n"
    )
    
    scExp <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = datamatrix), colData = annot)
    
    
    
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
    } else { # Row content is genes
      if (ref %in% c("hg38", "mm10"))
      {
        message(
          paste0(
            "ChromSCape::create_scExp - Selecting ",
            ref,
            " genes coordinates from Gencode."
          )
        )
        eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
        reference_annotation <-
          eval(parse(text = paste0("", ref, ".GeneTSS")))
      } else {
        stop("ChromSCape::create_scExp - For genes as rows matrices, please select
           a valid reference genome (hg38/mm10)")
      }
      
      genes = intersect(reference_annotation$gene, rownames(scExp))
      print(head(genes,15))
      print(head(reference_annotation$gene,15))
      print(head(rownames(scExp),15))
      if(length(genes) < 10){
        stop("ChromSCape::create_scExp - Error - The rownames of your matrix are
          recognized as being genes, but are not present in the reference genome
               you specified.")
      }

      chr = as.character(reference_annotation$chr[match(
         rownames(scExp),reference_annotation$gene)])
      start = as.character(reference_annotation$start[match(
         rownames(scExp),reference_annotation$gene)])
      end = as.character(reference_annotation$end[match(
         rownames(scExp),reference_annotation$gene)])
      
      rownames(scExp) = paste0(chr,"_",start,"_",end)
      
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
      scExp <- scExp[!duplicated(rownames(scExp)),] 
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
    
    suppressWarnings({scExp <- scater::calculateQCMetrics(scExp)})
    
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
filter_scExp =  function (scExp, min_cov_cell = 1600, quant_removal = 95, percentMin = 1, 
                          bin_min_count = 2, verbose = T) 
{
    stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(min_cov_cell), 
              is.numeric(quant_removal), is.numeric(percentMin), is.numeric(bin_min_count), 
              verbose %in% c(F, T))
    if (is.null(scExp)) {
        warning("ChromSCape::filter_scExp -\n                Please specify a SingleCellExperiment")
    }
    cellCounts <- Matrix::colSums(SingleCellExperiment::counts(scExp))
    thresh <- stats::quantile(cellCounts, probs = seq(0, 1, 0.01))
    sel1000 <- (cellCounts > 1000 & cellCounts <= thresh[quant_removal + 
                                                             1])
    sel <- (cellCounts > min_cov_cell & cellCounts <= thresh[quant_removal + 
                                                                 1])
    if (verbose) {
        cat("ChromSCape::filter_scExp -", length(which(sel)), 
            "cells pass the threshold of", min_cov_cell, "minimum reads and are lower than the ", 
            quant_removal, "th centile of library size ~=", round(thresh[quant_removal + 
                                                                             1]), "reads.\n")
    }
    rm(cellCounts)
    gc()
    bina_counts <- SingleCellExperiment::counts(scExp)[, sel1000]
    gc()
    sel_above_2 <- (bina_counts >= bin_min_count)
    gc()
    bina_counts = Matrix::sparseMatrix(i = 1,j = 1,x = 1,dims = dim(bina_counts))
    bina_counts[1,1] = 0
    gc()
    bina_counts[sel_above_2] <- 1
    gc()
    nCells_in_feature <- Matrix::rowSums(bina_counts)
    gc()
    fixedFeature <- which(nCells_in_feature > ((percentMin/100) * 
                                                   (ncol(bina_counts))))
    if (verbose) {
        cat("ChromSCape::filter_scExp -", length(fixedFeature), 
            "features pass the threshold of", percentMin, " % of total cells 'ON', representing a minimum of", 
            round((percentMin/100) * (ncol(bina_counts))), "cells.\n")
    }
    scExp <- scExp[, sel]
    scExp <- scExp[fixedFeature, ]
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
            suppressWarnings({ovrlps <-
                as.data.frame(GenomicRanges::findOverlaps(regions, excl_gr))[,
                                                                             1]})
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
            suppressWarnings({ovrlps <-
                GenomicRanges::intersect(features, features_to_exclude)})
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
    
    feature_ranges <- SummarizedExperiment::rowRanges(scExp)[,1]
    suppressWarnings({
      hits <-
        GenomicRanges::distanceToNearest(
            feature_ranges,
            reference_annotation,
            ignore.strand = T,
            select = "all"
        )
      })
    
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
#' @importFrom batchelor fastMNN
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
        rotation = NULL
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
                out <- pca_irlba_for_sparseMatrix(Matrix::t(mat), n)
                pca = out[["PCA"]]
                rotation = out[["rotation"]]
            } else
            {
                pca <- stats::prcomp(Matrix::t(mat),
                                     center = T,
                                     scale. = F)
                rotation = pca$rotation
                pca <- pca$x[, 1:n]
            }
            for (i in 1:length(b_names))
            {
                b_name <- b_names[i]
                batches[[i]] <- as.matrix(pca[scExp$batch_name == b_name,])
                adj_annot <- rbind(adj_annot,
                                   SummarizedExperiment::colData(scExp)[scExp$batch_name ==
                                                                            b_name,])
            }
          
            mnn.out <-
                do.call(batchelor::fastMNN, c(
                    batches,
                    list( 
                        k = 25,
                        d = 50,
                        ndist = 3,
                        pc.input = T,
                        auto.order = T,
                        cos.norm = F
                    )
                ))
            
            pca <- mnn.out$corrected
            SummarizedExperiment::colData(scExp) <- adj_annot
        } else
        {
            scExp$batch_id <- "batch_1"
            if (class(mat) %in% c("dgCMatrix", "dgTMatrix"))
            {
                out <- pca_irlba_for_sparseMatrix(Matrix::t(mat), n)
                pca = out[["PCA"]]
                rotation = out[["rotation"]]
            } else
            {
                pca <- stats::prcomp(Matrix::t(mat),
                                     center = T,
                                     scale. = F)
                rotation = pca$rotation
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
        if(!is.null(rotation)) {
          scExp@metadata$rotation = rotation 
          rownames(scExp@metadata$rotation) = rownames(scExp)
        }
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
    
    x.means <- Matrix::colMeans(x)
    svd.0 <- irlba::irlba(x, center = x.means, nv = n_comp)
    x. <- sweep(x, 2, x.means, "-")
    pca <- x. %*% svd.0$v
    
    return(list("PCA"=pca,"rotation" = svd.0$v))
}

#' Table of cells
#'
#' @param scExp A SingleCellExperiment object.
#'
#' @export
#' @return A formatted kable in HTML.
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr bind_rows tibble left_join
#' @importFrom kableExtra kable kable_styling group_rows
#'
num_cell_scExp <- function(annot)
{
  stopifnot(!is.null(annot))
  
  table <- as.data.frame(table(annot$sample_id))
  
  colnames(table) <- c("Sample", "#Cells")
  rownames(table) <- NULL

  table[, 1] <- as.character(table[, 1])
  table[nrow(table)+1,]= c("",sum(table[,2]))
  table %>% kableExtra::kable(escape = F, align = "c") %>%
    kableExtra::kable_styling(c("striped",
                                "condensed"), full_width = T) %>%
    kableExtra::group_rows("Total cell count",
                           dim(table)[1], dim(table)[1])
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

#' Subsample scExp
#'
#' Randomly sample x cells from each sample in a SingleCellExperiment to return
#' a subsampled SingleCellExperiment with all samples having maximum n cells. If
#' n is higher than the number of cell in a sample, this sample will not be
#' subsampled.
#'
#' @param scExp A SingleCellExperiment
#' @param n_cells An integer number of cells to subsample for each sample [500]
#'
#' @return
#' @export
#'
#' @importFrom SingleCellExperiment colData
#' @examples
subsample_scExp <- function(scExp, n_cells = 500) {
    stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(n_cells))
    
    annot = as.data.frame(SingleCellExperiment::colData(scExp))
    counts = SingleCellExperiment::counts(scExp)
    samples = as.character(unique(annot$sample_id))
    counts. = NULL
    annot. = NULL
    
    for (samp in samples) {
        cells = as.character(annot$cell_id[which(annot$sample_id == samp)])
        cells = sample(cells, min(n_cells, length(cells)), replace = F)
        if (is.null(counts.))
            counts. = counts[, cells]
        else
            counts. = Matrix::cbind2(counts., counts[, cells])
        if (is.null(annot.))
            annot. = annot[cells, ]
        else
            annot. = Matrix::rbind2(annot., annot[cells, ])
    }
    cat("ChromSCape::subsample_scExp -")
    ord = order(colnames(counts.))
    ord2 = order(rownames(annot.))
    counts. = counts.[, ord]
    annot. = annot.[ord2, ]
    scExp. = create_scExp(
        counts.,
        annot.,
        remove_zero_cells = F,
        remove_zero_features = F,
        remove_non_canonical = F,
        remove_chr_M = F
    )
    
    return(scExp.)
}

