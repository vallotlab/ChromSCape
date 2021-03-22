#' Concatenate single-cell BED into clusters
#'
#' @param affectation Annotation data.frame containing cluster information
#' @param files_list Named list of scBED file paths to concatenate. List Names 
#' must match affectation$sample_id and basenames must match
#'  affectation$barcode.
#' @param odir Output directory to write concatenate pseudo-bulk BEDs.
#'
#' @return Merge single-cell BED files into cluster BED files. Ungzip file if 
#' BED is gzipped.
#'
concatenate_scBed_into_clusters <- function(affectation, files_list, odir){
    unlink(file.path(odir, "C*.bed"))
    gzipped = grepl(".gz",files_list[[1]][1])
    suffix = ""
    if(gzipped) suffix = ".gz"
    print(head(files_list[[1]]))
    for(sample in names(files_list)){
        cat("Doing ", sample, "\n")
        message("ChromSCape:::concatenate_scBed_into_clusters - concatenating ",
                "files for ", sample, " files...")
        file_pool = files_list[[sample]]
        cat("Length file pool ", length(file_pool), "\n")
        for (class in levels(factor(affectation$cell_cluster)))
        {
            message("ChromSCape:::concatenate_scBed_into_clusters - concatenating ", class,"...")
            class_barcodes <- affectation$barcode[which(affectation$cell_cluster == as.character(class) &
                                                            affectation$sample_id == sample)]
            cat("Length class_barcodes", length(class_barcodes), "\n")
            files_class = file_pool[grep(paste(class_barcodes, collapse="|"),file_pool, perl = T)]
            print(head(files_class))
            cat("Length files_class", length(files_class), "\n")
            if(length(class_barcodes>0)){
                for(file in files_class){
                    command = paste0("cat '", file, "' >> '",
                                     file.path(odir, paste0(class, ".bed",suffix,"'")))
                    system(command)
                }
            }
        }
    }
    for(class in levels(factor(affectation$cell_cluster))){
        if(gzipped) {
            command = paste0("gzip -cd '", file.path(
                odir,paste0(class,".bed", suffix)),
                "' > '", file.path(odir,paste0(class, ".bed'")))
            system(command)
        }
    }
}


#' Smooth a vector of values with nb_bins left and righ values
#'
#' @param bin_score A numeric vector of values to be smoothed
#' @param nb_bins Number of values to take left and right
#'  
#'  @importFrom BiocParallel bpvec
#' @return A smooth vector of the same size 
#' 
smoothBin <- function(bin_score, nb_bins = 10){
    bin_original = bin_score
    fun <-function(v){
        v_or = v
        start = nb_bins + 1
        end = length(v) - nb_bins -1
        for(i in start:end){
            v[i] = mean(v_or[(i-nb_bins):(i+nb_bins)])
        }
        return(v)
    }
    bin_score = unlist(BiocParallel::bpvec(bin_score, fun))
    return(bin_score)
}

#' rawfile_ToBigWig : reads in BAM file and write out BigWig coverage file, 
#' normalized and smoothed
#'
#' @param filename Path to the BAM file (with index) or BED file
#' @param BigWig_filename Path to write the output BigWig file
#' @param format File format, either "BAM" or "BED"
#' @param bin_width Bin size for coverage
#' @param n_smoothBin Number of bins for smoothing values
#' @param ref Reference genome.
#' @param read_size Length of the reads.
#'
#' @importFrom rtracklayer export.bw
#' @importFrom GenomicRanges tileGenome width seqnames
#' @return Writes a BigWig file as output
#' 
rawfile_ToBigWig <- function(filename, BigWig_filename, format = "BAM",
                             bin_width = 150, n_smoothBin = 5, ref = "hg38",
                             read_size = 101){
    canonical_chr <- eval(parse(text = paste0("ChromSCape::",
                                              ref, ".chromosomes")))
    canonical_chr$start = 1
    canonical_chr <- as(canonical_chr, "GRanges")
    
    bins <- unlist(GenomicRanges::tileGenome(
        setNames(
            GenomicRanges::width(canonical_chr),
            GenomicRanges::seqnames(canonical_chr)
        ),
        tilewidth = bin_width
    ))
    gc()
    
    bins <- count_coverage(filename, format, bins, canonical_chr, n_smoothBin,
                           ref, read_size)
    ## export as bigWig
    rtracklayer::export.bw(bins, BigWig_filename)
}

#' Create a smoothed and normalized coverage track from a BAM file and
#' given a bin GenomicRanges object (same as deepTools bamCoverage)
#'
#' Normalization is CPM, smoothing is done by averaging on n_smoothBin regions
#' left and right of any given region.
#' 
#' @param filename Path towards the BAM to create coverage from
#' @param format File format, either "BAM" or "BED"
#' @param bins A GenomicRanges object of binned genome
#' @param canonical_chr GenomicRanges of the chromosomes to read the BAM file.
#' @param n_smoothBin Number of bins left and right to smooth the signal.
#' @param ref Genomic reference
#' @param read_size Length of the reads
#'   
#' @importFrom Rsamtools ScanBamParam scanBam scanBamWhat
#' @importFrom GenomicRanges tileGenome width seqnames findOverlaps
#' GRanges
#' @importFrom IRanges IRanges
#' 
#' 
#'  
#' @return A binned GenomicRanges that can be readily exported into bigwig file.
#'
count_coverage <-function(filename, format = "BAM", bins, canonical_chr,
                          n_smoothBin = 5, ref = "hg38", read_size = 101)
{
    if(format == "BAM"){
        param = Rsamtools::ScanBamParam(which = canonical_chr,
                                        what = Rsamtools::scanBamWhat()[c(3,5)])
        gr <- Rsamtools::scanBam(file = filename, param = param)
        chr = unlist(lapply(1:length(canonical_chr), function(i){gr[[i]]$rname}))
        start = unlist(lapply(1:length(canonical_chr), function(i){gr[[i]]$pos}))
        end = start + 101
        gr = GenomicRanges::GRanges(seqnames = chr,
                                    ranges = IRanges::IRanges("start" = start,
                                                              "end" = end))
    } else {
        gr <- rtracklayer::import(filename)
    }
    
    hits <- GenomicRanges::findOverlaps(
        bins, gr, minoverlap = 1, ignore.strand = TRUE)
    hits = as.matrix(hits)
    hits_agg = table(hits[,1])
    bins$score = 0
    bins$score[as.numeric(names(hits_agg))] = hits_agg
    rm(hits_agg)
    gc()
    bins$score = 10^6* bins$score / length(gr)
    bins$score = smoothBin(bins$score, n_smoothBin)
    bins = bins[-which(bins$score==0)]
    gc()
    return(bins)
}
