
#' Generate cell cluster pseudo-bulk coverage tracks
#' 
#' @description  Generate cell cluster pseudo-bulk coverage tracks. First, scBED
#' files are concatenated into cell clusters contained in the 'cell_cluster' 
#' column of your SingleCellExperiment object. To do so, for each sample in the 
#' given list, the barcodes of each cluster are grepped and BED files are 
#' merged into pseudo-bulk of clusters (C1,C2...). Two cells from different can
#' have the same barcode ID as cell affectation is done sample by sample.
#'  Then coverage of pseudo-bulk BED files is calculated by averaging &
#'   smoothing reads on small genomic window (150bp per default). The pseudo 
#'   bulk BED and BigWigs coverage tracks are writtend to the output directory.
#'   This functionality is not available on Windows as it uses the 'cat' and 
#'   'gzip' utilities from Unix OS.
#' 
#' @param scExp_cf A SingleCellExperiment with cluster selected.
#' (see \code{\link{choose_cluster_scExp}}). It is recommended having a minimum
#' of ~100 cells per cluster in order to obtain smooth tracks.
#' @param input Either a named list of character vector of path towards 
#' single-cell BED files or a sparse raw matrix of small bins (<<500bp). If 
#' a named list specifying scBEDn the names MUST correspond to the 'sample_id' 
#' column in your SingleCellExperiment object. The single-cell BED files names MUST 
#' match the  barcode names in your SingleCellExperiment (column 'barcode'). The
#' scBED files can be gzipped or not. 
#' @param odir The output directory to write the cumulative BED and BigWig
#'  files.
#' @param format  File format, either "BAM" or "BED"
#' @param ref_genome The genome of reference, used to constrain to canonical 
#' chromosomes. Either 'hg38' or 'mm10'. 'hg38' per default. 
#' @param bin_width The width of the bin to create the coverage track. The 
#' smaller the greater the resolution & runtime. Default to 150.
#' @param n_smoothBin Number of bins left & right to average ('smooth') the 
#' signal on. Default to 5.
#' @param read_size The estimated size of reads. Default to 101.
#' @param progress A Progress object for Shiny. Default to NULL.
#'
#' @return Generate coverage tracks (.bigwig) for each cluster in the 
#' SingleCellExperiment ("cell_cluster" column).
#' @export
#'
#' @examples \dontrun{
#' data(scExp)
#' input_files_coverage = list(
#'   "scChIP_Jurkat_K4me3" = paste0("/path/to/",scExp$barcode[1:51],".bed"),
#'   "scChIP_Ramos_K4me3" = paste0("/path/to/",scExp$barcode[52:106],".bed")
#' )
#' generate_coverage_tracks(scExp, input_files_coverage, "/path/to/output",
#' ref_genome = "hg38")
#' }
generate_coverage_tracks <- function(scExp_cf, input, odir, 
                                     format = "scBED",
                                     ref_genome = c("hg38","mm10")[1],
                                     bin_width = 150,  n_smoothBin = 5,
                                     read_size = 101,
                                     progress = NULL){
    stopifnot(is(scExp_cf,"SingleCellExperiment"), is.character(format),
              dir.exists(odir), ref_genome %in% c("hg38","mm10"), 
              is.numeric(bin_width), is.numeric(n_smoothBin),
              is.numeric(read_size))
    
    if(format == "scBED") stopifnot(length(intersect(scExp_cf$sample_id,names(input))) > 0)
    
    nclust = length(unique(scExp_cf$cell_cluster))
    if (!is.null(progress)) progress$set(message=paste0('Generating coverage tracks for k= ', nclust,' clusters ...'), value = 0.1)
    
    affectation = SingleCellExperiment::colData(scExp_cf)
    affectation$cell_cluster = as.factor(affectation$cell_cluster)
    
    if (!is.null(progress)) progress$set(detail = "Generating pseudo-bulk...", value = 0.2)
    if(format == "scBED") {
        concatenate_scBed_into_clusters(affectation, input, odir)
    } 
    
    if (!is.null(progress)) progress$set(detail = "Creating coverage files for each cluster...", value = 0.40)
    
    if(format == "raw_mat"){
        original_bins = rownames(input)
        original_bins =  strsplit(original_bins, "_", fixed = T)
        original_bins_chr = as.character(lapply(original_bins, function(x) x[1]))
        original_bins_start = as.numeric(lapply(original_bins, function(x) x[2]))
        original_bins_end = as.numeric(lapply(original_bins, function(x) x[3]))
        gc()
        original_bins = GenomicRanges::GRanges(
            seqnames = original_bins_chr, ranges = IRanges::IRanges(original_bins_start, original_bins_end))
        input = input[,match(colnames(scExp_cf), colnames(input))]
    }
    suffix = ".bed"
    n = 0
    for(class in unique(scExp_cf$cell_cluster)) {
        n = n + 1
        if (!is.null(progress)) progress$set(detail = paste0("Coverage ",class,"..."),
                                             value = 0.4 + n * (0.6/length(unique(scExp_cf$cell_cluster))) )
        out_bw = file.path(odir, paste0(class,".bw"))
        if(format == "scBED") {
            input = file.path(odir, paste0(class,suffix))
            rawfile_ToBigWig(input, out_bw, "BED", bin_width = bin_width,
                             n_smoothBin = n_smoothBin,  ref = ref_genome,
                             read_size = read_size)
        } else{
            input. = input[,which(scExp_cf$cell_cluster %in% class)]
            rawfile_ToBigWig(input = input., BigWig_filename = out_bw,
                             format = "raw_mat", bin_width = bin_width,
                             n_smoothBin = n_smoothBin,  ref = ref_genome,
                             read_size = read_size,
                             original_bins = original_bins)
        }
        }
    if (!is.null(progress)) progress$set(detail = "Done !", value = 0.95)
}

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
    gzipped = grepl(".gz", files_list[[1]][1])
    suffix = ""
    if(gzipped) suffix = ".gz"
    for(sample in names(files_list)){
        message("ChromSCape:::concatenate_scBed_into_clusters - concatenating ",
                "files for ", sample, " files...")
        file_pool = files_list[[sample]]
        for (class in levels(factor(affectation$cell_cluster)))
        {
            message("ChromSCape:::concatenate_scBed_into_clusters - concatenating ", class,"...")
            class_barcodes <- affectation$barcode[which(affectation$cell_cluster == as.character(class) &
                                                            affectation$sample_id == sample)]
   
            files_class = c()
            if(length(class_barcodes)>500){
                
                while(length(class_barcodes)>500){
                    class_barcodes_chunk = class_barcodes[1:500]
                    files_class = c(files_class,
                                    file_pool[grep(paste(class_barcodes_chunk, collapse="|"),
                                                   file_pool,useBytes = TRUE, perl = TRUE)])
                    class_barcodes = class_barcodes[-c(1:500)]
                }
                files_class = c(files_class,
                                file_pool[grep(paste(class_barcodes, collapse="|"),
                                               file_pool,useBytes = TRUE, perl = TRUE)])
                
            } else {
                files_class = file_pool[grep(paste(class_barcodes, collapse="|"),
                                             file_pool,useBytes = TRUE, perl = TRUE)]
            }
            if(length(class_barcodes)>0){
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
#' @param input Either a named list of character vector of path towards 
#' single-cell BED files or a sparse raw matrix of small bins (<<500bp). If 
#' a named list specifying scBEDn the names MUST correspond to the 'sample_id' 
#' column in your SingleCellExperiment object. The single-cell BED files names MUST 
#' match the  barcode names in your SingleCellExperiment (column 'barcode'). The
#' scBED files can be gzipped or not.
#' @param BigWig_filename Path to write the output BigWig file
#' @param format File format, either "BAM" or "BED"
#' @param bin_width Bin size for coverage
#' @param n_smoothBin Number of bins for smoothing values
#' @param ref Reference genome.
#' @param read_size Length of the reads.
#' @param original_bins Original bins GenomicRanges in case the format is raw
#' matrix.
#'
#' @importFrom rtracklayer export.bw
#' @importFrom GenomicRanges tileGenome width seqnames
#' @return Writes a BigWig file as output
#' 
rawfile_ToBigWig <- function(input, BigWig_filename, format = "BAM",
                             bin_width = 150, n_smoothBin = 5, ref = "hg38",
                             read_size = 101, original_bins = NULL){
    bins = NULL
    canonical_chr <- eval(parse(text = paste0("data(",ref, ".chromosomes)")))
    canonical_chr$start = 1
    canonical_chr <- GenomicRanges::GRanges(seqnames = canonical_chr$chr,
                                            ranges = IRanges::IRanges(
                                              start = canonical_chr$start,
                                              end = canonical_chr$end
                                              ))
    GenomeInfoDb::seqlengths(canonical_chr) = end(canonical_chr)
    if(format != "raw_mat") {
        message("ChromSCape:::rawfile_ToBigWig - generating bigwig for  ",
             basename(filename), " file...")

        
        bins <- unlist(GenomicRanges::tileGenome(
            setNames(
                GenomicRanges::width(canonical_chr),
                GenomicRanges::seqnames(canonical_chr)
            ),
            tilewidth = bin_width
        ))
        gc()
    } else if(format == "raw_mat")  {
        stopifnot(!is.null(original_bins), is(original_bins, "GRanges"))
           message("ChromSCape:::rawfile_ToBigWig - generating bigwig from the ",
                                     "raw matrix...")
    } else{
        message("ChromSCape:::rawfile_ToBigWig - format must be either 'raw_mat'",
                ", 'BAM' or 'BED'. Returning.")
        return()
    }
    
    bins <- count_coverage(input, format, bins, canonical_chr, n_smoothBin,
                           ref, read_size, original_bins)
    ## export as bigWig
    GenomeInfoDb::seqlengths(bins) = GenomeInfoDb::seqlengths(canonical_chr)[
      match(names(GenomeInfoDb::seqlengths(bins)),
            names(GenomeInfoDb::seqlengths(canonical_chr)))]
    rtracklayer::export.bw(bins, BigWig_filename)
}

#' Create a smoothed and normalized coverage track from a BAM file and
#' given a bin GenomicRanges object (same as deepTools bamCoverage)
#'
#' Normalization is CPM, smoothing is done by averaging on n_smoothBin regions
#' left and right of any given region.
#' 
#' @param input Either a named list of character vector of path towards 
#' single-cell BED files or a sparse raw matrix of small bins (<<500bp). If 
#' a named list specifying scBEDn the names MUST correspond to the 'sample_id' 
#' column in your SingleCellExperiment object. The single-cell BED files names MUST 
#' match the  barcode names in your SingleCellExperiment (column 'barcode'). The
#' scBED files can be gzipped or not.
#' @param format File format, either "BAM" or "BED"
#' @param bins A GenomicRanges object of binned genome
#' @param canonical_chr GenomicRanges of the chromosomes to read the BAM file.
#' @param n_smoothBin Number of bins left and right to smooth the signal.
#' @param ref Genomic reference
#' @param read_size Length of the reads
#' @param original_bins Original bins GenomicRanges in case the format is raw
#' 
#' matrix.
#' @importFrom Rsamtools ScanBamParam scanBam scanBamWhat
#' @importFrom GenomicRanges tileGenome width seqnames findOverlaps
#' GRanges
#' @importFrom IRanges IRanges
#' 
#' 
#'  
#' @return A binned GenomicRanges that can be readily exported into bigwig file.
#'
count_coverage <-function(input, format = "BAM", bins, canonical_chr,
                          n_smoothBin = 5, ref = "hg38", read_size = 101,
                          original_bins = NULL)
{
    if(format == "BAM"){
        param = Rsamtools::ScanBamParam(which = canonical_chr,
                                        what = Rsamtools::scanBamWhat()[c(3,5)])
        gr <- Rsamtools::scanBam(file = input, param = param)
        chr = unlist(lapply(1:length(canonical_chr), function(i){gr[[i]]$rname}))
        start = unlist(lapply(1:length(canonical_chr), function(i){gr[[i]]$pos}))
        end = start + 101
        gr = GenomicRanges::GRanges(seqnames = chr,
                                    ranges = IRanges::IRanges("start" = start,
                                                              "end" = end))
    } else if(format == "BED") {
        gr <- rtracklayer::import(input)
    } else if(format == "raw_mat"){
        bins <- original_bins
    } else {
        message("ChromSCape::count_coverage - format must be either 'BAM', 'BED',",
             " or 'raw_mat'. Returning.")
        return()
    }
    
    if(format != "raw_mat"){
        hits <- GenomicRanges::findOverlaps(
            bins, gr, minoverlap = 1, ignore.strand = TRUE)
        hits = as.matrix(hits)
        hits_agg = table(hits[,1])
        bins$score = 0
        bins$score[as.numeric(names(hits_agg))] = hits_agg
        rm(hits_agg)
        gc()
    } else{
        bins$score = rowSums(input)
    }
   
    bins$score = 10^6* bins$score / length(bins)
    bins$score = smoothBin(bins$score, n_smoothBin)
    bins = bins[-which(bins$score==0)]
    gc()
    return(bins)
}
