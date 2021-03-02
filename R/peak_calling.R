## Authors : Pacôme Prompsy, Celine Vallot

#' Peak calling on cell clusters
#' 
#' This functions does peak calling on each cell population in order to refine
#' gene annotation for large bins. For instance, a 50000bp bin might containt
#' the TSS of several genes, while in reality only one or two of these genes are
#' overlapping the signal (peak). To do so, first in-silico cell sorting is
#' applied based on previously defined clusters contained in the
#' SingleCellExperiment. Taking BAM files of each sample as input, samtools
#' pools then splits reads from each cell barcode into 1 BAM file per cell
#' cluster (pseudo-bulk). Then MACS2 calls peaks on each cluster. The peaks are
#' aggregated and merged if closer to a certain distance
#' defined by user (10000bp). Then, 
#' 
#' This function takes as input a SingleCellExperiment, that must contain a
#' 'cell_cluster' column in it's colData, an output directory where to store
#' temporary files, the list of BAM files corresponding to each sample and
#' containing the cell barcode information as a tag (for instance tag CB:Z:xxx,
#' XB:Z:xxx or else...), the p.value used by MACS2 to distinguish significant
#' peaks, the reference genome (either hg38 or mm10), the maximal merging
#' distance in bp and a data.frame containing gene TSS genomic cooridnates of
#' corresponding genome (if set to NULL, will automatically load geneTSS). The
#' output is a SingleCellExperiment with GRanges object containing ranges of
#' each merged peaks that falls within genomic bins of the SingleCellExperiment,
#' saving the bin range as additional column (window_chr, window_start,
#' window_end), as well as the closests genes and their distance relative to the
#' peak. The peaks may be present in several rows if multiple genes are close /
#' overlap to the peaks.
#'
#' Note that the user must have MACS2 installed and available in the
#' PATH. Users can open command terminal and type 'which
#' macs2' to verify the availability of these programs. Will only work on unix
#' operating system. Check operating system with 'print(.Platform)'.
#' 
#'
#' @param scExp A SingleCellExperiment object
#' @param odir Output directory where to write temporary files and each
#'   cluster's BAM file
#' @param input A character vector of file paths to each sample's BAM file,
#'   containing cell barcode information as tags. BAM files can be paired-end or
#'   single-end.
#' @param p.value a p-value to use for MACS2 to determine significant peaks.
#'   (0.05)
#' @param ref A reference genome, either hg38 or mm10. ('hg38')
#' @param peak_distance_to_merge Maximal distance to merge peaks together after
#'   peak calling , in bp. (10000)
#' @param geneTSS_annotation A data.frame annotation of genes TSS. If NULL will
#'   automatically load Gencode list of genes fro specified reference genome.
#' @param run_coverage Create coverage tracks (.bw) for each cluster ? 
#'
#' @return A SingleCellExperiment with refinded annotation
#' @export
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom GenomicRanges GRanges pintersect ranges
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom IRanges findOverlapPairs subsetByOverlaps
#' @importFrom utils read.table
#' @importFrom Rsamtools mergeBam indexBam
#' 
#' @examples 
#' \dontrun{
#' data("scExp")
#' subset_bam_call_peaks(scExp, "path/to/out/", list("sample1" = 
#'    "path/to/BAM/sample1.bam", "sample2" = "path/to/BAM/sample2.bam"),
#'    p.value = 0.05, ref = "hg38", peak_distance_to_merge = 10000, 
#'    geneTSS_annotation = NULL)
#'}
subset_bam_call_peaks <- function(scExp, odir, input, format = "BAM", p.value = 0.05,
                                ref = "hg38", peak_distance_to_merge = 10000,
                                geneTSS_annotation = NULL, run_coverage = TRUE)
    {
    stopifnot(is(scExp, "SingleCellExperiment"), dir.exists(odir),
            is.numeric(p.value), is.character(ref),
            is.numeric(peak_distance_to_merge))
    if (!ref %in% c("hg38", "mm10")) 
        stop("ChromSCape::subset_bam_call_peaks - 
            Reference genome (ref) must be 'hg38' or 'mm10'.")
    if (!"barcode" %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::subset_bam_call_peaks - 
            Colnames of cell annotation must contain 'barcode'.")
    if (!"cell_cluster" %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::subset_bam_call_peaks -
            Colnames of cell annotation must contain 'cell_cluster'.")
    
    print("Subset Bam Call peaks")
    print(format)
    if (is.null(geneTSS_annotation))
    {
        message("ChromSCape::gene_set_enrichment_analysis_scExp - ",
                    "Selecting ", 
            ref, " genes from Gencode.")
        eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
        geneTSS_annotation = 
            eval(parse(text = paste0("", ref, ".GeneTSS")))
        start = geneTSS_annotation$start
        geneTSS_annotation$start = ifelse(geneTSS_annotation$strand == "+",
                                          geneTSS_annotation$start,
                                          geneTSS_annotation$end)
        geneTSS_annotation$end = ifelse(geneTSS_annotation$strand == "+",
                                        start +1,
                                        geneTSS_annotation$end +1)
        geneTSS_annotation$strand = NULL
        geneTSS_annotation = as(geneTSS_annotation,"GRanges")
    } else geneTSS_annotation = as(geneTSS_annotation, "GRanges")
    if(format == "BAM"){
        if (length(input) > 1)
        {
            merged = file.path(odir, "merged.bam")
            Rsamtools::mergeBam(input, destination = merged_bam,
                                overwrite = TRUE, indexDestination = TRUE)
            
        } else merged = input[1]
        if(!file.exists(paste0(merged,".bai"))) 
            Rsamtools::indexBam(merged, overwrite= TRUE)
    } 
    print("Writting barcode files...")
    affectation = SingleCellExperiment::colData(scExp)
    affectation$cell_cluster = as.factor(affectation$cell_cluster)
    
    p.value = paste0(" -p ", p.value, " ")  # format for system call to macs2
    
    if(format == "BAM") separate_BAM_into_clusters(affectation, odir, merged) 
    else concatenate_scBed_into_clusters(affectation, input, odir, ref)
    
    merged_peaks <- call_macs2_merge_peaks(affectation, odir, p.value, 
                                           format, ref, peak_distance_to_merge)
    scExp@metadata$refined_annotation <- 
        annotation_from_merged_peaks(scExp, odir,
                                     merged_peaks, geneTSS_annotation)
    print("Finished creating BAM and annot")
    cat("Running coverage ? ", run_coverage)
    if(run_coverage){
        suffix = ifelse(format=="BAM", ".bam", ".bed")
        for(class in unique(scExp$cell_cluster)) {
            input_file = file.path(odir, paste0(class,suffix))
            out_bw = file.path(odir, paste0(class,".bw"))
            rawfile_ToBigWig(input_file, out_bw, format, bin_width = 150,
                             n_smoothBin = 5,  ref = ref, read_size = 101)
        }
    }
    return(scExp)
}

#' Concatenate single-cell BED into clusters
#'
#' @param affectation Annotation data.frame containing cluster information
#' @param files Input files to concatenate (must match affectation$barcode)
#' @param odir Output directory
#'
#' @return Merge single-cell BED files into cluster BED files. Ungzip file if 
#' BED is gzipped.
#'
concatenate_scBed_into_clusters <- function(affectation, files, odir,
                                            ref){
    unlink(file.path(odir, "*.bed"))
    # canonical_chr <- eval(parse(text = paste0("ChromSCape::",
    #                                           ref, ".chromosomes")))
    # canonical_chr$start = 1
    # canonical_chr <- as(canonical_chr, "GRanges")
    gzipped = grepl(".gz",files[1])
    suffix = ""
    if(gzipped) suffix = ".gz"
    print(files)
    for (class in levels(factor(affectation$cell_cluster)))
    {
        message("ChromSCape:::concatenate_scBed_into_clusters - generating ",
                class, " BED file...")
        class_barcodes <- affectation$barcode[which(affectation$cell_cluster == 
                                                        as.character(class))]
        files_class = files[grep(paste(class_barcodes, collapse="|"),files)]
        for(file in files_class){
            command = paste0("cat '", file, "' >> '",
                             file.path(odir, paste0(class, ".bed",suffix,"'")))
            system(command)
        }
        if(gzipped) {
            command = paste0("gzip -cd '", file.path(
                odir,paste0(class,".bed", suffix)),
                "' > '", file.path(odir,paste0(class, ".bed'")))
            system(command)
        }
        # bed = rtracklayer::import(file.path(odir, paste0(class, ".bed")))
        # bed= bed[GenomicRanges::findOverlaps(
        #     canonical_chr,bed,minoverlap = 1)@to]
        # rtracklayer::export.bed(bed, file.path(odir, paste0(class, ".bed")))
    }
}

#' Calling MACS2 peak caller and merging resulting peaks
#'
#' @param affectation Annotation data.frame with cell cluster and cell id
#'   information
#' @param odir Output directory to write MACS2 output
#' @param p.value P value to detect peaks, passed to MACS2
#' @param format File format, either "BAM" or "BED"
#' @param ref Reference genome
#' @param peak_distance_to_merge Distance to merge peaks
#'
#' @return A list of merged GRanges peaks
#'
#' 
call_macs2_merge_peaks <- function(affectation, odir, p.value, format = "BAM",
                                ref, peak_distance_to_merge){
    suffix = ".bam"
    if(format != "BAM") suffix = ".bed"
    merged_peaks=list()
    for (class in levels(factor(affectation$cell_cluster))){
        
        macs2_options = paste0(" -f ",format," --nomodel --extsize 300 ")
        
        command = paste0("macs2 callpeak ", p.value,
                    macs2_options, " --keep-dup all --broad -t '", 
                    file.path(odir, paste0(class, suffix)),
                    "' --outdir '", odir, "' --name ", class)
        print(command)
        system(command)
        merged_peaks[[class]] <-
            merge_MACS2_peaks(odir, class, peak_distance_to_merge, ref)
    }
    unlink(file.path(odir, "*.xls"))
    unlink(file.path(odir, "*.gappedPeak"))
    unlink(file.path(odir, "*_model.r"))
    unlink(file.path(odir, "header.sam"))
    unlink(file.path(odir, 'merged.bam'))
    return(merged_peaks)
}

#' Merge peak files from MACS2 peak caller
#'
#' @param odir Output directory
#' @param class Cell cluster
#' @param peak_distance_to_merge Maximum distance to merge two peaks 
#' @param ref Reference genome
#'
#' @return Peaks as GRanges
#'
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges GRanges reduce ranges
merge_MACS2_peaks <- function(odir,class,peak_distance_to_merge,ref){
    merged_peaks = list()
    peaks = read.table(
        file = file.path(odir, paste0(class, "_peaks.broadPeak")),
        colClasses = c("character", "integer", "integer", rep("NULL", 6)))
    colnames(peaks) = c("chr", "start", "end")
    peaks = GenomicRanges::GRanges(peaks)
    peaks = peaks[which(width(GenomicRanges::ranges(peaks)) >= 500), ]
    peaks = GenomicRanges::reduce(peaks, min.gapwidth = peak_distance_to_merge, 
                                ignore.strand = TRUE)
    
    ref_chromosomes = GenomicRanges::GRanges(
        eval(parse(text = paste0(ref, ".chromosomes"))))
    peaks = IRanges::subsetByOverlaps(peaks, ref_chromosomes, 
                                                    ignore.strand = TRUE)
    return(peaks)
}


#' Find nearest peaks of each gene and return refined annotation
#'
#' @param scExp A SingleCellExperiment object 
#' @param odir An output directory where to write the mergedpeaks BED file
#' @param merged_peaks A list of GRanges object containing the merged peaks
#' @param geneTSS_annotation A GRanges object with reference genes
#'
#' @return A data.frame with refined annotation
#' 
#' @importFrom GenomicRanges start end union distanceToNearest 
#' pintersect
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom IRanges findOverlapPairs
#' @importFrom rtracklayer export.bed
annotation_from_merged_peaks <- function(scExp, odir,
                                        merged_peaks,
                                        geneTSS_annotation){
    segmentation = SummarizedExperiment::rowRanges(scExp)
    S4Vectors::mcols(segmentation) = NULL
    segmentation$window_ID = paste(
        as.character(segmentation@seqnames), GenomicRanges::start(segmentation),
        GenomicRanges::end(segmentation), sep = "_")
    
    merged_peak = merged_peaks[[1]]
    for (i in 2:length(merged_peaks))
    {
        merged_peak = suppressWarnings(
            GenomicRanges::union(merged_peak, merged_peaks[[i]], 
                                ignore.strand = TRUE))
    }
    rtracklayer::export.bed(merged_peak, file.path(odir, "merged_peaks.bed"))
    pairs <- IRanges::findOverlapPairs(
        segmentation, merged_peak, ignore.strand = TRUE)
    refined_annotation = GenomicRanges::pintersect(pairs, ignore.strand = TRUE)
    S4Vectors::mcols(refined_annotation)$hit = NULL
    
    hits_genes = GenomicRanges::distanceToNearest(
        refined_annotation, geneTSS_annotation, 
        ignore.strand = TRUE, select = "all")
    
    refined_annotation = refined_annotation[S4Vectors::queryHits(hits_genes)]
    refined_annotation$Gene = as.character(
        geneTSS_annotation$gene[S4Vectors::subjectHits(hits_genes)])
    refined_annotation$distance = hits_genes@elementMetadata$distance
    
    refined_annotation$peak_ID = paste(
        as.character(refined_annotation@seqnames), 
        GenomicRanges::start(refined_annotation),
        GenomicRanges::end(refined_annotation), sep = "_")
    return(refined_annotation)
}

#' Separate BAM files into cell cluster BAM files
#'
#' @param affectation An annotation data.frame containing cell_id and 
#' cell_cluster columns
#' @param odir A valid output directory path
#' @param merged_bam A list of merged bam file paths
#'  
#'  @importFrom Rsamtools filterBam ScanBamParam
#' @return Create one BAM per cluster from one BAM per condition
separate_BAM_into_clusters <- function(affectation, odir, merged_bam){
    
    for (class in levels(factor(affectation$cell_cluster)))
    {
        message("ChromSCape:::separate_BAM_into_clusters - generating ",class,
                " BAM file...")
        class_barcodes <- affectation$barcode[which(affectation$cell_cluster == 
                                            as.character(class))]
        filt <- Rsamtools::ScanBamParam(tag=c("XB"),
                                        tagFilter=list(XB=class_barcodes))
        
        Rsamtools::filterBam(merged_bam,file.path(odir,paste0(class,".bam")),
                    indexDestination = TRUE, overwrite=TRUE, param=filt)
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
