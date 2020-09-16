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
#' Note that the user must have samtools & MACS2 installed and available in the
#' PATH. Users can open command terminal and type 'which samtools' & 'which
#' macs2' to verify the availability of these programs. Will only work on unix
#' operating system. Check operating system with 'print(.Platform)'.
#' 
#'
#' @param scExp A SingleCellExperiment object
#' @param odir Output directory where to write temporary files and each
#'   cluster's BAM file
#' @param inputBam A character vector of file paths to each sample's BAM file,
#'   containing cell barcode information as tags. BAM files can be paired-end or
#'   single-end.
#' @param p.value a p-value to use for MACS2 to determine significant peaks.
#'   (0.05)
#' @param ref A reference genome, either hg38 or mm10. ('hg38')
#' @param peak_distance_to_merge Maximal distance to merge peaks together after
#'   peak calling , in bp. (10000)
#' @param geneTSS_annotation A data.frame annotation of genes TSS. If NULL will
#'   automatically load Gencode list of genes fro specified reference genome.
#'
#' @return A SingleCellExperiment with refinded annotation
#' @export
#'
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom GenomicRanges GRanges pintersect ranges
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom IRanges findOverlapPairs subsetByOverlaps
#' @importFrom utils read.table
subset_bam_call_peaks <- function(scExp, odir, inputBam, p.value = 0.05, ref = "hg38", 
    peak_distance_to_merge = 10000, geneTSS_annotation = NULL)
    {
    stopifnot(is(scExp, "SingleCellExperiment"), dir.exists(odir), is.numeric(p.value), 
        is.character(ref), is.numeric(peak_distance_to_merge))
    if (!ref %in% c("hg38", "mm10")) 
        stop("ChromSCape::subset_bam_call_peaks - Reference genome (ref) must be 'hg38' or 'mm10'.")
    if (!"barcode" %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::subset_bam_call_peaks - Colnames of cell annotation must contain 'barcode'.")
    if (!"cell_cluster" %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::subset_bam_call_peaks - Colnames of cell annotation must contain 'cell_cluster'.")
    
    if (is.null(geneTSS_annotation))
    {
        message(paste0("ChromSCape::gene_set_enrichment_analysis_scExp - Selecting ", 
            ref, " genes from Gencode."))
        eval(parse(text = paste0("data(", ref, ".GeneTSS)")))
        geneTSS_annotation = as(eval(parse(text = paste0("", ref, ".GeneTSS"))), 
            "GRanges")
    } else geneTSS_annotation = as(geneTSS_annotation, "GRanges")
    if (length(inputBam) > 1)
    {
        write(inputBam, file = file.path(odir, "bam_list.txt"))
        system(paste0("samtools merge -@ 4 -f -h ",
                      inputBam[1], " -b ", file.path(odir, "bam_list.txt"),
                      " ", file.path(odir, "merged.bam")))
        merged_bam = file.path(odir, "merged.bam")
    } else merged_bam = inputBam[1]

    print("Writting barcode files...")
    affectation = SingleCellExperiment::colData(scExp)
    affectation$cell_cluster = as.factor(affectation$cell_cluster)
    
    p.value = paste0(" -p ", p.value, " ")  # format for system call to macs2
    
    separate_BAM_into_clusters(affectation, odir, merged_bam)
    merged_peaks <- call_macs2_merge_peaks(affectation,odir,p.value,ref,
                                           peak_distance_to_merge)
    scExp@metadata$refined_annotation <- 
        annotation_from_merged_peaks(scExp, merged_peaks, geneTSS_annotation)
    return(scExp)
}


#' Calling MACS2 peak caller and merging resulting peaks
#'
#' @param affectation Annotation data.frame with cell cluster and cell id
#'   information
#' @param odir Output directory to write MACS2 output
#' @param p.value P value to detect peaks, passed to MACS2
#' @param ref Reference genome
#' @param peak_distance_to_merge Distance to merge peaks
#'
#' @return A list of merged GRanges peaks
#'
#' 
call_macs2_merge_peaks <- function(affectation,odir,p.value,
                                   ref,peak_distance_to_merge){
    for (class in levels(factor(affectation$cell_cluster)))
    {
        command <- paste0(
            "samtools flagstat ",
            file.path(odir, paste0(class, ".bam")),
            " | grep \"properly paired\" | sed \"s/.*properly paired",
            " (//g\" | cut -f1 -d\" \" | sed \"s/\\%//g\"")
        print(command)
        percent_properlyPaired = as.double(system(command,intern = TRUE))
        if (!is.na(percent_properlyPaired) & percent_properlyPaired > 50)
        {
            print("Bam files provided have more than 50% of their reads properly
                  mapped, running macs2 with default parameters...")
            macs2_options = ""
        } else
        {
            print("Bam files provided have less than 50% of their reads properly
                  mapped, flagging all reads as single end running macs2 with
                  --nomodel --extsize 300...")
            macs2_options = " --nomodel --extsize 300 "
            print("Transforming bam so that all mapped reads are considered
                  as single-end...")
            system(paste0("samtools view -H ",
                          file.path(odir, paste0(class, ".bam")), 
                          " > ", file.path(odir, "header.sam")))
            system(paste0("samtools view -F4 ",
                          file.path(odir, paste0(class, ".bam")), 
                          " | awk -v OFS=\"\t\" \"{\\$2=0; print \\$0}\" >> ",
                          file.path(odir, "header.sam")))
            system(paste0("samtools view -b ",
                          file.path(odir, "header.sam"), " > ", 
                          file.path(odir, paste0(class, ".bam"))))
        }
        command = paste0("macs2 callpeak ", p.value,
                     macs2_options, " --keep-dup all --broad -t ", 
                     file.path(odir, paste0(class, ".bam")),
                     " --outdir ", odir, " --name ", class)
        system(command)
        merged_peaks <- merge_MACS2_peaks(odir,class,peak_distance_to_merge,ref)
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
#' @return A list of merged peaks
#'
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges GRanges reduce ranges
merge_MACS2_peaks <- function(odir,class,peak_distance_to_merge,
                              ref){
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
    peaks = suppressWarnings(IRanges::subsetByOverlaps(peaks, ref_chromosomes, 
                                                       ignore.strand = TRUE))
    merged_peaks[[class]] = peaks
    return(merged_peaks)
}


#' Find nearest peaks of each gene and return refined annotation
#'
#' @param scExp A SingleCellExperiment object 
#' @param merged_peaks A list of GRanges object containing the merged peaks
#' @param geneTSS_annotation A GRanges object with reference genes
#'
#' @return A data.frame with refined annotation
#' 
#' @importFrom GenomicRanges start end union distanceToNearest 
#' pintersect
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom IRanges findOverlapPairs
annotation_from_merged_peaks <- function(scExp,
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
separate_BAM_into_clusters <- function(affectation, odir, merged_bam){
    for (class in levels(factor(affectation$cell_cluster)))
    {
        vec <- affectation$barcode[which(affectation$cell_cluster == 
                                             as.character(class))]
        write(as.vector(vec), 
              file = file.path(odir, paste0(class, ".barcode_class")))
    }
    write(levels(factor(affectation$cell_cluster)),
          file = file.path(odir, "barcodes.barcode_class"))
    
    print("Using Samtools to extract reads from each cell...")
    system(paste0("samtools view -H ", merged_bam, " > ",
                  file.path(odir, "header.sam")))  # keeping header
    
    system(paste0(
        "for i in $(cat ",
        file.path(odir, "barcodes.barcode_class"), "); do samtools view -h ", 
        merged_bam, " | fgrep -w -f ", file.path(odir, "/$i.barcode_class"),
        " > ",file.path(odir, "$i.sam"), ";done")
    )    
    
    system(paste0(
        "for i in $(cat ",
        file.path(odir, "barcodes.barcode_class"), "); do cat ", 
        file.path(odir, "header.sam"), " ", file.path(odir, "$i.sam"),
        " | samtools view -b - > ", file.path(odir, "$i.bam"), " ; done"))
}
