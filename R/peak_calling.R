## Authors : Pacôme Prompsy, Celine Vallot

#' Peak calling on cell clusters
#' 
#' This functions does peak calling on each cell population in order to refine gene annotation
#' for large bins. For instance, a 50000bp bin might containt the TSS of several genes, while 
#' in reality only one or two of these genes are overlapping the signal (peak). To do so, first
#' in-silico cell sorting is applied based on previously defined clusters contained in the SingleCellExperiment.
#' Taking BAM files of each sample as input, samtools pools then splits reads from 
#' each cell barcode into 1 BAM file per cell cluster (pseudo-bulk). Then MACS2 calls 
#' peaks on each cluster. The peaks are aggregated and merged if closer to a certain distance 
#' defined by user [default to 10000bp]. Then, 
#' 
#' This function takes as input a SingleCellExperiment, that must contain 
#' a 'cell_cluster' column in it's colData, an output directory where to 
#' store temporary files, the list of BAM files corresponding to each sample and 
#' containing the cell barcode information as a tag (for instance tag CB:Z:xxx, XB:Z:xxx or
#'  else...), the p.value used by MACS2 to distinguish significant peaks, the reference
#'  genome (either hg38 or mm10), the maximal merging distance in bp and a data.frame containing
#'  gene TSS genomic cooridnates of corresponding genome (if set to NULL, will automatically load
#'  geneTSS). 
#'  The output is a SingleCellExperiment with GRanges object containing ranges of each merged peaks that falls within
#'  genomic bins of the SingleCellExperiment, saving the bin range as additional column (window_chr,
#'  window_start, window_end), as well as the closests genes and their distance relative to the peak.
#'  The peaks may be present in several rows if multiple genes are close / overlap to the peaks.
#'  
#'  Note that the user must have samtools & MACS2 installed and available in the PATH.
#'  Users can open command terminal and type 'which samtools' & 'which macs2' to verify
#'  the availability of these programs. Will only work on unix operating system. Check 
#'  operating system with 'print(.Platform[1])'.
#' 
#'
#' @param scExp 
#' @param odir Output directory where to write temporary files and each cluster's BAM file
#' @param inputBam A character vector of file paths to each sample's BAM file, containing cell
#' barcode information as tags. BAM files can be paired-end or single-end.
#' @param p.value a p-value to use for MACS2 to determine significant peaks. [0.05]
#' @param ref A reference genome, either hg38 or mm10. ['hg38']
#' @param peak_distance_to_merge Maximal distance to merge peaks together after peak calling
#' , in bp. [10000]
#' @param geneTSS_annotation A data.frame annotation of genes TSS. If NULL will automatically load 
#' Gencode list of genes fro specified reference genome.
#'
#' @return A SingleCellExperiment with refinded annotation
#' @export
#'
#' @examples
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom GenomicRanges GRanges start end union pintersect distanceToNearest ranges
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom IRanges findOverlapPairs subsetByOverlaps
#' @importFrom S4Vectors mcols
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
    } else
    {
        geneTSS_annotation = as(geneTSS_annotation, "GRanges")
    }
    print("Subsetting BAM files...")
    if (length(inputBam) > 1)
    {
        write(inputBam, file = file.path(odir, "bam_list.txt"))
        system(paste0("samtools merge -@ 4 -f -h ", inputBam[1], " -b ", file.path(odir, 
            "bam_list.txt"), " ", file.path(odir, "merged.bam")))
        merged_bam = file.path(odir, "merged.bam")
    } else
    {
        merged_bam = inputBam[1]
    }
    
    print("Writting barcode files...")
    affectation = SingleCellExperiment::colData(scExp)
    affectation$cell_cluster = as.factor(affectation$cell_cluster)
    
    p.value = paste0(" -p ", p.value, " ")  # format for system call to macs2
    
    for (class in levels(factor(affectation$cell_cluster)))
    {
        write(as.vector(affectation$barcode[which(affectation$cell_cluster == 
            as.character(class))]), file = file.path(odir, paste0(class, ".barcode_class")))
    }
    write(levels(factor(affectation$cell_cluster)), file = file.path(odir, "barcodes.barcode_class"))
    
    print("Using Samtools to extract reads from each cell...")
    system(paste0("samtools view -H ", merged_bam, " > ", file.path(odir, "header.sam")))  # keeping header
    
    system(paste0("for i in $(cat ", file.path(odir, "barcodes.barcode_class"), "); do samtools view -h ", 
        merged_bam, " | fgrep -w -f ", file.path(odir, "/$i.barcode_class"), " > ", 
        file.path(odir, "$i.sam"), ";done"))
    # grepping the cell barcode, as a word, pattern from X.barcode_class file.  BAM
    # file must contain cell barcode as a word (for instance tag CB:Z:, XB:Z: or
    # else...)
    
    # Reconvert to bam
    system(paste0("for i in $(cat ", file.path(odir, "barcodes.barcode_class"), "); do cat ", 
        file.path(odir, "header.sam"), " ", file.path(odir, "$i.sam"), " | samtools view -b - > ", 
        file.path(odir, "$i.bam"), " ; done"))
    
    # BamCoverage system(paste0('for i in $(cat ',
    # file.path(odir,'barcodes.barcode_class'), '); do samtools index ',
    # file.path(odir,'$i.bam'), '; done')) system(paste0('for i in $(cat ',
    # file.path(odir,'barcodes.barcode_class'), '); do bamCoverage --bam ',
    # file.path(odir,'$i.bam'), ' --outFileName ', file.path(odir,'$i.bw'), '
    # --binSize 50 --smoothLength 500 --extendReads 150 --ignoreForNormalization chrX
    # --numberOfProcessors 4 --normalizeUsing RPKM; done'))
    
    # system(paste0("rm ", file.path(odir, "*.barcode_class"), " ", file.path(odir, 
    #     "*.sam")))
    
    # Peak calling with macs2
    print("Using MACS2 to call peaks in each cluster...")
    merged_peaks = list()
    ref_chromosomes = GRanges(eval(parse(text = paste0(ref, ".chromosomes"))))
    # Count properly paired mapped reads
    for (class in levels(factor(affectation$cell_cluster)))
    {
        
        print(paste0("samtools flagstat ", file.path(odir, paste0(class, ".bam")), 
            " | grep \"properly paired\" | sed \"s/.*properly paired (//g\" | cut -f1 -d\" \" | sed \"s/\\%//g\""))
        percent_properlyPaired = as.double(system(paste0("samtools flagstat ", file.path(odir, 
            paste0(class, ".bam")), " | grep \"properly paired\" | sed \"s/.*properly paired (//g\" | cut -f1 -d\" \" | sed \"s/\\%//g\""), 
            intern = TRUE))
        
        # If there are enough paired and mapped reads -> use the model
        if (!is.na(percent_properlyPaired) & percent_properlyPaired > 50)
        {
            print("Bam files provided have more than 50% of their reads properly mapped, running macs2 with default parameters...")
            macs2_options = ""
        } else
        {
            print("Bam files provided have less than 50% of their reads properly mapped, flagging all reads as single end running macs2 with --nomodel --extsize 300...")
            macs2_options = " --nomodel --extsize 300 "
            
            # Transform the bam file so that all mapped reads are considered as single-end
            print("Transforming bam so that all mapped reads are considered as single-end...")
            system(paste0("samtools view -H ", file.path(odir, paste0(class, ".bam")), 
                " > ", file.path(odir, "header.sam")))
            system(paste0("samtools view -F4 ", file.path(odir, paste0(class, ".bam")), 
                " | awk -v OFS=\"\t\" \"{\\$2=0; print \\$0}\" >> ", file.path(odir, 
                  "header.sam")))
            system(paste0("samtools view -b ", file.path(odir, "header.sam"), " > ", 
                file.path(odir, paste0(class, ".bam"))))
            
        }
        # MACS2 CALL
        print(paste0("macs2 callpeak ", p.value, macs2_options, " --keep-dup all --broad -t ", 
            file.path(odir, paste0(class, ".bam")), " --outdir ", odir, " --name ", 
            class))
        system(paste0("macs2 callpeak ", p.value, macs2_options, " --keep-dup all --broad -t ", 
            file.path(odir, paste0(class, ".bam")), " --outdir ", odir, " --name ", 
            class))
        
        # Formatting and adding peaks to list
        peaks = read.table(file = file.path(odir, paste0(class, "_peaks.broadPeak")), 
            colClasses = c("character", "integer", "integer", rep("NULL", 6)))
        colnames(peaks) = c("chr", "start", "end")
        peaks = GenomicRanges::GRanges(peaks)
        peaks = peaks[which(width(GenomicRanges::ranges(peaks)) >= 500), ]
        peaks = GenomicRanges::reduce(peaks, min.gapwidth = peak_distance_to_merge, 
            ignore.strand = TRUE)
        peaks = suppressWarnings(IRanges::subsetByOverlaps(peaks, ref_chromosomes, 
            ignore.strand = TRUE))
        merged_peaks[[class]] = peaks
    }
    
    # Clean up files unlink(file.path(odir, 'bam_list.txt'))
    # unlink(file.path(odir, "*.xls"))
    # unlink(file.path(odir, "*.gappedPeak"))
    # unlink(file.path(odir, "*_model.r"))
    # unlink(file.path(odir, "header.sam"))
    # unlink(file.path(odir, 'merged.bam'))
    # 
    # call makePeakAnnot file
    print("Merging BAM files together...")
    
    segmentation = SummarizedExperiment::rowRanges(scExp)
    S4Vectors::mcols(segmentation) = NULL
    segmentation$window_ID = paste(as.character(segmentation@seqnames), GenomicRanges::start(segmentation), 
        GenomicRanges::end(segmentation), sep = "_")
    
    merged_peak = merged_peaks[[1]]
    for (i in 2:length(merged_peaks))
    {
        merged_peak = suppressWarnings(GenomicRanges::union(merged_peak, merged_peaks[[i]], 
            ignore.strand = TRUE))
    }
    
    pairs <- IRanges::findOverlapPairs(segmentation, merged_peak, ignore.strand = TRUE)
    refined_annotation = GenomicRanges::pintersect(pairs, ignore.strand = TRUE)
    S4Vectors::mcols(refined_annotation)$hit = NULL
    
    hits_genes = GenomicRanges::distanceToNearest(refined_annotation, geneTSS_annotation, 
        ignore.strand = TRUE, select = "all")
    
    refined_annotation = refined_annotation[S4Vectors::queryHits(hits_genes)]
    refined_annotation$Gene = as.character(geneTSS_annotation$gene[S4Vectors::subjectHits(hits_genes)])
    refined_annotation$distance = hits_genes@elementMetadata$distance
    
    refined_annotation$peak_ID = paste(as.character(refined_annotation@seqnames), 
        GenomicRanges::start(refined_annotation), GenomicRanges::end(refined_annotation), 
        sep = "_")
    scExp@metadata$refined_annotation = refined_annotation
    
    return(scExp)
}
