
#odir must be a full path
#inputBam must be a vector
#stat.value must be either '-p <thresh>' or '-q <thresh>'

subsetBamAndCallPeaks <- function(affectation, annotFeat, odir, inputBam, stat.value=" -p 0.05 ", anno_id,peak_distance_to_merge=20000){
  withProgress(message='Preparing peak data...', value=0, {
    
    incProgress(amount=0.1, detail=paste("merging BAM files"))
    if(length(inputBam) > 1) {
      write(inputBam, file=file.path(odir, "bam_list.txt"))
      system(paste0('samtools merge -@ 4 -f -h ', inputBam[1],' -b ', file.path(odir, "bam_list.txt"), ' ', file.path(odir, 'merged.bam')))
      merged_bam=file.path(odir, 'merged.bam')
    }else{
      merged_bam=inputBam[1]
    }

    for(class in levels(factor(affectation$ChromatinGroup))){
      write(as.vector(affectation$barcode[which(affectation$ChromatinGroup == as.character(class))]), file=file.path(odir, paste0(class, ".barcode_class")))
    }
    write(levels(factor(affectation$ChromatinGroup)),file = file.path(odir, "barcodes.barcode_class"))

    incProgress(amount=0.2, detail=paste("subsetting BAM files"))
    system(paste0('samtools view -H ', merged_bam, ' > ', file.path(odir, 'header.sam')))
    system(paste0('for i in $(cat ', file.path(odir, 'barcodes.barcode_class'), '); do samtools view -h ',merged_bam,' | fgrep -w -f ', file.path(odir,'/$i.barcode_class'), ' > ', file.path(odir,'$i.sam'), ';done'))

    #Reconvert to bam
    system(paste0('for i in $(cat ', file.path(odir,'barcodes.barcode_class'), '); do cat ', file.path(odir, 'header.sam'), ' ', file.path(odir, '$i.sam'), ' | samtools view -b - > ', file.path(odir, '$i.bam'), ' ; done'))

    #BamCoverage
    # system(paste0('for i in $(cat ', file.path(odir,'barcodes.barcode_class'), '); do samtools index ', file.path(odir,'$i.bam'), '; done'))
    # system(paste0('for i in $(cat ', file.path(odir,'barcodes.barcode_class'), '); do bamCoverage --bam ', file.path(odir,'$i.bam'), ' --outFileName ', file.path(odir,'$i.bw'), ' --binSize 50 --smoothLength 500 --extendReads 150 --ignoreForNormalization chrX --numberOfProcessors 4 --normalizeUsing RPKM; done'))

    system(paste0('rm ', file.path(odir,'*.barcode_class'), ' ', file.path(odir,'*.sam')))

    #Peak calling with macs2

    #Count properly paired mapped reads
    for(class in levels(factor(affectation$ChromatinGroup))){
      incProgress(amount=(0.4/length(levels(factor(affectation$ChromatinGroup)))), detail=paste("calling peaks for cluster", class))

      print(paste0('samtools flagstat ',file.path(odir, paste0(class,".bam")), ' | grep "properly paired" | sed "s/.*properly paired (//g" | cut -f1 -d" " | sed "s/\\%//g"'))
      percent_properlyPaired = as.double(system(paste0('samtools flagstat ',file.path(odir, paste0(class,".bam")),
                                                       ' | grep "properly paired" | sed "s/.*properly paired (//g" | cut -f1 -d" " | sed "s/\\%//g"'),intern = T))

      #If there are enough paired and mapped reads -> use the model
        if(percent_properlyPaired> 50){
        print("Bam files provided have more than 50% of their reads properly mapped, running macs2 with default parameters...")
        macs2_options=""
      }
      #If there are not enough properly paired reads use --nomodel --extsize 300 & flag all reads as single end
      else{
        print("Bam files provided have less than 50% of their reads properly mapped, flagging all reads as single end running macs2 with --nomodel --extsize 300...")
        macs2_options=" --nomodel --extsize 300 "

        #Transform the bam file so that all mapped reads are considered as single-end

          system(paste0('samtools view -H ', file.path(odir, paste0(class,".bam")),' > ', file.path(odir, 'header.sam')))
          system(paste0('samtools view -F4 ', file.path(odir, paste0(class,".bam")),' | awk -v OFS="\t" "{\\$2=0; print \\$0}" >> ', file.path(odir, 'header.sam')))
          system(paste0('samtools view -b ', file.path(odir, 'header.sam'),' > ', file.path(odir, paste0(class,".bam"))))

      }
      system(paste0('macs2 callpeak ', stat.value, macs2_options,' --keep-dup all --broad -t ', file.path(odir, paste0(class,".bam")), " --outdir ", odir," --name ",class))
      system(paste0("awk -v OFS='\t' '{if(($3-$2) >= 500){print $0}}' ", file.path(odir, paste0(class, '_peaks.broadPeak')), ' > ', file.path(odir, paste0(class, '.bed'))))
      system(paste0('bedtools merge -delim "\t" -d ',peak_distance_to_merge,' -i ', file.path(odir, paste0(class, '.bed')), ' > ', file.path(odir, paste0(class, '_merged.bed'))))
      
    }

    #Extract data for peak model plots
    # for(class in levels(factor(affectation$ChromatinGroup))){
    #   for(dat in c("p", "m", "xcorr", "ycorr")){
    #     system(paste(file.path("Modules", "get_pc_plotData.sh"), class, dat, odir, .Platform$file.sep))
    #   }
    # }
    
    #Clean up files
    unlink(file.path(odir, "bam_list.txt"))
    # unlink(file.path(odir, "*.bam"))
    # unlink(file.path(odir, "*.bam.bai"))
    unlink(file.path(odir, "*.xls"))
    unlink(file.path(odir, "*.gappedPeak"))
    unlink(file.path(odir, "*_model.r"))
    unlink(file.path(odir, "header.sam"))

    #call makePeakAnnot file
    incProgress(amount=0.1, detail=paste("annotating peaks"))
    mergeBams <- paste(sapply(levels(factor(affectation$ChromatinGroup)), function(x){ file.path(odir, paste0(x, "_merged.bed")) }), collapse = ';')
    mergeBams <- paste0('"', mergeBams, '"')

    write.table(annotFeat[,2:4],file.path(odir,"segmentation_file.bed"),quote = F,sep = "\t",row.names = F,col.names = F)
    system(paste0('bedtools sort -i ',file.path(odir,"segmentation_file.bed"),' > ',file.path(odir,"segmentation_file.sorted.bed")))
    unlink(file.path(odir, "segmentation_file.bed"))

    system(paste("bash", file.path("Modules", "makePeakAnnot.sh"), mergeBams, file.path("annotation", anno_id, "chrom.sizes.bed"), file.path(odir,"segmentation_file.sorted.bed"), file.path("annotation", anno_id, "Gencode_TSS_pc_lincRNA_antisense.bed"), paste0('"', odir, .Platform$file.sep, '"')))
    
    
    incProgress(amount=0.2, detail=paste("finished"))
  })
}
