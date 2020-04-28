context("testing peak calling")

load("../../Data/ChromSCape_Data/Reproducibility_new/datasets/HBCx_95_mm10_002/cor_filtered_data/HBCx_95_mm10_002_1600_1_95_uncorrected_99_1.RData")
scExp = data$scExp_cf
odir = "../../Data/ChromSCape_Data/Reproducibility_new/datasets/HBCx_95_mm10_002/peaks/HBCx_95_mm10_002_1600_1_95_uncorrected_99_1_k2/"
inputBam = c("/media/pacome/LaCie/InstitutCurie/Documents/Data/results/Matrices_mm10/GEO/HBCx_95_H3K27me3_flagged_rmPCR_RT_rmDup.bam",
               "/media/pacome/LaCie/InstitutCurie/Documents/Data/results/Matrices_mm10/GEO/HBCx_95_CapaR_H3K27me3_flagged_rmPCR_RT_rmDup.bam")
stat.value = " -p 0.05 "
ref = "mm10"
peak_distance_to_merge = 5000
geneTSS_annotation = NULL

segmentation = as(read.table("../../Data/ChromSCape_Data/MM468/old/datasets/MM468_exp3/peaks/MM468_exp3_2000_1_95_uncorrected_99_1_k3/segmentation_file.sorted.bed",
                ,col.names = c("chr","start","end")),"GRanges")

merged_beds = list()
ref_chromosomes = GRanges(eval(parse(text = paste0(ref,".chromosomes"))))
for(i in c("C1","C2","C3")){
  bed = read.table(paste0("../../Data/ChromSCape_Data/MM468/old/datasets/MM468_exp3/peaks/MM468_exp3_2000_1_95_uncorrected_99_1_k3/",i,"_peaks.broadPeak"),
                   colClasses = c("character","integer","integer", rep("NULL",6)), sep ="\t")
  
  
  colnames(bed) = c("chr","start","end")
  bed = GRanges(bed)
  bed = bed[which(width(ranges(bed)) >= 500),]
  bed = reduce(bed,min.gapwidth = 5000)
  bed = subsetByOverlaps(bed,ref_chromosomes, ignore.strand = TRUE)
  
  merged_beds[[i]] = bed 
  
}
