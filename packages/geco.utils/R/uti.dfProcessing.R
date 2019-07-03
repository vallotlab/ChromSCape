# Rd
# description >> Function to annotate positions in a data frame (dfPos) using segments in another data frame (dfSegm)
# argument
# item >> dfPos >> Data frame with positions to annotate
# item >> dfPos.chrom.col >> Chromosome column in dfPos
# item >> dfPos.pos.col >> Position column in dfPos
# item >> dfSegm >> Data frame with segments to use for annotating dfPos
# item >> dfSegm.chrom.col >> Chromosome column in dfSegm
# item >> dfSegm.start.col >> Start position column in dfSegm
# item >> dfSegm.end.col >> End position column in dfSegm
# item >> colsToAdd >> Names of columns in dfSegm that should be used to annotate dfPos
# item >> namesColsToAdd >> Column names to give in the columns added to dfPos
# item >> multseg >> How should the function behave when a single position belongs to several segmengs: NA returns NA; "first" returns only the first segment "all" returns all segments separated by a comma
# value >> annotated dfPos
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.dfPosXSegm <- function(dfPos=NULL,
							dfPos.chrom.col="chrom",
							dfPos.pos.col="pos",
							dfSegm=NULL,
							dfSegm.chrom.col="chrom",
							dfSegm.start.col="start",
							dfSegm.end.col="end",
							colsToAdd=NULL,
							namesColsToAdd=NULL,
							multseg=c(NA,"first","all")[1]
						   )
{
	for(col in namesColsToAdd)	dfPos[,col] <- NA
	dfSegm. <- split(dfSegm,dfSegm[,dfSegm.chrom.col])
	dfPos. <- split(dfPos,dfPos[,dfPos.chrom.col])
	if(is.na(multseg)){
		for(chr in intersect(names(dfSegm.),names(dfPos.)))
		{
			print(chr)
			ind <- unlist(sapply(dfPos.[[chr]][,dfPos.pos.col],function(pos){
				tmp <- which(dfSegm.[[chr]][,dfSegm.start.col] <= pos & dfSegm.[[chr]][,dfSegm.end.col] >= pos)
				if(length(tmp)==1)	tmp else{NA}
			}))
			dfPos.[[chr]][,namesColsToAdd] <- dfSegm.[[chr]][ind,colsToAdd]
		}
	}else{
		if(multseg=="first"){
			for(chr in intersect(names(dfSegm.),names(dfPos.)))
			{
				print(chr)
				ind <- unlist(sapply(dfPos.[[chr]][,dfPos.pos.col],function(pos){
					tmp <- which(dfSegm.[[chr]][,dfSegm.start.col] <= pos & dfSegm.[[chr]][,dfSegm.end.col] >= pos)[1]
				}))
				dfPos.[[chr]][,namesColsToAdd] <- dfSegm.[[chr]][ind,colsToAdd]
			}
		}	
		if(multseg=="all"){
			for(chr in intersect(names(dfSegm.),names(dfPos.)))
			{
				print(chr)
				ind <- sapply(dfPos.[[chr]][,dfPos.pos.col],function(pos){
					tmp <- which(dfSegm.[[chr]][,dfSegm.start.col] <= pos & dfSegm.[[chr]][,dfSegm.end.col] >= pos)
					tmp
				})
				for(j in 1:length(colsToAdd)){
					dfPos.[[chr]][,namesColsToAdd[j]] <- unlist(lapply(ind,function(z)	paste(dfSegm.[[chr]][z,colsToAdd[j]],collapse=",")))
					dfPos.[[chr]][which(sapply(ind,length)==0),namesColsToAdd[j]] <- NA
				}
			}
		}
	}
	dfPos <- unsplit(dfPos.,dfPos[,dfPos.chrom.col])
	dfPos
}

