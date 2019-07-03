# Rd
# description >> Function to export a matrix of numerical values in a file readable withIGV
# argument
# item >> mat >> Matrix of numerical values to export
# item >> feature >> Feature annotation table
# item >> chrom.col >> Name of chromosome column in feature table
# item >> start.col >> Name of start position column in feature table
# item >> end.col >> Name of end position column in feature table
# item >> probe.col >> Name of probe ID column in feature table
# item >> graphType >> Type of graph (bar, points or heatmap)
# item >> color >> RGB color for positive values
# item >> altColor >> RGB color for negative values
# item >> igvfile >> File to which the output should be written
# value >> No value returned
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.exportMatrixIGV <- function(mat=NA,
							    feature=NA,
							    chrom.col=NA,
							    start.col=NA,
							    end.col=NA,
							    probe.col=NA,
							    graphType="points",
							    color="0,150,50",
							    altColor="0,0,0",
							    igvfile=NA
							    )
{		
		write(paste("#track graphType=", graphType," color=",color," altColor=",altColor,sep=""),igvfile)
		igv <- data.frame(Chromosome=feature[,chrom.col],Start=feature[,start.col],End=feature[,end.col],Feature=feature[,probe.col],mat,check.names=F)
		mem <- colnames(igv)
		igv <- geco.genomOrder(igv,chrom="Chromosome",pos="Start")
		write.table(igv[,mem],igvfile,sep="\t",row.names=F,quote=F,append=T)
}							  




# Rd
# description >> Function to export sample annotations in IGV format
# argument
# item >> annot >> Sample annotation table
# item >> sample.col >> Name of column containing sample IDs in annot table
# item >> igvfile >> File to which the output should be written
# value >> No value returned
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.exportAnnotIGV <- function(annot=NA, 
							    sample.col="ID", 
							    igvfile=NA 
							   )
{
		if(!sample.col %in% colnames(annot)){print(paste(sample.col,"is not in annotation file."))}else{
			if(colnames(annot)[1]!=sample.col){annot <- annot[,c(sample.col,setdiff(colnames(annot),sample.col))]}
		}
		write.table(annot,igvfile,sep="\t",row.names=F,quote=F)
}							   




# Rd
# description >> Function to calculate the mean DNA methylation level
# argument
# item >> bval >> Beta value matrix
# item >> feature >> Feature annotation table
# value >> Data frame giving the mean beta value for each region of each gene
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.BvalSummaryByGene <- function(bval=NA, 
								   feature=NA
								  )
{
	if(!identical(rownames(bval),rownames(feature))){print("ERROR: Rownames of bval and feature must be identical.");break}
	
	allgenes <- unique(unlist(strsplit(feature$UCSC_REFGENE_NAME,split=";")))
	allgroups <- unique(unlist(strsplit(feature$UCSC_REFGENE_GROUP,split=";")))
	Gene <- Group <- Nprobe <- c();Bval <- matrix(nrow=0,ncol=ncol(bval),dimnames=list(c(),colnames(bval)))
	for(ge in allgenes)
	{
		feature. <- feature[grep(paste(paste("^",ge,"$",sep=""),paste("^",ge,";",sep=""),paste(";",ge,"$",sep=""),paste(";",ge,";",sep=""),sep="|"),feature$UCSC_REFGENE_NAME),]
		mygeneGROUP <- sapply(1:nrow(feature.),function(i){
			ungene <- unlist(strsplit(feature.[i,"UCSC_REFGENE_NAME"],";"))
			ungroup <- unlist(strsplit(feature.[i,"UCSC_REFGENE_GROUP"],";"))
			paste(ungroup[which(ungene==ge)],collapse=";")
		})
		for(gp in allgroups)
		{
			ind <- rownames(feature.)[grep(gp,mygeneGROUP)]
			Gene <- c(Gene,ge);Group <- c(Group,gp);Nprobe <- c(Nprobe,length(ind))
			if(length(ind)){Bval <- rbind(Bval,apply(matrix(bval[ind,],nrow=length(ind)),2,mean,na.rm=T))}
		}
	}
	data.frame(Gene,Group,Nprobe,Bval)
}




