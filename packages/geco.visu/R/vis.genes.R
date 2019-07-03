# Rd
# description >> Function to plot genes in the UCSC style for figures of peak regions for instance.
# argument
# item >> refgene >> Table of Refseq genes (for example "UCSC.hg19.refGene.RData")
# item >> genes >> Character vector indicating which genes should be represented.
# item >> boundaries >> Boundaries of the region to represent
# item >> alternate >> Should genes be represented in 2 vertical layers (for dense regions) ?
# item >> mycolfill >> Color inside the rectangeles representing exons.
# item >> mycolborder >> Color of the border of rectangles.
# item >> myylim >> Parameter ylim for the plot
# value >> None.
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.plotGenes <- function(refgene=NULL,
						  genes=NULL,
						  boundaries=NULL,
						  alternate=0,
						  mycolfill="lightblue",
						  mycolborder="blue",
						  myylim=NULL)
{
	
	if(is.null(myylim)){myylim <- c(0,1)}
	plot(-10,xlim=boundaries,ylim=myylim,axes=F,xlab="",ylab="")
	
	genes <- genes[order(refgene$cdsStart[match(genes,refgene$name2)])]
	
	for(g in genes)
	{
		refseq=refgene[which(refgene$name2==g),]
		exonstart=as.numeric(unlist(strsplit(paste(refseq$exonStarts,collapse=""),split=",")))
		exonend=as.numeric(unlist(strsplit(paste(refseq$exonEnds,collapse=""),split=",")))
		S=exonstart[1];E=exonend[1]
		if(length(exonstart) > 1)
		{
			for(i in 2:length(exonstart))
			{
  				if(!any(S==exonstart[i] & E==exonend[i])){S=c(S,exonstart[i]);E=c(E,exonend[i])}
			}
		}
		unikex=data.frame(S,E);unikex=unikex[sort(unikex$S,index.return=T)$ix,]
		for(i in 1:nrow(unikex))
		{
  			rect(unikex$S[i],0.4+alternate,unikex$E[i],0.6+alternate,col=mycolfill,border=mycolborder,lwd=2)
		}
		if(nrow(unikex) > 1)
		{
			for(i in 1:(nrow(unikex)-1))
			{
  				if(unikex$S[i]!=unikex$S[i+1]){segments(unikex$E[i],0.5+alternate,unikex$S[i+1],0.5+alternate,col= mycolborder,lwd=2)}
			}
		}
		if(alternate){alternate <- -alternate}
	}
}



