# Rd
# description >> Function to convert a table of annotations to a table of colors
# argument
# item >> annotT >> Annotation table
# item >> annotS >> selection of interest in Annotation Table
# item >> missing >> Fields that should be considered as missing data
# item >> anotype >> Vector indicating the type of each column ("categ","binary","quantit"). If not provided, categories will be determined automatically
# item >> maxnumcateg >> Maximum number of numerical values to consider a numeric column as categorical
# item >> categcol >> Colors for categorical variables (optional)
# item >> quantitcol >> Colors for quantitive variables (optional)
# item >> plotLegend >> Should the legends be plotted in a file?
# item >> plotLegendFile >> Name of the file where the legends should be plotted
# value >> Annotation matrix with color codes
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.annotToCol <- function(annotS=NULL,
							annotT=NULL, #table d'annotation
							missing=c("",NA),
							anotype=NULL,
							maxnumcateg=2,
							categCol=NULL,
							quantitCol=NULL,
							plotLegend=T,
							plotLegendFile=NULL
							)
{
    # Case with only one annotCol
    if(is.null(ncol(annotS))){
    		annotS<-data.frame(annotS)
    		colnames(annotS)=annotCol
    		rownames(annotS)=rownames(annotT)
    		}
    # Set all missing values to NA
    for(j in 1:ncol(annotS))	annotS[which(annotS[,j] %in% missing),j] <- NA
    
    # Get the type of each annotation column
    if(is.null(anotype)){
	    anotype <- rep("categ",ncol(annotS));names(anotype) <- colnames(annotS)
    	classes <- sapply(1:ncol(annotS),function(j) class(annotS[,j]))
	    nmodal <- sapply(1:ncol(annotS),function(j) length(unique(setdiff(annotS[,j],NA))))
    	anotype[which(classes %in% c("integer","numeric") & nmodal > maxnumcateg)] <- "quantit"
		anotype[which(nmodal==2)] <- "binary"
	}
	
	# Convert annotations to colors
	anocol <- annotS
	
	if(plotLegend)	pdf(plotLegendFile)

	# if(is.null(categCol))	categCol <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#bd18ea", "#2ef4ca", "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106", "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776")
	if(is.null(categCol))	categCol <- c("royalblue", "peru", "red", "palegreen4", "skyblue", "sienna2", "slateblue3", "pink2", "slategray", "black", "orange", "turquoise4", "yellow3", "slategray2", "orchid1", "palegreen2", "tan1", "red4", "palevioletred1", "orangered", "snow4", "purple", "sienna4", "turquoise1")


k <- 1
	for(j in which(anotype=="categ")){
		tmp <- as.factor(anocol[,j])
		classes <- as.character(levels(tmp))
		ncat <- length(levels(tmp))
		if(k+ncat > length(categCol))	categCol <- c(categCol,categCol)
		levels(tmp) <- categCol[k:(k+ncat-1)]
		fill <- as.character(levels(tmp))
		anocol[,j] <- as.character(tmp)
		k <- k+ncat
		if(plotLegend){
			par(mar=c(0,0,0,0))
			plot(-10,axes=F,xlim=c(0,5),ylim=c(0,5),xlab="",ylab="")
			legend(1,5,legend=classes,fill=fill,title=colnames(anocol)[j],xjust=0.5,yjust=1)
		}
	}
	
	memcol <- c()
	for(j in which(anotype=="binary")){
		new <- setdiff(anocol[,j],c(NA,memcol))
		if(length(new)==2){memcol <- c(memcol,c("dodgerblue4","firebrick"));names(memcol)[(length(memcol)-1):length(memcol)] <- sort(new)}
		if(length(new)==1){memcol <- c(memcol,setdiff(c("dodgerblue4","firebrick"),memcol[setdiff(anocol[,j],c(NA,new))]));names(memcol)[length(memcol)] <- new}
		# anocol[,j] <- memcol[anocol[,j]] # initial line

					## Modification as bug. 2nd memcol was getting NA then white color
					anocol[,j] <-  as.character(anocol[,j])
					for (z in 1:length(memcol)){
					anocol[which(anocol[,j]==names(memcol)[z]),j] <- memcol[z]
					}
					##########
		
		if(plotLegend){
			par(mar=c(0,0,0,0))
			plot(-10,axes=F,xlim=c(0,5),ylim=c(0,5),xlab="",ylab="")
			classes <- intersect(names(memcol),annotS[,j]);fill <- memcol[classes]
			legend(1,5,legend=classes,fill=fill,title=colnames(anocol)[j],xjust=0.5,yjust=1)
		}		
	}
	
	if(is.null(quantitCol))	quantitCol <- c("orangered1","darkgreen","darkblue","darkred","darkgoldenrod4","darkorchid4","darkolivegreen4","darkorange4","darkslategray")
	k <- 1
	for(j in which(anotype=="quantit")){
		colrange <- colorRampPalette(c("white",quantitCol[k]))(100)
		anocol[,j] <- colrange[round(geco.changeRange(anocol[,j],newmin=1,newmax=100))]
		if(k < length(quantitCol)){k <- k+1}else{k <- 1}
		if(plotLegend){
			par(mar=c(8,2,5,1))
			lims <- seq(-1,1,length.out=200)
			image(matrix(lims,nc=1),col= colrange,axes=F,xlab=colnames(anocol)[j])
		}
	}

	if(plotLegend)	dev.off()

	for(j in 1:ncol(anocol))	anocol[which(is.na(anocol[,j])),j] <- "white"
	as.matrix(anocol)
}



# Rd
# description >> Function to make a visual representation of a matrix containing color codes
# argument
# item >> matcol >> Color matrix
# item >> strat >> Column of matrix used to draw split lines
# item >> xlab >> Size of x axis labels (if 0 no labels are drawn)
# item >> ylab >> Size of y axis labels (if 0 no labels are drawn)
# item >> drawLines >> One of "none" (no line drawn), "h" (horizontal lines), "v" (vertical lines) or "b" (both)
# value >> None
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.imageCol <- function(matcol=NULL,
						  strat=NULL,
						  xlab.cex=0.5,
						  ylab.cex=0.5,
						  drawLines=c("none","h","v","b")[1],
						  ...
						  )
{
    if(is.null(ncol(matcol))){
    		matcol<-data.frame(matcol)
    		colnames(matcol)=colnames(anocol)
    		}
    	matcol <- matcol[,ncol(matcol):1]
     if(is.null(ncol(matcol))){
    		matcol<-data.frame(matcol)
    		colnames(matcol)=colnames(anocol)
    		}
    csc <- matcol
    csc.colors <- matrix()
    csc.names <- names(table(csc))
    csc.i <- 1
    for(csc.name in csc.names){
       csc.colors[csc.i] <- csc.name
       csc[csc == csc.name] <- csc.i
       csc.i <- csc.i + 1
    }
  
    if(dim(csc)[2]==1){
    		csc<-matrix(as.numeric(unlist(csc)), nrow = dim(csc)[1])
    		}	else {
    		csc <- matrix(as.numeric(csc), nrow = dim(csc)[1])
    		}
    
    image(csc, col = as.vector(csc.colors), axes = FALSE, ...)
    if(xlab.cex!=0){
        axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1),colnames(matcol),las = 2,tick = FALSE,cex.axis=xlab.cex, ...)
    }
    if(ylab.cex!=0){
        axis(3, 0:(dim(csc)[1] - 1)/(dim(csc)[1] - 1),rownames(matcol),las = 2,tick = FALSE,cex.axis=ylab.cex, ...)
    }
    if(drawLines %in% c("h","b"))	abline(h=-0.5:(dim(csc)[2] - 1)/(dim(csc)[2] - 1));box()
    if(drawLines %in% c("v","b"))	abline(v=0.5:(dim(csc)[1] - 1)/(dim(csc)[1] - 1));box()
    if(!is.null(strat)){
    	z <- factor(matcol[,strat]);levels(z) <- 1:length(levels(z))
    	z <- geco.vectorToSegments(as.numeric(z))
    	abline(v=geco.changeRange(c(0.5,z$Ind_K+0.5)/max(z$Ind_K),newmin=par()$usr[1],newmax=par()$usr[2]),lwd=2,lty=2)
    }	
}       




# Rd
# description >> This function tests the association between all colunms in an annotation table vs 2 or more groups
# argument
# item >> groups >> Vector giving the group of each sample (ordered like annot rows)
# item >> annotT >> Annotation table
# item >> annotS >> selection of interest in Annotation Table
# item >> groupOrder >> Ordering of the groups for the outputs (optional)
# item >> orientated >> Should tests for quantitative variables be orientated?
# value >> List containing a pval table (each annot vs all groups), a pvalgp matrix (each annot vs each group) and a pvalgpcol matrix (equivalent to pvalgp with a color code)
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.groupsVSannot <- function(groups=NULL,
							   annotS=NULL,
							   annotT=NULL,
							   groupOrder=NA,
#							   maxnumcateg=10,
							   orientated=F
							   )
{
    # Case with only one annotCol
    if(is.null(ncol(annotS))){
    		annotS<-data.frame(annotS)
    		colnames(annotS)=annotCol
    		rownames(annotS)=rownames(annotT)
    		}
    # Get the type of each annotation column
    anotype <- rep("categ",ncol(annotS));names(anotype) <- colnames(annotS)
    classes <- sapply(1:ncol(annotS),function(j) class(annotS[,j]))
	nmodal <- sapply(1:ncol(annotS),function(j) length(unique(setdiff(annotS[,j],NA))))
    anotype[which(classes %in% c("integer","numeric"))] <- "quantit" # & nmodal > maxnumcateg
    
    # Compute one pvalue per annot
    pval <- rep(NA,ncol(annotS));names(pval) <- colnames(annotS)
    test <- rep("Chi2",ncol(annotS));names(pval) <- colnames(annotS)
	for(j in which(anotype=="categ")){try(pval[j] <- chisq.test(groups,annotS[,j])$p.value)}
    ind <- which(anotype=="quantit")    
    if(length(unique(groups)) > 2){
    	test[ind] <- "ANOVA"
    	for(j in ind){try(pval[j] <- anova(lm(annotS[,j] ~ groups))$Pr[1])}
    }else{
    	test[ind] <- "Wilcoxon"
    	for(j in ind){try(pval[j] <- wilcox.test(annotS[,j] ~ groups)$p.value)}
    }
    pval <- data.frame(anotype,test,pval)	
    
    # Compute pvalues for each group separately
    if(is.na(groupOrder)){groupOrder <- names(table(groups))}
    if(orientated){myalt <- "greater"}else{myalt="two.sided"}
    pvalgp <- matrix(NA,nrow=length(groupOrder),ncol=ncol(annotS))
    rownames(pvalgp) <- groupOrder;colnames(pvalgp) <- colnames(annotS)
    for(g in groupOrder)
    {
       ind1 <- which(groups==g);ind0 <- which(groups!=g)
       for(j in which(anotype=="categ")){
            if(length(setdiff(unique(annotS[, j]),NA))==2){try(pvalgp[g,j] <- fisher.test(groups==g,annotS[,j],alternative="greater")$p.value,silent=TRUE)}else{
                 try(pvalgp[g,j] <- chisq.test(groups==g,annotS[,j])$p.value,silent=TRUE)}
       }
       for(j in which(anotype=="quantit")){try(pvalgp[g,j] <- wilcox.test(annotS[,j] ~ (groups==g),alternative=myalt)$p.value)}
    }
    
    pvalgpcol <- pvalgp
    pvalgpcol[which(pvalgp > 0.1)] <- "white"
    pvalgpcol[which(pvalgp <= 0.1 & pvalgp > 0.05)] <- heat.colors(10)[10]
    pvalgpcol[which(pvalgp <= 0.05 & pvalgp > 0.01)] <- heat.colors(10)[6]
    pvalgpcol[which(pvalgp <= 0.01 & pvalgp > 0.001)] <- heat.colors(10)[5]
    pvalgpcol[which(pvalgp <= 0.001 & pvalgp > 1e-4)] <- heat.colors(10)[4]
    pvalgpcol[which(pvalgp <= 1e-4 & pvalgp > 1e-6)] <- heat.colors(10)[3]
    pvalgpcol[which(pvalgp <= 1e-6)] <- heat.colors(10)[1]
    pvalgpcol[which(is.na(pvalgp))] <- "lightgrey"
    
    list(pval=pval,pvalgp=pvalgp,pvalgpcol=pvalgpcol)
}



# Rd
# description >> This function tests the association between all colunms in an annotation table vs one or more numerical vectors
# argument
# item >> num >> List of numerical vectors to be tested
# item >> annotT >> Annotation table
# item >> annotS >> selection of interest in Annotation Table
# value >> pval table (each annot vs all vectors)
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.numVSannot <- function(num=NULL,
							annotT=NULL,
							annotS=NULL,
							plotPval=TRUE
							)
{
    # Case with only one annotCol
    if(is.null(ncol(annotS))){
    		annotS<-data.frame(annotS)
    		colnames(annotS)=annotCol
    		rownames(annotS)=rownames(annotT)
    		}
    # Get the type of each annotation column
    anotype <- rep("categ",ncol(annotS));names(anotype) <- colnames(annotS)
    classes <- sapply(1:ncol(annotS),function(j) class(annotS[,j]))
	nmodal <- sapply(1:ncol(annotS),function(j) length(unique(setdiff(annotS[,j],NA))))
    anotype[which(classes %in% c("integer","numeric"))] <- "quantit" # & nmodal > maxnumcateg
    
    # Compute one pvalue per annot
    pval <- matrix(NA,nrow=ncol(annotS),ncol=length(num),dimnames=list(colnames(annotS),names(num)))
    test <- rep(NA,ncol(annotS));names(test) <- colnames(annotS)
  	ind <- which(anotype=="categ" & nmodal==2)
  	test[ind] <- "Wilcoxon"
  	for(j in ind){try(pval[j,] <- unlist(lapply(num,function(z)	wilcox.test(z ~ annotS[,j])$p.value)))}
  	ind <- which(anotype=="categ" & nmodal > 2)
  	test[ind] <- "ANOVA"
  	for(j in ind){try(pval[j,] <- unlist(lapply(num,function(z)	anova(lm(z~ annotS[,j]))$Pr[1])))}
  	ind <- which(anotype=="quantit")
  	test[ind] <- "Pearson"
  	for(j in ind){try(pval[j,] <- unlist(lapply(num,function(z)	cor.test(z,annotS[,j])$p.value)))}

	if(plotPval)	geco.imagePval(as.matrix(t(pval)))
    
    data.frame(anotype,test,pval)	
}



# Rd
# description >> Function to produce a graphical representation of a pvalue table
# argument
# item >> pval >> P-value table
# item >> ... >> required
# value >> None
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.imagePval <- function(pval=NULL,continuous=F,maxlogP=NULL,...)
{
	if(continuous){
		mycolramp <- colorRampPalette(c("white","lightgoldenrodyellow","orange","red","darkred","purple"))(maxlogP*10+1)
		logp <- as.matrix(-log10(pval))
		logp[which(logp > maxlogP)] <- maxlogP
		logpcol <- logp
		for(j in 1:ncol(logp))	logpcol[,j] <- mycolramp[round(logp[,j]*10)+1]
		logpcol[which(is.na(logp))] <- "lightgrey"
		geco.imageCol(logpcol,...)
	}else{
		pvalcol <- as.matrix(pval)
    	pvalcol[which(pval > 0.1)] <- "white"
    	pvalcol[which(pval <= 0.1 & pval > 0.05)] <- heat.colors(10)[10]
    	pvalcol[which(pval <= 0.05 & pval > 0.01)] <- heat.colors(10)[6]
    	pvalcol[which(pval <= 0.01 & pval > 0.001)] <- heat.colors(10)[5]
    	pvalcol[which(pval <= 0.001 & pval > 1e-4)] <- heat.colors(10)[4]
    	pvalcol[which(pval <= 1e-4 & pval > 1e-6)] <- heat.colors(10)[3]
    	pvalcol[which(pval <= 1e-6)] <- heat.colors(10)[1]
    	pvalcol[which(is.na(pval))] <- "lightgrey"
		geco.imageCol(pvalcol,...)	
	}	
}









































