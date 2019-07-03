# Rd
# description >> convert chromosomes X, Y and MT into numerics 23, 24 and 25 
# argument
# item >> d >> data.frame which columns include genomic position information 
# item >> chrom >> the column in \code{d} giving the chromosome position
# value >> \code{d} with added \code{chromNum} column
# author >> Eric Letouze
# keyword >> methods
# end
geco.chromString2num <- function (d, 
                                  chrom = "chrom",
                          		  chromNum  = "chrNum"
                                 )
{ 
    d[, chromNum] <- as.character(sub("chr","",d[, chrom]))
   
    maxchr <- suppressWarnings( max(as.numeric(d[,chromNum]),na.rm=TRUE) )
    
    w <- which(d[, chromNum] == "X")
    if (length(w) > 0)d[w, chromNum] <- maxchr+1
    
    w <- which(d[, chromNum] == "Y")
    if (length(w) > 0)d[w, chromNum] <- maxchr+2
    
    w <- which(d[, chromNum] %in% c("M","MT"))
    if (length(w) > 0)d[w, chromNum] <- maxchr+3
    
    w <- which(!(d[, chromNum] %in% as.character(1:25)))
    if (length(w) > 0) d[w, chromNum] <- NA   
    
    d[, chromNum] <- as.numeric(d[, chromNum])
    d
}


# Rd
# description >> required
# argument
# item >> d >> a \code{data.frame}
# item >> chrom >> chromosome column name
# item >> pos >> sequence location column name
# item >> startPos >> the column in \code{d} giving the start base pair position
# item >> endPos >> the column in \code{d} giving the end base pair position
# item >> chromNum >> the name of added column whith chromosome from 1 to 24
# value >> \code{d} ordered genomically
# author >> Eric Letouze
# keyword >> methods
# end
geco.genomOrder <- function(d,
                            chrom     = "chrom",
                            pos       = "pos",
                            startPos  = "start",
                            endPos    = "end" ,
                            chromNum  = "chrNum"
                            ) 
{
	d <- factochar(d)
	if(!chromNum %in% names(d))  d <- geco.chromString2num(d,chrom=chrom)

    if (is.null(pos)) {
        pos <- "meanPos"
        if (options("verbose")$verbose) 
           warning("missing value for parameter pos -> using default value (meanPos).")
    }

	if(!pos %in% colnames(d))
	{
        if(all(c(startPos,endPos) %in% colnames(d))){d$meanPos <- rowMeans(d[,c(startPos,endPos)],na.rm=T)}else{
        	stop("Please set valid pos or startPos and endPos values.")
		}
		pos <- "meanPos"
	}
  	ord <- order(d[,chromNum], as.numeric(d[,pos]), na.last = TRUE)
	d <- d[ord,]

	d
}



# Rd
# description >> return pangenomic absolute position, convenient for plots
# argument
# item >> x >> data.frame which columns include genomic position information 
# item >> chrom >> the column in \code{x} giving the chromosome position (in 1:22,X,Y)
# item >> pos >> the column in \code{x} giving the base pair position   (optional)
# item >> startPos >> the column in \code{x} giving the start base pair position
# item >> endPos >> the column in \code{x} giving the end base pair position
# item >> cyto >> cytoband object (optional)
# item >> absPos >> name of added column giving absolute pangenomic order
# item >> chromNum >> the column in \code{x} giving the chromosome position with X = 23 & Y =24
# value >> \code{x} in genormic order, with absolute position (base pair) added
# author >> Eric Letouze
# keyword >> graphs
# end
geco.pangenomCoord <- function( x,
                            chrom="chrom",
                            pos="meanPos",
                            startPos="start",   
                            endPos="end",                         
                            cyto=NULL,
                            absPos="absPos",
                            chromNum="chrNum"
                            )
{
        if(is.null(cyto))  stop("Cytoband object required")
     
         x <- geco.genomOrder(d=x,
                              chrom     = chrom,
                              pos       = pos,
                              startPos  = startPos,
                              endPos    = endPos ,
                              chromNum  = chromNum)   

         if(absPos %in% colnames(x))  return(x)

         if(!endPos %in% colnames(x)){
               endPos <- "endPos"
               if(pos %in% colnames(x)){
                     x[,endPos] <- x[,pos]
                     if(options("verbose")$verbose) warning("uncorrect value for parameter endPos -> using pos value.")
               }else{
                     stop("missing value for pos or endPos.")
               }
         }

         if(is.null(absPos)){         
             absPos <- "absPos"
             if(options("verbose")$verbose) warning("missing value for parameter absPos -> using default value.")
         }
         
         if(!absPos %in% names(x)){
             x$absPos <- as.numeric(x[,pos])
             if(pos %in% colnames(x)){
                   x[,absPos] <- as.numeric(x[,pos])
                   if(options("verbose")$verbose) warning("uncorrect value for parameter absPos -> using pos value.")
             }else{
                  stop("missing value for absPos or pos.")
             }
         }
         
         if(!is.null(cyto)){                        
             if(!"ChrNumeric" %in% names(cyto) ) cyto <-  geco.chromString2num(cyto,chrom="Chromosome",chromNum="ChrNumeric")
             
             allchrnum <- sort(unique(c(cyto$ChrNumeric,x[,chromNum])))
                 
             maxByChr <- sapply(allchrnum,function(z){
                                       m <- 0
                                       w <- which(cyto$"ChrNumeric"==z)
                                       if(length(w)>0) m <- max(cyto[w,"End"],na.rm=TRUE) 
                                       w <- which(x[,chromNum]==z)
                                       if(length(w)>0){
                                        if( which.max(c(m,max(x[w,endPos])))==2 )
                                          warning("In geco.pangenomCoord, some positions in your data are superior to max cytoband in cyto.\ncyto object used does not seem to be the right cytoband object version.\n")
                                        m <- max(c(m,max(x[w,endPos],na.rm=TRUE)))
                                       }
                                       m})

             maxByChr <- cumsum(as.numeric(maxByChr))

         }else{
             maxByChr <- cumsum(as.numeric(sapply(split(x[,endPos],x[,chromNum]),max)))
         }

         
         for(i in setdiff(unique(x[,chromNum]),1)){
              w <- which(x[,chromNum] == i)
              x[w,absPos] <-  x[w,absPos] + maxByChr[i-1]     
         }
         w <- which(is.na(x[,chromNum]))
         if(length(w)>0) x[w,absPos] <- NA
         x
}



# Rd
# description >> pangenomic plot
# argument
# item >> d >> data.frame which columns include genomic position information and the wanted y axis
# item >> ycol >> the column in \code{d} to be used as y axis
# item >> chrom >> the column in  \code{d} giving the chromosome position
# item >> pos >> the column in  \code{d} giving the base pair position   (optional)
# item >> startPos >> the column in  \code{d} giving the start base pair position
# item >> endPos >> the column in  \code{d} giving the end base pair position
# item >> colorscol >> (optional) column in \code{d} giving the color for each point from ycol or a unique color name or number (default black)
# item >> cyto >> cytoband object (optional)
# item >> absPos >> the column in  \code{d} giving the absolute base pair position (will be calculated if not present in \code{d})
# item >> rectangleOption >> boolean : if \code{TRUE} rectangles are plotted instead of points (default) - NB : for points, several 'looks' can be obtained (ex using parameter type , pch,...)
# item >> plotCentro >> boolean : should centromeric delimitation be plotted ?
# item >> plotCentro.color >> color of the centromeric delimitation
# item >> plotChrom >>  boolean : should chromosome delimitation be plotted ?
# item >> plotnew >>  boolean : should the a new graph be plotted (default=\code{TRUE}), or should something be added to the current graph (-> \code{plotnew=FALSE})
# item >> plotyaxis >>  boolean : should an Y axis be plotted ? (considered only if \code{plotnew = TRUE})
# item >> yaxmark >>  y axis 'at' parameter (considered only if \code{plotyaxis = TRUE})
# item >> yaxlab >>  y axis 'labels' parameter (considered only if \code{plotyaxis = TRUE})
# item >> chromToPlot >> a vector of the chromosomes to be plotted (default c(1:22,"X","Y"))
# item >> ... >> required
# value >> \code{d} in genormic order, with absolute position (base pair) added
# author >> Eric Letouze
# keyword >> graphs
# end
geco.pangenomPlot <- function(d=NULL,
                            ycol=NULL, 
                            chrom="chrom",
                            pos="meanPos",
                            startPos="start",
                            endPos="end",
                            colorscol=NULL,
                            cyto=NULL,
                            absPos="absPos",
                            rectangleOption=FALSE,
                            plotCentro = TRUE,
                            plotCentro.color = "lightgrey",
                            plotChrom =TRUE,
                            plotnew=TRUE,
                            plotyaxis=TRUE,
                            plotxaxis=TRUE,
                            yaxmark=NULL,
                            yaxlab=NULL,
                            chromToPlot=c(1:22,"X","Y"),
                            xlim=NULL,
                            CexAxis=0.8,
                            ...)
{                            
       if(! ycol %in% names(d))  stop("ycol ", ycol, " not found in d.\n")
       if(is.null(cyto))  stop("Cytoband object required")
 
       d <- geco.pangenomCoord(d,chrom=chrom,pos=pos,startPos=startPos,endPos=endPos,cyto=cyto)
       y <- d[,ycol]                                        

       if(is.null(colorscol) ){colorY <-"black"}else{colorY <- d[,colorscol]}

       if(!is.null(chromToPlot))	cyto <- cyto[which(cyto$Chromosome %in% chromToPlot),] # added to rm X and Y chromosome from graphics
       cyto <- geco.pangenomCoord (cyto, chrom="ChrNumeric", pos="End", startPos="Start", endPos="End", cyto=cyto, absPos="absPos") 
       if(is.null(xlim)) {xlim <- c(0,max(cyto$absPos,na.rm=T))}
       delta <- diff(xlim)*0.02;xlim <- xlim+c(-delta,delta)
       delta <- diff(xlim)*0.04/0.92;xlim <- xlim+c(delta,-delta)
      
       if(plotnew){  
            plot( d[,absPos],y,axes=FALSE,col=colorY,xlim=xlim, ...)

            box()

            if(plotyaxis){
                okk <- FALSE
                if(is.null(yaxmark)){
                    yaxmark <- pretty(y,n=10)
                    okk <- TRUE
                }
                if(is.null(yaxlab)) yaxlab <- yaxmark
                axis(2, #1=below, 2=left, 3=above and 4=right
                    at = yaxmark,
                    labels = yaxlab, cex.axis=CexAxis, las=2)
            }
       
       }else{
            points( d[,absPos],y, col=colorY, ...)               
       }
       
       wat <- NULL
       if(!is.null(cyto)){
           w <- which(cyto$"Centro"==1)
           wat <- sapply(split(cyto[w,"absPos"],cyto[w,"ChrNumeric"]),mean)
           if(length(w) == 0) wat <- sapply(split(cyto[,"absPos"],cyto[,"ChrNumeric"]),mean)

           if(plotCentro & length(w) != 0)abline(v=wat,lty=3,col=plotCentro.color)
      
           if(plotChrom ){

                if(plotxaxis){
                axis(1, #1=below, 2=left, 3=above and 4=right
                    at = wat,
                    labels = c(1:22,"X","Y")[c(1:22,"X","Y") %in% unique(cyto$Chromosome)],#c(1:22,"X","Y"),
                    cex.axis=CexAxis, las=2)
                }
                
                wat <-  sapply(split(cyto[,"absPos"],cyto$"ChrNumeric"),max)          
                abline(v=wat[-length(wat)],lty=1)
                #if(length(w) == 0) abline(h=0,lty=1)

           } 
       }            
       d                            
}

# Rd
# description >> chromosome-wide plot 
# argument
# item >> d >> data.frame which columns include genomic position information and the wanted \code{y} axis
# item >> ycol >> the column in \code{d} to be used as \code{y} axis
# item >> ytransform >> a transform to be applied to \code{y} in the plot
# item >> reverseYaxis >> boolean (default=\code{FALSE}) : use \code{TRUE} to reverse the y axis
# item >> colorscol >> column in \code{d} giving the color  (optional)
# item >> colortransform >> required
# item >> chromToPlot >>  chromosome for which data are to be plotted
# item >> chrom >> the column in  \code{d} giving the chromosome position
# item >> pos >> the column in  \code{d} giving the base pair position   (optional)
# item >> startPos >> the column in  \code{d} giving the start base pair position
# item >> endPos >> the column in  \code{d} giving the end base pair position
# item >> rectangleOption >> boolean: if \code{TRUE} rectangles are plotted instead of points (default) - NB: for points, several 'views' can be obtained (ex using parameters \code{type} , \code{pch},...)
# item >> segmentOption >> boolean: if \code{TRUE} segments are plotted instead of points (default) - NB: for points, several 'views' can be obtained (ex using parameters \code{type} , \code{pch},...)
# item >> cyto >> cytoband object (optional)
# item >> textInfoCol >> (optional) column in \code{d} giving information to be plotted
# item >> textInfoFilterCol >> (optional) column in \code{d} (\code{TRUE},\code{FALSE}) ou \code{(0,1)} with 1/\code{TRUE} indicates values from column \code{textInfoCol} to be plotted
# item >> textInfoColTransform >> (optional) a function to be applied to values in column \code{textInfoCol} 
# item >> plotCentro >> boolean : should centromeric delimitation be plotted ?
# item >> plotCentro.color >> required
# item >> plotCyto >> boolean: should cytoband delimitation be plotted ?
# item >> plotnew >>  boolean: should the a new graph be plotted (default=\code{TRUE}), or should something be added to the current graph (-> \code{plotnew=FALSE})
# item >> plotyaxis >>  boolean: should an Y axis be plotted ? (only if \code{plotnew = TRUE})
# item >> yaxmark >>  y axis 'at' parameter (only if \code{plotyaxis = TRUE})
# item >> yaxlab >>  y axis 'labels' parameter (only if \code{plotyaxis = TRUE})
# item >> col0 >> required
# item >> ... >> required
# author >> Eric Letouze
# keyword >> graphs
# end
geco.chromPlot <- function( d=NULL,
                            ycol=NULL,                            
                            ytransform= NULL,
                            reverseYaxis =FALSE,
                            colorscol=NULL,
                            colortransform=NULL,
                            chromToPlot=c(1:22,"X","Y")[1],
                            chrom="chrom",
                            pos="meanPos",
                            startPos="start",
                            endPos="end",
                            rectangleOption=FALSE,
                            segmentOption=FALSE,
                            cyto=NULL,
                            textInfoCol=NULL,                            
                            textInfoFilterCol=NULL,   
                            textInfoColTransform=NULL,
                            plotCentro = TRUE,
                            plotCentro.color = "blue",
                            plotCyto = TRUE,
                            plotnew=TRUE,
                            plotyaxis=TRUE,
                            yaxmark=NULL,
                            yaxlab=NULL,
                            col0 = 1,
                            xlim=NULL,
                            ylim=NULL,
                            ...)  {
                           
       wc <- which(d[,chrom]==chromToPlot)
       
       if(length(wc)==0)return(NULL)
       
       d <- d[wc,]
       
       if(is.null(pos)){
            d$meanPos <- rowMeans(d[,c(startPos,endPos)])
            pos <- "meanPos"
       }   
   
       d <- d[order(d[,pos]),]
  
       
       y <- d[,ycol]                                        
       if(!is.null(ytransform)) y <- ytransform(y)
       
       if(is.null(colortransform)) colortransform <- function(z){z}

        if(is.null(col0)) col0 <- 1
       if(is.null(colorscol) ){
           colorY <- col0
       }else{
          if( ! colorscol %in% names(d) )
             colorY <- col0
          else
           colorY <- colortransform(d[,colorscol])
       }

       if(!is.null(cyto)){
              cyto <- geco.pangenomCoord ( cyto,
                                        chrom="Chromosome",
                                        startPos="Start",
                                        endPos="End",
                                        cyto=cyto,
                                        absPos="absPos"  ) 
       }
       
       if(plotnew){
        
            if(!is.null(cyto)){
               if(is.null(xlim)){ 
                     xlim <- range(unlist(cyto[which(cyto$"Chromosome"==chromToPlot),c("Start","End")]))     
       				 delta <- diff(xlim)*0.02;xlim <- xlim+c(-delta,delta)
       				 delta <- diff(xlim)*0.04/0.92;xlim <- xlim+c(delta,-delta)
               }                      
            }else{
               if(is.null(xlim)) {
                     xlim <- c(0,max(d[,endPos]))
       				 delta <- diff(xlim)*0.02;xlim <- xlim+c(-delta,delta)
       				 delta <- diff(xlim)*0.04/0.92;xlim <- xlim+c(delta,-delta)
               }
            }
      
            
            if(is.null(ylim)) ylim <- range(ifelse(reverseYaxis,-1,1)*y)
            
            
            plot( xlim,ylim,axes=FALSE,col="white",xlim=xlim,ylim=ylim,...)
            
            if(rectangleOption){
                abline(h=0)
                if(reverseYaxis){
                    rect( d[,startPos],-y,d[,endPos],rep(0,nrow(d)),col=colorY, ...)
                }else{
                    rect( d[,startPos],y,d[,endPos],rep(0,nrow(d)),col=colorY, ...)
                }       
            }else{
                if(segmentOption){
                        if(reverseYaxis){
                            rect( d[,startPos],-y,d[,endPos],-y,border=colorY, ...)
                        }else{
                            rect( d[,startPos],y,d[,endPos],y,border=colorY, ...)
                        }      
                }else{
                        if(reverseYaxis){
                           points( d[,pos],-y,col=colorY, ...)
                        }else{
                           points( d[,pos],y,col=colorY, ...)
                        }
                }
            }
            box()           
            if(plotyaxis){
                okk <- FALSE
                if(is.null(yaxmark)){
                    yaxmark <- pretty(y,n=10)
                    okk <- TRUE
                }
                if(is.null(yaxlab)) yaxlab <- yaxmark
                
                if(reverseYaxis & okk) yaxmark <- -yaxmark
                
                axis(2, #1=below, 2=left, 3=above and 4=right
                    at = yaxmark,
                    labels = yaxlab,las=2) 
            }
       }else{
            if(rectangleOption){
                if(reverseYaxis){
                    rect( d[,startPos],-y,d[,endPos],rep(0,nrow(d)),col=colorY, ...)
                }else{
                    rect( d[,startPos],y,d[,endPos],rep(0,nrow(d)),col=colorY, ...)
                }     
            }else{
                if(segmentOption){
                    if(reverseYaxis){
                        rect( d[,startPos],-y,d[,endPos],-y,border=colorY, ...)
                    }else{
                        rect( d[,startPos],y,d[,endPos],y,border=colorY, ...)
                    }  
                }else{
                    if(reverseYaxis){
                       points( d[,pos],-y,col=colorY, ...)
                    }else{
                       points( d[,pos],y,col=colorY, ...)
                    }    
                }           
            }
       }
       
       if(!is.null(cyto)){
             
             if(plotCentro){
                 w <- which(cyto$"Centro"==1 & cyto$"Chromosome"==chromToPlot)
                 wat <- mean(as.numeric(apply(cyto[w,c("Start","End")],1,mean)))      
                 abline(v=wat,lty=1,col=plotCentro.color)
             }
             
             if(plotCyto){
                  w <- which(cyto$"Chromosome"==chromToPlot)
                  wat <- as.numeric(apply(cyto[w,c("Start","End")],1,mean))
                  axis(1, #1=below, 2=left, 3=above and 4=right
                      at = wat,
                      labels = cyto[w,"Band"],las=3)
                  if(plotCentro) w <- which(cyto$"Chromosome"==chromToPlot & cyto$"Centro"!=1)    
                  abline(v=cyto[w,"End"],lty=3,col="black")
             }

            
       }
       if(!is.null(textInfoCol) )
          if(textInfoCol %in% names(d)){
            w <- 1:nrow(d)
            if(!is.null(textInfoFilterCol))
              if( textInfoFilterCol %in% names(d)) w <- which(as.logical(d[,textInfoFilterCol]))
            if(length(w)>0){
                    labs <- as.character(d[w,textInfoCol])
                    if(!is.null(textInfoColTransform)) labs <-  textInfoColTransform(labs)
            
                    if(reverseYaxis){
                       text( d[w,pos],-y[w],labs, ...)
                    }else{
                       text( d[w,pos],y[w],labs, ...)
                    }
            
            }
       }
                          
}  




# Rd
# description >> pangenomic and chromosome-wide plot
# argument
# item >> d >> data.frame which columns include genomic position information and the wanted y axis
# item >> ycol >> the column in \code{d} to be used as y axis
# item >> ytransform >> a transform to be applied to y in the plot
# item >> reverseYaxis >> boolean (default=\code{FALSE}) : use \code{TRUE} to reverse the y axis
# item >> colorscol >> column in \code{d} giving the color  (optional)
# item >> colortransform >> required
# item >> chrom >> the column in  \code{d} giving the chromosome position
# item >> pos >> the column in \code{d} giving the base pair position   (optional)
# item >> startPos >> the column in \code{d} giving the start base pair position
# item >> endPos >> the column in \code{d} giving the end base pair position
# item >> rectangleOption >> boolean : if \code{TRUE} rectangles are plotted instead of points (default) - NB : for points, several 'looks' can be obtained (ex using parameter \code{type} , \code{pch},...)
# item >> cyto >> cytoband object (optional)
# item >> absPos >> the column in  \code{d} giving the absolute base pair position (will be calculated if not present in d)
# item >> textInfoCol >> (optional) column in \code{d} giving information to be plotted
# item >> textInfoFilterCol >> (optional) column in \code{d} (\code{TRUE},\code{FALSE}) or \code{(0,1)} with 1/\code{TRUE} indicates values from column \code{textInfoCol} to be plotted
# item >> textInfoColTransform >> (optional) a function to be applied to values in column \code{textInfoCol} 
# item >> plotCentro >> boolean : should centromeric delimitation be plotted ?
# item >> plotChrom >>  boolean : should chromosome delimitation be plotted ?
# item >> plotCyto >> boolean : should cytoband delimitation be plotted ?
# item >> plotnew >>  boolean : should the a new graph be plotted (default=\code{TRUE}), or should something be added to the current graph (-> \code{plotnew=FALSE})
# item >> plotyaxis >>  boolean : should an Y axis be plotted ? (considered only if \code{plotnew = TRUE})
# item >> yaxmark >>  y axis 'at' parameter (considered only if \code{plotyaxis = TRUE})
# item >> yaxlab >>  y axis 'labels' parameter (considered only if \code{plotyaxis = TRUE})
# item >> col0 >> required
# item >> ... >> required
# item >> pages >>  a vector of values in \code{c("all",1:22,"X","Y")} - for each value a graph will be plotted, all stand for pangenomic view, the rest is for chromosome views
# value >> \code{d} in genormic order, with absolute position (base pair) added
# author >> Eric Letouze
# keyword >> graphs
# end
geco.pangenomPlotMore <- function(  d=NULL,
                            ycol=NULL, 
                            ytransform= NULL,
                            reverseYaxis =FALSE,
                            colorscol=NULL,
                            colortransform=NULL,
                            chrom="chrom",
                            pos="meanPos",
                            startPos="start",
                            endPos="end",
                            rectangleOption=FALSE,
                            cyto=geco.cytoband(),
                            absPos="absPos",
                            textInfoCol=NULL,
                            textInfoFilterCol=NULL,
                            textInfoColTransform=NULL,
                            plotCentro = TRUE,
                            plotChrom =TRUE,
                            plotCyto = TRUE,
                            plotnew=TRUE,
                            plotyaxis=TRUE,
                            yaxmark=NULL,
                            yaxlab=NULL,
                            col0 = 1,
                            pages =c("all",1:22,"X","Y"),
							CexAxis=0.8,
                            ...)  {
                            
	   if( ! ycol %in% names(d) ) stop("Column ycol ", ycol, " not found in d.\n")
          

                            
       d <- geco.pangenomCoord(d, 
                             chrom = chrom, 
                             pos = pos, 
                             startPos = startPos, 
                             endPos = endPos, 
                             cyto = cyto, 
                             absPos = absPos) 
        

       if("all" %in% pages){
            tmp <- geco.pangenomPlot(d=d,
                            ycol=ycol, 
                            ytransform= ytransform,
                            reverseYaxis =reverseYaxis,
                            colorscol=colorscol,
                            colortransform=colortransform,
                            chrom=chrom,
                            pos=pos,
                            startPos=startPos,
                            endPos=endPos,
                            rectangleOption=rectangleOption,
                            cyto=cyto,
                            absPos=absPos,
                            plotCentro =plotCentro,
                            plotChrom =plotChrom,
                            plotnew=plotnew,
                            plotyaxis=plotyaxis,
                            yaxmark=yaxmark,
                            yaxlab=yaxlab,
                            col0 = col0,
							CexAxis=CexAxis,
                            ...)
               
       }
       
       pages <- setdiff(pages,"all")
       for(pa in pages){
          tmp <- geco.chromPlot( d=d,
                                ycol=ycol, 
                                ytransform= ytransform,
                                reverseYaxis =reverseYaxis,  
                                colorscol=colorscol,  
                                colortransform=colortransform,                            
                                chromToPlot=pa,                                
                                chrom=chrom,
                                pos=pos,
                                startPos=startPos,
                                endPos=endPos,
                                rectangleOption=rectangleOption,
                                cyto=cyto, 
                                textInfoCol=textInfoCol,
                                textInfoFilterCol=textInfoFilterCol,
                                textInfoColTransform=textInfoColTransform,                              
                                plotCentro =plotCentro,
                                plotCyto = plotCyto,                                
                                plotnew=plotnew,
                                plotyaxis=plotyaxis,
                                yaxmark=yaxmark,
                                yaxlab=yaxlab,
                                col0 = col0,
                                ...)
         mtext(paste("chr",pa),line=.5)

       }  
       d                           

} 
                           
                            
   

