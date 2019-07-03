# Rd
# description >> Function that returns the R object stored in filename
# argument
# item >> filename >> File to load
# value >> R object stored in file
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.load <- function (filename) # File to be loaded
{
    if (file.exists(filename)) 
        return(eval(parse(text = load(filename))))
    cat(paste("error - function geco.load : file ", filename, 
        " not found"))
    NULL
}

# Rd
# description >> Function that estimates the number of test under H1 hypothesis using Storey method
# argument
# item >> pv >> Vector of pvalues
# item >> lambda >> Lambda parameter for Storey method
# value >> Proportion of tests under H1 hypothesis
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.H1proportion <- function(pv=NA,
				              lambda = 0.5
				             ) 
{
    pi1 = 1 - mean(pv > lambda, na.rm = TRUE)/(1 - lambda)
    if (pi1 < 0) {
        warning(paste("estimated pi1 =",round(pi1, digit = 4),"set to 0"))
        pi1 = 0
    }
    if (pi1 > 1) {
        warning(paste("estimated pi1 =",round(pi1, digit = 4),"set to 1"))
        pi1 = 1
    }
    return(pi1)
}

# Rd
# description >> Function to convert factors factors in a data frame to characters
# argument
# item >> d >> Data frame to be converted
# value >> Data frame with factors converted to characters
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
factochar <- function(d) {
  for(i in 1:ncol(d)) if(is.factor(d[,i])) d[,i] <- as.character(d[,i])
  d
}

# Rd
# description >> Function to convert factors in a data frame to the appropriate format (character, numeric...)
# argument
# item >> d >> Data frame to be converted
# item >> ncmax >> Maximum number of characters for numeric fields
# value >> required
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
factoall <- function (d, ncmax = 10) 
{
    n <- ncol(d)
    for (i in 1:n) {
        if (is.factor(d[, i])) {
            d[, i] <- as.character(d[, i])
            na <- which(is.na(d[, i]))
            num <- suppressWarnings(as.numeric(d[, i]))
            nanum <- which(is.na(num))
            if (length(nanum) == length(na)) {
#                int <- suppressWarnings(as.integer(d[, i]))
#                naint <- which(is.na(int))
#                nc <- nchar(num)
#                if (length(naint) == length(nanum) & all(nc < ncmax)) {
#                  d[, i] <- int
#                }
#                else {
                  d[, i] <- num
#                }
            }
        }
    }
    d
}

# Rd
# description >> Function to convert a numeric vector in a series of segments
# argument
# item >> v >> Numeric vector to convert into a matrix of segments
# item >> replaceNA >> Numeric value that temporarily replaces NAs in the function. Should be different from any value in v.
# item >> d >> Data frame with other characteristics to add to the result. Each column should have a unique value for each continuous segment in v.
# value >> Data frame giving the coordinates of segment boundaries and values within each segment from v, plus additional characteristics if d was supplied.
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.vectorToSegments <- function(v,
                                 replaceNA=999,
                                 d=NULL)
{
       if(!is.null(replaceNA)){
          wna <- which(is.na(v))
          if(length(wna)>0) v[wna] <- replaceNA
       }
       w <- which(diff(v)!=0)
       if(length(w)==0){res <- data.frame(Ind=1,Ind_K=length(v),Len=length(v),value=v[1])}else{
            res <- data.frame(Ind=c(1,w+1),Ind_K=c(w,length(v)),Len=c(w,length(v))-c(1,w+1)+1,value=v[c(w,length(v))])
       }
       if(!is.null(replaceNA)) {
            wna <- which(res$value==replaceNA)
            if(length(wna)>0) res[wna,"value"] <- NA
       }
       if(!is.null(d)) res <- cbind(res[,1:3],d[res[,1],])
       res
}

# Rd
# description >> Discretize a continous variable by specified \code{lim} cut-off(s)
# argument
# item >> x >>  a vector
# item >> lim >>  cut-off(s)
# item >> quant >>  if TRUE (default FALSE) the \code{x} is discretize by quantile and \code{lim} is considered as cut-off(s) for quantile, ie 0<lim<1
# item >> addlevels >>  add character levels indicating the cut-offs (i.e. for un cut-off iqq levels=c("<iqq",">=iqq") )
# value >> a vector of integers
# author >> Eric Letouze
# keyword >> utilities
# end
geco.discretize <- function(x, lim, quant = FALSE, addlevels = FALSE){

  lim <- sort(lim)

  if( quant & ( any(lim>1) | any(lim<0) ) )
    stop("lim must be [0;1] as quant=T\n")

  res <- rep(NA, length(x))

  if(quant) lim  <- quantile(x, probs=lim, na.rm=TRUE)
  n <- length(lim)
  for(i in n:1) res[which(x<lim[i])] <- i
  res[which(x>=lim[n])] <- n+1
  
  if (addlevels) {
    res <- as.factor(res)

    if(quant)
      lim <- gsub(" ", "", prettyNum(lim, format="g", digits=1))

    if(length(lim)==1)
      lev <- c(paste("<", lim[1], sep = ""), paste(">=", lim[1], sep = ""))
    else
      lev <- c(paste("<", lim[1], sep = ""), paste("[",lim[-length(lim)],";",lim[-1],"[",sep=""), paste(">=", lim[length(lim)], sep = ""))
    levels(res) <- lev[as.numeric(levels(res))]
  }
  res
}

# Rd
# description >> wrapping of the function 'density' adding down and top points to the result 
# argument
# item >> x >> a numeric vector
# item >> doplot >> boolean
# item >> pc >> ...
# value >> a list with 2 elements : down and top points
# author >> Eric Letouze
# keyword >> methods
# end
geco.density <- function(x,doplot=FALSE,pc=.05,...){
      dx <- density(x,na.rm=TRUE,...)
      ymax <- diff(range(dx$y))
      n <- length(dx$y)
      croissance <- as.numeric((dx$y[-1] - dx$y[-n])>0)
      wB <- 1+which(diff(croissance)== 1)
      wH <- 1+which(diff(croissance)== -1)
      pointsBas   <- dx$x[wB]
      pointsHauts <- dx$x[intersect(wH,which(dx$y>pc*ymax))]
      
      if(length(pointsHauts)>0)pointsBas   <- sapply(split(pointsBas,geco.discretize(pointsBas,pointsHauts)),median)

      if(doplot){
         plot(dx,...)
         abline(v=pointsBas,lty=3,col="blue")
         abline(v=pointsHauts,lty=3,col="red")
      }
      L <- c(dx,list("down"=pointsBas,"top"=pointsHauts))
      attr(L,"class") <- "density"
      L
}

# Rd
# description >> wrapping of the function 'geco.density' to control the maximum number of peaks
# argument
# item >> x >> a numeric vector
# item >> doplot >> boolean
# item >> percentHighestPeak >> ...
# item >> maxNbPeaks >> ...
# item >> minDeltaBetweenPeaks >> ...
# item >> deltaApproach >> ...
# item >> ... >> additionnal paramters to be passed to plot
# value >> a list with 3 elements : top (peaks) and down (inter-peaks) points  , number of peaks
# author >> Eric Letouze
# keyword >> methods
# end  
geco.peaks <- function(x,
                      percentHighestPeak=.2,
                      maxNbPeaks=NULL,
                      minDeltaBetweenPeaks=.03,
                      deltaApproach=1,
                      doplot=FALSE,
                      ...){

        if(is.na( minDeltaBetweenPeaks ))minDeltaBetweenPeaks <- NULL
        v <- geco.density(x,pc =percentHighestPeak)
        m <- v$y[which(v$x %in% v$top)]
        w <- which(m/max(m) > percentHighestPeak)

        xw <- v$top[w]

        which.eq <- function(z){ order(abs(v$x-z))[1]}

       # umin[which(apply(t(sapply(umin,function(z)c( max(which(umax <=z )),min(which(umax>=z))))),1,function(z)!any(is.infinite(z))))]

      if(deltaApproach==1 & !is.null( minDeltaBetweenPeaks ) & length(xw)>1){
             for(i in 1:(length(xw)-1)){
                if(xw[i+1]-xw[i] <  minDeltaBetweenPeaks) {
                    xw[i:(i+1)] <- xw[c(i:(i+1))][which.max(m[i:(i+1)])]
                }
             }

             xw <- unique(xw)

        }


        if(!is.null(maxNbPeaks)){
               if(length(xw)> maxNbPeaks)
                   for(i in 1:(length(xw)-maxNbPeaks)){
                          oneamong <- xw[c(0,1)+which.min(diff(xw))]
                          w1 <-  which.eq(oneamong[1])
                          w2 <-  which.eq(oneamong[2])
                          out <- oneamong[which.min(v$y[c(w1,w2)])]
                          xw <- setdiff(xw,out)

                   }
        }


        umin <- NULL
        if(length(xw)>1){
            for(i in 1:(length(xw)-1)){
                 possiblemin <- v$down[which(v$down >= xw[i] & v$down <= xw[i+1])]
                 if(length(possiblemin)>1){
                      w<- sapply(possiblemin,which.eq)
                      possiblemin <- possiblemin[which.min(v$y[w])]
                 }
                 umin <- c(umin,possiblemin)
            }

        }


      #  if(doplot){
#            plot(v,...)
#            abline(v=xw,col="red")
#            if(length(umin)>0)abline(v=umin,col="green",lty=3)
#        }

        L <- list("x abciss big peaks"=xw,
             "x abciss inter-peaks"=umin,
             "nb big peaks"=length(xw))
        
        peaks <- L[[1]]
        interpeaks <- L[[2]]     
        temp <- as.data.frame(t(sapply(peaks,function(pic) {
                         wleft <- which(interpeaks < pic)
                         if(length(wleft)>0){
                            wleft <- max(interpeaks[wleft])
                          }else{
                            wleft <- min(x,na.rm=TRUE)
                          }
                          wright <- which(interpeaks > pic)
                          if(length(wright)>0){
                            wright <- min(interpeaks[wright])
                          }else{
                            wright <- max(x,na.rm=TRUE)
                          } 
                          c(pic,wleft,wright,length(which(x>=wleft & x<=wright)))
                          }) ))
        names(temp) <- c("peak","left born","right born","size")                           
        L$"peaks x size" <- temp   
        
        
        if(deltaApproach==2 &!is.null( minDeltaBetweenPeaks ) & L$'nb big peaks'>1){
            Lini <- L
            names(Lini) <- paste("initial",names(L))

            while(any( diff(L$"x abciss big peaks") < minDeltaBetweenPeaks) & L$"nb big peaks" >1){

                  w <- which.min(abs(diff(L$"x abciss big peaks")) )


                  lb <- min(L$"peaks x size"[w:(w+1),"left born"])
                  rb <- max(L$"peaks x size"[w:(w+1),"right born"])
                  si <- sum(L$"peaks x size"[w:(w+1),"size"])
                  pe <- median(x[which(x>=lb & x<=rb)])
                  
                  L$"peaks x size"[w,] <- c(pe,lb,rb,si)
                  L$"peaks x size" <- L$"peaks x size"[-(w+1),]
                  
                  L$"nb big peaks" <- L$"nb big peaks" - 1
                  L$"x abciss big peaks"[w] <- pe
                  L$"x abciss big peaks" <- L$"x abciss big peaks"[-(w+1)]
                  L$"x abciss inter-peaks" <- L$"x abciss inter-peaks"[-(w+1)]
            }

            L <- c(Lini,L)
        }
        
        if(doplot){
            plot(v,...)
            abline(v=L$"x abciss big peaks",col="red")
            if(length(L$"x abciss inter-peaks")>0)abline(v=L$"x abciss inter-peaks",col="green",lty=3)
        }
        
        L     
}                

# Rd
# description >> Function that linearly changes the range of a numeric vector
# argument
# item >> v >> Vector to be transformed
# item >> newmin >> New min value
# item >> newmax >> New max value
# value >> Vector linearly transformed with the new range
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.changeRange <- function (v, newmin = 1, newmax = 10) 
{
    oldmin <- min(v, na.rm = TRUE)
    oldmax <- max(v, na.rm = TRUE)
    newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
}












