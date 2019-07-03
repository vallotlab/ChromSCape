my.heatmap.plus <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
    distfun = dist, hclustfun = hclust, reorderfun = function(d, 
        w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, 
        "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
    margins = c(5,5,5,5), ColSideColors, RowSideColors, cexRow = 0.2 + 
        1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
    verbose = getOption("verbose"), colAnnot.height=0.2, rowAnnot.width=0.2, cex.annot=1, ...) 
{
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("'x' must have at least 2 rows and 2 columns")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) 
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm) 
                x
            else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv) 
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
        if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
        if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
        4)
    if (!missing(ColSideColors)) {
        if (!is.matrix(ColSideColors)) 
            stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || dim(ColSideColors)[1] != 
            nc) 
            stop("'ColSideColors' dim()[2] must be of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1], colAnnot.height, lhei[2])
    }
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)) 
            stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || dim(RowSideColors)[1] != 
            nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
            1), lmat[, 2] + 1)
        lwid <- c(lwid[1], rowAnnot.width, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, 
            "; lmat=\n")
        print(lmat)
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], margins[2], 0, 0.5))
        rsc = RowSideColors[rowInd, ]
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
            rsc.colors[rsc.i] = rsc.name
            rsc[rsc == rsc.name] = rsc.i
            rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(colnames(RowSideColors)) > 0) {
            axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), 
                las = 2, tick = FALSE, cex.axis=cex.annot)
        }
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, margins[2], 0.5, margins[4]))
        csc = ColSideColors[colInd, ]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
            csc.colors[csc.i] = csc.name
            csc[csc == csc.name] = csc.i
            csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
            axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), colnames(ColSideColors), 
                las = 2, tick = FALSE, cex.axis=cex.annot)
        }
    }
    par(mar = c(margins[1],margins[2],0,margins[4]))
    if (!symm || scale != "none") {
        x <- t(x)
    }
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 0))
    if (doRdend) 
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
    if (doCdend) 
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main)) 
        frame()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
        doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}


geco.HeatmapPlot <- function(x=NULL,
							 hmColors=NULL,
							 hmRowNames=TRUE,
							 hmRowNames.cex=1,
							 hmColNames=TRUE,
							 hmColNames.cex=1
							 )
{
	image(t(x),axes=FALSE,xlab="",ylab="",col=hmColors)
	box()
	if(hmRowNames){
		axis(3,at=seq(0,1,length.out=ncol(x)),labels=colnames(x),las=2,cex.axis=hmRowNames.cex)
	}
	if(hmColNames){
		axis(4,at=seq(0,1,length.out=nrow(x)),labels=rownames(x),las=1,cex.axis=hmNames.cex)
	}
}							


geco.AnnotHeatmapPlot <- function(x=NULL,
							       hmColors=NULL,
							       anocol=NULL,
							       xpos=c(0.2,0.8),
							       ypos=c(0.1,0.75,0.75,0.95),
							 	   xlab.cex=0.8,
							       hmRowNames=TRUE,
	   							   hmRowNames.cex=0.5
							 )
{
	#layout(matrix(1:3,3),heights=lhei)
	par(fig=c(xpos[1],xpos[2],ypos[3],ypos[4]), new=FALSE, mar=rep(0,4))
	geco.imageCol(anocol,xlab.cex=xlab.cex,ylab.cex=0)
	par(fig=c(xpos[1],xpos[2],ypos[1],ypos[2]), new=TRUE, mar=rep(0,4))
	image(t(x),axes=FALSE,xlab="",ylab="",col=hmColors)
	box()
	if(hmRowNames){
		axis(4,at=seq(0,1,length.out=nrow(x)),labels=rownames(x),las=1,cex.axis=hmRowNames.cex)
	}
}							


geco.hclustAnnotPlot <- function(hc=NULL,
							      hmColors=NULL,
							 	  anocol=NULL,
							 	  xpos=c(0.1,0.9,0.114,0.885),
							 	  ypos=c(0.1,0.4,0.4,0.95),
							 	  dendro.cex=1,
							 	  xlab.cex=0.8
							 )
{
	#layout(matrix(1:3,3),heights=lhei)
	par(fig=c(xpos[1],xpos[2],ypos[3],ypos[4]), new=FALSE, mar=c(0,0,1.5,0))#c(0.1,0.9,0.3,1)
	plot(hc,main="Hierarchical clustering", xlab="", sub="", las=2,cex=dendro.cex,cex.axis=dendro.cex)
	par(fig=c(xpos[3],xpos[4],ypos[1],ypos[2]), new=TRUE, mar=rep(0,4))
	geco.imageCol(anocol,xlab.cex=xlab.cex,ylab.cex=0)# [hc$order,]
}							


geco.hclustAnnotHeatmapPlot <- function(x=NULL,
							 hc=NULL,
							 hmColors=NULL,
							 anocol=NULL,
							 xpos=c(0.1,0.9,0.114,0.885),
							 ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
							 dendro.cex=1,
							 xlab.cex=0.8,
							 hmRowNames=FALSE,
							 hmRowNames.cex=0.5
							 )
{
	#layout(matrix(1:3,3),heights=lhei)
	par(fig=c(xpos[1],xpos[2],ypos[5],ypos[6]), new=FALSE, mar=c(0,0,1.5,0))#c(0.1,0.9,0.3,1)
	plot(hc,main="Hierarchical clustering", xlab="", sub="", las=2,cex=dendro.cex,cex.axis=dendro.cex)
	par(fig=c(xpos[3],xpos[4],ypos[3],ypos[4]), new=TRUE, mar=rep(0,4))
	geco.imageCol(anocol,xlab.cex=xlab.cex,ylab.cex=0)# [hc$order,]
	par(fig=c(xpos[3],xpos[4],ypos[1],ypos[2]), new=TRUE, mar=rep(0,4))
	image(t(x),axes=FALSE,xlab="",ylab="",col=hmColors)
	box()
	if(hmRowNames){
		axis(4,at=seq(0,1,length.out=nrow(x)),labels=rownames(x),las=1,cex.axis=hmRowNames.cex)
	}
}							


