#' distPearson
#'
#' @param m 
#'
#' @return
#' @importFrom stats cor as.dist
#'
distPearson <- function(m)
{
    stats::as.dist(1 - stats::cor(t(m), method = "pearson"))
}

#' CompareWilcox
#'
#' @param dataMat 
#' @param annot 
#' @param ref 
#' @param groups 
#' @param featureTab 
#' @param block 
#'
#' @return
#' @importFrom stats p.adjust wilcox.test
#' @export
#'
CompareWilcox <- function(dataMat = NULL, annot = NULL, ref = NULL, groups = NULL, 
    featureTab = NULL, block = NULL)
    {
    res = featureTab
    # res = res[order(res$ID), ] dataMat = dataMat[order(row.names(dataMat)), ]
    for (k in 1:length(groups))
    {
        if (length(ref) == 1)
        {
            refsamp <- ref[[1]]
        } else
        {
            refsamp <- ref[[k]]
        }
        gpsamp <- groups[[k]]
        
        annot. <- annot[c(refsamp, gpsamp), c("barcode", "cell_id", "batch_id")]
        annot.$Condition <- c(rep("ref", length(refsamp)), rep("gpsamp", length(gpsamp)))
        
        cells_cluster = data.frame(cell_id = as.character(colnames(dataMat)))
        cells_cluster = dplyr::left_join(cells_cluster, annot., by = "cell_id")
        
        if (!is.null(block))
        {
            log2Datamat = log2(dataMat + 1)
            testWilc <- scran::pairwiseWilcox(x = log2Datamat, clusters = cells_cluster$Condition, 
                block = cells_cluster$batch_id, direction = "any")
            pval.gpsamp <- testWilc$statistics[[1]]$p.value
        } else
        {
            testWilc <- apply(dataMat, 1, function(x) stats::wilcox.test(as.numeric(x[as.character(refsamp)]), 
                as.numeric(x[as.character(gpsamp)])))
            pval.gpsamp <- unlist(lapply(testWilc, function(x) x$p.value))
        }
        
        qval.gpsamp <- stats::p.adjust(pval.gpsamp, method = "BH")
        Count.gpsamp <- apply(dataMat, 1, function(x) mean(x[as.character(gpsamp)]))
        cdiff.gpsamp <- apply(dataMat, 1, function(x) log(mean(x[as.character(gpsamp)])/mean(x[as.character(refsamp)]), 
            2))
        
        Rank.gpsamp <- rank(qval.gpsamp)  # This is different from the rank used in the Wilcox.test !! 
        
        res <- data.frame(res, Rank.gpsamp, Count.gpsamp, cdiff.gpsamp, pval.gpsamp, 
            qval.gpsamp)
        
        colnames(res) <- sub("ref", names(ref)[min(c(k, length(ref)))], sub("gpsamp", 
            names(groups)[k], colnames(res)))
    }
    res
}




#' Creates a summary table with the number of genes under- or overexpressed in each group and outputs several graphical representations
#'
#' @param ref List containing one or more vectors of reference samples. Name of the vectors will be used in the results table. The length of this list should be 1 or the same length as the groups list
#' @param dataMat reads matrix
#' @param annot selected annotation of interest
#' @param featureTab Feature annotations to be added to the results table
#' @param norm_method 
#' @param groups List containing the IDs of groups to be compared with the reference samples. Names of the vectors will be used in the results table
#'
#' @return Results table
#'
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT
#' @author Eric Letouze & Celine Vallot
#' @examples
CompareedgeRGLM <- function(dataMat=NULL,
                                 annot=NULL,
                                 ref=NULL,
                                 groups=NULL,
                                 featureTab=NULL,
                                 norm_method="RLE"
)
{
    res <- featureTab
    for(k in 1:length(groups))
    {
        print(paste("Comparing",names(ref)[min(c(k,length(ref)))],"versus",names(groups)[k]))
        if(length(ref)==1){refsamp <- ref[[1]]}else{refsamp <- ref[[k]]}
        gpsamp <- groups[[k]]
        annot. <- annot[c(refsamp,gpsamp),1:2];annot.$Condition <- c(rep("ref",length(refsamp)),rep("gpsamp",length(gpsamp)))
        mat.<-dataMat [,c(as.character(refsamp),as.character(gpsamp))]
        edgeRgroup <- c(rep(1,length(refsamp)),rep(2,length(gpsamp)))	
        y <- edgeR::DGEList(counts=mat.,group=edgeRgroup)
        y <- edgeR::calcNormFactors(y,method=norm_method)
        design <- stats::model.matrix(~edgeRgroup)
        y <- edgeR::estimateDisp(y, design)		
        fit <- edgeR::glmFit(y, design)
        comp <- edgeR::glmLRT(fit, coef=2)
        pvals <- comp$table
        colnames(pvals) <- c("log2FC.gpsamp","logCPM.gpsamp","LR.gpsamp","pval.gpsamp")
        pvals$qval.gpsamp <- stats::p.adjust(pvals$pval.gpsamp, method = "BH")
        res <- data.frame(res,pvals[,c("log2FC.gpsamp","logCPM.gpsamp","pval.gpsamp","qval.gpsamp")])
        colnames(res) <- sub("ref",names(ref)[min(c(k,length(ref)))],sub("gpsamp",names(groups)[k],colnames(res)))		
        
    }
    
    res
}	



#' changeRange
#'
#' @param v 
#' @param newmin 
#' @param newmax 
#'
#' @return
#'
changeRange <- function(v, newmin = 1, newmax = 10)
{
    oldmin <- min(v, na.rm = TRUE)
    oldmax <- max(v, na.rm = TRUE)
    newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
}

#' H1proportion
#'
#' @param pv 
#' @param lambda 
#'
#' @return
#' @export
#'
#' @examples
H1proportion <- function(pv = NA, lambda = 0.5)
{
    pi1 = 1 - mean(pv > lambda, na.rm = TRUE)/(1 - lambda)
    if (pi1 < 0)
    {
        warning(paste("estimated pi1 =", round(pi1, digits = 4), "set to 0"))
        pi1 = 0
    }
    if (pi1 > 1)
    {
        warning(paste("estimated pi1 =", round(pi1, digits = 4), "set to 1"))
        pi1 = 1
    }
    return(pi1)
}

#' enrichmentTest
#'
#' @param gene.sets 
#' @param mylist 
#' @param possibleIds 
#' @param sep 
#' @param silent 
#'
#' @importFrom stats phyper
#' @return
enrichmentTest <- function(gene.sets, mylist, possibleIds, sep = ";", silent = F)
{
    possibleIds <- unique(possibleIds)
    mylist <- unique(mylist)
    gene.sets <- lapply(gene.sets, unique)
    nids <- length(possibleIds)
    gene.sets <- lapply(gene.sets, function(x) intersect(x, possibleIds))
    nref <- sapply(gene.sets, length)
    if (all(nref == 0)) 
        stop("Error: no intersection between gene sets and possible IDs.")
    if (any(nref == 0)) 
        print("Warning: some of the gene sets have no intersection with possibleIds")
    if (!all(mylist %in% possibleIds)) 
        stop("Error: some genes in mylist are not in possibleIds")
    if (!silent) 
        cat(paste("NB : enrichment tests are based on", nids, "distinct ids.\n"))
    gene.sets <- gene.sets[nref > 0]
    n <- length(mylist)
    fun <- function(x)
    {
        y <- intersect(x, mylist)
        nx <- length(x)
        ny <- length(y)
        pval <- stats::phyper(ny - 1, nx, nids - nx, n, lower.tail = F)
        c(nx, ny, pval, paste(y, collapse = sep))
    }
    tmp <- as.data.frame(t(sapply(gene.sets, fun)))
    rownames(tmp) <- names(gene.sets)
    for (i in 1:3) tmp[, i] <- as.numeric(as.character(tmp[, i]))
    tmp <- data.frame(tmp[, 1:3], p.adjust(tmp[, 3], method = "BH"), tmp[, 4])
    names(tmp) <- c("Nb_of_genes", "Nb_of_deregulated_genes", "p-value", "q-value", 
        "Deregulated_genes")
    tmp
}

#' hclustAnnotHeatmapPlot
#'
#' @param x 
#' @param hc 
#' @param hmColors 
#' @param anocol 
#' @param xpos 
#' @param ypos 
#' @param dendro.cex 
#' @param xlab.cex 
#' @param hmRowNames 
#' @param hmRowNames.cex 
#'
#' @return
#'
hclustAnnotHeatmapPlot <- function(x = NULL, hc = NULL, hmColors = NULL, anocol = NULL, 
    xpos = c(0.1, 0.9, 0.114, 0.885), ypos = c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95), dendro.cex = 1, 
    xlab.cex = 0.8, hmRowNames = FALSE, hmRowNames.cex = 0.5)
    {
    # layout(matrix(1:3,3),heights=lhei)
    par(fig = c(xpos[1], xpos[2], ypos[5], ypos[6]), new = FALSE, mar = c(0, 0, 1.5, 
        0))  #c(0.1,0.9,0.3,1)
    plot(hc, main = "Hierarchical clustering", xlab = "", sub = "", las = 2, cex = dendro.cex, 
        cex.axis = dendro.cex)
    par(fig = c(xpos[3], xpos[4], ypos[3], ypos[4]), new = TRUE, mar = rep(0, 4))
    imageCol(anocol, xlab.cex = xlab.cex, ylab.cex = 0)  # [hc$order,]
    par(fig = c(xpos[3], xpos[4], ypos[1], ypos[2]), new = TRUE, mar = rep(0, 4))
    image(t(x), axes = FALSE, xlab = "", ylab = "", col = hmColors)
    graphics::box()
    if (hmRowNames)
    {
        axis(4, at = seq(0, 1, length.out = nrow(x)), labels = rownames(x), las = 1, 
            cex.axis = hmRowNames.cex)
    }
}

#' imageCol
#'
#' @param matcol 
#' @param strat 
#' @param xlab.cex 
#' @param ylab.cex 
#' @param drawLines 
#' @param ... 
#'
#' @return
#'
imageCol <- function(matcol = NULL, strat = NULL, xlab.cex = 0.5, ylab.cex = 0.5, 
    drawLines = c("none", "h", "v", "b")[1], ...)
    {
    if (is.null(ncol(matcol)))
    {
        matcol <- data.frame(matcol)
        colnames(matcol) = colnames(anocol)
    }
    matcol <- matcol[, ncol(matcol):1]
    if (is.null(ncol(matcol)))
    {
        matcol <- data.frame(matcol)
        colnames(matcol) = colnames(anocol)
    }
    csc <- matcol
    csc.colors <- matrix()
    csc.names <- names(table(csc))
    csc.i <- 1
    for (csc.name in csc.names)
    {
        csc.colors[csc.i] <- csc.name
        csc[csc == csc.name] <- csc.i
        csc.i <- csc.i + 1
    }
    
    if (dim(csc)[2] == 1)
    {
        csc <- matrix(as.numeric(unlist(csc)), nrow = dim(csc)[1])
    } else
    {
        csc <- matrix(as.numeric(csc), nrow = dim(csc)[1])
    }
    
    image(csc, col = as.vector(csc.colors), axes = FALSE, ...)
    if (xlab.cex != 0)
    {
        axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), colnames(matcol), las = 2, 
            tick = FALSE, cex.axis = xlab.cex, ...)
    }
    if (ylab.cex != 0)
    {
        axis(3, 0:(dim(csc)[1] - 1)/(dim(csc)[1] - 1), rownames(matcol), las = 2, 
            tick = FALSE, cex.axis = ylab.cex, ...)
    }
    if (drawLines %in% c("h", "b")) 
        abline(h = -0.5:(dim(csc)[2] - 1)/(dim(csc)[2] - 1))
    graphics::box()
    if (drawLines %in% c("v", "b")) 
        abline(v = 0.5:(dim(csc)[1] - 1)/(dim(csc)[1] - 1))
    graphics::box()
    if (!is.null(strat))
    {
        z <- factor(matcol[, strat])
        levels(z) <- 1:length(levels(z))
        z <- vectorToSegments(as.numeric(z))
        abline(v = changeRange(c(0.5, z$Ind_K + 0.5)/max(z$Ind_K), newmin = par()$usr[1], 
            newmax = par()$usr[2]), lwd = 2, lty = 2)
    }
}


#' annotToCol2
#'
#' @param annotS 
#' @param annotT 
#' @param missing 
#' @param anotype 
#' @param maxnumcateg 
#' @param categCol 
#' @param quantitCol 
#' @param plotLegend 
#' @param plotLegendFile 
#'
#' @return
#'
annotToCol2 <- function(annotS = NULL, annotT = NULL, missing = c("", NA), anotype = NULL, 
    maxnumcateg = 2, categCol = NULL, quantitCol = NULL, plotLegend = T, plotLegendFile = NULL)
    {
    if (is.null(ncol(annotS)))
    {
        annotS <- data.frame(annotS)
        colnames(annotS) = annotCol
        rownames(annotS) = rownames(annotT)
    }
    for (j in 1:ncol(annotS)) annotS[which(annotS[, j] %in% missing), j] <- NA
    if (is.null(anotype))
    {
        anotype <- rep("categ", ncol(annotS))
        names(anotype) <- colnames(annotS)
        classes <- sapply(1:ncol(annotS), function(j) class(annotS[, j]))
        nmodal <- sapply(1:ncol(annotS), function(j) length(unique(setdiff(annotS[, 
            j], NA))))
        anotype[which(classes %in% c("integer", "numeric") & nmodal > maxnumcateg)] <- "quantit"
        anotype[which(nmodal == 2)] <- "binary"
    }
    anocol <- annotS
    if (plotLegend) 
        pdf(plotLegendFile)
    if (is.null(categCol)) 
        categCol <- c("#FFEB3B", "#FFC107", "#9CCC65", "#FF9800", "#FF5722", "#009688", 
                      "gray88", "gray74", "#474448", "orangered4", "orchid", "palegreen2", 
                      "orchid4", "red4", "peru", "orangered", "palevioletred4", "purple", "sienna4", 
                      "turquoise1")
        # categCol <- c("#4285F4", "#DB4437", "#F4B400", "#0F9D58", "slategray", "black", 
        #     "orange", "turquoise4", "yellow3", "orangered4", "orchid", "palegreen2", 
        #     "orchid4", "red4", "peru", "orangered", "palevioletred4", "purple", "sienna4", 
        #     "turquoise1")
    # categCol <- c('royalblue', 'palevioletred1', 'red', 'palegreen4', 'skyblue',
    # 'sienna2', 'slateblue3', 'pink2', 'slategray', 'black', 'orange', 'turquoise4',
    # 'yellow3', 'orangered4', 'orchid', 'palegreen2', 'orchid4', 'red4', 'peru',
    # 'orangered', 'palevioletred4', 'purple', 'sienna4', 'turquoise1')
    k <- 1
    for (j in which(anotype == "categ"))
    {
        tmp <- as.factor(anocol[, j])
        classes <- as.character(levels(tmp))
        ncat <- length(levels(tmp))
        if (k + ncat > length(categCol)) 
            categCol <- c(categCol, categCol)
        levels(tmp) <- categCol[k:(k + ncat - 1)]
        fill <- as.character(levels(tmp))
        anocol[, j] <- as.character(tmp)
        k <- k + ncat
        if (plotLegend)
        {
            par(mar = c(0, 0, 0, 0))
            plot(-10, axes = F, xlim = c(0, 5), ylim = c(0, 5), xlab = "", ylab = "")
            legend(1, 5, legend = classes, fill = fill, title = colnames(anocol)[j], 
                xjust = 0.5, yjust = 1)
        }
    }
    memcol <- c()
    for (j in which(anotype == "binary"))
    {
        new <- setdiff(anocol[, j], c(NA, memcol))
        if (length(new) == 2)
        {
            memcol <- c(memcol, c("aquamarine", "plum1"))
            names(memcol)[(length(memcol) - 1):length(memcol)] <- sort(new)
        }
        if (length(new) == 1)
        {
            memcol <- c(memcol, setdiff(c("dodgerblue4", "firebrick"), memcol[setdiff(anocol[, 
                j], c(NA, new))]))
            names(memcol)[length(memcol)] <- new
        }
        anocol[, j] <- as.character(anocol[, j])
        for (z in 1:length(memcol))
        {
            anocol[which(anocol[, j] == names(memcol)[z]), j] <- memcol[z]
        }
        if (plotLegend)
        {
            par(mar = c(0, 0, 0, 0))
            plot(-10, axes = F, xlim = c(0, 5), ylim = c(0, 5), xlab = "", ylab = "")
            classes <- intersect(names(memcol), annotS[, j])
            fill <- memcol[classes]
            legend(1, 5, legend = classes, fill = fill, title = colnames(anocol)[j], 
                xjust = 0.5, yjust = 1)
        }
    }
    if (is.null(quantitCol)) 
        quantitCol <- c("darkgreen", "darkblue", "darkgoldenrod4", "darkorchid4", 
            "darkolivegreen4", "darkorange4", "darkslategray")
    k <- 1
    for (j in which(anotype == "quantit"))
    {
        colrange <- matlab.like(100)
        anocol[, j] <- colrange[round(changeRange(anocol[, j], newmin = 1, newmax = 100))]
        if (k < length(quantitCol))
        {
            k <- k + 1
        } else
        {
            k <- 1
        }
        if (plotLegend)
        {
            par(mar = c(8, 2, 5, 1))
            lims <- seq(-1, 1, length.out = 200)
            image(matrix(lims, ncol = 1), col = colrange, axes = F, xlab = colnames(anocol)[j])
        }
    }
    if (plotLegend) 
        dev.off()
    for (j in 1:ncol(anocol)) anocol[which(is.na(anocol[, j])), j] <- "white"
    as.matrix(anocol)
}

#' Title
#'
#' @param mat 
#' @param margin 
#' @param groups 
#' @param method 
#'
#' @return
#'
#' @examples
groupMat <- function(mat = NA, margin = 1, groups = NA, method = "mean")
{
    if (!method %in% c("mean", "median"))
    {
        print("Method must be mean or median")
        return()
    }
    if (!margin %in% 1:2)
    {
        print("Margin must be 1 or 2")
        return()
    }
    for (i in 1:length(groups))
    {
        if (margin == 1)
        {
            if (length(groups[[i]]) == 1)
            {
                v <- mat[, groups[[i]]]
            } else
            {
                if (method == "mean")
                {
                  v <- apply(mat[, groups[[i]]], margin, mean)
                } else
                {
                  v <- apply(mat[, groups[[i]]], margin, median)
                }
            }
            if (i == 1)
            {
                res <- matrix(v, ncol = 1)
            } else
            {
                res <- cbind(res, v)
            }
        } else
        {
            if (length(groups[[i]]) == 1)
            {
                v <- mat[groups[[i]], ]
            } else
            {
                if (method == "mean")
                {
                  v <- apply(mat[groups[[i]], ], margin, mean)
                } else
                {
                  v <- apply(mat[groups[[i]], ], margin, median)
                }
            }
            if (i == 1)
            {
                res <- matrix(v, nrow = 1)
            } else
            {
                res <- rbind(res, v)
            }
        }
    }
    if (margin == 1)
    {
        rownames(res) <- rownames(mat)
        colnames(res) <- names(groups)
    } else
    {
        rownames(res) <- names(groups)
        colnames(res) <- colnames(mat)
    }
    res
}


#' Sparse covariance calculation
#' Directly copied from Ron, Joris Meys & mike
#' SO https://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r 
#'
#' @param x 
#'
#' @return A covariance matrix
#'
#' @examples
sparse_cov <- function(x){
    n <- nrow(x)
    cMeans <- Matrix::colMeans(x)
    covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
    return(covmat)
}
