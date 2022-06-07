
#' Get experiment names from a SingleCellExperiment
#'
#' @param scExp A SingleCellExperiment with named mainExp and altExps.
#'
#' @return Character vector of unique experiment names
#' @importFrom SingleCellExperiment mainExpName altExpNames
#' @export
#'
#' @examples
#' data(scExp)
#' getExperimentNames(scExp)
#' 
getExperimentNames<- function(scExp){
    return(unique(c(SingleCellExperiment::mainExpName(scExp),
                    SingleCellExperiment::altExpNames(scExp))))
}

#' Get Main experiment of a SingleCellExperiment
#'
#' @param scExp A SingleCellExperiment with named mainExp and altExps.
#'
#' @return The swapped SingleCellExperiment towards "main" experiment
#' @export
#' @importFrom SingleCellExperiment mainExpName altExpNames
#' @examples
#' data(scExp)
#' getMainExperiment(scExp)
#' 
getMainExperiment  <- function(scExp){
    if(SingleCellExperiment::mainExpName(scExp) == "main"){
        return(scExp) 
    } else if("main" %in% SingleCellExperiment::altExpNames(scExp)){
        return(swapAltExp_sameColData(scExp, "main"))
    } else{
        return(scExp)
    }
}


#' Swap main & alternative Experiments, with fixed colData
#'
#' @param scExp A SingleCellExperiment
#' @param alt Name of the alternative experiment
#'
#' @return A swapped SingleCellExperiment with the exact same colData.
#' @export
#' @importFrom SingleCellExperiment mainExpName altExpNames swapAltExp
#' @importFrom SummarizedExperiment colData
#' @examples
#' data(scExp) 
#' swapAltExp_sameColData(scExp, "peaks")
#' 
swapAltExp_sameColData <- function(scExp, alt){
    if(alt %in% SingleCellExperiment::altExpNames(scExp)){ 
        cd = SummarizedExperiment::colData(scExp)
        SummarizedExperiment::colData(scExp) = NULL
        scExp = SingleCellExperiment::swapAltExp(
            scExp, alt, SingleCellExperiment::mainExpName(scExp),
            withColData = FALSE)
        SummarizedExperiment::colData(scExp) = cd
    } else {
        message("ChromSCape::swapAltExp_sameColData - ", alt, " not in ",
                "altExpNames(scExp), returning unchanged scExp.")
    }
    return(scExp)
}

#' distPearson
#'
#' @param m A matrix
#'
#' @return A dist object
#' @importFrom stats cor as.dist
#'
distPearson <- function(m)
{
    stats::as.dist(1 - stats::cor(t(m), method = "pearson"))
}

#' Find comparable variable scExp
#'
#' @param scExp A SingleCellExperiment
#' @param allExp A logical indicating wether alternative experiments comparable
#' variables should also be fetch.
#'
#' @return A character vector with the comparable variable names
#'
comparable_variables <- function(scExp, allExp = TRUE)
{
    comparable = c()
    exps = SingleCellExperiment::mainExpName(scExp)
    if(allExp) exps = c(exps, SingleCellExperiment::altExpNames(scExp))
    
    for(exp in exps){
        if(exp == exps[1]) scExp. = scExp else
            scExp. = SingleCellExperiment::swapAltExp(scExp, exp)
        annot = SingleCellExperiment::colData(scExp.)
        comparable = c(comparable,
                       names(which(unlist(lapply(annot, class)) %in% c("character", "factor") &
                                       unlist(lapply(lapply(annot, table), length)) > 1 &
                                       unlist(lapply(lapply(annot, table), length)) < 100 &
                                       !grepl("_color", colnames(annot))))
        )
    }
    
    return(unique(comparable))
}



#' CompareWilcox
#'
#' @param dataMat A raw count matrix
#' @param annot A cell annotation data.frame
#' @param ref_group List with cells in reference group(s)
#' @param groups List with cells in group(s) to test
#' @param featureTab data.frame with feature annotation
#' @param block Use a blocking factor to conteract batch effect ? 
#' @param BPPARAM BPPARAM object for multiprocessing. See
#'  \link[BiocParallel]{bpparam} for more informations. Will take the default
#'  BPPARAM set in your R session.
#'
#' @return A dataframe containing the foldchange and p.value of each feature
#' @importFrom stats p.adjust wilcox.test
#' @importFrom DelayedArray blockApply colAutoGrid
#' @importFrom scran pairwiseWilcox 
#' @importFrom matrixTests col_wilcoxon_twosample 
#' @export
#' @author Eric Letouze & Celine Vallot & Pacome Prompsy
#' 
#' @examples
#' data("scExp")
#' scExp_cf = correlation_and_hierarchical_clust_scExp(scExp)
#' scExp_cf = choose_cluster_scExp(scExp_cf,nclust=2,consensus=FALSE)
#' featureTab = as.data.frame(SummarizedExperiment::rowRanges(scExp_cf))
#' rownames(featureTab) = featureTab$ID
#' ref_group = list("C1"=scExp_cf$cell_id[which(scExp_cf$cell_cluster=="C1")])
#' groups = list("C2"=scExp_cf$cell_id[which(scExp_cf$cell_cluster=="C2")])
#' myres = CompareWilcox(as.matrix(SingleCellExperiment::normcounts(scExp_cf)),
#' annot=as.data.frame(SingleCellExperiment::colData(scExp_cf)),
#'    ref_group=ref_group,groups=groups, featureTab=featureTab)
#' 
CompareWilcox <- function(dataMat = NULL, annot = NULL, ref_group = NULL,
                        groups = NULL, featureTab = NULL, block = NULL,
                        BPPARAM = BiocParallel::bpparam())
{
    res = featureTab
    mat = Matrix::Matrix((dataMat > 0) + 0, sparse = TRUE)
    
    for (k in seq_along(groups)){
    if (length(ref_group) == 1) refsamp <- ref_group[[1]] else
        refsamp <- ref_group[[k]]
    
    gpsamp <- groups[[k]]
    refidx = match(refsamp,colnames(dataMat))
    gpidx = match(gpsamp,colnames(dataMat))
    
    cluster_bin_mat = mat[, gpidx]
    reference_bin_mat = mat[, refidx]
    
    group_sum = Matrix::rowSums(cluster_bin_mat)
    group_activation = group_sum / ncol(cluster_bin_mat)
    
    reference_sum = Matrix::rowSums(reference_bin_mat) 
    reference_activation = reference_sum / ncol(reference_bin_mat)
      
    if (!is.null(block))
    {
        testWilc <- scran::pairwiseWilcox(
            dataMat[,c(refidx, gpidx)],
            groups = as.factor(c(rep(1, length(refidx)),
                                 rep(2, length(gpidx)))),
            block = annot.$batch_id[c(refidx, gpidx)],
            BPPARAM = BPPARAM)
        
        pval.gpsamp <- testWilc$statistics[[1]]$p.value
    } else
    {
        dataMat = t(dataMat)
        
        system.time({
            testWilc =  DelayedArray::blockApply(
                dataMat,
                grid = DelayedArray::colAutoGrid(dataMat, ncol = 1000),
                BPPARAM = BPPARAM,
                function(X, gpidx, refidx){
                    matrixTests::col_wilcoxon_twosample(X[gpidx,], X[refidx,])
                }, gpidx, refidx)
        })
        testWilc = do.call("rbind", testWilc)
        pval.gpsamp <- testWilc$pvalue
        dataMat = t(dataMat)
    }
    
    qval.gpsamp <- stats::p.adjust(pval.gpsamp, method = "BH")
    Count.gpsamp <- rowMeans(dataMat[,gpidx])
    Count.refsamp <- rowMeans(dataMat[,refidx])
    logFC.gpsamp <- log(Count.gpsamp/Count.refsamp,2)
    
    if(length(which(logFC.gpsamp == Inf))>0)
      logFC.gpsamp[which(logFC.gpsamp == Inf)] = 
      max(logFC.gpsamp[is.finite(logFC.gpsamp)])
    if(length(which(logFC.gpsamp == -Inf))>0) 
      logFC.gpsamp[which(logFC.gpsamp == -Inf)] = 
      min(logFC.gpsamp[is.finite(logFC.gpsamp)])
    
    res <- cbind(res, data.frame(logFC.gpsamp = logFC.gpsamp,
                                 qval.gpsamp = qval.gpsamp,
                                 group_activation.gpsamp = group_activation,
                                 reference_activation.gpsamp = reference_activation)
    )
    colnames(res) = gsub("gpsamp", names(groups)[k], colnames(res))
    }
    return(res)
}


#' Creates a summary table with the number of genes under- or overexpressed in
#' each group and outputs several graphical representations
#'
#' @param ref_group List containing one or more vectors of reference samples.
#'   Name of the vectors will be used in the results table. The length of this
#'   list should be 1 or the same length as the groups list
#' @param dataMat reads matrix
#' @param annot selected annotation of interest
#' @param featureTab Feature annotations to be added to the results table
#' @param norm_method Which method to use for normalizing ('upperquantile')
#' @param groups List containing the IDs of groups to be compared with the
#'   reference samples. Names of the vectors will be used in the results table
#'
#' @return A dataframe containing the foldchange and p.value of each feature
#'
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT
#' @author Eric Letouze & Celine Vallot
#' @export
#' @examples
#' data("scExp")
#' scExp_cf = correlation_and_hierarchical_clust_scExp(scExp)
#' scExp_cf = choose_cluster_scExp(scExp_cf,nclust=2,consensus=FALSE)
#' featureTab = as.data.frame(SummarizedExperiment::rowRanges(scExp_cf))
#' rownames(featureTab) = featureTab$ID
#' ref_group = list("C1"=scExp_cf$cell_id[which(scExp_cf$cell_cluster=="C1")])
#' groups = list("C2"=scExp_cf$cell_id[which(scExp_cf$cell_cluster=="C2")])
#' myres = CompareedgeRGLM(as.matrix(SingleCellExperiment::counts(scExp_cf)),
#' annot=as.data.frame(SingleCellExperiment::colData(scExp_cf)),
#'    ref_group=ref_group,groups=groups, featureTab=featureTab)
#' 
CompareedgeRGLM <- function(
    dataMat=NULL, annot=NULL, ref_group=NULL, groups=NULL, featureTab=NULL, 
    norm_method="TMMwsp"){
    res <- featureTab
    mat = Matrix::Matrix((dataMat > 0) + 0, sparse = TRUE)
    
    for(k in seq_along(groups)){
        print(
            paste("Comparing",names(ref_group)[min(c(k,length(ref_group)))],
                        "versus",names(groups)[k]))
        if(length(ref_group)==1){
          refsamp <- ref_group[[1]]
        }else {refsamp <- ref_group[[k]]}
        gpsamp <- groups[[k]]
        annot. <- annot[c(refsamp,gpsamp),seq_len(2)]
        annot.$Condition <- c(
            rep("ref",length(refsamp)),rep("gpsamp",length(gpsamp)))
        mat. <- dataMat[,c(as.character(refsamp), as.character(gpsamp))]
        
        cluster_bin_mat = mat[, as.character(gpsamp)]
        reference_bin_mat = mat[, as.character(refsamp)]
        
        group_sum = Matrix::rowSums(cluster_bin_mat)
        group_activation = group_sum / ncol(cluster_bin_mat)
        
        reference_sum = Matrix::rowSums(reference_bin_mat) 
        reference_activation = reference_sum / ncol(reference_bin_mat)
        
        edgeRgroup <- c(rep(1,length(refsamp)),rep(2,length(gpsamp)))
        y <- edgeR::DGEList(counts=mat.,group=edgeRgroup)
        y <- edgeR::calcNormFactors(y,method=norm_method)
        design <- stats::model.matrix(~edgeRgroup)
        y <- edgeR::estimateDisp(y, design)
        fit <- edgeR::glmFit(y, design)
        comp <- edgeR::glmLRT(fit, coef=2)
        pvals <- comp$table
        colnames(pvals) <- c(
            "logFC.gpsamp","logCPM.gpsamp","LR.gpsamp","pval.gpsamp")
        pvals$qval.gpsamp <- stats::p.adjust(pvals$pval.gpsamp, method = "BH")
        res <- data.frame(res,pvals[,c(
            "logFC.gpsamp","qval.gpsamp")])
        colnames(res) <-sub("gpsamp", names(groups)[k],colnames(res))
        res[,paste0("group_activation.", names(groups)[k])] = group_activation
        res[,paste0("reference_activation.", names(groups)[k])] = reference_activation
    }
    
    res
}

#' changeRange
#'
#' @param v A numeric vector
#' @param newmin New min
#' @param newmax New max
#'
#' @return A matrix with values scaled between newmin and newmax
changeRange <- function(v, newmin = 1, newmax = 10)
{
    oldmin <- min(v, na.rm = TRUE)
    oldmax <- max(v, na.rm = TRUE)
    newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
}

#' H1proportion
#'
#' @param pv P.value vector
#' @param lambda Lambda value
#'
#' @return H1 proportion value
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
#' @param gene.sets A list of reference gene sets
#' @param mylist A list of genes to test
#' @param possibleIds All existing genes
#' @param sep Separator used to collapse genes
#' @param silent Silent mode ?
#'
#' @importFrom stats phyper
#' @return A dataframe with the gene sets and their enrichment p.value
#' 
enrichmentTest <- function(
    gene.sets, mylist, possibleIds, sep = ";", silent = FALSE)
{
    possibleIds <- unique(possibleIds)
    mylist = intersect(mylist, possibleIds)
    nids <- length(possibleIds)
    possibleIds
    nref <- as.numeric(lapply(gene.sets, length))
    if (all(nref == 0)) 
        stop("Error: no intersection between gene sets and possible IDs.")

    if (!silent) 
        cat("NB : enrichment tests are based on", nids, "distinct ids.\n")
    gene.sets <- gene.sets[nref > 0]
    n <- length(mylist)
    fun <- function(x)
    {
        y <- intersect(x, mylist)
        nx <- length(x)
        ny <- length(y)
        pval <- stats::phyper(ny - 1, nx, nids - nx, n, lower.tail = FALSE)
        c(nx, ny, pval, paste(y, collapse = sep))
    }
    tmp <- as.data.frame(t(as.matrix(
        vapply(gene.sets, fun, FUN.VALUE = c(
            "Nb_of_genes" = 0, "Nb_of_deregulated_genes" = 0,
            "p-value" =0, "Deregulated_genes" = ""
        ))
    )))
    rownames(tmp) <- names(gene.sets)
    for (i in seq_len(3)) tmp[, i] <- as.numeric(
        as.character(tmp[, i]))
    tmp <- data.frame(
        tmp[, seq_len(3)], p.adjust(tmp[, 3], method = "BH"), tmp[, 4])
    names(tmp) <- c("Nb_of_genes", "Nb_of_deregulated_genes",
                                    "p-value", "q-value", "Deregulated_genes")
    tmp
}

#' hclustAnnotHeatmapPlot
#'
#' @param x A correlation matrix
#' @param hc An hclust object
#' @param hmColors A color palette
#' @param anocol A matrix of colors
#' @param xpos Xpos
#' @param ypos Ypos
#' @param dendro.cex Size of denro names 
#' @param xlab.cex Size of x label
#' @param hmRowNames Write rownames ?
#' @param hmRowNames.cex Size of rownames ?
#'
#' @importFrom graphics par plot image axis
#' 
#' @return A heatmap
#'
hclustAnnotHeatmapPlot <- function(
    x = NULL, hc = NULL, hmColors = NULL, anocol = NULL, 
    xpos = c(0.1, 0.9, 0.114, 0.885), ypos = c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95),
    dendro.cex = 1, xlab.cex = 0.8, hmRowNames = FALSE, hmRowNames.cex = 0.5)
{
    par(fig = c(
        xpos[1], xpos[2], ypos[5], ypos[6]), new = FALSE,
        mar = c(0, 0, 1.5, 0))
    plot(hc, main = "Hierarchical clustering", xlab = "", sub = "", las = 2,
            cex = dendro.cex, cex.axis = dendro.cex)
    par(fig = c(xpos[3], xpos[4], ypos[3], ypos[4]), new = TRUE,
        mar = rep(0, 4))
    imageCol(anocol, xlab.cex = xlab.cex, ylab.cex = 0)
    par(fig = c(xpos[3], xpos[4], ypos[1], ypos[2]), new = TRUE,
        mar = rep(0, 4))
    image(t(x), axes = FALSE, xlab = "", ylab = "", col = hmColors)
    graphics::box()
    if (hmRowNames)
    {
        axis(4, at = seq(0, 1, length.out = nrow(x)), labels = rownames(x),
            las = 1, cex.axis = hmRowNames.cex)
    }
}

#' imageCol
#'
#' @param matcol A matrix of colors
#' @param strat Strat
#' @param xlab.cex X label size
#' @param ylab.cex Y label size
#' @param drawLines Draw lines ?
#' @param ... Additional parameters
#'
#' @importFrom graphics image axis abline par
#' @return A rectangular image
#'
imageCol <- function(
    matcol = NULL, strat = NULL, xlab.cex = 0.5,
    ylab.cex = 0.5, drawLines = c("none", "h", "v", "b")[1], ...)
{
    stopifnot(!is.null(ncol(matcol)))
    matcol <- matcol[, ncol(matcol):1,drop=FALSE]
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
    if (dim(csc)[2] == 1){
        csc <- matrix(as.numeric(unlist(csc)), nrow = dim(csc)[1])
    } else{
        csc <- matrix(as.numeric(csc), nrow = dim(csc)[1])
    }
    image(csc, col = as.vector(csc.colors), axes = FALSE, ...)
    if (xlab.cex != 0){
        if(ncol(csc)>1) axis(
            2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1),colnames(matcol), las = 2,
            tick = FALSE, cex.axis = xlab.cex, ...)
        if(ncol(csc)==1) axis(2, 0, colnames(matcol), las = 2, tick = FALSE,
                            cex.axis = xlab.cex, ...)
    }
    if (ylab.cex != 0){
        if(nrow(csc)>1) axis(3, 0:(dim(csc)[1] - 1)/(dim(csc)[1] - 1),
                            rownames(matcol), las = 2, tick = FALSE,
                            cex.axis = ylab.cex, ...)
        if(nrow(csc)==1)        axis(3, 0, rownames(matcol), las = 2, 
                                    tick = FALSE, cex.axis = ylab.cex, ...)
    }
    if (drawLines %in% c("h", "b")) 
        abline(h = -0.5:(dim(csc)[2] - 1)/(dim(csc)[2] - 1))
    graphics::box()
    if (drawLines %in% c("v", "b")) 
        abline(v = 0.5:(dim(csc)[1] - 1)/(dim(csc)[1] - 1))
    graphics::box()
}


#' annotToCol2
#'
#' @param annotS A color matrix
#' @param annotT A color matrix
#' @param missing Convert missing to NA
#' @param anotype Annotation type
#' @param maxnumcateg Maximum number of categories
#' @param categCol Categorical columns 
#' @param quantitCol Quantitative columns
#' @param plotLegend Plot legend ?
#' @param plotLegendFile Which file to plot legend ?
#'
#' @importFrom graphics par plot legend image 
#' @return A matrix of continuous or discrete colors
#' @export
#' @examples
#' data("scExp")
#' annotToCol2(SingleCellExperiment::colData(scExp), plotLegend = FALSE)
#' 
annotToCol2 <- function(
    annotS = NULL, annotT = NULL, missing = c("", NA), anotype = NULL, 
    maxnumcateg = 2, categCol = NULL, quantitCol = NULL, plotLegend = TRUE,
    plotLegendFile = NULL)
{
    if (is.null(ncol(annotS))){
        annotS <- data.frame(annotS)
        colnames(annotS) = colnames(annotT)
        rownames(annotS) = rownames(annotT)
    }
    for (j in seq_len(ncol(annotS)) ) {
    annotS[which(annotS[, j] %in% missing), j] <- NA}
    if (is.null(anotype)){
        anotype <- rep("categ", ncol(annotS))
        names(anotype) <- colnames(annotS)
        classes <- as.character(
        lapply(seq_len(ncol(annotS)), function(j) class(annotS[, j])))
        nmodal <- as.numeric(lapply(seq_len(ncol(annotS)), function(j) length(
            unique(setdiff(annotS[,j], NA)))))
        anotype[which(classes %in% c(
            "integer", "numeric") & nmodal > maxnumcateg)] <- "quantit"
        anotype[which(nmodal == 2)] <- "binary"
    }
    anocol <- annotS
    if (plotLegend) pdf(plotLegendFile)
    anocol <- anocol_categorical(anocol, categCol, anotype, plotLegend, annotS)
    anocol <- anocol_binary(anocol,anotype, plotLegend, annotS)
    if (is.null(quantitCol)) 
        quantitCol <- c(
            "darkgreen", "darkblue", "darkgoldenrod4", "darkorchid4", 
            "darkolivegreen4", "darkorange4", "darkslategray")
    k <- 1
    for (j in which(anotype == "quantit")){
        colrange <- matlab.like(100)
        anocol[, j] <- colrange[round(changeRange(anocol[, j],
                                                newmin = 1, newmax = 100))]
        if (k < length(quantitCol)) k <- k + 1 else k <- 1
        if (plotLegend){
            par(mar = c(8, 2, 5, 1))
            lims <- seq(-1, 1, length.out = 200)
            image(matrix(lims, ncol = 1), col = colrange,
                axes = FALSE, xlab = colnames(anocol)[j])
        }
    }
    if (plotLegend) dev.off()
    for (j in seq_len(ncol(anocol))) 
    anocol[which(is.na(anocol[, j])), j] <- "white"
    as.matrix(anocol)
}

#' Helper binary column for anocol function
#'
#' @param anocol The color feature matrix
#' @param anotype The feature types
#' @param plotLegend Plot legend ?
#' @param annotS A color matrix
#' @param categCol Colors for categorical features
#'
#' @return A color matrix similar to anocol with binrary columns colored
#'
anocol_categorical <- function(anocol, categCol, anotype, plotLegend, annotS){
    
    anocol. <- anocol
    if (is.null(categCol)) 
    categCol <- c(
        "#4285F4", "#DB4437", "#F4B400", "#0F9D58", "slategray", "black", 
        "orange", "turquoise4", "yellow3", "orangered4", "orchid",
        "palegreen2", "orchid4", "red4", "peru", "orangered",
        "palevioletred4", "purple", "sienna4", "turquoise1")
k <- 1
for (j in which(anotype == "categ")){
    tmp <- as.factor(anocol.[, j])
    classes <- as.character(levels(tmp))
    ncat <- length(levels(tmp))
    if(ncat > length(categCol)){
        color = grDevices::colors()[grep('gr(a|e)y',
                                        grDevices::colors(), invert = TRUE)]
        categCol <- sample(color, 200)
    }
    if (k + ncat > length(categCol)) 
        categCol <- c(categCol, categCol)
    levels(tmp) <- categCol[k:(k + ncat - 1)]
    fill <- as.character(levels(tmp))
    anocol.[, j] <- as.character(tmp)
    k <- k + ncat
    if (plotLegend){
        par(mar = c(0, 0, 0, 0))
        plot(-10, axes = FALSE, xlim = c(0, 5), ylim = c(0, 5), xlab = "",
                ylab = "")
        legend(1, 5, legend = classes, fill = fill,
                    title = colnames(anocol.)[j], xjust = 0.5, yjust = 1)
    }
}
return(anocol.)
}

#' Helper binary column for anocol function
#'
#' @param anocol The color feature matrix
#' @param anotype The feature types
#' @param plotLegend Plot legend ?
#' @param annotS A color matrix
#'
#' @return A color matrix similar to anocol with binrary columns colored
#'
anocol_binary <- function(anocol,anotype, plotLegend, annotS){
    memcol <- c()  
    anocol. <- anocol
    for (j in which(anotype == "binary"))
    {
        new <- setdiff(anocol.[, j], c(NA, memcol))
        if (length(new) == 2){
        memcol <- c(memcol, c("aquamarine", "plum1"))
        names(memcol)[(length(memcol) - 1):length(memcol)] <- sort(new)
        }
        if (length(new) == 1){
        memcol <- c(memcol, setdiff(c("dodgerblue4", "firebrick"),
                                    memcol[setdiff(anocol.[,j], c(NA,new))]))
        names(memcol)[length(memcol)] <- new
        }
        anocol.[, j] <- as.character(anocol.[, j])
        for (z in seq_along(memcol)){
            anocol.[which(anocol.[, j] == names(memcol)[z]), j] <- memcol[z]
        }
        if (plotLegend){
        par(mar = c(0, 0, 0, 0))
        plot(-10, axes = FALSE, xlim = c(0, 5), ylim = c(0, 5), xlab = "",
            ylab = "")
        classes <- intersect(names(memcol), annotS[, j])
        fill <- memcol[classes]
        legend(1, 5, legend = classes, fill = fill,
                title = colnames(anocol.)[j], xjust = 0.5, yjust = 1)
        }
    }
    return(anocol.)
}

#' groupMat
#'
#' @param mat A matrix
#' @param margin By row or columns ?
#' @param groups Groups
#' @param method Method to group
#'
#' @return A grouped matrix
#'
groupMat <- function(mat = NA, margin = 1, groups = NA, method = "mean")
{
    if (!method %in% c("mean", "median")){
        print("Method must be mean or median")
        return() }
    if (!margin %in% seq_len(2)){
        print("Margin must be 1 or 2")
        return() }
    for (i in seq_along(groups)){
        if (margin == 1){
            if (length(groups[[i]]) == 1){
                v <- mat[, groups[[i]]]
            } else{
                if (method == "mean") v <- apply(mat[, groups[[i]]], margin,
                                                mean)
                else v <- apply(mat[, groups[[i]]], margin, median)
            }
            if (i == 1) res <- matrix(v, ncol = 1)
            else res <- cbind(res, v)
        } else{
            if (length(groups[[i]]) == 1)
            {
                v <- mat[groups[[i]], ]
            } else{
                if (method == "mean"){
                    v <- apply(mat[groups[[i]], ], margin, mean)
                } else{
                    v <- apply(mat[groups[[i]], ], margin, median)
                }
            }
            if (i == 1){
                res <- matrix(v, nrow = 1)
            } else{
                res <- rbind(res, v)
            }
        }
    }
    if (margin == 1){
        rownames(res) <- rownames(mat)
        colnames(res) <- names(groups)
    } else{
        rownames(res) <- names(groups)
        colnames(res) <- colnames(mat)
    }
    res
}