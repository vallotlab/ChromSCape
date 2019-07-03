# Rd
# description >> This function tests for gene set enrichements in a gene list using the hypergeometric test
# argument
# item >> gene.sets >> List of gene sets to be tested
# item >> mylist >> List of selected genes
# item >> possibleIds >> List of genes that could have possibly been selected (e.g. genes present on the array)
# item >> sep >> required
# item >> silent >> required
# value >> required
# author >> Eric Letouze
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
geco.enrichmentTest <- function (gene.sets, mylist, possibleIds, sep = ";", silent = F) 
{
    possibleIds <- unique(possibleIds)
    mylist <- unique(mylist)
    gene.sets <- lapply(gene.sets, unique)
    nids <- length(possibleIds)
    gene.sets <- lapply(gene.sets, function(x) intersect(x, possibleIds))
    nref <- sapply(gene.sets, length)
    if (all(nref == 0)) stop("Error: no intersection between gene sets and possible IDs.")
    if (any(nref == 0)) print("Warning: some of the gene sets have no intersection with possibleIds")
    if (!all(mylist %in% possibleIds)) stop("Error: some genes in mylist are not in possibleIds")
    if (!silent) cat(paste("NB : enrichment tests are based on", nids, "distinct ids.\n"))
    gene.sets <- gene.sets[nref > 0]
    n <- length(mylist)
    fun <- function(x) {
        y <- intersect(x, mylist)
        nx <- length(x)
        ny <- length(y)
        pval <- phyper(ny - 1, nx, nids - nx, n, lower.tail = F)
        c(nx, ny, pval,paste(y, collapse = sep))
    }
    tmp <- as.data.frame(t(sapply(gene.sets, fun)))
    rownames(tmp) <- names(gene.sets)
    for (i in 1:3) tmp[,i] <- as.numeric(as.character(tmp[,i]))
    tmp <- data.frame(tmp[,1:3],p.adjust(tmp[,3],method="BH"),tmp[,4])
    names(tmp) <- c("Nb_of_genes","Nb_of_deregulated_genes","p-value","q-value","Deregulated_genes")
    tmp
}

