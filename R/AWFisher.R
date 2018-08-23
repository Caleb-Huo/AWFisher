##' R package for fast computing for adaptively weighted fisher's method
##'
##' fast computing for adaptively weighted fisher's method
##' @title AWFisher
##' @param p.values Input G by K p-value matrix. Each row represent a gene and each column represent a study. Note that K has to be >=2 and <=100.
##' @return A list consisting of AWFisher pvalues and AWweight.
##' \item{pvalues}{AWFisher pvalues.}
##' \item{weights}{G by K binary weight matrix W. $W_{gk} = 1$ represents for gene $g$, study $k$ contributes to the meta-analysis result. $W_{gk} = 0$ otherwise.}
##' @useDynLib AWFisher
##' @author Caleb
##' @export
##' @examples
##' K <- 40
##' G <- 10000
##' p.values = matrix(rbeta(K*G, 1,1), ncol=K)
##' res = AWFisher.pvalue(p.values)
##' hist(res$pvalues, breaks=40)
##' table(rowSums(res$weights))
##' pvalues=res$pvalues[order(res$pvalues)]
##' plot(-log10((1:NROW(pvalues))/(1+NROW(pvalues))), -log10(pvalues),xlab="theoretical quantile", ylab="observed quantile")
##' lines(c(0,100), c(0,100), col=2)

AWFisher.pvalue <- function(p.values) {
    if(NCOL(p.values) == 1) p.values= t(p.values)

    n = NCOL(p.values)
	if(n<2 | n>100){
		stop("number of studies K has to be >0 2 and <= 100.")
	}

    out <- .C("AWpvalue", best_stat=rep(0, NROW(p.values)), sum_weight=as.integer(rep(0, NROW(p.values))), weights=as.integer(rep(0, length(p.values))),
                    pval=t(p.values), nrow=as.integer(NROW(p.values)), ncol=as.integer(NCOL(p.values)))
    bestStat <- out $ best_stat
    weights <- t(matrix(out $ weights, ncol=NROW(p.values)))
    list(pvalues = aw.fisher.stat(bestStat, n), weights = weights)

}