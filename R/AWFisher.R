##' R package for fast computing for adaptively weighted fisher's method
##'
##' fast computing for adaptively weighted fisher's method
##' @title AW
##' @param p.values 1
##' @param method 1
##' @param log 1
##' @return 1
##' @author Caleb
##' @export
##' @examples
##' n=154
##' p.values = matrix(rbeta(n*10000, 1,1), ncol=n)
##' res = AWFisher.pvalue(p.values)
##' hist(res$pvalues, breaks=40)
##' table(rowSums(res$weights))
##' pvalues=res$pvalues[order(res$pvalues)]
##' plot(-log10((1:NROW(pvalues))/(1+NROW(pvalues))), -log10(pvalues),xlab="theoretical quantile", ylab="observed quantile")
##' lines(c(0,100), c(0,100), col=2)

AWFisher.pvalue <- function(p.values) {
    if(NCOL(p.values) == 1) p.values= t(p.values)

    n = NCOL(p.values)

    out <- .C("AWpvalue", best_stat=rep(0, NROW(p.values)), sum_weight=as.integer(rep(0, NROW(p.values))), weights=as.integer(rep(0, length(p.values))),
                    pval=t(p.values), nrow=as.integer(NROW(p.values)), ncol=as.integer(NCOL(p.values)))
    bestStat <- out $ best_stat
    weights <- t(matrix(out $ weights, ncol=NROW(p.values)))
    list(pvalues = aw.fisher.stat(bestStat, n), weights = weights)

}
