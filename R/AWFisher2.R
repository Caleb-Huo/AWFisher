##' R package for fast computing for adaptively weighted fisher's method
##'
##' fast computing for adaptively weighted fisher's method
##' @title AWFisher
##' @param p.values Input G by K p-value matrix. Each row represent a gene and each column represent a study.
##' @return A list consisting of AWFisher pvalues and AWweight.
##' \item{pvalues}{AWFisher pvalues.}
##' \item{weights}{G by K binary weight matrix W. $W_{gk} = 1$ represents for gene $g$, study $k$ contributes to the meta-analysis result. $W_{gk} = 0$ otherwise.}
##' @author Caleb
##' @export
##' @examples
##' K <- 154
##' G <- 10000
##' p.values = matrix(rbeta(K*G, 1,1), ncol=K)
##' res = AWFisher.pvalue(p.values)
##' hist(res$pvalues, breaks=40)
##' table(rowSums(res$weights))
##' pvalues=res$pvalues[order(res$pvalues)]
##' plot(-log10((1:NROW(pvalues))/(1+NROW(pvalues))), -log10(pvalues),xlab="theoretical quantile", ylab="observed quantile")
##' lines(c(0,100), c(0,100), col=2)

AWFisher.pvalue2 <- function(p.values) {
    if(NCOL(p.values) == 1) p.values= t(p.values)

    k = NCOL(p.values)
	stopifnot(k>=2)
	n <- NROW(p.values)
	
	awStat.all <- p.values
    orderMatrix = t(apply(p.values, 1, order))
    for(i in 1:n) {
        awStat.all[i,] =  - pchisq(-2*cumsum(log(p.values)[i,orderMatrix[i,]]), 2*(1:k), lower.tail=F, log=T) 		
    }
	bestStat <- apply(awStat.all,1,max)
	weights <- ifelse(orderMatrix - apply(awStat.all,1,which.max) <= 0, 1, 0)
	
    list(pvalues = aw.fisher.stat(bestStat, k), weights = weights)

}
