##' How to struct hierarchical tree stucture from a dissimilarity matrix
##'
##' How to struct hierarchical tree stucture from a dissimilarity matrix
##' @title hierarchical tree clustering
##' @param dissimilarity Input dissimilarity matrix
##' @param method hierarchical clustering method.
##' @return a dendrogram. Note the the dendrogram is re-ordered by row mean of the dissimilarity matrix.
##' @author Caleb
##' @export
##' @examples
##' N0 = 10
##' G <- 1000
##' GDEp <- 50
##' GDEn <- 50
##' K = 4
##'
##' studies <- NULL
##' set.seed(15213)
##' for(k in 1:K){
##' 	astudy <- matrix(rnorm(N0*2*G),nrow=G,ncol=N0*2)
##' 	ControlLabel <- 1:N0
##' 	caseLabel <- (N0 + 1):(2*N0)
##'
##' 	astudy[1:GDEp,caseLabel] <- astudy[1:GDEp,caseLabel] + 2
##' 	astudy[1:GDEp + GDEn,caseLabel] <- astudy[1:GDEp + GDEn,caseLabel] - 2
##'
##' 	alabel = c(rep(0,length(ControlLabel)),rep(1,length(caseLabel)))
##'
##' 	studies[[k]] <- list(data=astudy, label=alabel)
##' }
##'
##'
##' result <- biomarkerCategorization(studies,function_limma,B=100,DEindex=NULL,seed = 15213)
##' dissimilarity <- result$dissimilarity
##' atree <- hclustTree(dissimilarity)
##' plot(atree)

hclustTree <- function(dissimilarity,method='ward.D'){
	reorderfun=function(d, w) reorder(d, w)
	hclustfun = function(x) hclust(x,method = 'ward.D')
	Rowv <- rowMeans(dissimilarity)
	distr <- dist(dissimilarity)
	hcr <- hclustfun(distr)
	ddr <- as.dendrogram(hcr)
	ddr <- reorderfun(ddr, Rowv)
	#rowInd <- order.dendrogram(ddr)
	return(ddr)
}
