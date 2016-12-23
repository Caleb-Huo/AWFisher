##' Visualize Tree Structure
##'
##' Visualize Tree Structure
##' @title Visualize Tree Structure
##' @param ddr Input dendrogram
##' @param clusterMembership clusterMembership result from previous cut tree
##' @param ... Other parameters for plot
##' @return NULL
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
##' ddr <- hclustTree(dissimilarity)
##' clusterMembership <- cutTree(ddr)
##' visualizeTree(ddr, clusterMembership)
##'
##'
visualizeTree <- function(ddr, clusterMembership, ...){
	k <- length(unique(clusterMembership))
	ddr <- color_branches(ddr, k = k)
	ddr <- color_labels(ddr, k = k)
	plot(ddr,...)
}

