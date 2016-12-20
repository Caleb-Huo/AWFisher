##' Visualize dissimilarity heatmap
##'
##' Visualize dissimilarity heatmap
##' @title Visualize dissimilarity heatmap
##' @param dissimilarity Input dissimilarity matrix
##' @param ddr Input dendrogram
##' @param clusterMembership clusterMembership result from previous cut tree
##' @param ... Other parameters for heatmap.2
##' @return heatmap.2 object
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
##' visualizeDissimilarity(dissimilarity, ddr, clusterMembership)
##'
##'
visualizeDissimilarity <- function(dissimilarity, ddr, clusterMembership,...){
	reorderLabel <- function(alabel,aorder){
	  uniLab = unique(alabel)
	  resLabel = NULL
	  for(i in 1:length(uniLab)){
	    resLabel[uniLab[i]==alabel]=aorder[i]
	  }
	  return(resLabel)
	}

	anOrder <- order.dendrogram(ddr)
	dissimilarity_order <- dissimilarity[anOrder,anOrder]

	clusterMembership_order <- clusterMembership[order.dendrogram(ddr)]

	colorBar <- NULL
	membershipTable <- table(clusterMembership_order)
	n <- length(membershipTable)
	allColor <- rainbow_fun(n)
	coll0 <- allColor[clusterMembership_order]
	coll <- reorderLabel(coll0, allColor)
	#allColor[n] <- '#777777'

	my_palette <- colorRampPalette(c("yellow2", "black", "skyblue"))(n = 300)

	heatmap.2(dissimilarity_order,ColSideColors = coll,col=my_palette, trace="none",Rowv=NA,Colv=NA,key=FALSE,...)

}

