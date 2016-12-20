##' visualizeVariability
##'
##' visualizeVariability
##' @title visualizeVariability
##' @param biomarkerCategorization biomarkerCategorization object from biomarkerCategorization function
##' @param clusterMembership clustering membership
##' @param groups re-order the module order.
##' @param ... extra parameters
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
##' visualizeVariability(result, clusterMembership, labRow='',labCol='',na.color=par("bg"),lwid=c(0.1,4), lhei=c(0.1,15),margins = c(1, 1))
##'
visualizeVariability <- function(biomarkerCategorization, clusterMembership, groups=NULL,...){
	DEindex <- biomarkerCategorization$DEindex

	if(is.null(groups)){
		groups <- 1:length(unique(clusterMembership))
	}

	rowDataVariability <- NULL
	for(i in groups){
		## visualize the weight with direction
		geneOrder <- order(rowSums(biomarkerCategorization$varibility[DEindex,])[clusterMembership==i])
		acolor <- rainbow_fun(max(clusterMembership))[i]

		## varibility score
		avaribility <- biomarkerCategorization$varibility[DEindex,][clusterMembership==i,][geneOrder,]
		## sum of varibility score
		aRowVaribility <- rowMeans(avaribility)
		aRowVariability <- cbind(avaribility,aRowVaribility)
		rowDataVariability <- rbind(rowDataVariability,aRowVariability,matrix(NA,nrow=10,ncol=ncol(aRowVariability)))
	}

	B = 100
	heatmap.2(rowDataVariability,col=gray.colors(B), trace="none",Rowv=NA,Colv=NA,key=FALSE, ...)

}


