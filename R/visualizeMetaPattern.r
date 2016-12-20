##' visualize MetaPattern
##'
##'visualize MetaPattern
##' @title visualize MetaPattern
##' @param studies Input list of studies. See biomarkerCategorization function for detail information.
##' @param biomarkerCategorization biomarkerCategorization result
##' @param clusterMembership clustering membership from cutTree funciton
##' @param groups relabel the gene modules.
##' @param ... Extra paramenters
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
##' visualizeMetaPattern(studies, result, clusterMembership, labRow='',labCol='',na.color=par("bg"),lwid=c(0.1,4), lhei=c(0.1,15),margins = c(1, 1))

visualizeMetaPattern <- function(studies, biomarkerCategorization, clusterMembership, groups=NULL,...){
	rainbow_fun <- function(n, c=90, l=50, ...) {
	if(requireNamespace("colorspace")) {
	colorspace::rainbow_hcl(n, c = c, l = l, ...)		
		} else {
			rainbow(n, ...)
		}
	}
	
	DEindex <- biomarkerCategorization$DEindex
	DEdata <- lapply(studies,function(x) {x$data <- x$data[DEindex,]; return(x)})

	if(is.null(groups)){
		groups <- 1:length(unique(clusterMembership))
	}
	allData <- NULL

	for(i in groups){
		cat('membership: ',i,'\n')

		## gene order based on variability score
		geneOrder <- order(rowSums(biomarkerCategorization$varibility[DEindex,])[clusterMembership==i])
		acolor <- rainbow_fun(max(clusterMembership))[i]

		rowData <- NULL
		rowLabel <- NULL

		## visualize the heatmap
		for(j in 1:length(DEdata)){
			adata <- DEdata[[j]]$data[clusterMembership==i,][geneOrder,]
			bdata <- t(scale(t(adata)))
			alabel <- DEdata[[j]]$label + 1
			rowData <- cbind(rowData,bdata,NA)
			rowLabel <- c(rowLabel, alabel, NA)
		}
		allData <- rbind(allData,rowData, matrix(NA,nrow=10,ncol=ncol(rowData)))
	}

	colColor <- character(length(rowLabel))
	colColor[is.na(rowLabel)] <- 'white'
	colColor[rowLabel==1] <- 'black'
	colColor[rowLabel==2] <- 'orange'

	heatmap.2(allData,ColSideColors = colColor,col="greenred", trace="none",Rowv=NA,Colv=NA,key=FALSE,...)
}

