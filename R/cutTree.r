##' Cut a dendrogram
##'
##' Automatically cut a dendrogram. The tree cut is starting is using fix height. Each time lower the height until two new clusters occur. Define the new clusters as active clusters. Stop until the new active clusters are not generated from the old active clusters.
##' @title Cut a dendrogram
##' @param ddr A dendrogram object.
##' @param K Maximum number of clusters.
##' @return clustering membership.
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
##' set.seed(15212)
##' for(k in 1:K){
##' 	astudy <- matrix(rnorm(N0*2*G),nrow=G,ncol=N0*2)
##' 	ControlLabel <- 1:N0
##' 	caseLabel <- (N0 + 1):(2*N0)
##'
##' 	astudy[1:GDEp,caseLabel] <- astudy[1:GDEp,caseLabel] + 2
##' 	astudy[1:GDEp + GDEn,caseLabel] <- astudy[1:GDEp,caseLabel] - 2
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
##' table(clusterMembership)

cutTree <- function(ddr, K=10){
	clusterMembership0 <- cutree(ddr,k=1:K) # it now works on a dendrogram
	k <- 2
	activeMembership <- unique(clusterMembership0[,k])
	while(k<K){
		print(k)
		logic1 <- clusterMembership0[,k] == activeMembership[1]
		logic2 <- clusterMembership0[,k] == activeMembership[2]
		nextClusterMembership <- clusterMembership0[,k+1]
		unique1 <- unique(nextClusterMembership[logic1])
		unique2 <- unique(nextClusterMembership[logic2])
		if(length(unique1)==2){
			activeMembership <- unique1
		} else if(length(unique2)==2) {
			activeMembership <- unique2
		} else {
			break
		}
		k <- k + 1
	}
	return(clusterMembership0[,k])
}
