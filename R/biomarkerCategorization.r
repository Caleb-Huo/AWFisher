##' biomarker categrorization
##'
##' biomarker categrorization via boostrap AW weight.
##' @title biomarker categrorization
##' @param studies a list of K studies. 
##' Each element (kth study) of the list is another list 
##' consisting gene expression matrix and and label information. 
##' @param afunction A function for DE analysis. 
##' Options can be function_limma or function_edgeR. 
##' Default option is function_limma. 
##' However, use could define their own function. 
##' The input of afunction should be list(data, label) 
##' which is consistent with one element of the studies list/argument. 
##' The return of afunction should be 
##' list(pvalue=apvalue, effectSize=aeffectsize)
##' @param B number of permutation should be used. B=1000 is suggested.
##' @param DEindex If NULL, 
##' BH method will be applied to p-values and FDR 0.05 will be used. 
##' User could specify a logical vector as DEindex.
##' @param fdr Default is 0.05. 
##' The co-membership matrix calculation will base on 
##' genes with this specified fdr.
##' @param silence If TRUE, will print out the bootstrapping procedure.
##' @return A list consisting of biomarker categrorization result.
##' \item{varibility}{Varibility index for all genes}
##' \item{dissimilarity}{Dissimilarity matrix of genes of DEindex==TRUE}
##' \item{DEindex}{DEindex for Dissimilarity}
##' @author Zhiguang Huo
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
##' for(k in seq_len(K)){
##'     astudy <- matrix(rnorm(N0*2*G),nrow=G,ncol=N0*2)
##'     ControlLabel <- seq_len(N0)
##'     caseLabel <- (N0 + 1):(2*N0)
##' 
##'     astudy[1:GDEp,caseLabel] <- astudy[1:GDEp,caseLabel] + 2
##'     astudy[1:GDEp + GDEn,caseLabel] <- astudy[1:GDEp + GDEn,caseLabel] - 2
##' 
##'     alabel = c(rep(0,length(ControlLabel)),rep(1,length(caseLabel)))
##' 
##'     studies[[k]] <- list(data=astudy, label=alabel)
##' }
##' 
##' 
##' result <- biomarkerCategorization(studies,function_limma,B=100,DEindex=NULL)
##' sum(result$DEindex)
##' head(result$varibility)
##' print(result$dissimilarity[1:4,1:4])

biomarkerCategorization <- 
  function(studies, afunction, B = 10, 
           DEindex = NULL, fdr = NULL, silence = FALSE) {
    if (B <= 2) stop("B has be to greater than 2!")
    if (is.null(fdr)) fdr = 0.05
    if (is.null(DEindex)) {
        cat("generate DE index since it is NULL", "\n")
        cat("based on AW fdr ", fdr, "\n")
        res.obs <- getPvalueAll(studies, afunction)
        pval.obs <- res.obs$p.matrix
        awres.obs <- AWFisher_pvalue(pval.obs)
        awres.obs$es <- res.obs$es.matrix
        fdr.obs <- p.adjust(awres.obs$pvalues, "BH")
        DEindex <- fdr.obs <= fdr
    }
    Tscore <- 0
    es.null <- NULL
    pval.null <- NULL
    weight.null <- NULL
    selfDistDirection <- matrix(0, nrow = sum(DEindex), ncol = sum(DEindex))
    if (!silence) 
        cat("calculating permutated score, b = 1,2,..., B (= ", 
            B, 
            ")  [one \".\" per sample]:\n", sep = "")
    for (b in seq_len(B)) {
        cat(".", if (b%%50 == 0) 
            paste(b, "\n"))
        studies.b <- permuteLabels(studies)
        res.null <- getPvalueAll(studies.b, afunction)
        pval.null <- res.null$p.matrix
        es.null <- res.null$es.matrix
        awres.null <- AWFisher_pvalue(pval.null)
        weight.null[[b]] <- awres.null$weights
        aDEweightDirection <- (awres.null$weights * sign(es.null))[DEindex, ]
        maxDistDirection <- dist(aDEweightDirection, "maximum")
        bselfDistDirection <- ifelse(as.matrix(maxDistDirection) == 0, 1, 0)
        selfDistDirection <- selfDistDirection + as.matrix(bselfDistDirection)/B
    }  ## b for end of B
    cat("\n", "calculating variability index", "\n")
    aveScore <- Reduce("+", weight.null)/B
    for (b in seq_len(B)) {
        Tscore = Tscore + (weight.null[[b]] - aveScore)^2/B
    }
    result <- list(varibility = Tscore * 4, 
                   dissimilarity = selfDistDirection, 
                   AWres = awres.obs, 
                   DEindex = DEindex)
    result
}

