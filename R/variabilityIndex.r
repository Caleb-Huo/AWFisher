##' Variability Index
##'
##' Variability Index via boostrap AW weight.
##' @title Variability Index
##' @param studies a list of K studies. 
##' Each element (kth study) of the list is another list consisting 
##' gene expression matrix and and label information. 
##' @param afunction A function for DE analysis. 
##' Options can be function_limma or function_edgeR. 
##' Default option is function_limma. 
##' However, use could define their own function. 
##' The input of afunction should be list(data, label) 
##' which is consistent with one element of the studies list/argument. 
##' The return of afunction should be list(pvalue=apvalue, effectSize=aeffectsize)
##' @param B number of permutation should be used. B=1000 is suggested.
##' @param silence If TRUE, will print out the bootstrapping procedure.
##' @return A list consisting of biomarker categrorization result.
##' \item{varibility}{Varibility index for all genes}
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
##' for(k in 1:K){
##'     astudy <- matrix(rnorm(N0*2*G),nrow=G,ncol=N0*2)
##'     ControlLabel <- 1:N0
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
##' result <- variabilityIndex(studies,function_limma,B=100)
##' head(result)

variabilityIndex <- function(studies, afunction, B = 10, silence = FALSE) {
    Tscore <- 0
    pval.null <- NULL
    weight.null <- NULL
    
    if (B <= 2) {
        stop("B has be to greater than 2!")
    }
    
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
        
        awres.null <- AWFisher_pvalue(pval.null)
        weight.null[[b]] <- awres.null$weights
    }  ## b for end of B
    cat("\n", "calculating variability index", "\n")
    
    aveScore <- Reduce("+", weight.null)/B
    for (b in seq_len(B)) {
        Tscore = Tscore + (weight.null[[b]] - aveScore)^2/B
    }
    
    varibility <- Tscore * 4
    varibility
}

