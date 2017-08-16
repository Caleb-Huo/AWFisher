##' use limma function to get pvalue
##'
##' use limma function to get pvalue
##' @title use limma function to get pvalue
##' @param astudy A list contains a data matrix and a vector of group label
##' @return A list of pvalue and effect size
##' @author Caleb
##' @export
##' @examples
##'
##' N0 = 10
##' G <- 1000
##' GDEp <- 50
##' GDEn <- 50
##'
##' set.seed(15213)
##'
##' astudy <- matrix(rpois(N0*2*G,10),nrow=G,ncol=N0*2)
##' ControlLabel <- 1:N0
##' caseLabel <- (N0 + 1):(2*N0)
##'
##' astudy[1:GDEp,caseLabel] <- astudy[1:GDEp,caseLabel] + 2
##' astudy[1:GDEp + GDEn,caseLabel] <- astudy[1:GDEp,caseLabel] - 2
##'
##' alabel <- c(rep(0,length(ControlLabel)),rep(1,length(caseLabel)))
##' Study <- list(data=astudy, label=alabel)
##'
##' result <- function_edgeR(Study)
##' fdr <- p.adjust(result$pvalue)
##' sum(fdr<=0.05)
##'

function_edgeR <- function(astudy){
	alabel <- astudy$label
	aData <- astudy$data
	
	group<-factor(alabel)
	d=DGEList(counts=aData,group=group)
	d <- calcNormFactors(d) # normalization part
	design=model.matrix(~group)
	d=estimateGLMCommonDisp(d)
	d=estimateGLMTrendedDisp(d)
	d=estimateGLMTagwiseDisp(d)
	fit <- glmFit(d, design)
	lrt <- glmLRT(fit, coef=2)
	apvalue=lrt$table[,4]
	aeffectsize <- fit$coefficient[,2]
	
	return(list(pvalue=apvalue, effectSize=aeffectsize))
}

