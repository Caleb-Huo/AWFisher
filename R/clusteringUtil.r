permuteLabels <- function(studies,seed=15213){
	set.seed(seed)
	for(k in 1:length(studies)){
		studies[[k]]$label <- sample(studies[[k]]$label)
	}
	studies
}


function_limma <- function(astudy){
	alabel <- astudy$label
	aData <- astudy$data
	
	adesign = cbind(rep(1,length(alabel)),alabel)
	afit <- lmFit(aData,adesign)
	afit <- eBayes(afit)		

	aeffectsize <- afit$coefficients[,2]
	apvalue <- afit$p.value[,2]
	return(list(pvalue=apvalue, effectSize=aeffectsize))
}

getPvalueSingle <- function(astudy, afunction){
	## input a study and label as a list
	## return p value and effect size
	return(afunction(astudy))
}


getPvalueAll <- function(studies, afunction){
	K <- length(studies)
	
	p.matrix <- NULL
	es.matrix <- NULL
	
	for(k in 1:K){
		astudy <- studies[[k]]
		ares <- afunction(astudy)
		p.matrix <- cbind(p.matrix, ares$pvalue)
		es.matrix <- cbind(es.matrix, ares$effectSize)		
	}
	return(list(p.matrix=p.matrix, es.matrix=es.matrix))
}

