permuteX <- function(dataExp, labels) {
    permx <- list()
    for (i in seq_along(dataExp)) {
        X <- matrix(, nrow = nrow(dataExp[[i]]), ncol(dataExp[[i]]))
        aData = dataExp[[i]]
        alabel = labels[[i]]
        
        for (j in seq_along(alabel)) {
            aalabel <- alabel[[j]]
            for (k in aalabel) {
                X[, k] <- aData[, sample(aalabel, 1)]
            }
        }
        permx[[i]] = X
    }
    permx
}

permuteLabels <- function(studies) {
    for (k in seq_along(studies)) {
        adata <- studies[[k]]$data
        alabel <- studies[[k]]$label
        uniqueLabels <- unique(alabel)
        for (j in seq_along(uniqueLabels)) {
            aalabel <- which(uniqueLabels[j] == alabel)
            adata[, aalabel] <- adata[, sample(aalabel, replace = TRUE)]
            studies[[k]]$data <- adata
        }
    }
    studies
}


getPvalueSingle <- function(astudy, afunction) {
    ## input a study and label as a list return p value and effect size
    return(afunction(astudy))
}


getPvalueAll <- function(studies, afunction) {
    p.matrix <- NULL
    es.matrix <- NULL
    
    for (k in seq_along(studies)) {
        astudy <- studies[[k]]
        ares <- afunction(astudy)
        p.matrix <- cbind(p.matrix, ares$pvalue)
        es.matrix <- cbind(es.matrix, ares$effectSize)
    }
    return(list(p.matrix = p.matrix, es.matrix = es.matrix))
}

