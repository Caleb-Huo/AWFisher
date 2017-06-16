##' @export
AWFisher_permutation <- function(dataExp, labels, B=1000, seed=15213){

	getAWstat <- function(avec){
		asort <- sort(avec)
		alog <- -2 * log(asort)
		cumlog <- cumsum(alog)
		aAWstat <- min(pchisq(cumlog, 2*seq_along(cumlog), lower.tail = FALSE))
		return(aAWstat)
	}

	## observed AW stat
	Z <- myLimma(dataExp,labels)
	p.matrix <- 2*pnorm(-abs(Z))
	AWstat <- apply(p.matrix,1,getAWstat)

	AWstat_null0 <- matrix(,nrow=length(AWstat),ncol=B)

	## observed AW stat
	for(b in 1:B){
		set.seed(seed+b)
		labels_b <- lapply(labels,function(x){
			pool <- c(x[[1]],x[[2]])
			x[[1]] <- sample(pool, length(x[[1]]))
			x[[2]] <- setdiff(pool, x[[1]])
			x
		})
		sample(labels)
		Z_b <- myLimma(dataExp,labels_b)
		p.matrix_b <- 2*pnorm(-abs(Z_b))
		AWstat_null0[,b] <- apply(p.matrix_b,1,getAWstat)		
		AWstat_null <- c(AWstat, AWstat_null0)
	}

	pval <- rank(AWstat_null)[seq_along(AWstat)]/length(AWstat_null)
	
    return(pval)
}

