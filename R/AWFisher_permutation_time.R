##' @export
AWFisher_permutation_time <- function(dataExp, labels, B=1000, seed=15213){

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
	
	time_start <- Sys.time()
	
	time_limma <- time_start - time_start
	time_AWstat <- time_start - time_start
	time_saveNull <- time_start - time_start
	time_rank <- time_start - time_start
	

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
		
		time_start <- Sys.time()
		Z_b <- myLimma(dataExp,labels_b)
		time_end <- Sys.time()
		time_limma <- time_limma + time_end - time_start
		
		p.matrix_b <- 2*pnorm(-abs(Z_b))

		time_start <- Sys.time()
		AWstat_null0[,b] <- apply(p.matrix_b,1,getAWstat)		
		time_end <- Sys.time()
		time_AWstat <- time_AWstat + time_end - time_start

		time_start <- Sys.time()
		AWstat_null <- c(AWstat, AWstat_null0)
		time_end <- Sys.time()
		time_saveNull <- time_saveNull + time_end - time_start
	}

	time_start <- Sys.time()
	pval <- rank(AWstat_null)[seq_along(AWstat)]/length(AWstat_null)
	time_end <- Sys.time()
	time_rank <- time_rank + time_end - time_start
	
	time <- list(time_limma=time_limma, time_AWstat=time_AWstat, time_saveNull=time_saveNull, time_rank=time_rank)
    return(time)
}

