##' @export
AWFisher_close2 <- function(x){
  if(NCOL(x) == 1) x= t(x)
	AW.Fisher <- function(t){
	  ### this function is used to compute the p-values ########
	  ### t: vector of test statistics
	  t1 <- t[,1]
	  t2 <- t[,2]
 
	  temp = exp(-t1/2)
	  p <- temp*(2-temp)
	  select <- t1 >= t2/2
	  p0 <- 2*exp(-t1/2) + ( - t2/2 + t1-1)*exp(-t2/2)
	  p[select] <- p0[select]
	  return(p)
	}

  #### this function is used to return the p-value and weights #######
  a <- -2*log(x)
  T2 <- rowSums(a)
  
  p1 <- pmin(x[,1],x[,2])
  p2 <- pchisq(T2,4,lower.tail=FALSE)
  
  p_aw <- pmin(p1,p2)
  
  t <- cbind(qchisq(p_aw,2,lower.tail=FALSE),qchisq(p_aw,4,lower.tail=FALSE))
  p <- AW.Fisher(t)

  return(p)
}

if(F){
  ## test example
  pvalues <- matrix(c(0.1,0.01,0.2,0.02,0.3,0.03),nrow=3)
  AWFisher_close2(pvalues)
  
}
