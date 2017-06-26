WD <- "/mnt/glusterfs/zhh18/AW/AWFisher/importanceSampling"

setwd(WD)

library(foreach)
library(doParallel)


args = commandArgs(trailingOnly = TRUE)
## args: N, K, sigma1, seed
args_num <- as.numeric(args)

k <- args_num[1]
n1 <- args_num[2]
n0 <- args_num[3]

cat("k: ", k, "\n")
cat("n1: ", n1, "\n")
cat("n0: ", n0, "\n")

if(F){
	n0 <- 1e3
	n1 <- 1e6
	k <- 2
}


pTargetList0 <- c(seq(0.99,0.01, by=-0.01))
pTargetList1 <- c(0.1^{3:100})

pTargetList <- c(pTargetList0, pTargetList1)


medianReturn <- function(a, n, k, seed = 15213) {
  set.seed(seed)
  pmatrix <- matrix(rbeta(n*k,a,1),n,k)
  
  pCumLog <- pmatrix
  
  for(j in 1:n){
    pCumLog[j,] <- cumsum(log(sort(pmatrix[j,])))
  }
  
  logweights <-  - k*log(a) - (a-1) * pCumLog[,k]
  statNew = -pchisq(-2*pCumLog[,1], 2, lower.tail=F, log=T)
  for(i in 2:k) 
    statNew = pmax(statNew, -pchisq(-2*pCumLog[,i], 2*i, lower.tail=F, log=T))

  ordStat = order(-statNew)
  logProbNew = -log(cumsum(exp(logweights[ordStat]))/(n+1))
  
  logProbNew[logProbNew>800] = 800
  
  median(logProbNew)
}

	
cl<-makeCluster(28)
registerDoParallel(cl)

result <- foreach(p = 1:length(pTargetList)) %dopar% {	
	apTarget <- pTargetList[p]

	setwd(WD)
	dyn.load('importanceSampling.so')
	kFolder <- paste0("k",k)
	system(paste("mkdir -p", kFolder))
	setwd(kFolder)
		
	## problematic
    a <- uniroot(function(a, n, k, target) medianReturn(a, n, k) - -log(target), 
            lower=1e-10,upper=1e3, n=n0, k = k, target = apTarget)$root		

	awStats <- .C('importanceSampling_R',a=as.double(a),n=as.integer(n1),k=as.integer(k),
				pTarget=as.double(apTarget),J=as.integer(length(apTarget)),awStat_vec=as.double(apTarget))$awStat_vec
				
	end <- Sys.time()
	timeDiff <- end - start
	results <- list(apTarget=apTarget, awStats=awStats, k=k, n1=n1, n0=n0, time=timeDiff)
	filename <- paste0("awStat_k",k,"_apTarget_",signif(apTarget[1],1),".rdata")
	if(apTarget>=0.01)
		filename <- paste0("awStat_k",k,"_apTarget_",apTarget[1],".rdata")
	save(results, file=filename)	
}

stopCluster(cl)


#g++ -I/usr/include/R -DNDEBUG   -std=c++11   -fpic  -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c importanceSampling.cpp -o importanceSampling.o
#R CMD SHLIB importanceSampling.cpp
