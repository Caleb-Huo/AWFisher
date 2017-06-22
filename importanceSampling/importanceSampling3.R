WD <- "/mnt/glusterfs/zhh18/AW/AWFisher/importanceSampling"

setwd(WD)



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
	n1 <- 1e8
	k <- 2
}


pTargetList0 <- c(seq(0.99,0.01, by=-0.01))
pTargetList1 <- c(0.1^{2:100})
pTargetList <- c(pTargetList0, pTargetList1)

awStats0 <- numeric(length(pTargetList0))
awStats1 <- numeric(length(pTargetList1))

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

	
## g++ -I/usr/lib64/R/include -DNDEBUG   -std=c++0x   -fpic  -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c importanceSampling.cpp -o importanceSampling.o
## R CMD SHLIB importanceSampling.cpp


dyn.load('importanceSampling.so')


start <- Sys.time()
cat("k =", k)

awStats0 <- .C('importanceSampling_R',a=as.double(1),n=as.integer(n1),k=as.integer(k),
			pTarget=as.double(pTargetList0),J=as.integer(length(pTargetList0)),awStat_vec=as.double(pTargetList0))$awStat_vec
#
end <- Sys.time()
timeDiff <- end - start

timeDiff


for(i in 1:length(pTargetList1)){
  cat(".")
  apTarget <- pTargetList1[i]
  a <- uniroot(function(a, n, k, target) medianReturn(a, n, k) - -log(target), 
          lower=1e-10,upper=1, n=n0, k = k, target = apTarget)$root

  awStats1[i] <-  .C('importanceSampling_R',a=as.double(a),n=as.integer(n1),k=as.integer(k),
			pTarget=as.double(apTarget),J=as.integer(length(apTarget)),awStat_vec=as.double(apTarget))$awStat_vec
  
}
awStats <- c(awStats0, awStats1)

cat("\n")
end <- Sys.time()
timeDiff <- end - start


results <- list(pTargetList=pTargetList, awStats=awStats, k=k, n1=n1, n0=n0, time=timeDiff)

fileName <- paste0("awStat_k",k,".rdata")

save(results, file=fileName)



if(F){
	WD <- "/mnt/glusterfs/zhh18/AW/AWFisher/importanceSampling"
	setwd(WD)

	n0 <- 1e3
	n1 <- 1e6
	k <- 2

	pTargetList0 <- c(seq(0.99,0.01, by=-0.01))
	pTargetList1 <- c(0.1^{2:100})
	pTargetList <- c(pTargetList0, pTargetList1)

	awStats0 <- numeric(length(pTargetList0))
	awStats1 <- numeric(length(pTargetList1))
	
	dyn.load('importanceSampling.so')

	awStats0 <- .C('importanceSampling_R',a=as.double(1.0),n=as.integer(1000000),k=as.integer(k),
				pTarget=as.double(pTargetList0[6]),J=as.integer(length(pTargetList0[6])),awStat_vec=as.double(pTargetList0[6]))$awStat_vec

	
}



#g++ -I/usr/include/R -DNDEBUG   -std=c++11   -fpic  -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c importanceSampling.cpp -o importanceSampling.o
#R CMD SHLIB importanceSampling.cpp
