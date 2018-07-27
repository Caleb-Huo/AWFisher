WD0 <- "/data/home/zhuo/research/Tseng/AW/sysdata"
WD <- "/data/home/zhuo/research/Tseng/AW/sysdata/data"
dir.create(WD)
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
	k <- 5
}


pTargetList0 <- c(seq(0.99,0.01, by=-0.01))
pTargetList1 <- c(0.1^{3:100})

pTargetList <- list(pTargetList0)
for(i in 1:length(pTargetList1)){
	pTargetList[[i+1]] <- pTargetList1[i]
}

set.seed(15213)
pmatrix0 <- matrix(runif(n0*k),n0,k)


medianReturn <- function(a, n, k, pmatrix0) {
  pmatrix <- pmatrix0^(1/a)
  
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

	
cl<-makeCluster(33)
registerDoParallel(cl)

result <- foreach(p = 1:length(pTargetList)) %dopar% {	
	apTarget <- pTargetList[[p]]

	setwd(WD0)
	dyn.load('importanceSampling.so')

	setwd(WD)
	kFolder <- paste0("k",k)
	system(paste("mkdir -p", kFolder))
	setwd(kFolder)
		
	
	start <- Sys.time()
	if(p==1){
		a <- 1
	} else {
	    a <- uniroot(function(a, n, k, target) medianReturn(a, n, k, pmatrix0) - -log(target), 
	            lower=1e-10,upper=1, n=n0, k = k, target = apTarget)$root		
	}
	awStats <- .C('importanceSampling_R',a=as.double(a),n=as.integer(n1),k=as.integer(k),
				pTarget=as.double(apTarget),J=as.integer(length(apTarget)),awStat_vec=as.double(apTarget))$awStat_vec				
	end <- Sys.time()
	timeDiff <- end - start
	results <- list(apTarget=apTarget, awStats=awStats, k=k, n1=n1, n0=n0, time=timeDiff)
	filename <- paste0("awStat_k",k,"_apTarget_",signif(apTarget[1],1),".rdata")
	save(results, file=filename)	
}

stopCluster(cl)



#g++ -I/usr/include/R -DNDEBUG   -std=c++11   -fpic  -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c importanceSampling_debugChuck.cpp -o importanceSampling_debugChuck.o
#R CMD SHLIB importanceSampling_debugChuck.cpp



# cd /data/home/zhuo/research/Tseng/AW/sysdata
# g++ -I/usr/share/R/include -DNDEBUG   -std=c++11   -fpic  -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c importanceSampling.cpp -o importanceSampling.o
# R CMD SHLIB importanceSampling.cpp

