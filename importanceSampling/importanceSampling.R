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
	n1 <- 1e3 
	k <- 2
}


pTargetList <- c(seq(0.95,0.05, by=-0.05), 0.1^{2:10}, 1e-12,1e-15,1e-20,1e-30,1e-50,1e-80,1e-200)

awStats <- numeric(length(pTargetList))

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

importantSampling <- function(a, n, k, pTarget, seed = 15213) {
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
  statOrdered <- statNew[ordStat]
  
  ## may be wrong
  statOrdered[which(logProbNew <= -log(pTarget))[1]]
}

start <- Sys.time()
cat("k =", k)
for(i in 1:length(pTargetList)){
  cat(".")
  apTarget <- pTargetList[i]
  if(apTarget >= 0.1){
    a <- 1
  } else {
    a <- uniroot(function(a, n, k, target) medianReturn(a, n, k) - -log(target), 
            lower=1e-10,upper=1, n=n, k = k, target = apTarget)$root
  }
  awStats[i] <- importantSampling(a, n=n1, k=k, pTarget= apTarget)
}
cat("\n")
end <- Sys.time()
timeDiff <- end - start


results <- list(pTargetList=pTargetList, awStats=awStats, k=k, n1=n1, n0=n0, time=timeDiff)

fileName <- paste0("awStat_k",k,".rdata")

save(results, file=fileName)

