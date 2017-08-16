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
  resultStat <- numeric(length(pTarget))
  for(i in 1:length(pTarget)){
  	resultStat[i] <- statOrdered[which(logProbNew <= -log(pTarget[i]))[1]]
  }
  resultStat
}

start <- Sys.time()
cat("k =", k)

awStats0 <- importantSampling(a=1, n=n1, k=k, pTarget= pTargetList0)


for(i in 1:length(pTargetList1)){
  cat(".")
  apTarget <- pTargetList1[i]
  if(apTarget >= 0.01){
    a <- 1
  } else {
    a <- uniroot(function(a, n, k, target) medianReturn(a, n, k) - -log(target), 
            lower=1e-10,upper=1, n=n0, k = k, target = apTarget)$root
  }
  awStats1[i] <- importantSampling(a, n=n1, k=k, pTarget= apTarget)
}
awStats <- c(awStats0, awStats1)

cat("\n")
end <- Sys.time()
timeDiff <- end - start


results <- list(pTargetList=pTargetList, awStats=awStats, k=k, n1=n1, n0=n0, time=timeDiff)

fileName <- paste0("awStat_k",k,".rdata")

save(results, file=fileName)

