WD0 <- "/data/home/zhuo/research/Tseng/AW/sysdata/"
WD <- "/data/home/zhuo/research/Tseng/AW/sysdata/nn"
dir.create(WD)
setwd(WD)

library(foreach)
library(doParallel)

sysdata <- get(load("sysdata.rda"))
rownames(sysdata$original)[sysdata$original[,1] == 0.99]


if(T){
	n0 <- sysdata$n0
	n1 <- sysdata$n1
	n1 <- 1e6
	ks <- as.numeric(rownames(sysdata$original)[sysdata$original[,1] == 0.99])
}


pTargetList0 <- c(seq(0.99,0.01, by=-0.01))


	
cl<-makeCluster(length(ks))
registerDoParallel(cl)

result <- foreach(k = ks) %dopar% {	
	apTarget <- pTargetList0

	set.seed(15213)
	pmatrix0 <- matrix(runif(n0*k),n0,k)

	setwd(WD0)
	dyn.load('importanceSampling.so')

	setwd(WD)
	kFolder <- paste0("k",k)
	system(paste("mkdir -p", kFolder))
	setwd(kFolder)
	filename <- paste0("awStat_k",k,"_apTarget_",signif(apTarget[1],1),".rdata")

	start <- Sys.time()
	a <- 1
	awStats <- .C('importanceSampling_R',a=as.double(a),n=as.integer(n1),k=as.integer(k),
				pTarget=as.double(apTarget),J=as.integer(length(apTarget)),awStat_vec=as.double(apTarget))$awStat_vec				
	end <- Sys.time()
	timeDiff <- end - start
	results <- list(apTarget=apTarget, awStats=awStats, k=k, n1=n1, n0=n0, time=timeDiff)
	save(results, file=filename)
}

stopCluster(cl)



#g++ -I/usr/include/R -DNDEBUG   -std=c++11   -fpic  -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c importanceSampling_debugChuck.cpp -o importanceSampling_debugChuck.o
#R CMD SHLIB importanceSampling_debugChuck.cpp



# cd /data/home/zhuo/research/Tseng/AW/sysdata
# g++ -I/usr/share/R/include -DNDEBUG   -std=c++11   -fpic  -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c importanceSampling.cpp -o importanceSampling.o
# R CMD SHLIB importanceSampling.cpp

