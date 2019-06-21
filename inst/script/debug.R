rm(list=ls())

WD <- "/Users/caleb/Box Sync/Shaowu/AW/AWFisher/importanceSampling"

compile <- "g++ -I/Library/Frameworks/R.framework/Resources/include  -DNDEBUG   -std=c++0x  -fpic  -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security -D_FORTIFY_SOURCE=2 -g  -c importanceSampling.cpp -o importanceSampling.o"
build <- "R CMD SHLIB importanceSampling.cpp"

setwd(WD)
system(compile)
system(build)

dyn.load('importanceSampling.so')

n0 <- 1e3
n1 <- 1e6
k <- 2

pTargetList0 <- c(seq(0.99,0.01, by=-0.01))
pTargetList1 <- c(0.1^{2:100})
pTargetList <- c(pTargetList0, pTargetList1)

awStats0 <- numeric(length(pTargetList0))
awStats1 <- numeric(length(pTargetList1))

dyn.load('importanceSampling.so')

awStats0 <- .C('importanceSampling_R',a=as.double(1),n=as.integer(1000),k=as.integer(k),
               pTarget=as.double(pTargetList0[6]),J=as.integer(length(pTargetList0[6])),awStat_vec=as.double(pTargetList0[6]))$awStat_vec


a <- as.double(1)
n <- as.integer(10000000)
k <- as.integer(2)
pTarget <- as.double(c(seq(0.99, 0.90, -0.01)))
J <- as.integer(length(pTarget))
awStat_vec <- pTarget

obj <- .C('importanceSampling_R',a=a,n=n,k=k,pTarget=pTarget,J=J,awStat_vec=awStat_vec)
obj
