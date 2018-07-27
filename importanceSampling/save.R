#WD <- "/mnt/glusterfs/zhh18/AW/AWFisher/importanceSampling_tmp"
WD <- "/data/home/zhuo/research/Tseng/AW/sysdata"

kRange <- c(seq(2,117,1), 120, 180, 300, 500, 1000)

logPTarget <- NULL
original <- NULL
time <- NULL

n1 <- NULL
n0 <- NULL

for(k in kRange){
	awStatFolder_k <- paste0('k',k)
	setwd(WD)
	setwd(awStatFolder_k)
	awStatFiles_p <- dir(pattern="awStat_")

	aPTarget <- NULL
	aAWstat <- NULL
	time0 <- Sys.time()
	timeDiff <- time0 - time0
	
	for(aawStatFiles_p in awStatFiles_p){
		results <- get(load(aawStatFiles_p))
		aPTarget <- c(aPTarget, results$apTarget)
		aAWstat <- c(aAWstat, results$awStats)
		timeDiff <- timeDiff + results$time
		if(is.null(n1)){
			n1 <- results$n1
		} else {
			stopifnot(all(n1 == results$n1))			
		}
		if(is.null(n0)){
			n0 <- results$n0
		} else {
			stopifnot(all(n0 == results$n0))			
		}
	}
	
	aP_order <- order(-aPTarget)
	bPTarget <- aPTarget[aP_order]
	bAWstat <- aAWstat[aP_order]
	
	if(is.null(logPTarget)){
		logPTarget <- -log(bPTarget)
	} else {
		stopifnot(all(logPTarget == -log(bPTarget)))
	}
	
	original <- rbind(original, bAWstat)
	time <- c(time, as.double(timeDiff,units="hours"))
	
}

sum(time)

rownames(original) <- kRange

setwd(WD)
pdf("timeVsK.pdf")
plot(kRange, time)
dev.off()

sysdata = list(logPTarget=logPTarget, original=original, nList=kRange, n1=n1, n0=n0)
save(sysdata, file ="sysdata.rda", compress="xz")

