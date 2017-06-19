WD <- "/mnt/glusterfs/zhh18/AW/AWFisher/importanceSampling"

setwd(WD)

awStatfiles0 <- dir(pattern="awStat_")
awStatfiles1 <- gsub("awStat_k","",awStatfiles0)
awStatfiles2 <- gsub(".rdata","",awStatfiles1)

ks <- sort(as.numeric(awStatfiles2))

logPTarget <- NULL
original <- NULL
time <- NULL

for(i in 1:length(ks)){
	k <- ks[i]
	fileName <- paste0("awStat_k",k,".rdata")
	results <- get(load(fileName))
	if(is.null(logPTarget)){
		logPTarget <- -log(results$pTargetList)
	} else {
		stopifnot(all(logPTarget == -log(results$pTargetList)))
	}
	original <- rbind(original, results$awStats)
	time <- c(time, results$time)
}

pdf("timeVsK.pdf")
plot(ks, time)
dev.off()

sysdata = list(logPTarget=logPTarget, original=original, nList=ks)
save(sysdata, file ="sysdata.rda", compress="xz")

