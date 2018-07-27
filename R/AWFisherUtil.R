aw.fisher.stat <- function(pstat, n) {
    numN = NCOL(sysdata$original)
    quant = rep(0, numN)
    for(i in 1:numN) {
        f = splinefun(c(1, sysdata[["nList"]]), c(sysdata[["logPTarget"]][i], sysdata$original[,i]))
        quant[i] = f(n)
    }

    ##### Estimating ###########
    f = splinefun(c(0,quant), c(0,sysdata[["logPTarget"]]), method="monoH.FC")
    exp(-f(pstat))
}

