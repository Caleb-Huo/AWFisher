aw.fisher.stat <- function(pstat, n) {
    index = match(n, sysdata[["nList"]])
    
    ##### Estimating ###########
    f = splinefun(c(0, sysdata$original[index, ]), 
        c(0, sysdata[["logPTarget"]]), method = "monoH.FC")
    exp(-f(pstat))
}
