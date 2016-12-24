##' This function could perform the tight clustering algorithm using dissimilarity matrix as input.
##'
##' Tight clustering method is a resampling-evaluated clustering method that aims to directly identify tight clusters in a high-dimensional complex data set and allow a set of scattered objects without being clustered. The method was originally developed for gene cluster analysis in microarray data but can be applied in any complex data. The most important parameter is k.min. A large k.min results in smaller and tighter clusters. Normally k.min>=target+5 is suggested. All other parameters do not affect the quality of final clustering results too much and are suggested to remain unchanged. This packages is directly inherited from R tightClust package
##' @title Tight Clustering
##' @param adist Input dissimilarity matrix.
##' @param target The total number of clusters that the user aims to find.
##' @param k.min The starting point of k0. See 'Details' for more information.
##' @param alpha The threshold of comembership index. Default value is suggested to be used.
##' @param beta The threshold of clusters stably found in consecutive k0. Default value is suggested to be used.
##' @param top.can The number of top (size) candidate clusters for a specific k0. Default value is suggested to be used.
##' @param seq.num The number of subsequent k0 that finds the tight cluster. Default value is suggested to be used.
##' @param resamp.num Total number of resampling to obtain comembership matrix. Default value is suggested to be used.
##' @param samp.p Percentage of subsamples selected. Default value is suggested to be used.
##' @param nstart Number of different random inital for K-means. Default value is suggested to be used.
##' @param remain.p Stop searching when the percentage of remaining points <= remain.p. Default value is suggested to be used.
##' @param k.stop Stop decreasing k0 when k0<=k.stop. Default value is suggested to be used.
##' @param random.seed If random.seed is NULL no random seed will be set. If random.seed is a number, it will be used as the random seed. This parameters should be used to get the same result for different runs.
##' @return Returned value is a "tight.clust" object (list). The first element is the original data matrix. The second element is a vector representing the cluster identity (-1: scattered gene set; 1: the first cluster; 2: the second cluster; ...). The third element is a vector of the size of each tight cluster.
##' @author Chi Song <song.1188@osu.edu>
##' @export
##' @examples
##' N0 = 10
##' G <- 1000
##' GDEp <- 50
##' GDEn <- 50
##' K = 4
##'
##' studies <- NULL
##' set.seed(15213)
##' for(k in 1:K){
##' 	astudy <- matrix(rnorm(N0*2*G),nrow=G,ncol=N0*2)
##' 	ControlLabel <- 1:N0
##' 	caseLabel <- (N0 + 1):(2*N0)
##'
##' 	astudy[1:GDEp,caseLabel] <- astudy[1:GDEp,caseLabel] + 2
##' 	astudy[1:GDEp + GDEn,caseLabel] <- astudy[1:GDEp + GDEn,caseLabel] - 2
##'
##' 	alabel = c(rep(0,length(ControlLabel)),rep(1,length(caseLabel)))
##'
##' 	studies[[k]] <- list(data=astudy, label=alabel)
##' }
##'
##'
##' result <- biomarkerCategorization(studies,function_limma,B=100,DEindex=NULL,seed = 15213)
##' dissimilarity <- result$dissimilarity
##' tightClustResult <- tightClustPam(dissimilarity, target=2, k.min=5)
tightClustPam <- function (adist, target, k.min, alpha = 0.1, beta = 0.6, top.can = 7,
    seq.num = 2, resamp.num = 10, samp.p = 0.7, nstart = 1, remain.p = 0.1,
    k.stop = 5, random.seed = NULL)
{
    find.candidates <- function(adist, k, alpha = 0.1, top.can = 7,
        resamp.num = 10, samp.p = 0.7, nstart = 1) {
        pam.classify <- function(adist, k, p = 0.7, ...) {
            size <- dim(adist)[1]
            partial.size <- round(p * size)
            partial <- sample(rownames(adist), partial.size)
            centers <-  pam(adist, k = k, diss=TRUE)$id.med
            return(suppressWarnings(pam(adist, k=k, medoids = centers, diss=TRUE)$clustering))
        }
        find.candidates.one <- function(x) {
            tmp <- apply(x == 1, 1, sum)
            return(which(x[, which(tmp == max(tmp))[1]] == 1))
        }
        extend.candidate <- function(D, can, alpha = 0.1) {
            can.ex <- which(apply(as.matrix(D[, can] >= 1 - alpha),
                1, all))
            D.temp <- D[can.ex, can.ex]
            if (!is.matrix(D.temp)) {
                D.temp <- as.matrix(D.temp)
                colnames(D.temp) <- names(can.ex)
            }
            D.bad <- apply(as.matrix(D.temp < 1 - alpha), 1,
                sum)
            while (sum(D.bad) > 0) {
                index <- which(D.bad == max(D.bad))[1]
                D.temp <- D.temp[-index, -index]
                D.bad <- apply(as.matrix(D.temp < 1 - alpha),
                  1, sum)
            }
            return(can.ex[colnames(D.temp)])
        }
        N <- dim(adist)[1]
        Dbar <- matrix(0, N, N)
        for (i in 1:resamp.num) {
            cl <- pam.classify(adist, k, p = samp.p)
            D <- outer(cl, cl, function(a, b) a == b)
            Dbar = Dbar + D
        }
        Dbar = Dbar/resamp.num
        colnames(Dbar) <- 1:N
        rownames(Dbar) <- 1:N
        i = 1
        D.temp <- Dbar
        res <- list()
        while (i <= top.can * 2 && dim(D.temp)[1] > 0) {
            candidate.one <- find.candidates.one(D.temp)
            candidate <- extend.candidate(D.temp, candidate.one,
                alpha = alpha)
            D.temp <- D.temp[-candidate, -candidate]
            res[[i]] <- names(candidate)
            mode(res[[i]]) <- "numeric"
            i = i + 1
        }
        res <- res[order(unlist(lapply(res, length)), decreasing = TRUE)][1:top.can]
        return(res)
    }
    if (!is.null(random.seed))
        set.seed(random.seed)
    original.data <- adist
    k.max <- k.min + 10
    id <- rownames(adist)
    N <- dim(adist)[1]
    write(paste("Number of points:", N, "\tDimension:", dim(adist)[2],
        "\n"), "")
    rownames(adist) <- 1:N
    index.m <- as.matrix(expand.grid(lapply(1:seq.num, function(x) 1:top.can)))
    remain <- N
    nfound <- 0
    found <- TRUE
    k0 <- k.min
    k <- k0
    candidates <- list()
    tclust <- list()
    while (nfound < target && remain/N >= remain.p && (found ||
        k <= k.max)) {
        if (found) {
            write(paste("Looking for tight cluster", nfound +
                1, "..."), "")
            k <- k0
            for (i in 1:seq.num) {
                write(paste("k =", k + i - 1), "")
                candidates[[i]] <- find.candidates(adist, k + i -
                  1, alpha = alpha, top.can = top.can, resamp.num = resamp.num,
                  samp.p = samp.p)
            }
        }
        else {
            candidates <- candidates[-1]
            candidates[[seq.num]] <- find.candidates(adist, k + seq.num -
                1, alpha = alpha, top.can = top.can, resamp.num = resamp.num,
                samp.p = samp.p, nstart = nstart)
        }
        calc.beta <- function(y) {
            temp <- lapply(1:seq.num, function(z) candidates[[z]][[y[z]]])
            i.temp <- temp[[1]]
#            u.temp <- temp[[i]]
			u.temp <- temp[[1]]
            for (j in 2:seq.num) {
                i.temp <- intersect(i.temp, temp[[j]])
                u.temp <- union(u.temp, temp[[j]])
            }
            return(length(i.temp)/length(u.temp))
        }
        beta.temp <- unlist(apply(index.m, 1, calc.beta))
        if (any(beta.temp >= beta)) {
            found = TRUE
            nfound = nfound + 1
            write(paste(nfound, "tight cluster(s) found!"), "")
            if (k0 > k.stop)
                k0 = k0 - 1
            found.temp <- candidates[[seq.num]][[index.m[which(beta.temp >=
                beta)[1], seq.num]]]
            tclust[[nfound]] <- rownames(adist)[found.temp]
            mode(tclust[[nfound]]) <- "numeric"
            adist <- adist[-found.temp, -found.temp]
            remain <- remain - length(tclust[[nfound]])
            write(paste("Cluster size:", length(tclust[[nfound]]),
                "\tRemaining number of points:", remain, "\n"),
                "")
        }
        else {
            found = FALSE
            k = k + 1
        }
    }
    clust.id <- rep(-1, N)
    size <- unlist(lapply(tclust, length))
    for (i in 1:length(tclust)) clust.id[tclust[[i]]] <- i
    res <- list(data = original.data, cluster = clust.id, size = size)
    class(res) <- "tight.clust"
    return(res)
}
