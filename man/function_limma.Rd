% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_limma.r
\name{function_limma}
\alias{function_limma}
\title{use limma function to get pvalue}
\usage{
function_limma(astudy)
}
\arguments{
\item{astudy}{A list contains a data matrix and a vector of group label}
}
\value{
A list of pvalue and effect size
}
\description{
use limma function to get pvalue
}
\details{
use limma function to get pvalue
}
\examples{

N0 = 10
G <- 1000
GDEp <- 50
GDEn <- 50

set.seed(15213)

astudy <- matrix(rnorm(N0*2*G),nrow=G,ncol=N0*2)
ControlLabel <- 1:N0
caseLabel <- (N0 + 1):(2*N0)

astudy[1:GDEp,caseLabel] <- astudy[1:GDEp,caseLabel] + 2
astudy[1:GDEp + GDEn,caseLabel] <- astudy[1:GDEp,caseLabel] - 2

alabel <- c(rep(0,length(ControlLabel)),rep(1,length(caseLabel)))
Study <- list(data=astudy, label=alabel)

result <- function_limma(Study)
fdr <- p.adjust(result$pvalue)
sum(fdr<=0.05)

}
\author{
Zhiguang Huo
}
