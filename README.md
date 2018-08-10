# AWFisher
Github repository for adaptively weighted fisher's method (AW-Fisher)


## Install This Package from github
First you need R `devtools` package installed.

* In R console

```{R}
library(devtools)
install_github("Caleb-Huo/AWFisher") 
```

## Short tutorial

```{R}
library(AWFisher)

K <- 50
G <- 10000

set.seed(15213)
p.values = matrix(runif(K*G), ncol=K)
res = AWFisher.pvalue(p.values)

hist(res$pvalues, breaks=40)

ks<-ks.test(res$pvalues, "punif", min=0, max=1, alternative = "two.sided"); ## KS test to test if the AW p-values are uniformly distributed under the null
ks

```


