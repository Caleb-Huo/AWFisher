# AWFisher
Github repository for adaptively weighted fisher's method (AW-Fisher)


## Install This Package from github
* In R console

```{R}
library(devtools)
install_github("Caleb-Huo/AWFisher") 
```

## Citation

* P-value evaluation, variability index and biomarker categorization for adaptively weighted Fisher's meta-analysis method in omics applications. Zhiguang Huo, Shaowu Tang, Yongseok Park, George Tseng. (2018+, under review)

* Arxiv preprint can be find here: https://arxiv.org/abs/1708.05084

## Short tutorial

* Currently, only K=2, 3, ..., 100 (number of studies) are allowed in the R package.

```{R}
library(AWFisher)

K <- 50 ## combining K studies
G <- 10000 ## simulate G genes

set.seed(15213)
p.values = matrix(runif(K*G), ncol=K)
res = AWFisher.pvalue(p.values)

hist(res$pvalues, breaks=40)

ks<-ks.test(res$pvalues, "punif", min=0, max=1, alternative = "two.sided"); ## KS test to test if the AW p-values are uniformly distributed under the null
ks

```


