# AWFisher
Github repository for adaptively weighted fisher's method (AW-Fisher)


## Install This Package from github
* In R console

```{R}
library(devtools)
install_github("Caleb-Huo/AWFisher") 
```

## Citation

* Huo, Z., Tang, S., Park, Y. and Tseng, G., 2020. P-value evaluation, variability index and biomarker categorization for adaptively weighted Fisher’s meta-analysis method in omics applications. *Bioinformatics*, 36(2), pp.524-532.

* The manuscript can be found here: https://www.ncbi.nlm.nih.gov/pubmed/31359040

## Full tutorial

* Including a real data example using mouse metabolism data of three tissues
* Perform transcriptomic meta-analysis and differential expression pattern detection

http://htmlpreview.github.io/?https://github.com/Caleb-Huo/AWFisher/blob/master/vignettes/AWFisher.html


## Short tutorial

* This short tutorial is about how to perform meta-analysis combining p-values from multiple studies.
* Currently, only K=2, 3, ..., 100 (number of studies) are allowed in the R package.

```{R}
library(AWFisher)

K <- 50 ## combining K studies
G <- 10000 ## simulate G genes

set.seed(15213)
p.values = matrix(runif(K*G), ncol=K)
res = AWFisher_pvalue(p.values)

hist(res$pvalues, breaks=40)

ks<-ks.test(res$pvalues, "punif", min=0, max=1, alternative = "two.sided"); ## KS test to test if the AW p-values are uniformly distributed under the null
ks

```


