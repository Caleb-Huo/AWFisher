# AWFisher
Github repository for adaptively weighted fisher's method (AW-Fisher)


## Required Package
* If for AW-Fisher p-value calculation, not declared yet
* If for metaPattern, limma, edgeR, gplots, tightClust

## Install This Package from github
First you need R `devtools` package installed.
* In command line:
```

* In R console
```R
library(devtools)
install_github("Caleb-Huo/AWFisher", build_vignettes=TRUE) ## smaller package, only for AW-Fisher p-value calculation
install_github("Caleb-Huo/AWFisher", ref="metaPattern", build_vignettes=TRUE) ## smaller package, only for AW-Fisher p-value calculation
```


* Or install from a released package.
    - First, download the latest released package 
[https://github.com/Caleb-Huo/MetaSparseKmeans/releases](https://github.com/Caleb-Huo/MetaSparseKmeans/releases)
    - Install this downloaded package directly.
```R
install.packages("MetaSparseKmeans_0.0.3.tar.gz",repos=NULL,type="source")
```

## How to use this R package:

* After installing this package from github, In R console:
```R
browseVignettes("AWFisher")
```