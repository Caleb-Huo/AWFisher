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
install_github("Caleb-Huo/AWFisher", ref="metaPattern", build_vignettes=TRUE) ## For metaPattern
```


* Or install from a released package.
    - First, download the latest released package (for p-value calculation)
[https://github.com/Caleb-Huo/AWFisher/releases/tag/p1.0.0](https://github.com/Caleb-Huo/AWFisher/releases/tag/p1.0.0)
    - Or, download the latest released package (for metaPattern)
[https://github.com/Caleb-Huo/AWFisher/releases/tag/m1.0.0](https://github.com/Caleb-Huo/AWFisher/releases/tag/m1.0.0)
    - Install this downloaded package directly.
```R
install.packages("AWFisher-p1.0.0.tar.gz",repos=NULL,type="source") ## smaller package, only for AW-Fisher p-value calculation
install.packages("AWFisher-m1.0.0.tar.gz",repos=NULL,type="source")  ## For metaPattern
```

## How to use this R package (Under development):

* After installing this package from github, In R console:
```R
browseVignettes("AWFisher")
```