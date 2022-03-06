# NxtIRFcore: Command line NxtIRF for Bioconductor

NxtIRF quantifies Intron Retention and Alternative Splicing from BAM files using the IRFinder engine. Features interactive visualisation including RNA-seq coverage plots normalised by condition at the splice junction level.

NxtIRF release version: http://www.bioconductor.org/packages/release/bioc/html/NxtIRFcore.html

NxtIRF devel version: http://www.bioconductor.org/packages/devel/bioc/html/NxtIRFcore.html

## Installation

### On devel R (>= 4.2.0)

* Requires Bioconductor devel version:

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::valid()              # checks for out of date packages

BiocManager::install("NxtIRFcore")
```

* Development version from Github (NxtIRFcore and its dependency NxtIRFdata):
```
library("devtools")
install_github("alexchwong/NxtIRFdata")
install_github("alexchwong/NxtIRFcore", dependencies=TRUE, build_vignettes=TRUE)
```

### On current version of R == 4.1.x

* Requires Bioconductor release version:

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::valid()              # checks for out of date packages

BiocManager::install("NxtIRFcore")
```

### On R == 4.0.x

* As NxtIRFcore's vignette relies on a demo dataset that is deposited in AnnotationHub for Bioconductor 3.14, the vignette cannot be installed in R 4.0. The following should still work for Bioconductor 3.12 or 3.13. A warning that NxtIRF on old versions of Bioconductor have not been extensively tested, and unknown bugs may exist.

```
library("devtools")
install_github("alexchwong/NxtIRFdata")
install_github("alexchwong/NxtIRFcore", dependencies=TRUE, build_vignettes=FALSE)
```

## Vignettes

* From within R, after installing NxtIRFcore:

```
browseVignettes("NxtIRFcore")
```

* You can also view the online version at Bioconductor here: http://www.bioconductor.org/packages/release/bioc/vignettes/NxtIRFcore/inst/doc/NxtIRF.html
