# NxtIRFcore: Command line NxtIRF for Bioconductor

NxtIRF quantifies Intron Retention and Alternative Splicing from BAM files using the IRFinder engine. Features interactive visualisation including RNA-seq coverage plots normalised by condition at the splice junction level.

NxtIRF release version: http://www.bioconductor.org/packages/release/bioc/html/NxtIRFcore.html

NxtIRF devel version: http://www.bioconductor.org/packages/devel/bioc/html/NxtIRFcore.html

## Installation

### On current R (>= 4.2.0)

* Requires Bioconductor devel version:

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::valid()              # checks for out of date packages
```

* Development version from Github (NxtIRFcore and its dependency NxtIRFdata):
```
library("devtools")
install_github("alexchwong/NxtIRFdata")
install_github("alexchwong/NxtIRFcore", dependencies=TRUE, build_vignettes=TRUE)
```

## Vignettes

```
browseVignettes("NxtIRFcore")
```
