# NxtIRFcore: Command line NxtIRF for Bioconductor

NxtIRF quantifies Intron Retention and Alternative Splicing from BAM files using the IRFinder engine. Features interactive visualisation including RNA-seq coverage plots normalised by condition at the splice junction level.

## Installation

### On current R (>= 4.0.0)
* Development version from Github:
```
library("devtools")
install_github("alexchwong/NxtIRFcore", dependencies=TRUE, build_vignettes=TRUE)
```

## Vignettes

```
browseVignettes("NxtIRFcore")
```
