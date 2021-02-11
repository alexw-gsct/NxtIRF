# NxtIRF
NxtIRF quantifies Intron Retention and Alternative Splicing from BAM files using the IRFinder engine. Features interactive visualisation including RNA-seq coverage plots normalised by condition at the splice junction level.

## Installation

### On current R (>= 4.0.0)
* Development version from Github:
```
library("devtools")
install_github("alexw-gsct/NxtIRF", ref = "Thesis_branch", dependencies=TRUE)
```

## Vignettes

```
browseVignettes("NxtIRF")
```

## Saving interactive plots from within the NxtIRF shinydashboard app
Please note that plotly/orca needs to be installed to use this functionality. See https://github.com/plotly/orca#installation for detailed instructions on how to install plotly/orca on your system.
