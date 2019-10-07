# scCancer

The `scCancer` package focuses on processing and analyzing droplet-based scRNA-seq data for cancer research. Except routine data processing steps, this package takes several special considerations for cancer-specific features. 

The `scCancer` pipeline consists of two parts: `scStatistics` and `scAnnotation`.
* The `scStatistics` performs basic statistical analysis of raw data and quality control.
* The `scAnnotation` performs functional data analyses and visualizations, such as low dimensional representation, clustering, cell type classification, malignancy estimation, cellular phenotype scoring, gene signature analysis, etc.

After these analyses, user-friendly graphic reports will be generated.

## Systems Requirement



## Installation

The `scCancer` package can be installed via

```R
library(devtools)
devtools::install_github("wguo-research/scCancer")
```

**Hint:**  If you haven't installed the package `devtools`, please run `install.packages("devtools")` firstly.

## Usage

The vignette of `scCancer` can be found in the project [wiki](https://github.com/wguo-research/scCancer/wiki/scCancer-vignettes).


## Citation
Please use the following citation:

