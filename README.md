# scCancer

The `scCancer` package focuses on processing and analyzing droplet-based scRNA-seq data for cancer research. Except routine data processing steps, this package takes several special considerations for cancer-specific features. 

The `scCancer` pipeline consists of two parts: `scStatistics` and `scAnnotation`.
* The `scStatistics` implements quality control for cells and genes.
* The `scAnnotation` implements normal analysis and annotate the cells' low-dimension coordinates, cluster, cell type, cell cycle, malignancy, gene set signatures, expression program, and so on. 

After these analyses, user-friendly graphic reports will be generated.


## Installation

The `scCancer` package can be installed via

    library(devtools)
    devtools::install_github("wguo-research/scCancer")


## Usage

The vignette of `scCancer` can be found in the project [wiki](https://github.com/wguo-research/scCancer/wiki/scCancer-vignettes).


## Citation
Please use the following citation:

