# scCancer – A package for automated processing of single cell RNA-seq data in cancer

The cowplot package provides various features that help with creating publication-quality figures, such as a set of themes, functions to align plots and arrange them into complex compound figures, and functions that make it easy to annotate plots and or mix plots with images. The package was originally written for internal use in the Wilke lab, hence the name (Claus O. Wilke's plot package). It has also been used extensively in the book  [Fundamentals of Data Visualization.](https://www.amazon.com/gp/product/1492031089)

The `scCancer` package focuses on processing and analyzing droplet-based scRNA-seq data for cancer research. Except routine data processing steps, this package takes several special considerations for cancer-specific features. 

The `scCancer` pipeline consists of two parts: `scStatistics` and `scAnnotation`.
* The `scStatistics` implements quality control for cells and genes.
* The `scAnnotation` implements normal analysis and annotate the cells' low-dimension coordinates, cluster, cell type, cell cycle, malignancy, gene set signatures, expression program, and so on. 

After these analyses, a user-friendly graphic report will be generated.


# Installation

The scCancer package can be installed via

    library(devtools)
    devtools::install_github("wguo-research/scCancer")


# Usage

The vignette of scCancer can be found here.


# Citation
Please use the following citation:

