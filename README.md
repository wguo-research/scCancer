# scCancer

## Introduction

The `scCancer` package focuses on processing and analyzing droplet-based scRNA-seq data for cancer research. Except basic data processing steps, this package takes several special considerations for cancer-specific features.

The workflow of  `scCancer` mainly consists of two parts: `scStatistics` and `scAnnotation`.
* The `scStatistics` performs basic statistical analysis of raw data and quality control.
* The `scAnnotation` performs functional data analyses and visualizations, such as low dimensional representation, clustering, cell type classification, malignancy estimation, cellular phenotype scoring, gene signature analysis, etc.

After these analyses, user-friendly graphic reports will be generated.

![scCancer-workflow](http://lifeome.net/software/sccancer/scCancer-workflow.png)

## System Requirements

* **Memery**:  >= 16G  (for a data with ~10000 cells)
* **R version**: >= 3.4.0 

## Installation

Please install or update the package `devtools` by running `install.packages("devtools")` firstly. 

Then the `scCancer` can be installed via

```R
library(devtools)
devtools::install_github("wguo-research/scCancer")
```

## Usage

The vignette of `scCancer` can be found in the project [wiki](https://github.com/wguo-research/scCancer/wiki/scCancer-vignettes).

We also provide a [kidney cancer example data](http://lifeome.net/software/sccancer/KC-example.zip), and here are the generated reports: [`report-scStat.html`](http://lifeome.net/software/sccancer/KC-example-results/report-scStat.html), [`report-scAnno.html`](http://lifeome.net/software/sccancer/KC-example-results/report-scAnno.html).


## Citation
Please use the following citation:

## License

GPL-3   &copy; [G-Lab](http://lifeome.net/glab/jgu/), [Tsinghua University](http://www.tsinghua.edu.cn)