# scCancer

## Introduction

The `scCancer` package focuses on processing and analyzing droplet-based scRNA-seq data for cancer research. Except basic data processing steps, this package takes several special considerations for cancer-specific features.

The workflow of  `scCancer` mainly consists of two parts: `scStatistics` and `scAnnotation`.
* The `scStatistics` performs basic statistical analysis of raw data and quality control.
* The `scAnnotation` performs functional data analyses and visualizations, such as low dimensional representation, clustering, cell type classification, malignancy estimation, cellular phenotype scoring, gene signature analysis, etc.

After these analyses, user-friendly graphic reports will be generated.

<img src="http://lifeome.net/software/sccancer/scCancer-workflow.png" width="70%" alt="scCancer-workflow" align=center>

## System Requirements

* **Memery**:  >= 32G  (for a data with ~10000 cells)
* **R version**: >= 3.5.0 

## Installation

Firstly, please install or update the package `devtools` by running 
```R
install.packages("devtools")
```

Then the `scCancer` can be installed via

```R
library(devtools)
devtools::install_github("wguo-research/scCancer")
```

## Usage

The vignette of `scCancer` can be found in the project [wiki]( https://github.com/wguo-research/scCancer/wiki) or vignette [page](http://lifeome.net/software/sccancer/scCancer-vignette.html).

We also provide an [example data](http://lifeome.net/software/sccancer/KC-example.tar.gz) of kidney cancer from 10X Genomics, and here are the generated reports: 

* [`report-scStat.html`](http://lifeome.net/software/sccancer/KC-example-results/report-scStat.html)
* [`report-scAnno.html`](http://lifeome.net/software/sccancer/KC-example-results/report-scAnno.html)


## Citation
Please use the following citation:

Wenbo Guo, Dongfang Wang, Shicheng Wang, Yiran Shan, Jin Gu. 2019. bioRxiv doi: https://doi.org/10.1101/800490

## License
GPL-3