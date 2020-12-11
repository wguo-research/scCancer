# scCancer

## Introduction

The `scCancer` package focuses on processing and analyzing droplet-based scRNA-seq data for cancer research. Except basic data processing steps, this package takes several special considerations for cancer-specific features.

The workflow of  `scCancer` mainly consists of three modules: `scStatistics`, `scAnnotation`, and `scCombination`.
* The `scStatistics` performs basic statistical analyses of raw data and quality control.
* The `scAnnotation` performs functional data analyses and visualizations, such as low dimensional representation, clustering, cell type classification, cell malignancy estimation, cellular phenotype analyses, gene signature analyses, cell-cell interaction analyses, etc.
* The `scCombination` perform multiple samples data integration, batch effect correction and analyses visualization.

After the computational analyses, detailed and graphical reports were generated in user-friendly HTML format.

<img src="http://lifeome.net/software/sccancer/scCancer-workflow.png" width="70%" alt="scCancer-workflow" align=center>

([Click to view larger workflow picture](http://lifeome.net/software/sccancer/scCancer-workflow.png))


## System Requirements
* R version: >= 3.5.0 (**suggest** to use R 3.5 or 3.6, **not 4.0**)
* **! ! ! Hint:  if your R version is >=4.0 under Windows system**, you may need update your Rtools to 4.0 from https://cran.r-project.org/bin/windows/Rtools/. So currently, we **don't** suggest to use R (>=4.0), if you are not familiar with R environment configuration.

## Current version

* scCancer 2.2.0 (update at 2020.12.10)
* (All version log](https://github.com/wguo-research/scCancer/wiki/Version-Log)

## Installation

The detailed installation instruction can be found in the project [wiki]( https://github.com/wguo-research/scCancer/wiki/2.-Installation).


## Usage

The vignette of `scCancer` can be found in the project [wiki]( https://github.com/wguo-research/scCancer/wiki).

* [Quick start](https://github.com/wguo-research/scCancer/wiki/3.-Quick-start)
* [Step by step introduction](https://github.com/wguo-research/scCancer/wiki/4.-Step-by-step-introduction)
* [Other personalized settings](https://github.com/wguo-research/scCancer/wiki/5.-Other-personalized-settings)

We provide an [example data](http://lifeome.net/software/sccancer/KC-example.tar.gz) of kidney cancer from 10X Genomics, and following are the generated HTML reports:

* [`report-scStat.html`](http://lifeome.net/software/sccancer/KC-example-report-scStat.html)
* [`report-scAnno.html`](http://lifeome.net/software/sccancer/KC-example-report-scAnno.html)

For multi-datasets, following is a generated HTML report for three kidney cancer samples integration analysis:

* [`report-scAnnoComb.html`](http://lifeome.net/software/sccancer/KC123-report-scAnnoComb.html)


## Citation
Please use the following citation:

Wenbo Guo, Dongfang Wang, Shicheng Wang, Yiran Shan, Changyi Liu, Jin Gu, scCancer: a package for automated processing of single-cell RNA-seq data in cancer, _Briefings in Bioinformatics_,  bbaa127, [https://doi.org/10.1093/bib/bbaa127](https://doi.org/10.1093/bib/bbaa127)

## License
GPL-3
