
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scflow

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/combiz/scflow.svg?branch=master)](https://travis-ci.org/combiz/scflow)
[![Codecov test
coverage](https://codecov.io/gh/combiz/scflow/branch/master/graph/badge.svg)](https://codecov.io/gh/combiz/scflow?branch=master)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/combiz/scflow?branch=master&svg=true)](https://ci.appveyor.com/project/combiz/scflow)
<!-- badges: end -->

The goal of scflow is to provide tools in R to build a complete analysis
workflow for single-cell/nuclei RNA sequencing data.

  - Quality control of gene-cell matrices
      - Filtering of matrices by counts and features
      - Filtering of mitochondrial genes and mitochondrial counts
        thresholding
      - Doublet and multiplet identification and removal with
        [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
      - Rich QC metrics annotation with
        [scater](https://github.com/davismcc/scater)
  - Dimensionality reduction and celltype identification
      - Louvain clustering, UMAP dimensionality reduction, and cluster
        marker gene identification with
        [monocle](https://github.com/cole-trapnell-lab/monocle-release)
      - Celltype annotation with
        [EWCE](https://github.com/NathanSkene/EWCE) and
        [Liger](https://github.com/MacoskoLab/liger)
      - Cluster identity mapping against the [Allen Human Brain
        Atlas](https://www.brain-map.org) and Mouse Nervous System Data
        from [Zeisel 2018](https://doi.org/10.1016/j.cell.2018.06.021)
  - Differential gene expression implementations
      - Zero-inflated regression model with
        [MAST](https://github.com/RGLab/MAST)
      - Random effects model with [Limma](https://github.com/cran/limma)
      - Negative binomial distribution pseudobulking model with
        [DESeq2](https://github.com/mikelove/DESeq2)
      - Pseudobulk generalized likelihood ratio tests with
        [EdgeR](https://github.com/StoreyLab/edge)
  - Pathway and functional category enrichment analysis
      - Interface to the Enrichr database with
        [EnrichR](https://github.com/cran/enrichR)
      - Interface to the WebGestalt tool with
        [WebGestaltR](http://www.webgestalt.org/)
  - Publication quality plots and analysis reports
      - QC plots and tabular metrics suitable for reports.
      - UMAP plots for cell features and gene expression.
      - Violin plots for gene expression.
      - Pathway and gene enrichment plots

The package functions are designed to interface neatly with
[NextFlow](https://www.nextflow.io/) for scalable and containerized
pipelines deployed locally, on high-performance computing clusters, or
in the cloud. An accompanying NextFlow pipeline is in the works - TBA.

## Installation

You can install the development version of scflow from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("combiz/scflow")
```

## Example

Remove this . This is a basic example which shows you how to solve a
common problem:

``` r
library(scflow)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
