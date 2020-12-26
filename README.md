<!-- README.md is generated from README.Rmd. Please edit that file -->
<img src=https://github.com/scFlow/blob/master/docs/logo/scflow.png width=300 />

scFlow
======

<!-- badges: start -->
[![Travis build
status](https://travis-ci.org/combiz/scFlow.svg?branch=master)](https://travis-ci.org/combiz/scFlow)
[![Codecov test
coverage](https://codecov.io/gh/combiz/scFlow/branch/master/graph/badge.svg)](https://codecov.io/gh/combiz/scFlow?branch=master)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/combiz/scFlow?branch=master&svg=true)](https://ci.appveyor.com/project/combiz/scFlow)
[![R build
status](https://github.com/combiz/scFlow/workflows/R-CMD-check/badge.svg)](https://github.com/combiz/scFlow/actions)
<!-- badges: end -->

The goal of scFlow is to provide tools in R to build a complete analysis
workflow for single-cell/nuclei RNA sequencing data.

-   Quality control of gene-cell matrices
    -   Filtering of matrices by counts and features
    -   Filtering of mitochondrial genes and mitochondrial counts
        thresholding
    -   Doublet and multiplet identification and removal with
        [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
    -   Rich QC metrics annotation with
        [scater](https://github.com/davismcc/scater)
-   Dimensionality reduction and celltype identification
    -   Louvain clustering, UMAP dimensionality reduction, and cluster
        marker gene identification with
        [monocle](https://github.com/cole-trapnell-lab/monocle-release)
    -   Celltype annotation with
        [EWCE](https://github.com/NathanSkene/EWCE) and
        [Liger](https://github.com/MacoskoLab/liger)
    -   Cluster identity mapping against the [Allen Human Brain
        Atlas](https://www.brain-map.org) and Mouse Nervous System Data
        from [Zeisel 2018](https://doi.org/10.1016/j.cell.2018.06.021)
-   Differential gene expression implementations
    -   Zero-inflated regression model with
        [MAST](https://github.com/RGLab/MAST)
    -   Random effects model with [Limma](https://github.com/cran/limma)
    -   Negative binomial distribution pseudobulking model with
        [DESeq2](https://github.com/mikelove/DESeq2)
    -   Pseudobulk generalized likelihood ratio tests with
        [EdgeR](https://github.com/StoreyLab/edge)
-   Pathway and functional category enrichment analysis
    -   Interface to the Enrichr database with
        [EnrichR](https://github.com/cran/enrichR)
    -   Interface to the WebGestalt tool with
        [WebGestaltR](http://www.webgestalt.org/)
-   Publication quality plots and analysis reports
    -   QC plots and tabular metrics suitable for reports.
    -   UMAP plots for cell features and gene expression.
    -   Violin plots for gene expression.
    -   Pathway and gene enrichment plots

The package functions are designed to interface neatly with
[NextFlow](https://www.nextflow.io/) for scalable and containerized
pipelines deployed locally, on high-performance computing clusters, or
in the cloud. See https://nf-co.re/scflow for the accompanying NextFlow pipeline.

Get Started and Documentation
-----------------------------

Our primary documentation is at https://combiz.github.io/scflow-manual/.

You can install the development version of scFlow from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("combiz/scFlow")
```

Support
-------

- Ask a question [on Stack Overflow with the `scFlow` tag](https://stackoverflow.com/questions/ask?tags=scFlow), we monitor this for new questions.
- Discuss on the [scFlow Slack team](https://ukdri.slack.com/archives/C01BDM7AJ6A).
- Open **bug reports** and **feature requests** (not questions) on [GitHub issues](https://github.com/combiz/scFlow/issues).

How to Contribute
-----------------

Check [CONTRIBUTING](https://github.com/microsoft/scFlow/blob/master/CONTRIBUTING.md) page.

You may need to install scFlow using a Personal Access Token
(Github-&gt;Settings-&gt;Developer Settings): -

``` r
# install.packages("devtools")
devtools::install_github("combiz/scFlow", auth_token = "YOURTOKEN")
```

Reference Papers
----------------

Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen. "[The nf-core framework for community-curated bioinformatics pipelines](https://www.nature.com/articles/s41587-020-0439-x)".  Nat Biotechnology (2020); doi: 10.1038/s41587-020-0439-x

**Note**: If you use scFlow in your GitHub projects, please add `scFlow` in the `requirements.txt`.

License
-------

This project is licensed under the terms of the GPLv3 license. See [LICENSE](https://github.com/combiz/scFlow/blob/master/LICENSE.md) for additional details.
