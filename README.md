<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" width=200 style="align:left;" />

# scFlow - Single-Cell Workflow

<!-- badges: start -->

[![R GitHub
Actions](https://github.com/combiz/scFlow/actions/workflows/r_package.yaml/badge.svg)](https://github.com/combiz/scFlow/actions/workflows/r_package.yaml)
[![R GitHub Actions
Dev](https://github.com/combiz/scFlow/actions/workflows/r_package2.yaml/badge.svg)](https://github.com/combiz/scFlow/actions/workflows/r_package2.yaml)
[![Documentation](https://img.shields.io/badge/docs-passing-brightgreen.svg?style=flat)](https://combiz.github.io/scFlow/)
[![Manual](https://img.shields.io/badge/manual-passing-brightgreen.svg?style=flat)](https://combiz.github.io/scflow-manual/)
[![contributions
welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/combiz/scFlow/issues)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-green.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

The scFlow R package provides the analytical back-end to the
[nf-core/scflow](https://nf-co.re/scflow) Nextflow pipeline for the
orchestration of automated, scalable, and reproducible single-cell RNA
sequencing analyses.

The scFlow R package is built to enable standardized workflows following
best practices on top of popular single-cell R packages, including
Seurat, Monocle, scater, emptyDrops, DoubletFinder, LIGER, and MAST (Hao
2021, Cao 2019, McCarthy 2017, Lun 2019, McGinnis 2019, Welch 2019).
Common analytical tasks required by users that involve multiple tools
can be performed in scFlow with a single command (i.e. a higher level of
abstraction). Key plots, tables, and other analysis outputs are
automatically generated, together with interactive HTML reports for each
stage of the analysis. These reports provide informative summaries of
analytical steps in ways that can highlight the impact of parameter
choices and guide their revision when needed.

The following example illustrates a complete sample quality-control with
default parameters using the scFlow R package, including ambient RNA
profiling, gene/cell annotation, thresholding, doublet/multiplet
removal, and generation of an interactive HTML report with key plots: -

    sce <- read_sparse_matrix(matrix_path) %>%
    generate_sce(metadata) %>%
    find_cells() %>%
    annotate_sce() %>%
    filter_sce() %>%
    find_singlets() %>%
    filter_sce() %>%
    report_qc_sce()

Overview of scFlow features

-   Quality control of gene-cell matrices
    -   Profiling of ambient RNA with
        [emptyDrops](https://github.com/MarioniLab/DropletUtils)
    -   Filtering of matrices by counts and features including optional
        adaptive filtering
    -   Filtering of mitochondrial and ribosomal genes and thresholding
        of counts
    -   Doublet and multiplet identification and removal with
        [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
    -   Rich QC metrics annotation with
        [scater](https://github.com/davismcc/scater)
-   Integration
    -   Calculation of latent-gene metafactors for integration
        [Liger](https://github.com/MacoskoLab/rliger)
    -   Evaluation of integration performance across major sources of
        experimental variance with
        [kBet](https://github.com/theislab/kBET)
    -   Latent factor contribution analysis with
        [UpSetR](https://github.com/hms-dbmi/UpSetR/) plots
-   Dimensionality reduction and Clustering
    -   Dimensionality reduction with tSNE and/or UMAP with/without
        integration
    -   Community detection with the Louvain or Leiden clustering
        algorithms using
        [monocle](https://github.com/cole-trapnell-lab/monocle-release)
-   Cell-type identification
    -   Cluster marker gene identification with
        [monocle](https://github.com/cole-trapnell-lab/monocle-release)
    -   Automated cell-type annotation with
        [EWCE](https://github.com/NathanSkene/EWCE)
    -   Cluster identity mapping against the [Allen Human Brain
        Atlas](https://www.brain-map.org)
-   Differential gene expression implementations
    -   Flexible pre-processing including gene expressivity filtering,
        and options for numeric covariate scaling/centering, matrix
        normalization, and pseudobulking
    -   Zero-inflated mixed-effect regression model with
        [MAST](https://github.com/RGLab/MAST)
    -   Random effects models with
        [Limma](https://github.com/cran/limma)
-   Pathway and functional category enrichment analysis
    -   Interface to the Enrichr database with
        [EnrichR](https://github.com/cran/enrichR)
    -   Interface to the WebGestalt tool with
        [WebGestaltR](http://www.webgestalt.org/)
-   Cell-type composition analysis
    -   Dirichlet modeling of differential cell-type composition using
        [DirichletReg](https://cran.r-project.org/web/packages/DirichletReg/index.html)
-   Publication quality plots and analysis reports
    -   Eight interactive HTML reports with key plots and algorithm
        performance metrics

## Get Started and Documentation

Our primary documentation is at
<https://combiz.github.io/scflow-manual/>.

You can install the development version of scFlow from GitHub with: -

``` r
# install.packages("devtools")
devtools::install_github("combiz/scFlow")
```

## Running RStudio Server with pre-installed scFlow via Docker

If using linux system, open a cmd terminal and run the following
commands:

``` r
docker pull nfancy/scflow

docker run --rm -d \
-v $HOME:/home/rstudio/home \
-e ROOT=true \
-e PASSWORD=password \
-p 8950:8787 \
nfancy/scflow
```

Then access the RStudio Server on your browser at IP:8950. The default
username is `rstudio` and password is `password`. To understand more on
how to run docker Rstudio follow the link
[here](https://davetang.org/muse/2021/04/24/running-rstudio-server-with-docker/).

## Support

-   Ask a question on Stack Overflow with the scFlow tag, we monitor
    this for new questions.
-   Discuss on the scFlow Slack team.
-   Open bug reports and feature requests (not questions) on GitHub
    issues.

## How to Contribute

Check the
[CONTRIBUTING](https://github.com/scFlow/blob/master/CONTRIBUTING.md)
page.

## Reference Papers

Combiz Khozoie, Nurun Fancy, Mahdi M. Marjaneh, Alan E. Murphy, Paul M.
Matthews, Nathan Skene. “scFlow: A Scalable and Reproducible Analysis
Pipeline for Single-Cell RNA Sequencing Data.” bioRxiv 2021 August 19.
doi: 10.22541/au.162912533.38489960/v1.

Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes
Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven
Nahnsen. “The nf-core framework for community-curated bioinformatics
pipelines”. Nat Biotechnology (2020); doi: 10.1038/s41587-020-0439-x

Note: If you use scFlow in your GitHub projects, please add scFlow in
the requirements.txt.
