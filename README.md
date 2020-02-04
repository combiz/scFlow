
<!-- README.md is generated from README.Rmd. Please edit that file -->
scflow
======

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/combiz/scflow.svg?branch=master)](https://travis-ci.org/combiz/scflow) [![Codecov test coverage](https://codecov.io/gh/combiz/scflow/branch/master/graph/badge.svg)](https://codecov.io/gh/combiz/scflow?branch=master) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/combiz/scflow?branch=master&svg=true)](https://ci.appveyor.com/project/combiz/scflow) <!-- badges: end -->

The goal of scflow is to provide tools in R to build a complete analysis workflow for single-cell/nuclei RNA sequencing data.

-   Quality control of gene-cell matrices
-   Filtering of matrices by counts and features
-   Filtering of mitochondrial genes and mitochondrial counts thresholding
-   Doublet and multiplet identification and removal with [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
-   Rich QC metrics annotation with [scater](https://github.com/davismcc/scater)

-   Dimensionality reduction and celltype identification
-   Louvain clustering, UMAP dimensionality reduction, and cluster marker gene identification with [monocle](https://github.com/cole-trapnell-lab/monocle-release)
-   Celltype annotation with [EWCE](https://github.com/NathanSkene/EWCE) and [Liger](https://github.com/MacoskoLab/liger)
-   Cluster identity mapping against the [Allen Human Brain Atlas](https://www.brain-map.org) and Mouse Nervous System Data from [Zeisel 2018](https://doi.org/10.1016/j.cell.2018.06.021)

-   Differential gene expression implementations
-   Zero-inflated regression model with [MAST](https://github.com/RGLab/MAST)
-   Random effects model with [Limma](https://github.com/cran/limma)
-   Negative binomial distribution pseudobulking model with [DESeq2](https://github.com/mikelove/DESeq2)
-   Pseudobulk generalized likelihood ratio tests with [EdgeR](https://github.com/StoreyLab/edge)

-   Pathway and functional category enrichment analysis
-   Interface to the Enrichr database with [EnrichR](https://github.com/cran/enrichR)
-   Interface to the WebGestalt tool with [WebGestaltR](http://www.webgestalt.org/)

-   Publication quality plots and analysis reports
-   QC plots and tabular metrics suitable for reports.
-   UMAP plots for cell features and gene expression.
-   Violin plots for gene expression.
-   Pathway and gene enrichment plots

The package functions are designed to interface neatly with [NextFlow](https://www.nextflow.io/) for scalable and containerized pipelines deployed locally, on high-performance computing clusters, or in the cloud. An accompanying NextFlow pipeline is in the works - TBA.

Installation
------------

You can install the development version of scflow from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("combiz/scflow")
```

Developers
----------

You may need to install scflow using a Personal Access Token (Github-&gt;Settings-&gt;Developer Settings): -

``` r
# install.packages("devtools")
devtools::install_github("combiz/scflow", auth_token = "YOURTOKEN")
```

Running *scflow* using an example dataset
-----------------------------------------

Preprocess the `scFlowExample` dataset following this [link](https://github.com/neurogenomics/scFlowExample).

Running *scflow*
----------------

The basic `scflow` workflow for sample QC begins with the import of the feature-barcode sparse matrix with `read_sparse_matrix`. The metadata for the sample is then imported from a sample sheet with `read_metadata`. A SingleCellExperiment object is created from the matrix and the metadata using `generate_sce` which is then annotated with both gene and cell-level data using `annotate_sce`. We then filter the SingleCellExperiment to select only cells and genes meeting our QC criteria using `filter_sce`. We can then optionally find singlets in our filtered SingleCellExperiment using `find_singlets` before filtering them out again with `filter_sce`. A complete QC report can then be generated using `report_qc_sce` before saving the filtered and quality-controlled SingleCellExperiment with `write_sce`.

### Step one - import the matrix and metadata

``` r
library(scflow)
mat <- read_sparse_matrix(matrix_fp)
```

Next we retrieve the metadata by pointing to a Sample Sheet and specifying a unique identifier (unique\_key) in a specific column (key\_colname):

`r   metadata <- read_metadata(   unique_key = "lodut",   key_colname = "manifest",   samplesheet_path = samplesheet_fp   )`
