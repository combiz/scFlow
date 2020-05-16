#install.packages(c("devtools", "roxygen2", "usethis", "spelling"))
library(devtools)
library(roxygen2)
library(usethis)
library(spelling)



# package setup following approx Emil Hvitfeldt

## Initial package creation steps
#create_package("~/Documents/scflow")
#options(usethis.full_name = "Combiz Khozoie")
#use_git()
#use_github()
#use_mit_license()
#use_readme_rmd()
#use_travis()
#use_data_raw()
#use_appveyor()
#use_coverage(type = c("codecov"))
#use_testthat()
#use_spell_check()
#use_news_md()
#use_pkgdown()

use_rmarkdown_template(
  template_name = "Quality Control",
  template_description = "Quality control report for a SingleCellExperiment",
  template_create_dir = FALSE
)

## Ongoing package dev
use_r("function_name")
use_r("retrieve_sample_metadata")

use_r("sce_to_seurat")

use_r("read_feature_barcode_matrix")
use_r("map_ensembl_gene_id")
use_r("generate_sce_from_mat_and_meta")
use_r("annotate_sce_genes")

use_r("annotate_sce")
use_r("generate_sce")
use_r("annotate_sce_cells")
use_r("filter_sce")
use_r("write_sce")
use_r("read_sce")

use_r("merge_sce")

use_r("citations")

use_r("report_qc_sce")

use_r("find_singlets_with_doubletfinder")

use_r("write_feature_barcode_matrix")

use_r(".sce_to_cds")
use_r("reduce_dims_sce")
use_r("cluster_sce")

use_r("diffexp_models")

use_r("map_with_ewce")
use_r("reduced_dimension_plots")

?use_r

use_test("function_name")


use_package("dplyr") # add to DESCRIPTION
use_package("rmarkdown")
use_package("DT")
use_package("stats")
use_package("RcppParallel")
use_package("monocle3")
use_package("plotly")
use_package("threejs")
use_package("assertthat")
use_package("EWCE")
use_package("plyr")
use_package("httr")
use_package("MAST")
use_package("bib2df")
use_package("limma")
use_package("methods")
use_package("paletteer")
use_package("scater")
use_package("leaflet")
use_package("liger")
use_package("RANN.L1")
use_package("ROntoTools")
use_package("BiocGenerics")
use_package("S4Vectors")
use_package("IRanges")
use_package("GenomicRanges")
use_package("DelayedArray")
use_package("gdtools")
use_package("WebGestaltR")
use_package("graphite")
use_package("uniftest")
use_package("liger")
use_package("batchelor")
use_package("DelayedMatrixStats")
use_package("limma")
use_package("Matrix.utils")
use_package("grDevices")
use_package("ggridges")
use_package("UpSetR")
use_package("nVennR")
use_package("kBET")

use_vignette("basic-qc", title = "Guided tutorial for sample quality control")

#use_vignette("How to do this cool analysis") # later

### BEFORE EVERY COMMIT
#Restart R Session Cmd+Shift+F10 (Ctrl+Shift+F10 for Windows)
#Document Package Cmd+Shift+D (Ctrl+Shift+D for Windows)
#Check Package Cmd+Shift+E (Ctrl+Shift+E for Windows)

## before every release
# knit the readme.Rmd <<--
# update the site
pkgdown::build_site()

use_version()
build()

#########################
use_package("purrr")
use_package("cli")
use_package("Matrix")
use_package("utils")
use_package("biomaRt")
use_package("stringr")
use_package("data.table")
use_package("SummarizedExperiment")
use_package("SingleCellExperiment")
use_package("vroom")
use_package("magrittr")
use_package("english")
use_package("Seurat")
use_package("DoubletFinder")
use_package("ggplot2")
use_package("rlang")
use_package("sctransform")
use_package("future")
use_package("future.apply")
