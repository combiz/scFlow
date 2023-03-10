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
use_r("volcano_plot")

use_r("map_with_ewce")
use_r("reduced_dimension_plots")

?use_r

use_test("function_name")

# add to DESCRIPTION

use_package("assertthat", min_version = "0.2.1")
#use_package("batchelor", min_version = "1.14.0")
use_package("bib2df")
use_package("BiocGenerics")
use_package("biomaRt", min_version = "2.54.0")
use_package("cli", min_version = "3.4.1")
use_package("cowplot", min_version = "1.1.1")
use_package("data.table", min_version = "1.14.6")
use_package("DelayedArray", min_version = "0.24.0")
use_package("DelayedMatrixStats", min_version = "1.20.0")
use_package("DirichletReg", min_version = "0.7.1")
use_package("DoubletFinder", min_version = "2.0.3")
use_package("dplyr", min_version = "1.1.0")
use_package("DropletUtils", min_version = "1.18.1")
use_package("DT", min_version = "0.26")
use_package("edgeR", min_version = "3.40.0")
use_package("english", min_version = "1.2.6")
use_package("enrichR", min_version = "3.1")
use_package("EWCE", min_version = "1.6.0")
use_package("forcats", min_version = "0.5.2")
use_package("formattable", min_version = "0.2.1")
use_package("future", min_version = "1.29.0")
use_package("future.apply", min_version = "1.10.0")
use_package("gdtools")
#use_package("GenomicRanges")
use_package("ggdendro", min_version = "0.1.23")
use_package("ggplot2", min_version = "3.4.0")
use_package("ggpubr", min_version = "0.5.0")
use_package("ggrepel", min_version = "0.9.2")
use_package("ggridges", min_version = "0.5.4")
use_package("httr", min_version = "1.4.4")
use_package("igraph", min_version = "1.3.5")
#use_package("IRanges")
use_package("kBET", min_version = "0.99.6")
use_package("leaflet", min_version = "2.1.1")
use_package("limma", min_version = "3.54.0")
use_package("lme4", min_version = "1.1.31")
use_package("magrittr", min_version = "2.0.3")
use_package("MAST", min_version = "1.24.0")
use_package("Matrix")
use_package("monocle3", min_version = "1.3.1")
use_package("paletteer", min_version = "1.5.0")
use_package("patchwork")
use_package("plotly")
use_package("plyr", min_version = "1.8.8")
use_package("preprocessCore")
use_package("prettydoc")
use_package("qs")
use_package("R.utils")
use_package("purrr", min_version = "0.3.5")
use_package("RANN", min_version = "2.6.1")
use_package("RcppParallel", min_version = "5.1.5")
use_package("rlang", min_version = "1.0.6")
use_package("rliger", min_version = "1.0.0")
use_package("rmarkdown", min_version = "2.18")
use_package("Rtsne", min_version = "0.16")
use_package("S4Vectors", min_version = "0.36.0")
use_package("scales")
use_package("scater", min_version = "1.26.1")
#use_package("sctransform")
use_package("Seurat", min_version = "4.3.0")
use_package("SingleCellExperiment", min_version = "1.20.0")
use_package("stringr", min_version = "1.4.1")
use_package("SummarizedExperiment", min_version = "1.28.0")
use_package("threejs", min_version = "0.3.3")
use_package("tibble", min_version = "3.1.8")
use_package("tidyr", min_version = "1.2.1")
use_package("tidyselect", min_version = "1.2.0")
use_package("UpSetR", min_version = "1.4.0")
use_package("utils", min_version = "4.2.1")
use_package("uwot", min_version = "0.1.14.9000")
#use_package("vroom", min_version = "1.6.0")
use_package("WebGestaltR", min_version = "0.4.4")
use_package("XML")



use_vignette("basic-qc", title = "Guided tutorial for sample quality control")

#use_vignette("How to do this cool analysis") # later
use_github_actions()
usethis::use_github_action("check-standard")
usethis::use_github_action("test-coverage")
usethis::use_github_action("lint")
usethis::use_github_actions_badge(name = "R-CMD-check", repo_spec = NULL)
usethis::use_github_action_check_standard()

pkgdown::init_site() #create favicons from pkg logo

### BEFORE EVERY COMMIT
#Restart R Session Cmd+Shift+F10 (Ctrl+Shift+F10 for Windows)
#Document Package Cmd+Shift+D (Ctrl+Shift+D for Windows)
#Check Package Cmd+Shift+E (Ctrl+Shift+E for Windows)

## before every release
# knit the readme.Rmd <<--
# update the site
use_version()
pkgdown::build_site()


build()

options(rgl.useNULL = TRUE)

#########################
