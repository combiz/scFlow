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

use_r("report_qc_sce")

use_r("find_singlets_with_doubletfinder")

use_r("write_feature_barcode_matrix")

use_test("function_name")

use_package("dplyr") # add to DESCRIPTION

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
