#!/usr/bin/env Rscript
#' Annotate integrated, dims reduced and clustered sce object
# Mahdi Moradi Marjaneh

# ____________________________________________________________________________
# Initialization ####

options(mc.cores = future::availableCores())

## ............................................................................
## Load packages ####
library(argparse)
library(scFlow)

## ............................................................................
## Parse command-line arguments ####

# create parser object
parser <- ArgumentParser()

# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")

required$add_argument(
  "--sce_path",
  help = "-path to the SingleCellExperiment",
  metavar = "dir",
  required = TRUE
)

required$add_argument(
  "--categorical_covariates",
  help = "-categorical covariates",
  metavar = "individual,diagnosis,region,sex",
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args ####

args <- parser$parse_args()
args <- purrr::map(args, function(x) {
if (length(x) == 1) {
if (toupper(x) == "TRUE") {
return(TRUE)
}
if (toupper(x) == "FALSE") {
return(FALSE)
}
if (toupper(x) == "NULL") {
return(NULL)
}
}
return(x)
})

## ............................................................................
## Annotate integrated sce ####

sce <- read_sce(args$sce_path)

sce <- annotate_integrated_sce(
sce,
categorical_covariates = args$categorical_covariates
)

dir.create(file.path(getwd(), "integration_report"))

report_integrated_sce(
  sce = sce,
  report_folder_path = file.path(getwd(), "integration_report"),
  report_file = "integrate_reduceDims_cluster_report_scflow",
)

print("Annotation complete, saving outputs..")

## ............................................................................
## Save Outputs ####

# Save SingleCellExperiment
write_sce(
sce = sce,
folder_path = file.path(getwd(), "integrated_sce")
)
