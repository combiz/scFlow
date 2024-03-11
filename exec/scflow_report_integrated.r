#!/usr/bin/env Rscript
# Reduce dimensions for a SCE
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = future::availableCores())

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)
library(parallel)
library(SingleCellExperiment) # due to monocle3 missing namespace::
library(knitr) # due to missing knitr:: namespace in the integrate report

##  ............................................................................
##  Parse command-line arguments                                            ####

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

required$add_argument(
  "--input_reduced_dim",
  help = "reduced dimension embedding to use for the integration report",
  metavar = "UMAP",
  required = TRUE
)

required$add_argument(
  "--reddimplot_pointsize",
  default = 0.1,
  type = "double",
  required = TRUE,
  help = "Point size for reduced dimension plots",
  metavar = "N"
)

required$add_argument(
  "--reddimplot_alpha",
  default = 0.2,
  type = "double",
  required = TRUE,
  help = "Alpha value for reduced dimension plots",
  metavar = "N"
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()
args$categorical_covariates <- strsplit(args$categorical_covariates, ",")[[1]]

options("scflow_reddimplot_pointsize" = args$reddimplot_pointsize)
options("scflow_reddimplot_alpha" = args$reddimplot_alpha)


##  ............................................................................
##  Start                                                                   ####

sce <- read_sce(args$sce_path, read_metadata = TRUE)

sce <- annotate_integrated_sce(
  sce,
  categorical_covariates = args$categorical_covariates,
  input_reduced_dim = args$input_reduced_dim
)

##  ............................................................................
##  Save Outputs                                                            ####

dir.create(file.path(getwd(), "integration_report"))

report_integrated_sce(
  sce = sce,
  report_folder_path = file.path(getwd(), "integration_report"),
  report_file = "integrate_report_scflow"
)

##  ............................................................................
##  Clean up                                                                ####

# Clear biomart cache
