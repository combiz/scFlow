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
  "--cluster_method",
  help = "method to use for clustering",
  metavar = "louvain",
  required = TRUE
)

required$add_argument(
  "--reduction_method",
  help = "reduced dimension embedding to use for clustering",
  metavar = "UMAP",
  required = TRUE
)

required$add_argument(
  "--res",
  type = "double",
  default = 0.00001,
  help = "clustering resolution",
  metavar = "N",
  required = TRUE
)

required$add_argument(
  "--k",
  type = "integer",
  default = 100,
  help = "the number of kNN",
  metavar = "N",
  required = TRUE
)

required$add_argument(
  "--louvain_iter",
  type = "integer",
  default = 1,
  help = "number of iterations used for Louvain clustering",
  metavar = "N",
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()

##  ............................................................................
##  Start                                                                   ####

sce <- read_sce(args$sce_path, read_metadata = TRUE)

sce <- cluster_sce(
  sce,
  cluster_method = args$cluster_method,
  reduction_method = args$reduction_method,
  res = args$res,
  k = args$k,
  louvain_iter = args$louvain_iter
)

##  ............................................................................
##  Save Outputs                                                            ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "clustered_sce"),
  write_metadata = TRUE
)

##  ............................................................................
##  Clean up                                                                ####

# Clear biomart cache
