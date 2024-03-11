#!/usr/bin/env Rscript
# Model relative celltype abundance with a Dirichlet model
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = future::availableCores())

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)
library(cli)

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")

required$add_argument(
  "--sce_path",
  help = "path to SingleCellExperiment directory",
  metavar = "/dir/sce/",
  required = TRUE
)

required$add_argument(
  "--unique_id_var",
  help = "unique sample variable name",
  metavar = "manifest",
  required = TRUE
)

required$add_argument(
  "--celltype_var",
  help = "celltype variable name",
  metavar = "cluster_celltype",
  required = TRUE
)

required$add_argument(
  "--dependent_var",
  help = "name of the dependent variable",
  metavar = "group",
  required = TRUE
)

required$add_argument(
  "--ref_class",
  help = "class of the reference variable withing dependent_var",
  metavar = "control",
  required = TRUE
)

required$add_argument(
  "--var_order",
  help = "factor order for dependent variable",
  metavar = "c",
  required = TRUE
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
args <- parser$parse_args()
args$var_order <- strsplit(args$var_order, ",")[[1]]
if (tolower(args$var_order) == "null") { args$var_order <- NULL }

#   ____________________________________________________________________________
#   Start                                                                   ####

sce <- read_sce(args$sce_path)

results <- model_celltype_freqs(
  sce,
  unique_id_var = args$unique_id_var,
  celltype_var = args$celltype_var,
  dependent_var = args$dependent_var,
  ref_class = args$ref_class,
  var_order = args$var_order
)

## ............................................................................
## Save Outputs ####

new_dirs <- c(
  "dirichlet_report"
)

#make dirs
purrr::walk(new_dirs, ~ dir.create(file.path(getwd(), .)))

report_celltype_model(
  results,
  report_folder_path = file.path(
    getwd(),
    "dirichlet_report"
  ),
  report_file = paste0(
    args$celltype_var,
    args$dependent_var,
    "dirichlet_report",
    sep = "_")
)
