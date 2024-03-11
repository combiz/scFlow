#!/usr/bin/env Rscript
# Map celltypes for a SCE
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = future::availableCores())

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)
library(parallel)

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
  "--ctd_folder",
  help = "path to a folder containing ewce ctd files",
  metavar = "foo/bar",
  required = TRUE
)

required$add_argument(
  "--clusters_colname",
  help = "the sce colData variable storing cluster numbers",
  metavar = "foo/bar",
  required = TRUE
)

required$add_argument(
  "--cells_to_sample",
  type = "integer",
  default = 10000,
  help = "the number of cells to sample with ewce",
  metavar = "N",
  required = TRUE
)

required$add_argument(
  "--species",
  help = "the biological species (e.g. mouse, human)",
  default = "human",
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

options("scflow_species" = args$species)
options("scflow_reddimplot_pointsize" = args$reddimplot_pointsize)
options("scflow_reddimplot_alpha" = args$reddimplot_alpha)

##  ............................................................................
##  Start                                                                   ####

cat(print(tempdir()))

sce <- read_sce(args$sce_path)

sce <- map_celltypes_sce(
  sce,
  ctd_folder = args$ctd_folder,
  clusters_colname = args$clusters_colname,
  cells_to_sample = args$cells_to_sample,
  species = args$species
)

##  ............................................................................
##  Save Outputs                                                            ####

write_celltype_mappings(sce, folder_path = getwd())

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "celltype_mapped_sce")
)

##  ............................................................................
##  Clean up                                                                ####

# Clear biomart cache
