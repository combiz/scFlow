#!/usr/bin/env Rscript
# Integrate multiple single cell datasets (samples)
# Mahdi Moradi Marjaneh

# ____________________________________________________________________________
# Initialization ####

options(mc.cores = future::availableCores())

## ............................................................................
## Load packages ####
library(argparse)
library(scFlow)
library(parallel)

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
  "--method",
  required = TRUE,
  help = "The integration method to use",
  metavar = "Liger"
)

required$add_argument(
  "--unique_id_var",
  required = TRUE,
  help = "Unique id variable",
  metavar = "manifest"
)

required$add_argument(
  "--take_gene_union",
  default = FALSE,
  required = TRUE,
  help = "Whether to fill out raw.data matrices with union of genes",
  metavar = "Boolean"
)

required$add_argument(
  "--remove_missing",
  default = TRUE,
  required = TRUE,
  help = "Remove non-expressive genes and cells",
  metavar = "Boolean"
)

required$add_argument(
  "--num_genes",
  default = 3000,
  type = "integer",
  required = TRUE,
  help = "Number of genes to find for each dataset",
  metavar = "N"
)

required$add_argument(
  "--combine",
  default = "union",
  required = TRUE,
  help = "How to combine variable genes across experiments",
  metavar = "union,intersect"
)

required$add_argument(
  "--keep_unique",
  default = FALSE,
  required = TRUE,
  help = "Keep genes that occur only in one dataset",
  metavar = "Boolean"
)

required$add_argument(
  "--capitalize",
  default = FALSE,
  required = TRUE,
  help = "Capitalize gene names to match homologous genes(i.e. across species)",
  metavar = "Boolean"
)

required$add_argument(
  "--use_cols",
  default = TRUE,
  required = TRUE,
  help = "Treat each column as a cell",
  metavar = "Boolean"
)

required$add_argument(
  "--k",
  default = 30,
  type = "integer",
  required = TRUE,
  help = "Inner dimension of factorization (number of factors)",
  metavar = "N"
)

required$add_argument(
  "--lambda",
  default = 5.0,
  type = "double",
  required = TRUE,
  help = "Regularization parameter",
  metavar = "N"
)

required$add_argument(
  "--thresh",
  default = 0.0001,
  type = "double",
  required = TRUE,
  help = "Convergence threshold.",
  metavar = "N"
)

required$add_argument(
  "--max_iters",
  default = 100,
  type = "integer",
  required = TRUE,
  help = "Maximum number of block coordinate descent iterations to perform",
  metavar = "N"
)

required$add_argument(
  "--nrep",
  default = 1,
  type = "integer",
  required = TRUE,
  help = "Number of restarts to perform",
  metavar = "N"
)

required$add_argument(
  "--rand_seed",
  default = 1,
  type = "integer",
  required = TRUE,
  help = "Random seed to allow reproducible results",
  metavar = "N"
)

required$add_argument(
  "--knn_k",
  default = 20,
  type = "integer",
  required = TRUE,
  help = "Number of nearest neighbors for within-dataset knn graph",
  metavar = "N"
)

required$add_argument(
  "--k2",
  default = 500,
  type = "integer",
  required = TRUE,
  help = "Horizon parameter for shared nearest factor graph",
  metavar = "N"
)

required$add_argument(
  "--prune_thresh",
  default = 0.2,
  type = "double",
  required = TRUE,
  help = "Minimum allowed edge weight. Any edges below this are removed",
  metavar = "N"
)

required$add_argument(
  "--ref_dataset",
  default = "",
  required = TRUE,
  help = "Name of dataset to use as a reference for normalization",
  metavar = "ref"
)

required$add_argument(
  "--min_cells",
  default = 2,
  type = "integer",
  required = TRUE,
  help = "Minimum number of cells to consider a cluster shared across datasets",
  metavar = "N"
)

required$add_argument(
  "--quantiles",
  default = 50,
  type = "integer",
  required = TRUE,
  help = "Number of quantiles to use for quantile normalization",
  metavar = "N"
)

required$add_argument(
  "--nstart",
  default = 10,
  type = "integer",
  required = TRUE,
  help = "Number of times to perform Louvain community detection",
  metavar = "N"
)

required$add_argument(
  "--resolution",
  default = 1,
  type = "double",
  required = TRUE,
  help = "Controls the number of communities detected",
  metavar = "N"
)

required$add_argument(
  "--dims_use",
  default = "null",
  required = TRUE,
  help = "Indices of factors to use for shared nearest factor determination",
  metavar = "Indices"
)

required$add_argument(
  "--dist_use",
  default = "CR",
  required = TRUE,
  help = "Distance metric to use in calculating nearest neighbors",
  metavar = "CR"
)

required$add_argument(
  "--center",
  default = FALSE,
  required = TRUE,
  help = "Centers the data when scaling factors",
  metavar = "Boolean"
)

required$add_argument(
  "--small_clust_thresh",
  default = 0,
  type = "double",
  required = TRUE,
  help = "Extracts small clusters loading highly on single factor",
  metavar = "N"
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
## Integrate sce ####

sce <- read_sce(args$sce_path)

sce <- integrate_sce(
  sce,
  method = args$method,
  unique_id_var = args$unique_id_var,
  take_gene_union = args$take_gene_union,
  remove.missing = args$remove_missing,
  make.sparse = T,
  num_genes = args$num_genes,
  combine = args$combine,
  keep_unique = args$keep_unique,
  capitalize = args$capitalize,
  use_cols = args$use_cols,
  k = args$k,
  lambda = args$lambda,
  thresh = args$thresh,
  max_iters = args$max_iters,
  nrep = args$nrep,
  H_init = NULL,
  W_init = NULL,
  V_init = NULL,
  rand_seed = args$rand_seed,
  knn_k = args$knn_k,
  k2 = args$k2,
  prune_thresh = args$prune_thresh,
  ref_dataset = args$ref_dataset,
  min_cells = args$min_cells,
  quantiles = args$quantiles,
  nstart = args$nstart,
  resolution = args$resolution,
  dims_use = args$dims_use,
  dist_use = args$dist_use,
  center = args$center,
  small_clust_thresh = args$small_clust_thresh,
  do_plot = FALSE,
  id_number = NULL,
  print_obj = FALSE,
  print_mod = FALSE,
  print_align_summary = FALSE
)

## ............................................................................
## Save Outputs ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "integrated_sce"),
  write_metadata = TRUE
)
