# QC demo
# 27/11/19

##  ............................................................................
##  Initialize                                                              ####
#options(mc.cores = parallel::detectCores())
#library(parallel)
library(scFlow)

# v2 chemistry
matpath <- "~/Documents/ms-sc/data/raw/testfbmatrix/outs/raw_feature_bc_matrix"
# v3 chemistry, enriched
#matpath <- "~/Documents/testmatrices/enriched"

ensembl_fp <- "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv"
samplesheet_fp <- "~/Documents/nf-sc/refs/SampleSheet.tsv"
ctd_fp <- "~/Documents/nf-sc/refs/ctd/"

##  ............................................................................
##  Start QC                                                                ####

x <- Sys.time()

mat <- read_sparse_matrix(matpath)

metadata <- read_metadata(
  unique_key = "hajov",
  key_colname = "manifest",
  samplesheet_path = samplesheet_fp)

sce <- generate_sce(mat, metadata)

sce <- find_cells(sce, lower = 100)

sce <- annotate_sce(sce, ensembl_mapping_file = ensembl_fp, min_library_size = 100)

sce <- filter_sce(sce)

sce <- find_singlets(sce, "doubletfinder", pK = 0.005, vars_to_regress_out = NULL)

sce <- filter_sce(sce)

report_qc_sce(sce)

totaltime <- Sys.time() - x

print(totaltime)

alarm()

##  ............................................................................
##  Bonus - cluster and identify celltypes                                  ####

x <- Sys.time()

sce <- reduce_dims_sce(sce, pca_dims = 5, reduction_methods = c("tSNE", "UMAP"))

sce <- cluster_sce(sce)

sce <- map_celltypes_sce(sce, ctd_folder = ctd_fp)

totaltime <- Sys.time() - x

print(totaltime)
alarm()
