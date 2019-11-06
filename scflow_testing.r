## SET THE ENVIRONMENT VAR BEFORE LOADING PARALLEL
#Sys.setenv("MC_CORES" = 10L)
#options("mc.cores")
#library(parallel)
#library(future)
library(scflow)
#library(devtools)

#plan(list(multiprocess, sequential))

#install_github("satijalab/seurat", ref = "develop")

matpath <- "~/Documents/ms-sc/data/raw/testfbmatrix/outs/raw_feature_bc_matrix"

#ensembl_tsv <- read.delim("~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv")

mat <- read_feature_barcode_matrix(matpath)

ss_classes <- c(
  batch = "factor",
  capdate = "factor",
  prepdate = "factor",
  seqdate = "factor",
  aplevel = "factor"
)

metadata <- retrieve_sample_metadata(
  unique_id = "MS542",
  id_colname = "individual",
  samplesheet_path = "~/Documents/ms-sc/refs/sample_metadata.tsv",
  colClasses = ss_classes
)

sce <- generate_sce(mat, metadata)
rm(mat, metadata)

sce <- annotate_sce(
  sce,
  ensembl_mapping_file = "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv"
)

x <- filter_sce(sce)
y <- find_singlets(x, "doubletfinder", pK = 0.005, vars_to_regress_out = NULL)

#sce <- filter_sce(
#  sce,
#  filter_genes = TRUE, filter_cells = TRUE, drop_unmapped = TRUE, drop_mito = TRUE, drop_ribo = FALSE)
#sce <- filter_sce(sce)

#sce <- find_singlets(sce, "doubletfinder", pK = 0.005, vars_to_regress_out = "pc_mito")
sce <- find_singlets(sce, "doubletfinder", pK = 0.005, vars_to_regress_out = NULL)
#sce <- find_singlets(sce, "doubletfinder")

sce <- report_qc_sce(sce)
saveRDS(sce@metadata, "../junk/")

#start_time <- Sys.time()
#metadata_tmp_path <- file.path(tempdir(), "metadata.rds")
#saveRDS(sce@metadata, metadata_tmp_path)
#end_time <- Sys.time()
#end_time - start_time
