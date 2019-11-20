## SET THE ENVIRONMENT VAR BEFORE LOADING PARALLEL
#Sys.setenv("MC_CORES" = 10L)
#options("mc.cores")
#library(parallel)
#library(future)
library(scflow)
#library(devtools)
x <- Sys.time()
#plan(list(multiprocess, sequential))

#install_github("satijalab/seurat", ref = "develop")

matpath <- "~/Documents/ms-sc/data/raw/testfbmatrix/outs/raw_feature_bc_matrix"
#matpath <- "~/Documents/testmatrices/enriched"

#ensembl_tsv <- read.delim("~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv")

mat <- read_sparse_matrix(matpath)

ss_classes <- c(
  batch = "factor",
  capdate = "factor",
  prepdate = "factor",
  seqdate = "factor",
  aplevel = "factor"
)

metadata <- retrieve_metadata(
  unique_id = "MS542",
  id_colname = "individual",
  samplesheet_path = "~/Documents/ms-sc/refs/sample_metadata.tsv",
  colClasses = ss_classes
)

sce <- generate_sce(mat, metadata)
sce <- .append_citation_sce(sce, key = c("doubletfinder", "seurat"))

rm(mat, metadata)

sce <- annotate_sce(
  sce,
  ensembl_mapping_file = "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv"
)

sce <- filter_sce(sce, filter_genes = TRUE, filter_cells = TRUE)
sce <- find_singlets(sce, "doubletfinder", pK = 0.005, vars_to_regress_out = NULL)
#sce <- find_singlets(sce, "doubletfinder", vars_to_regress_out = c("nCount_RNA", "pc_mito"))
sce <- filter_sce(sce)
report_qc_sce(sce)
#saveRDS(sce@metadata, "../junk/sce/metadata.rds")
totaltime <- Sys.time() - x
print(totaltime)

#start_time <- Sys.time()
#metadata_tmp_path <- file.path(tempdir(), "metadata.rds")
#saveRDS(sce@metadata, metadata_tmp_path)
#end_time <- Sys.time()
#end_time - start_time

install.packages("mvbutils")
library(mvbutils)

foodweb( where = "package:scflow",
        border = TRUE, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=1)
mtext("scflow functions")

saveRDS(sce@metadata$qc_plots$number_genes_vs_count_depth, "test.rds")
#BiocManager::install("Rgraphviz")
