# QC demo
# 27/11/19

##  ............................................................................
##  Initialize                                                              ####
options(mc.cores = future::availableCores() - 2)
library(parallel)
library(scFlow)

# v2 chemistry
matpath <- "~/Documents/ms-sc/data/raw/testfbmatrix/outs/raw_feature_bc_matrix"
#matpath <- "~/Documents/junk/MS535" #bad sample
#matpath <- "~/Documents/junk/MS461" #ok sample
matpath <- "~/Documents/junk/MS426/outs/raw_feature_bc_matrix" # emptydrops peculiar sample (large pass #)
#matpath <- "~/Documents/junk/MS523/outs/raw_feature_bc_matrix/"# emptydrops peculiar sample (all < retain fail)
# v3 chemistry, enriched
#matpath <- "~/Documents/testmatrices/enriched"

ensembl_fp <- "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv"
samplesheet_fp <- "~/Documents/junk/refs/SampleSheet.tsv"
ctd_fp <- "~/Documents/junk/refs/ctd/"

##  ............................................................................
##  Start QC                                                                ####

x <- Sys.time()

mat <- read_sparse_matrix(matpath)

metadata <- read_metadata(
  unique_key = "hajov",
  key_colname = "manifest",
  samplesheet_path = samplesheet_fp)

sce <- generate_sce(mat, metadata)

#sce <- find_cells(sce, lower = 125, retain = NULL)

sce <- find_cells(sce, lower = 100, retain = 12000, niters = 30000)
sce <- annotate_sce(sce, ensembl_mapping_file = ensembl_fp, min_library_size = 250, max_library_size = "adaptive", nmads = 4)
#

#sce <- annotate_sce(sce, ensembl_mapping_file = ensembl_fp, min_library_size = 100, max_library_size = "adaptive")
#sce <- annotate_sce(sce, ensembl_mapping_file = ensembl_fp, min_library_size = 300, max_library_size = 10000, max_features = 500)

sce <- filter_sce(sce)

qs::qsave(sce, "sce.qs")
sce <- qs::qread("sce.qs", nthreads = future::availableCores())

###

###

#df2 <- read.csv("~/Documents/scFlow/misc/emptydropssweeping2.csv")

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


## with piping


metadata <- read_metadata(
  unique_key = "hajov",
  key_colname = "manifest",
  samplesheet_path = samplesheet_fp)

sce <- read_sparse_matrix(matpath) %>%
  generate_sce(metadata) %>%
  find_cells() %>%
  annotate_sce() %>%
  filter_sce() %>%
  find_singlets() %>%
  filter_sce() %>%
  report_qc_sce()

write_sce()

sce <- merge_sce() %>%
  report_merged_sce() %>%
  integrate_sce() %>%
  reduce_dims_sce() %>%
  cluster_sce()

write_sce()

de_res <- perform_de(sce)
ip_res <- find_impacted_pathways(de_res)




