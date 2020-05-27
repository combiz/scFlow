##  ............................................................................
##  Initialize                                                              ####
options(mc.cores = future::availableCores())
library(parallel)
library(scFlow)
print(getwd())

# v2 chemistry
#matpath <- "~/Documents/ms-sc/data/raw/testfbmatrix/outs/raw_feature_bc_matrix"
#matpath <- "~/Documents/junk/MS535" #bad sample
matpath <- "~/Documents/junk/MS461" #ok sample
# v3 chemistry, enriched
#matpath <- "~/Documents/testmatrices/enriched"

ensembl_fp <- "~/Documents/nf-sc/src/ensembl-ids/ensembl_mappings.tsv"
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

#sce <- find_cells(sce, lower = 100, retain = 300)
#sce <- find_cells(sce, lower = 100, retain = 674)
# v2
df <- data.frame()
for (lower_param in seq(50, 150, 5)) {
  print(lower_param)
  sce <- find_cells(sce, lower = lower_param, retain = NULL, niters = 10000)
  df <- rbind(df, as.data.frame(sce@metadata$emptydrops_params))
}
write.table(df, "emptydropssweeping2.csv", sep = ",", row.names = F, col.names = T)
