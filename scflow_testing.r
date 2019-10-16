library(scflow)

ensembl_tsv <- read.delim("~/Documents/nf-sc/src/ensembl-ids/ensemble_ids_all.tsv")
mat <- read_feature_barcode_matrix("~/Documents/ms-sc/data/raw/testfbmatrix/outs/raw_feature_bc_matrix")


x <- map_en
