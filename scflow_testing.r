library(scflow)
library(utils)

matpath <- "~/Documents/ms-sc/data/raw/testfbmatrix/outs/raw_feature_bc_matrix"

ensembl_tsv <- read.delim("~/Documents/nf-sc/src/ensembl-ids/ensemble_ids_all.tsv")

mat <- read_feature_barcode_matrix("~/Documents/ms-sc/data/raw/testfbmatrix/outs/raw_feature_bc_matrix")


x <- rownames(mat)

x <- map_ensembl_gene_id(rownames(mat))

write.table(x, file = "ensembl_mappings.tsv",  quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


x <- retrieve_sample_metadata(unique_id = "MS542",
                         id_colname = "individual",
                         samplesheet_path = "~/Documents/ms-sc/refs/sample_metadata.tsv")
sample_metadata.tsv

