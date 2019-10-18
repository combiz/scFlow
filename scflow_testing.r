library(scflow)

matpath <- "~/Documents/ms-sc/data/raw/testfbmatrix/outs/raw_feature_bc_matrix"

#ensembl_tsv <- read.delim("~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv")

mat <- read_feature_barcode_matrix(matpath)

metadata <- retrieve_sample_metadata(unique_id = "MS542",
                                     id_colname = "individual",
                                     samplesheet_path = "~/Documents/ms-sc/refs/sample_metadata.tsv")

sce <- generate_sce(mat, metadata)

sce <- annotate_sce_genes(sce, "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv")

sce <- annotate_sce_cells(sce)



