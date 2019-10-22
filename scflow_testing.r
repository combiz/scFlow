library(scflow)

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

metadata <- retrieve_sample_metadata(unique_id = "MS542",
                                     id_colname = "individual",
                                     samplesheet_path = "~/Documents/ms-sc/refs/sample_metadata.tsv",
                                     colClasses = ss_classes)

sce <- generate_sce(mat, metadata)

sce <- annotate_sce(
  sce,
  ensembl_mapping_file = "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv"
)

# DO QC PLOTS AND TABLE HERE!

x <- filter_sce(
  sce,
  filter_genes = TRUE, filter_cells = TRUE, drop_unmapped = TRUE, drop_mito = TRUE, drop_ribo = FALSE)

sce <- find_singlets(sce, "doubletfinder")

#sce <- annotate_sce_genes(sce, "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv")

#sce <- annotate_sce_cells(sce)

sce <- filter_sce(
  sce,
  filter_genes = TRUE,
  filter_cells = TRUE,
  keep_mito = FALSE,
  keep_ribo = TRUE
)


table(x$doublet_finder_annotation)

df <- data.frame(reducedDim(x, "seurat_umap_by_individual"))
df$is_singlet <- x$is_singlet

ggplot(data = df)+
  geom_point(aes(x = UMAP_1, y = UMAP_2, colour = is_singlet))


write_sce(sce, file.path(getwd(), "junk"))

write_feature_barcode_matrix(SingleCellExperiment::counts(sce), file.path(getwd(), "junk"))


