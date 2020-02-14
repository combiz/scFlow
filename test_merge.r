library(scflow)

a <- "~/Documents/nf-sc/results/qc/sce/nivir_sce"
b <- "~/Documents/nf-sc/results/qc/sce/tavij_sce"

sce_l <- lapply(list(a, b), read_sce)
sce <- merge_sce(
  sce_l,
  ensembl_mapping_file = "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv")

x <- annotate_sce(
  sce,
  annotate_cells = FALSE,
  ensembl_mapping_file = "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv"
)

keep_idx <- !(
  startsWith(colnames(colData(sce)), "qc_metric") |
  startsWith(colnames(colData(sce)), "pc_")
  )

colData(sce) <- colData(sce)[, keep_idx]
x <- annotate_sce(sce, ensembl_mapping_file = "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv")
sce$
