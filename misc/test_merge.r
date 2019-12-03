library(scflow)

a <- "~/Documents/nf-sc/results/qc/sce/nivir_sce"
b <- "~/Documents/nf-sc/results/qc/sce/tavij_sce"

sce_l <- lapply(list(a, b), read_sce)
sce <- merge_sce(
  sce_l,
  ensembl_mapping_file = "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv")

class(SingleCellExperiment::counts(sce))

write_sce(sce = sce,
          folder_path = file.path(getwd(), "test_sce"))
library(SingleCellExperiment)
class(counts(sce))

