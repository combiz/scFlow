library(scFlow)
library(SingleCellExperiment)

#sce_l <- list.dirs("~/Documents/scFlowExample/tmp_Zeisel2015_scflow")[-1]
a <- "~/Documents/junk/scFlowExampleSCE/bavaf_sce/"
b <- "~/Documents/junk/scFlowExampleSCE/majos_sce/"
c <- "~/Documents/junk/scFlowExampleSCE/nihoj_sce/"
#a <- "~/Documents/nf-sc/results/qc/sce/nivir_sce"
#b <- "~/Documents/nf-sc/results/qc/sce/tavij_sce"
#c <- "~/Documents/nf-sc/results/qc/sce/gatud_sce/"
sce_l <- lapply(list(a, b, c), read_sce)

sce <- merge_sce(
  sce_l,
  ensembl_mapping_file = "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv")

sce <- annotate_merged_sce(sce, unique_id_var = "manifest")
report_merged_sce(sce)

devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")


plot_vars = c("total_features_by_counts",
              "total_counts", "pc_mito",
              "pc_ribo")
unique_id_var = "manifest"
facet_vars = NULL
outlier_vars = c("total_features_by_counts",
                 "total_counts")
outlier_mads = 3
