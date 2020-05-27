library(scFlow)
sce <- read_sce("~/Documents/junk/mergedsce")
sce_all <- sce
sce_subset <- sce[, sce$cluster_celltype == "Micro"]


sce_pb <- pseudobulk_sce(
  sce_subset,
  keep_vars = c("individual", "group", "sex", "age", "PMI", "RIN", "seqdate"),
  assay_name = "counts",
  celltype_var = "cluster_celltype",
  sample_var = "individual"
  )

de_results <- perform_de(
  sce_pb,
  confounding_vars = c("cngeneson", "sex", "age", "PMI", "RIN", "seqdate")
)

report_de(de_results[[1]])
