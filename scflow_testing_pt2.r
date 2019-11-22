library(scflow)
#sce <- read_sce("../junk/subsetsce")
#sce <- read_sce("../junk/mergedsce")
sce <- read_sce("../junk/endosce")
#sce <- reduce_dims_sce(sce, pca_dims = 5)
#x <- reduce_dims_sce(sce, pca_dims = 5)
#sce <- cluster_sce(sce)

#sce <- read_sce("../junk/subsetscewithanno/")


#ctd_path = "~/Documents/nf-sc/refs/ctd/"
#sce <- map_celltypes_sce(sce, ctd_folder = ctd_path)
#table(sce$cluster_celltype)


sce_all <- sce

sce <- sce_all[, sce$cluster_celltype == "Endo"]
write_sce(sce, "../junk/endosce")

de_results <- perform_de(sce)
sce
sce

de_method = "MASTZLM"
min_counts = 1
min_cells_pc = 0.10
rescale_numerics = TRUE
dependent_var = "group"
ref_class = "Control"
confounding_vars = c("individual", "cngeneson", "sex", "age", "PMI", "RIN", "seqdate", "pc_mito")
random_effects_var = NULL
fc_threshold = 1.1
pval_cutoff = 0.05

args <- list()
args$de_method = "MASTZLM"
args$min_counts = 1
args$min_cells_pc = 0.10
args$rescale_numerics = TRUE
args$dependent_var = "group"
args$ref_class = "Control"
args$confounding_vars = c("individual", "cngeneson", "sex", "age", "PMI", "RIN", "seqdate", "pc_mito")
args$random_effects_var = NULL
args$fc_threshold = 1.1
args$pval_cutoff = 0.05
args$re_vars = NULL
