library(scflow)
#sce <- read_sce("../junk/subsetsce")
sce <- reduce_dims_sce(sce, pca_dims = 5)
sce <- cluster_sce(sce)

sce <- read_sce("../junk/subsetscewithanno/")


ctd_path = "~/Documents/nf-sc/refs/ctd/"
sce <- map_celltypes_sce(sce, ctd_folder = ctd_path)
table(sce$cluster_celltype)

