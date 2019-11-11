library(scflow)
library(monocle3)
library(uwot)
sce <- read_sce("../junk/subsetsce")

sce <-reduce_dims_sce(sce, pca_dims = 5)

colData(cds) <- colData(cds)[, 1:10]
x <- cluster_cells(cds)

