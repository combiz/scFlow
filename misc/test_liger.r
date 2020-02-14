library(caret)
library(scflow)
#sce_all <- read_sce("~/Documents/Amy_Glia_Expt/sce") # not working
sce_all <- read_sce("~/Documents/nf-sc/results/celltype_mapped_sce/celltype_mapped_sce")
idx <- as.numeric(createDataPartition(sce_all$manifest, p = .03, list = FALSE)) # 15% subset
sce <- sce_all[, idx]

sce <- sce[, sce$manifest != "kurus"]

sce <- integrate_sce(sce, unique_id_var = "manifest")

# sce <- read_sce("~/Documents/liger_test")

sce <- reduce_dims_sce(sce, input_reduced_dim = c("PCA", "Liger"), unique_id_var = "manifest")

umap_res <- do.call(
  uwot::umap,
  c(
    list(X = reducedDim(sce, "Liger")),
    fargs[names(fargs) %in%
            names(as.list(args(uwot::umap)))]
  )
)
row.names(umap_res) <- colnames(x)
SingleCellExperiment::reducedDim(x, "UMAP_Liger") <- umap_res

plot_reduced_dim(sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger", alpha = 1, size = 1)
plot_reduced_dim(sce, feature_dim = "manifest", reduced_dim = "UMAP_PCA", alpha = 1, size = 1)

plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_PCA", alpha = 1, size = 1)
plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger", alpha = 1, size = 1)

plot_reduced_dim(sce, feature_dim = "manifest", reduced_dim = "tSNE_PCA", alpha = 1, size = 1)
plot_reduced_dim(sce, feature_dim = "manifest", reduced_dim = "tSNE_Liger", alpha = 1, size = 1)

sce_bck <- sce

sce <- cluster_sce(sce, )


reduction_methods = c("tSNE", "UMAP", "UMAP3D")

