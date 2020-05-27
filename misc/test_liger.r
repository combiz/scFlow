library(caret)
library(SingleCellExperiment)
library(scFlow)
#sce_all <- read_sce("~/Documents/Amy_Glia_Expt/sce") # not working
sce_all <- read_sce("~/Documents/nf-sc/results/celltype_mapped_sce/celltype_mapped_sce")
#idx <- as.numeric(caret::createDataPartition(sce_all$manifest, p = .03, list = FALSE)) # 15% subset
idx <- as.numeric(caret::createDataPartition(sce_all$individual, p = .03, list = FALSE)) # 15% subset
mini_sce <- sce_all[, idx]
sce <- sce_all[, idx]

sce <- sce[, sce$manifest != "kurus"]

sce <- integrate_sce(sce, unique_id_var = "manifest", k = 20)
sce <- integrate_sce(sce, unique_id_var = "individual", k = 20)

write_sce(sce, "~/Documents/junk/mini_sce_liger")

# sce <- read_sce("~/Documents/liger_test")

sce <- reduce_dims_sce(sce, input_reduced_dim = c("PCA", "Liger"), unique_id_var = "manifest")
sce <- reduce_dims_sce(sce, input_reduced_dim = c("PCA", "Liger"), unique_id_var = "individual")

sce <- cluster_sce(sce)

sce <- annotate_integrated_sce(sce, categorical_covariates = c("group", "sex"))
sce <- annotate_integrated_sce(sce, categorical_covariates = c("sex"))

report_integrated_sce(sce)

#write_sce(sce, "~/Documents/liger_test2")

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

sce <- cluster_sce(sce, reduction_method = "UMAP_Liger", resolution = 0.01)
plot_reduced_dim(sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", alpha = 1, size = 1)
plot_reduced_dim(sce, feature_dim = "clusters", reduced_dim = "UMAP_PCA", alpha = 1, size = 1)


plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_PCA", alpha = 1, size = 1)
plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger", alpha = 1, size = 1)

plot_reduced_dim(sce, feature_dim = "manifest", reduced_dim = "tSNE_PCA", alpha = 1, size = 1)
plot_reduced_dim(sce, feature_dim = "clusters", reduced_dim = "tSNE_Liger", alpha = 1, size = 1)

plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "tSNE_Liger", alpha = 1, size = 1)

sce_bck <- sce

sce <- cluster_sce(sce, reduction_method = "UMAP_Liger")
sce <- cluster_sce(sce, reduction_method = "tSNE_Liger", resolution = 0.01, k = 20)

sce <- map_celltypes_sce(sce, ctd_folder = ctd_fp)
reduction_methods = c("tSNE", "UMAP", "UMAP3D")

