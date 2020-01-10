library(scflow)

sce <- read_sce("~/Documents/nf-sc/results/celltype_mapped_sce/celltype_mapped_sce/")

unique_var <- "individual"

####





plot_umap_with_feature(x, feature_dim = "group", reduced_dim = "UMAP", size = 6, alpha = 1)
plot_umap_with_feature(x, feature_dim = "individual", reduced_dim = "UMAP", size = 6, alpha = 1, label_clusters = TRUE)
plot_umap_with_feature(x, feature_dim = "individual", reduced_dim = "PCA", size = 6, alpha = 1, label_clusters = TRUE)
x@metadata$reduced_dim_plots$umap3d_plot_ly


unique_id_var <- "individual"

plot_umap_with_feature(sce, reduced_dim = "PCA_PB", feature_dim = "individual", label_clusters = TRUE)
plot_umap_with_feature(sce, reduced_dim = "UMAP_PB", feature_dim = "individual", label_clusters = TRUE)
plot_umap_with_feature(sce, reduced_dim = "UMAP3D_PB", feature_dim = "individual", label_clusters = TRUE)
library(tidyr)
#mydata <- t(table(x$sex, x$group))
mydata <- as.data.frame(table(sex = x$sex))
library(plotly)
library(scales)

#plot_vars <- c("total_features_by_counts", "total_counts", "pc_mito", "pc_ribo")
plot_vars <- c("age")
plot_vars <- unique(
  c("total_features_by_counts", "total_counts", "pc_mito", "pc_ribo",
    plot_vars)
)
facet_vars <- c("group", "seqdate")
usv <- "manifest"
