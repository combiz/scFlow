library(scflow)

sce <- read_sce("~/Documents/nf-sc/results/celltype_mapped_sce/celltype_mapped_sce")

sce <- annotate_merged_sce(sce, facet_vars = c("group", "diagnosis", "seqdate"))

saveRDS(sce@metadata, file = "scemetadata_newest.rds")

unique_var <- "individual"

####





plot_reduced_dim(sce, feature_dim = "group", reduced_dim = "UMAP", size = 6, alpha = 1)
plot_reduced_dim(x, feature_dim = "individual", reduced_dim = "UMAP", size = 6, alpha = 1, label_clusters = TRUE)
plot_reduced_dim(x, feature_dim = "individual", reduced_dim = "PCA", size = 6, alpha = 1, label_clusters = TRUE)
x@metadata$reduced_dim_plots$umap3d_plot_ly


unique_id_var <- "individual"

plot_reduced_dim(sce, reduced_dim = "PCA_PB", feature_dim = "individual", label_clusters = TRUE)
plot_reduced_dim(sce, reduced_dim = "UMAP_PB", feature_dim = "individual", label_clusters = TRUE)
plot_reduced_dim(sce, reduced_dim = "UMAP3D_PB", feature_dim = "individual", label_clusters = TRUE)
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
