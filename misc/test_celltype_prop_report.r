
library(magrittr)
sce$age <- as.factor(sce$age)
sce$Cluster <- as.factor(sce$Cluster)
x <- annotate_celltype_plots(sce, cluster_var = "Cluster", unique_id_var = "individual", facet_var = c("age", "sex", "group"), metric_vars = c("pc_mito", "total_counts", "total_features_by_counts"))
