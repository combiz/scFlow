library(scFlow)
library(ggplot2)

##
sce <- read_sce("~/Documents/junk/MS_Custom_Mapped_SCE")

res <- model_celltype_freqs(sce)

butcher::weigh(res$dirichlet_plot)

res <- model_celltype_freqs(sce, var_order = c("Control", "Low", "High"), palette = c("#3B4992FF", "#E64B35FF","#EE0000FF"))
report_celltype_model(res)



sce <- read_sce("~/Documents/junk/celltype_mapped_sce/", read_metadata = FALSE)
plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger", label_clusters = TRUE, size = .5)
plot_reduced_dim(sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = TRUE, size = .5)

ctm <- read_celltype_mappings("~/Documents/junk/celltype_mappings.tsv")
sce$clusters <- as.numeric(as.character(sce$clusters))
sce <- map_custom_celltypes(sce, ctm, cols = "cluster_celltype")

sce$clusters <- as.character(as.numeric(sce$clusters))
plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger", size = 0.5, label_clusters = TRUE)
sce_all <- sce
sce <- sce[, sce$cluster_celltype != "Doublets"]
#write_sce(sce, "~/Documents/junk/MS_Custom_Mapped_SCE")



###

sce <- annotate_celltype_metrics(sce, facet_vars = c("manifest", "group", "sex", "diagnosis", "seqdate", "capdate", "prepdate"),input_reduced_dim = "UMAP_Liger")
sce <- annotate_celltype_metrics(sce, facet_vars = c("manifest", "group", "sex", "diagnosis", "seqdate", "capdate", "prepdate"),input_reduced_dim = "UMAP_Liger")

lapply(sce@metadata$celltype_annotations$reddim_plots$cluster_var, function(x) format(utils::object.size(x), units = "auto"))

format(utils::object.size(sce@metadata$celltype_annotations$reddim_plots$cluster_var), units = "auto")

res <- model_celltype_freqs(sce)
report_celltype_model(res)
#qs::qsave(sce@metadata$celltype_annotations$reddim_plots$cluster_var, "test.qs")

