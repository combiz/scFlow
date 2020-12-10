library(scFlow)
library(ggplot2)

##
#sce <- read_sce("~/Documents/junk/celltype_mapped_sce/")
#sce <- read_sce("~/Documents/junk/MS_Custom_Mapped_SCE")
sce <- read_sce("~/Documents/junk/ms26072020/")

## old vs new
old <- read.delim("~/Documents/ms-sc/data/tidy/sce-fullcoldata.tsv")
old$barcode <- paste(old$individual, old$group, old$barcode, sep = "_")
old$in_new <- old$barcode %in% sce$barcode
old_dt <- data.frame(barcode = as.character(old$barcode), cluster_celltype_old = as.character(old$cluster_celltype))
new_dt <- data.frame(barcode = as.character(sce$barcode), cluster_celltype_new = as.character(sce$cluster_celltype))
dt <- dplyr::left_join(new_dt, old_dt, by = "barcode")
dt$cluster_celltype_old[is.na(dt$cluster_celltype_old)] <- "NA"
sce$cluster_celltype_old <- dt$cluster_celltype_old

# how many cells from each celltype were dropped by the new qc etc
sce$cluster_celltype_old <- dt$cluster_celltype_old
df <- as.data.frame.matrix(table(old_dt$cluster_celltype_old, old$in_new))
df$cluster_celltype <- rownames(df)
df$retained_pc <- (df$`TRUE`/(df$`FALSE` + df$`TRUE`))
df <- df[order(df$retained_pc),]
df$cluster_celltype <- factor(df$cluster_celltype, levels = df$cluster_celltype)
ggplot(df, aes(x = cluster_celltype, y = retained_pc))+
  geom_bar(stat = "identity")+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    axis.text = element_text(size = 16)
  )
###
library(dplyr)
freqs <- as.data.frame.matrix(table(sce$clusters, sce$cluster_celltype_old))
freqs$clusters <- row.names(freqs)
freqs$old_name <- colnames(freqs)[apply(freqs,1,which.max)]
new_cluster_celltype <- freqs %>%
  select(clusters, old_name) %>%
  arrange(as.numeric(as.character(clusters)))

x <- map_custom_celltypes(sce, new_cluster_celltype)

write.table(new_cluster_celltype, "~/junk/new_celltype_mappings.tsv", row.names = F, col.names = T, sep = "\t")

plot_reduced_dim(x, feature_dim = "old_name", reduced_dim = "UMAP_Liger", size = .5, label_clusters = TRUE)

setdiff(unique(new_cluster_celltype$old_name), unique(sce$cluster_celltype_old))
###



plot_reduced_dim(sce, feature_dim = "cluster_celltype_old", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim(sce, feature_dim = "cluster_celltype_old", reduced_dim = "UMAP_Liger", size = .5, label_clusters = TRUE)
plot_reduced_dim(sce, feature_dim = "cluster_celltype_old", reduced_dim = "UMAP_Liger", size = .5, label_clusters = TRUE)
plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger", size = .5, label_clusters = TRUE)
plot_reduced_dim(sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", size = .5, label_clusters = TRUE)

plot_reduced_dim(sce, feature_dim = "cluster_celltype_old", reduced_dim = "UMAP_Liger", highlight_feature = "B-cells", size = .1)

plot_reduced_dim_gene(sce, gene = "SST", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "PVALB", reduced_dim = "UMAP_Liger", size = .5)

plot_reduced_dim_gene(sce, gene = "SV2C", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "VIP", reduced_dim = "UMAP_Liger", size = .5)

plot_reduced_dim_gene(sce, gene = "PLP1", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "CD74", reduced_dim = "UMAP_Liger", size = .5)

plot_reduced_dim_gene(sce, gene = "DRD1", reduced_dim = "UMAP_Liger", size = .5)

plot_reduced_dim_gene(sce, gene = "CUX1", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "PLXDC2", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "RBMS3", reduced_dim = "UMAP_Liger", size = .5) #stromal
plot_reduced_dim_gene(sce, gene = "B2M", reduced_dim = "UMAP_Liger", size = .5) #b
plot_reduced_dim_gene(sce, gene = "SSR4", reduced_dim = "UMAP_Liger", size = .5) #b
plot_reduced_dim_gene(sce, gene = "PVALB", reduced_dim = "UMAP_Liger", size = .5) #b
plot_reduced_dim_gene(sce, gene = "VIP", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "SV2C", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "SLC17A7", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "SYT1", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "GAD2", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "IL26", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "IL26", reduced_dim = "UMAP_Liger", size = .5)

plot_reduced_dim_gene(sce, gene = "ETV1", reduced_dim = "UMAP_Liger", size = .5)

plot_reduced_dim_gene(sce, gene = "CD133", reduced_dim = "UMAP_Liger", size = .5) # EPendymal
plot_reduced_dim_gene(sce, gene = "ETV4", reduced_dim = "UMAP_Liger", size = .5)
plot_reduced_dim_gene(sce, gene = "NRG1", reduced_dim = "UMAP_Liger", size = .5)

table(sce$cluster_celltype_old)


sce_old <- read_sce("~/Documents/ms-sc/data/tidy/sce/")

res <- model_celltype_freqs(sce)

butcher::weigh(res$dirichlet_plot)

res <- model_celltype_freqs(sce, var_order = c("Control", "Low", "High"), palette = c("#3B4992FF", "#E64B35FF","#EE0000FF"))
report_celltype_model(res)



sce <- read_sce("~/Documents/junk/celltype_mapped_sce/", read_metadata = FALSE)
plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger", label_clusters = TRUE, size = .5)
plot_reduced_dim(sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = TRUE, size = .5)

ctm <- read_celltype_mappings("~/junk/new_celltype_mappings.tsv")
#sce$clusters <- as.numeric(as.character(sce$clusters))
sce <- map_custom_celltypes(sce, ctm, cols = "new_cluster_celltype")
table(sce$cluster_celltype)
sce$cluster_celltype <- NULL
sce$cluster_celltype <- sce$new_cluster_celltype
sce$clusters <- as.character(as.numeric(sce$clusters))
plot_reduced_dim(sce, feature_dim = "new_cluster_celltype", reduced_dim = "UMAP_Liger", size = 0.5, label_clusters = TRUE)
sce_all <- sce
sce <- sce[, sce$cluster_celltype != "Doublets"]
write_sce(sce, "~/Documents/junk/MS_Custom_Mapped_SCE")

sce$cluster_celltype <- sce$new_cluster_celltype

###
sce$clusters <- as.character(as.numeric(sce$clusters))
sce <- annotate_celltype_metrics(sce, facet_vars = c("manifest", "group", "sex", "diagnosis", "seqdate", "capdate", "prepdate"),input_reduced_dim = "UMAP_Liger")


lapply(sce@metadata$celltype_annotations$reddim_plots$cluster_var, function(x) format(utils::object.size(x), units = "auto"))

format(utils::object.size(sce@metadata$celltype_annotations$reddim_plots$cluster_var), units = "auto")

res <- model_celltype_freqs(sce, unique_id_var = "manifest", celltype_var = "cluster_celltype", dependent_var = "group", ref_class = "Control", var_order = c("Control", "Low", "High"), palette = c("#3B4992FF", "#E64B35FF","#EE0000FF"))
res <- model_celltype_freqs(sce, unique_id_var = "manifest", celltype_var = "cluster_celltype", dependent_var = "diagnosis", ref_class = "Control", var_order = c("Control", "MS"), palette = c("#3B4992FF", "#EE0000FF"))
report_celltype_model(res)
#qs::qsave(sce@metadata$celltype_annotations$reddim_plots$cluster_var, "test.qs")



###############
ensembl_fp <- "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv"
samplesheet_fp <- "~/Documents/junk/refs/SampleSheet.tsv"
ctd_fp <- "~/Documents/junk/refs/ctd/"
sce_all <- sce
idx <- as.numeric(caret::createDataPartition(sce$individual, p = .10, list = FALSE)) # 15% subset
sce <- sce[, idx]

sce <- map_celltypes_sce(sce,ctd_folder = ctd_fp, clusters_colname = "clusters")
class(sce$clusters)
