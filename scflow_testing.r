library(scflow)

matpath <- "~/Documents/ms-sc/data/raw/testfbmatrix/outs/raw_feature_bc_matrix"

#ensembl_tsv <- read.delim("~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv")

mat <- read_feature_barcode_matrix(matpath)

ss_classes <- c(
  batch = "factor",
  capdate = "factor",
  prepdate = "factor",
  seqdate = "factor",
  aplevel = "factor"
)

metadata <- retrieve_sample_metadata(
  unique_id = "MS542",
  id_colname = "individual",
  samplesheet_path = "~/Documents/ms-sc/refs/sample_metadata.tsv",
  colClasses = ss_classes
)

sce <- generate_sce(mat, metadata)

#sce <- annotate_sce(sce)

sce <- annotate_sce(
  sce,
  ensembl_mapping_file = "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv"
)

sce <- .qc_plot_count_depth_distribution(sce)
sce <- .qc_plot_count_depth_histogram(sce)
sce <- .qc_plot_features_histogram(sce)
sce <- .qc_plot_features_vs_count_depth(sce)
sce <- .qc_plot_mito_fraction_histogram(sce)
sce <- .qc_plot_ribo_fraction_histogram(sce)

# DO QC PLOTS AND TABLE HERE!

sce <- filter_sce(
  sce,
  filter_genes = TRUE, filter_cells = TRUE, drop_unmapped = TRUE, drop_mito = TRUE, drop_ribo = FALSE)

sce <- find_singlets(sce, "doubletfinder")

sce <- sce[, sce$is_singlet == TRUE]

write_sce(sce[1:5000, ], "../junk/a")
write_sce(sce[3000:7000, ], "../junk/b")
write_sce(sce[5000:10000, ], "../junk/c")

####
fp_l <- c("../junk/a", "../junk/b", "../junk/c")

x <- merge_sce(
  fp_l,
  ensembl_mapping_file = "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv"
)

##################

library(tidyr)
library(dplyr)
library(purrr)

dt <- as.data.frame(sce@colData) %>%
  select(
    colnames(sce@colData)[startsWith(colnames(sce@colData), "qc_metric")],
    pc_mito,
    pc_ribo,
    barcode)

dt <- dt %>% tidyr::pivot_longer(c(contains("qc_metric_"), starts_with("pc_")))

ggplot(dt)+
  geom_tile(aes(x = barcode, y = name, fill = value)) +
  scale_fill_gradientn(colours = c("red", "blue"))+
  theme_void()

library(ggplot2)
ggplot(data = dt[dt$name == "qc_metric_pc_mito_ok",])+
  geom_tile(aes(x = barcode, y = name, fill = value)) +
  scale_fill_gradient(low="grey60", high="red") +
  labs(x="", y="") +
  theme_bw()

table(x$doublet_finder_annotation)

df <- data.frame(SingleCellExperiment::reducedDim(x, "seurat_umap_by_individual"))
df$is_singlet <- x$is_singlet

ggplot(data = df)+
  geom_point(aes(x = UMAP_1, y = UMAP_2, colour = is_singlet))

write_sce(sce, file.path(getwd(), "junk"))

write_feature_barcode_matrix(SingleCellExperiment::counts(sce), file.path(getwd(), "junk"))


library(ggplot2)
library(dplyr)

ggarrange(
  sce@metadata$qc_plots$count_depth_distribution,
  sce@metadata$qc_plots$count_depth_histogram,
  sce@metadata$qc_plots$number_genes_histogram,
  sce@metadata$qc_plots$number_genes_vs_count_depth,
  sce@metadata$qc_plots$mito_fraction_histogram,
  sce@metadata$qc_plots$ribo_fraction_histogram
)

