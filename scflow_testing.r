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

sce <- annotate_sce(
  sce,
  ensembl_mapping_file = "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv"
)

#sce <- filter_sce(
#  sce,
#  filter_genes = TRUE, filter_cells = TRUE, drop_unmapped = TRUE, drop_mito = TRUE, drop_ribo = FALSE)

sce <- find_singlets(sce, "doubletfinder", pK = 0.005)

sce <- report_qc_sce(sce)

start_time <- Sys.time()
metadata_tmp_path <- file.path(tempdir(), "metadata.rds")
saveRDS(sce@metadata, metadata_tmp_path)
end_time <- Sys.time()
end_time - start_time

# after dev
system.file(
  "rmarkdown/templates/quality-control/skeleton/skeleton.Rmd",
  package = "scflow")

rmarkdown::render(
  file.path(getwd(), "inst/rmarkdown/templates/quality-control/skeleton/skeleton.Rmd"),
  params = list(
    metadata_path = metadata_tmp_path
  ),
  output_dir = getwd(),
  output_file = "qc"
)



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

table(sce$is_singlet)



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

library(ggpubr)
x <- ggplot(data.frame(pc_mito = sce[, sce$total_counts > 10]$pc_mito)) + geom_jitter(aes(x = "", y = pc_mito), size = 0.1, alpha = 0.1) + theme_bw() +     theme(panel.grid.major = element_blank(),
                                                                                                                                                                  panel.grid.minor = element_blank(),
                                                                                                                                                                  panel.background = element_blank(), axis.title.x = element_blank(),
                                                                                                                                                                  axis.text.x = element_blank(),
                                                                                                                                                                  axis.ticks.x = element_blank())

#y <- ggplot(data.frame(pc_mito = sce[, sce$total_counts > 10]$pc_mito)) + geom_density(aes(pc_mito, stat(count), fill = "sample")) + theme_bw()


ggarrange(x  + coord_flip(), metadata$qc_plots$mito_fraction_histogram, ncol = 1, heights = c(1,2))


