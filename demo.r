# QC demo
# 27/11/19

##  ............................................................................
##  Initialize                                                              ####

library(scflow)

# v2 chemistry
matpath <- "~/Documents/ms-sc/data/raw/testfbmatrix/outs/raw_feature_bc_matrix"
# v3 chemistry, enriched
#matpath <- "~/Documents/testmatrices/enriched"

ensembl_fp <- "~/Documents/ms-sc/src/ensembl-ids/ensembl_mappings.tsv"
samplesheet_fp <- "~/Documents/nfl/refs/SampleSheet.tsv"
ctd_fp <- "~/Documents/nf-sc/refs/ctd/"

##  ............................................................................
##  Start QC                                                                ####

x <- Sys.time()

mat <- read_sparse_matrix(matpath)

metadata <- read_metadata(
  unique_key = "hajov",
  key_colname = "manifest",
  samplesheet_path = samplesheet_fp)

sce <- generate_sce(mat, metadata)

sce <- annotate_sce(sce, ensembl_mapping_file = ensembl_fp)

sce <- filter_sce(sce)

sce <- find_singlets(sce, "doubletfinder", pK = 0.005, vars_to_regress_out = NULL)

sce <- filter_sce(sce)

report_qc_sce(sce)

totaltime <- Sys.time() - x

print(totaltime)
alarm()

##  ............................................................................
##  Bonus - cluster and identify celltypes                                  ####

x <- Sys.time()

sce <- reduce_dims_sce(sce, pca_dims = 5)

sce <- cluster_sce(sce)

sce <- map_celltypes_sce(sce, ctd_folder = ctd_fp)

totaltime <- Sys.time() - x

print(totaltime)
alarm()

##

plot_umap_with_feature(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP", label_clusters = TRUE, size = 1)
plot_umap_with_feature(sce, feature_dim = "cluster_celltype", reduced_dim = "tSNE", label_clusters = TRUE, size = 1)
plot_umap_with_feature(sce, feature_dim = "cluster_celltype", reduced_dim = "PCA", label_clusters = TRUE, size = 1)



# 3D


umap_res <- SingleCellExperiment::reducedDim(sce, "UMAP3D")

pal <- colorFactor(palette = "Accent", domain = sce$cluster_celltype.y)
celltype_pal <- pal(sce$cluster_celltype.y)

threejs::scatterplot3js(
  x = umap_res[, 1],
  y = umap_res[, 2],
  z = umap_res[, 3],
  color = celltype_pal,
  axis = FALSE,
  num.ticks = NULL,
  stroke = "#4682b4",
  size = .01
)


plotly::plot_ly(
  x = umap_res[, 1],
  y = umap_res[, 2],
  z = umap_res[, 3],
  color = celltype_pal,
  name = sce$cluster_celltype.y,
  type = "scatter3d",
  mode = "markers",
  size = 0.1
)
