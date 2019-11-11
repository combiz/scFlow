################################################################################
#' Cluster SingleCellExperiment with monocle3::cluster_cells
#'
#' @param sce a SingleCellExperiment object
#' @param reduction_methods one or more of "PCA", "tSNE", "UMAP", "UMAP3D"
#' @param pca_dims the number of pca dimensions used
#' @param ... see uwot::umap for umap options
#'
#' @return sce a SingleCellExperiment object annotated with reducedDims
#'
#' @family clustering and dimensionality reduction
#' @import cli monocle3
#' @importFrom SingleCellExperiment reducedDim reducedDims
#' @importFrom assertthat assert_that
#'
#' @export

cluster_sce <- function(sce, ...) {

  fargs <- list(
    cluster_method = "louvain",
    reduction_method = "UMAP",
    res = 1e-5,
    k = 100,
    louvain_iter = 1,
    verbose = T
  )

  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  assertthat::assert_that(
    fargs$reduction_method %in% names(SingleCellExperiment::reducedDims(sce))
  )

  cds <- .sce_to_cds(sce)
  cds <- do.call(
    monocle3::cluster_cells,
    c(list(cds = cds), fargs)
  )
  sce$clusters <-
    as.factor(cds@clusters[[fargs$reduction_method]]$clusters)
  sce$partitions <-
    as.factor(cds@clusters[[fargs$reduction_method]]$partitions)

  sce@metadata$cluster_params <- fargs

  return(sce)

}
