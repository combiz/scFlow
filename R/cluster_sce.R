################################################################################
#' Cluster SingleCellExperiment with monocle3::cluster_cells
#'
#' @param sce a SingleCellExperiment object
#' @param ... see uwot::umap for umap options. Includes reduction_methods one 
#' or more of "PCA", "tSNE", "UMAP", "UMAP3D"
#'
#' @return sce a SingleCellExperiment object annotated with reducedDims
#'
#' @family clustering and dimensionality reduction
#' @importFrom SingleCellExperiment reducedDim reducedDims
#' @importFrom monocle3 cluster_cells
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

  rd_method_bck <- NULL

  assertthat::assert_that(
    fargs$reduction_method %in% names(SingleCellExperiment::reducedDims(sce))
  )

  cds <- .sce_to_cds(sce)

  if(!(fargs$reduction_method %in% c("UMAP", "tSNE", "PCA", "LSI",
                                     "Aligned"))) {
    SingleCellExperiment::reducedDim(cds, "PCA") <-
      SingleCellExperiment::reducedDim(cds, fargs$reduction_method)
    rd_method_bck <- fargs$reduction_method
    fargs$reduction_method <- "PCA"
  }

  cds <- do.call(
    monocle3::cluster_cells,
    c(list(cds = cds), fargs)
  )
  sce$clusters <-
    as.factor(cds@clusters[[fargs$reduction_method]]$clusters)
  sce$partitions <-
    as.factor(cds@clusters[[fargs$reduction_method]]$partitions)

  if (!is.null(rd_method_bck)) { fargs$reduction_method <- rd_method_bck }
  sce@metadata$cluster_params <- fargs

  return(sce)

}
