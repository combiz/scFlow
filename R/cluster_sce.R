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
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom monocle3 cluster_cells
#' @importFrom assertthat assert_that
#'
#' @export

cluster_sce <- function(sce, ...) {

  fargs <- list(
    cluster_method = "leiden",
    reduction_method = "UMAP_Liger",
    resolution = 1e-5,
    k = 100,
    louvain_iter = 1,
    verbose = T
  )

  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  if(grepl("UMAP|tSNE", fargs$reduction_method)) {
    fargs[["nn_control_metric"]] <- "euclidean"
  } else {
    fargs[["nn_control_metric"]] <- "cosine"
  }

  rd_method_bck <- NULL

  assertthat::assert_that(
    fargs$reduction_method %in% SingleCellExperiment::reducedDimNames(sce)
  )

  cds <- .sce_to_cds(sce)

  if(!(fargs$reduction_method %in% c("UMAP", "tSNE", "PCA", "LSI",
                                     "Aligned"))) {
    SingleCellExperiment::reducedDim(cds, "PCA") <-
      SingleCellExperiment::reducedDim(cds, fargs$reduction_method)
    rd_method_bck <- fargs$reduction_method
    fargs$reduction_method <- "PCA"
  }

  cds <- monocle3::cluster_cells(cds = cds,
                                 cluster_method = fargs$cluster_method,
                                 reduction_method = fargs$reduction_method,
                                 resolution = fargs$resolution,
                                 k = fargs$k,
                                 louvain_iter = fargs$louvain_iter,
                                 verbose = fargs$verbose,
                                 nn_control = list(metric = fargs$nn_control_metric)
                                 )

  sce$clusters <-
    as.factor(cds@clusters[[fargs$reduction_method]]$clusters)
  sce$partitions <-
    as.factor(cds@clusters[[fargs$reduction_method]]$partitions)

  if (!is.null(rd_method_bck)) { fargs$reduction_method <- rd_method_bck }
  sce@metadata$cluster_params <- fargs

  return(sce)

}
