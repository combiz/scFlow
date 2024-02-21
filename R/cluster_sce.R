################################################################################
#' Cluster SingleCellExperiment with monocle3::cluster_cells
#'
#' @param sce a SingleCellExperiment object
#' @param cluster_method Available methods include leiden and louvain.
#' @param reduction_method Input slot for clustering. Available options are PCA,
#' UMAP_Liger, tSNE_Liger.
#' @param resolution Clustering resolution. If NULL, clustering method will be
#' set to louvain.
#' @param k Integer number of nearest neighbors to use when creating the k
#' nearest neighbor graph for Louvain/Leiden clustering.
#' k is related to the resolution of the clustering result, a bigger k will
#' result in lower resolution and vice versa. Default is 50.
#' @param louvain_iter Integer number of iterations used for Louvain clustering.
#' The clustering result giving the largest modularity score will be used
#' as the final clustering result. Default is 1. Note that if num_iter
#' is greater than 1, the random_seed argument will be ignored
#' for the louvain method.
#' @param verbose A logic flag to determine whether or not we should print
#' the run details.
#' @param ... see monocle3::cluster_cells for more clustering options.
#'
#' @return sce a SingleCellExperiment object annotated with reducedDims
#'
#' @family clustering and dimensionality reduction
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom monocle3 cluster_cells
#' @importFrom assertthat assert_that
#'
#' @export

cluster_sce <- function(sce,
                        cluster_method = "louvain",
                        reduction_method = "UMAP_Liger",
                        resolution = NULL,
                        k = 50,
                        louvain_iter = 1,
                        verbose = T,
                        ...) {

  fargs <- list(
    cluster_method = cluster_method,
    reduction_method = reduction_method,
    resolution = resolution,
    k = k,
    louvain_iter = louvain_iter,
    verbose = verbose
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
