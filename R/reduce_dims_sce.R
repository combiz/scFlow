################################################################################
#' Calculate dimensionality reductions for a SingleCellExperiment object
#' or merged SingleCellExperiment objects using PCA, tSNE, UMAP, UMAP3D, Liger
#'
#' @param sce a SingleCellExperiment object or merged sce objects
#' @param reduction_methods one or more of "PCA","tSNE","UMAP","UMAP3D","Liger"
#' @param pca_dims the number of pca dimensions used
#' @param ... see uwot::umap for umap options
#'
#' @return sce SingleCellExperiment object annotated with reducedDims
#'
#' @family clustering and dimensionality reduction
#' @importFrom monocle3 preprocess_cds
#' @importFrom SingleCellExperiment reducedDim reducedDims
#' @importFrom RcppParallel defaultNumThreads
#' @importFrom purrr map_lgl
#' @importFrom threejs scatterplot3js
#' @importFrom plotly plot_ly
#'
#' @export

reduce_dims_sce <- function(sce,
                            reduction_methods = c(
                              "PCA", "tSNE", "UMAP", "UMAP3D", "Liger"
                            ),
                            vars_to_regress_out = NULL,
                            pca_dims = 20,
                            ...) {
  fargs <- list(
    n_neighbors = 30L,
    n_components = 2L,
    init = "pca",
    metric = "euclidean",
    n_epochs = NULL,
    learning_rate = 1.0,
    min_dist = 0.3,
    spread = 1.0,
    set_op_mix_ratio = 1.0,
    local_connectivity = 1L,
    repulsion_strength = 1,
    negative_sample_rate = 5,
    fast_sgd = FALSE,
    n_threads = max(1, RcppParallel::defaultNumThreads() - 2)
  )
  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  sce@metadata$reduced_dim_plots <- list()

  if (!all(purrr::map_lgl(
    reduction_methods,
    ~ . %in% c("PCA", "tSNE", "UMAP", "UMAP3D", "Liger")
  ))) {
    stop("reduction methods must be from: PCA, tSNE, UMAP, UMAP3D, Liger")
  }

  # Generate res mod formula for preprocess
  if (!is.null(vars_to_regress_out)) {
    res_mod_formula_str <- paste0(
      "~", paste0(vars_to_regress_out, collapse = "+")
    )
  }

  # Preprocess cds

  cds <- .sce_to_cds(sce)

  message(sprintf("Reducing Dimensions with PCA"))
  if (!is.null(vars_to_regress_out)) {
    message(sprintf("Regressing out: %s", res_mod_formula_str))
  }
  cds <- monocle3::preprocess_cds(
    cds,
    num_dim = pca_dims,
    residual_model_formula_str = res_mod_formula_str
  )


  SingleCellExperiment::reducedDim(sce, "PCA") <-
    SingleCellExperiment::reducedDim(cds, "PCA")

  for (reddim_method in reduction_methods) {
    message(sprintf("Reducing Dimensions with %s", reddim_method))

    # Reduce dimensions with tSNE
    if (reddim_method == "tSNE") {
      cds <- monocle3::reduce_dimension(
        cds,
        preprocess_method = "PCA",
        reduction_method = "tSNE"
      )
      SingleCellExperiment::reducedDim(sce, "tSNE") <-
        SingleCellExperiment::reducedDim(cds, "tSNE")
    }

    # Reduce dimensions with UMAP
    if (reddim_method %in% c("UMAP", "UMAP3D")) {
      mat <- as.matrix(SingleCellExperiment::reducedDim(cds, "PCA"))

      if (reddim_method == "UMAP3D") {
        fargs$n_components <- 3
      }

      umap_res <- do.call(
        uwot::umap,
        c(
          list(X = mat),
          fargs[names(fargs) %in%
            names(as.list(args(uwot::umap)))]
        )
      )

      row.names(umap_res) <- colnames(cds)
      SingleCellExperiment::reducedDim(sce, reddim_method) <- umap_res

      if (reddim_method == "UMAP3D") {
        sce@metadata$reduced_dim_plots$umap3d_plot_ly <-
          plotly::plot_ly(
            x = umap_res[, 1],
            y = umap_res[, 2],
            z = umap_res[, 3],
            type = "scatter3d",
            mode = "markers",
            size = 0.1
          )

        sce@metadata$reduced_dim_plots$umap3d <-
          threejs::scatterplot3js(
            x = umap_res[, 1],
            y = umap_res[, 2],
            z = umap_res[, 3],
            axis = FALSE,
            num.ticks = NULL,
            color = rep("#4682b4", dim(umap_res)[[1]]),
            stroke = "#4682b4",
            size = .01
          )
      }
    }
    # Reduce dimensions with Liger
    if (reddim_method == "Liger") {

      # Preprocess with Liger
      ligerex <- liger_preprocess(sce, ...)
      sce@metadata$liger_params$liger_preprocess <-
        ligerex@parameters$liger_params$liger_preprocess

      # Reduce dimensions with Liger
      ligerex <- liger_reduce_dims(ligerex, ...)
      sce@metadata$liger_params$liger_reduce_dims <-
        ligerex@parameters$liger_params$liger_reduce_dims

      SingleCellExperiment::reducedDim(sce, "Liger") <- ligerex@H.norm
    }
  }

  return(sce)
}
