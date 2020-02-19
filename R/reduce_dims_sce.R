################################################################################
#' Calculate dimensionality reductions for a SingleCellExperiment object
#' or merged SingleCellExperiment objects using tSNE, UMAP, UMAP3D
#'
#' @param sce a SingleCellExperiment object or merged sce objects
#' @param input_reduced_dim the input reducedDim
#' @param reduction_methods one or more of "tSNE","UMAP","UMAP3D"
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
#' @importFrom Rtsne Rtsne
#'
#' @export

reduce_dims_sce <- function(sce,
                            input_reduced_dim = c("PCA", "Liger"),
                            reduction_methods = c("tSNE", "UMAP", "UMAP3D"),
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
    n_threads = max(1, RcppParallel::defaultNumThreads() - 2),
    # Rtsne
    dims = 2,
    initial_dims = 50,
    perplexity = 30,
    theta = 0.5,
    check_duplicates = FALSE,
    partial_pca = FALSE,
    max_iter = 1000,
    is_distance = FALSE,
    Y_init = NULL,
    pca_center = TRUE,
    pca_scale = FALSE,
    normalize = TRUE,
    momentum = 0.5,
    final_momentum = 0.8,
    eta = 200,
    exaggeration_factor = 12,
    num_threads = max(1, RcppParallel::defaultNumThreads() - 2)
  )

  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  sce@metadata$reduced_dim_plots <- list()
  sce@metadata$reduced_dim_params <- list()

  if (!all(purrr::map_lgl(
    reduction_methods,
    ~ . %in% c("tSNE", "UMAP", "UMAP3D")
  ))) {
    stop("reduction methods must be from: tSNE, UMAP, UMAP3D")
  }

  before_rd_names <- SingleCellExperiment::reducedDimNames(sce)

  cli::cli_h1("Starting Dimensionality Reduction")

  # Generate res mod formula for preprocess
  if (!is.null(vars_to_regress_out)) {
    res_mod_formula_str <- paste0(
      "~", paste0(vars_to_regress_out, collapse = "+")
    )
  }

  # Preprocess cds

  cds <- .sce_to_cds(sce)

  cli::cli_h2("Computing Principal Component Analysis (PCA)")
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
  cli::cli_alert_success("{.strong PCA} was completed successfully")

  for (reddim_method in reduction_methods[!reduction_methods == "PCA"]) {

    header_msg <- dplyr::case_when(
      reddim_method == "tSNE" ~
        "Computing T-distributed Stochastic Neighbor Embedding (tSNE)",
      reddim_method == "UMAP" ~
        "Computing Uniform Manifold Approximation and Projection (UMAP)",
      reddim_method == "UMAP3D" ~
        "Computing 3D Uniform Manifold Approximation and Projection (UMAP3D)"
    )
    cli::cli_h2(header_msg)

    for (input_rd in input_reduced_dim) {

      # Reduce dimensions with tSNE
      if (reddim_method == "tSNE") {

        mat <- as.matrix(SingleCellExperiment::reducedDim(cds, input_rd))

        tsne_args <- fargs[names(fargs) %in%
                             names(as.list(args(Rtsne:::Rtsne.default)))]

        tsne_res <- do.call(
          Rtsne::Rtsne, c(list(X = mat, pca = FALSE), tsne_args)
          )

        tsne_data <- tsne_res$Y[, 1:fargs$dims]
        row.names(tsne_data) <- colnames(tsne_data)

        rd_name <- paste(reddim_method, input_rd, sep = "_")

        SingleCellExperiment::reducedDim(sce, rd_name) <- tsne_data
        sce@metadata$reduced_dim_params[[rd_name]] <- tsne_args

        cli::cli_alert_success(c(
          "{.strong {reddim_method}} was computed successfully ",
          "with {.strong {input_rd}} input"))
      }

      # Reduce dimensions with UMAP
      if (reddim_method %in% c("UMAP", "UMAP3D")) {

        mat <- as.matrix(SingleCellExperiment::reducedDim(cds, input_rd))

        umap_args <- fargs[names(fargs) %in%
                             names(as.list(args(uwot::umap)))]

        if (reddim_method == "UMAP3D") {
          umap_args$n_components <- 3
        }

        umap_res <- do.call(
          uwot::umap, c(list(X = mat), umap_args)
        )

        cli::cli_alert_success(c(
          "{.strong {reddim_method}} was computed successfully ",
          "with {.strong {input_rd}} input"))

        row.names(umap_res) <- colnames(cds)
        rd_name <- paste(reddim_method, input_rd, sep = "_")
        SingleCellExperiment::reducedDim(sce, rd_name) <- umap_res
        sce@metadata$reduced_dim_params[[rd_name]] <- umap_args

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
    }
  }

  new_reddims <- dplyr::setdiff(
    SingleCellExperiment::reducedDimNames(sce),
    before_rd_names
    )

  cli::cli_alert_success(
    c("Successfully added {.val {length(new_reddims)}} reducedDim ",
      "slots: {.var {paste0(new_reddims, collapse = \", \")}}"))

  return(sce)
}
