################################################################################
#' Dimensionality reduction using Liger factorization
#'
#' @param ligerex SingleCellExperiment object preprocessed for Liger
#'
#' Perform iNMF on scaled datasets:
#' @param k Inner dimension of factorization (number of factors).
#'   Set to k=30 as default.
#' @param lambda Regularization parameter. Larger values penalize
#'   dataset-specific effects more strongly (ie. alignment should increase
#'   as lambda increases). Set to lambda=5.0 as default.
#' @param thresh Convergence threshold. Convergence occurs
#'   when |obj0-obj|/(mean(obj0,obj)) < thresh (default 1e-4).
#' @param max_iters Maximum number of block coordinate descent iterations
#'   to perform (default 100).
#' @param nrep Number of restarts to perform (iNMF objective function
#'   is non-convex, so taking the best objective from multiple successive
#'   initializations is recommended). For easier reproducibility,
#'   this increments the random seed by 1 for each consecutive restart,
#'   so future factorizations of the same dataset can be run with
#'   one rep if necessary. (default 1)
#' @param h_init Initial values to use for H matrices. (default NULL)
#' @param w_init Initial values to use for W matrix (default NULL)
#' @param v_init Initial values to use for V matrices (default NULL)
#' @param rand_seed Random seed to allow reproducible results (default 1).
#' @param print_obj Print objective function values after convergence
#'   (default FALSE).
#'
#' Quantile align (normalize) factor loadings:
#' @param quantiles Number of quantiles to use for quantile normalization
#'   (default 50).
#' @param ref_dataset Name of dataset to use as a "reference"
#'   for normalization. By default, the dataset with the largest number
#'   of cells is used.
#' @param min_cells Minimum number of cells to consider a cluster shared
#'   across datasets (default 2)
#' @param knn_k Number of nearest neighbors for within-dataset knn graph
#'   (default 20).
#' @param center Centers the data when scaling factors (useful for less
#'   sparse modalities like methylation data). (default FALSE)
#' @param resolution Controls the number of communities detected.
#'   Higher resolution -> more communities. (default 1)
#' @param ... Additional arguments.
#'
#' @return liger object with H, H.norm, W, and V slots sets.
#'
#' @importFrom rliger optimizeALS quantile_norm
#' @importFrom cli cli_alert
#'
#' @family Data integration
#'
#' @export

liger_reduce_dims <- function(ligerex,
                              k = 30,
                              lambda = 5.0,
                              thresh = 1e-4,
                              max_iters = 100,
                              nrep = 1,
                              h_init = NULL,
                              w_init = NULL,
                              v_init = NULL,
                              rand_seed = 1,
                              print_obj = FALSE,
                              quantiles = 50,
                              ref_dataset = NULL,
                              min_cells = 2,
                              knn_k = 20,
                              center = FALSE,
                              resolution = 1,
                              ...) {
  fargs <- as.list(environment())
  fargs <- fargs[fargs = c(
    "k",
    "lambda",
    "thresh",
    "max_iters",
    "nrep",
    "h_init",
    "w_init",
    "v_init",
    "rand_seed",
    "print_obj",
    "quantiles",
    "ref_dataset",
    "min_cells",
    "knn_k",
    "center",
    "resolution"
  )]
  ligerex@parameters$liger_params$liger_reduce_dims <- fargs
  ### Factorization
  # Perform iNMF on scaled datasets
  cli::cli_alert(
    "Performing integrative non-negative matrix factorization (iNMF)")
  ligerex <- rliger::optimizeALS(ligerex,
    k = k, lambda = lambda, thresh = thresh,
    max.iters = max_iters, nrep = nrep,
    H.init = h_init, W.init = w_init,
    V.init = v_init, rand.seed = rand_seed,
    print.obj = print_obj
  )
  ### Quantile Alignment/Normalization

  # Quantile align (normalize) factor loadings
  cli::cli_alert("Normalizing factor loadings")
  ligerex <- rliger::quantile_norm(ligerex,
    quantiles = quantiles,
    ref_dataset = ref_dataset,
    min_cells = min_cells,
    knn_k = knn_k,
    dims.use = seq_len(ncol(ligerex@H[[1]])),
    do.center = center,
    max_sample = 1000,
    eps = 0.9,
    refine.knn = TRUE,
    rand.seed = rand_seed,
    resolution = resolution
  )
  return(ligerex)
}
