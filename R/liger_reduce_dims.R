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
#' @param max.iters Maximum number of block coordinate descent iterations
#'   to perform (default 100).
#' @param nrep Number of restarts to perform (iNMF objective function
#'   is non-convex, so taking the best objective from multiple successive
#'   initializations is recommended). For easier reproducibility,
#'   this increments the random seed by 1 for each consecutive restart,
#'   so future factorizations of the same dataset can be run with
#'   one rep if necessary. (default 1)
#' @param H.init Initial values to use for H matrices. (default NULL)
#' @param W.init Initial values to use for W matrix (default NULL)
#' @param V.init Initial values to use for V matrices (default NULL)
#' @param rand.seed Random seed to allow reproducible results (default 1).
#' @param print.obj Print objective function values after convergence
#'   (default FALSE).
#'
#' Quantile align (normalize) factor loadings:
#' @param knn_k Number of nearest neighbors for within-dataset knn graph
#'   (default 20).
#' @param k2 Horizon parameter for shared nearest factor graph.
#'   Distances to all but the k2 nearest neighbors are set to 0
#'   (cuts down on memory usage for very large graphs) (default 500).
#' @param prune.thresh Minimum allowed edge weight. Any edges below
#'   this are removed (given weight 0) (default 0.2).
#' @param ref_dataset Name of dataset to use as a "reference"
#'   for normalization. By default, the dataset with the largest number
#'   of cells is used.
#' @param min_cells Minimum number of cells to consider a cluster shared
#'   across datasets (default 2)
#' @param quantiles Number of quantiles to use for quantile normalization
#'   (default 50).
#' @param nstart Number of times to perform Louvain community detection with
#'   different random starts (default 10).
#' @param resolution Controls the number of communities detected.
#'   Higher resolution -> more communities. (default 1)
#' @param dims.use Indices of factors to use for shared nearest
#'   factor determination (default 1:ncol(H[[1]])).
#' @param dist.use Distance metric to use in calculating nearest
#'   neighbors (default "CR").
#' @param center Centers the data when scaling factors (useful for less
#'   sparse modalities like methylation data). (default FALSE)
#' @param small.clust.thresh Extracts small clusters loading highly
#'   on single factor with fewer cells than this before regular alignment
#'   (default 0 -- no small cluster extraction).
#' @param id.number Number to use for identifying edge file
#'   (when running in parallel) (generates random value by default).
#' @param print.mod Print modularity output from clustering algorithm
#'   (default FALSE).
#' @param print.align.summary Print summary of clusters which did not
#'   align normally (default FALSE).
#'
#' @return liger object with H, H.norm, W, and V slots sets.
#'
#' @importFrom liger optimizeALS
#' @importFrom liger quantileAlignSNF
#' @importFrom cli cli_h2 cli_alert
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
                              knn_k = 20,
                              k2 = 500,
                              ref_dataset = NULL,
                              resolution = 1,
                              dist_use = "CR",
                              prune_thresh = 0.2,
                              min_cells = 2,
                              quantiles = 50,
                              nstart = 10,
                              center = FALSE,
                              small_clust_thresh = 0,
                              id_number = NULL,
                              print_mod = FALSE,
                              print_align_summary = FALSE,
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
    "knn_k",
    "k2",
    "ref_dataset",
    "resolution",
    "dist_use",
    "prune_thresh",
    "min_cells",
    "quantiles",
    "nstart",
    "center",
    "small_clust_thresh",
    "id_number",
    "print_mod",
    "print_align_summary"
  )]
  
  ligerex@parameters$liger_params$liger_reduce_dims <- fargs
  
  ### Factorization
  
  # Perform iNMF on scaled datasets
  cli::cli_alert("Performing integrative non-negative matrix factorization (iNMF)")
  ligerex <- liger::optimizeALS(ligerex,
                                k = k, lambda = lambda, thresh = thresh,
                                max.iters = max_iters, nrep = nrep, H.init = h_init, W.init = w_init,
                                V.init = v_init, rand.seed = rand_seed, print.obj = print_obj
  )
  
  ### Quantile Alignment/Normalization
  
  # Quantile align (normalize) factor loadings
  cli::cli_alert("Normalizing factor loadings")
  ligerex <- liger::quantileAlignSNF(ligerex,
                                     knn_k = knn_k, k2 = k2, prune.thresh = prune_thresh,
                                     ref_dataset = ref_dataset, min_cells = min_cells, quantiles = quantiles,
                                     nstart = nstart, resolution = resolution,
                                     dims.use = seq_len(ncol(ligerex@H[[1]])),
                                     dist.use = dist_use, center = center,
                                     small.clust.thresh = small_clust_thresh, id.number = id_number,
                                     print.mod = print_mod, print.align.summary = print_align_summary
  )
  
  return(ligerex)
}