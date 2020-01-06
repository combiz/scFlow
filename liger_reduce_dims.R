################################################################################
#' Calculate dimensionality reductions for a SingleCellExperiment using Liger factorization
#'
#' @param sce a SingleCellExperiment object
#'
#' Make a Liger object:
#' @param raw.data List of expression matrices (gene by cell). Should be named by dataset.
#' @param make.sparse Whether to convert raw data into sparse matrices (default TRUE).
#' @param take.gene.union Whether to fill out raw.data matrices with union of genes across all
#'   datasets (filling in 0 for missing data) (requires make.sparse=T) (default FALSE).
#' @param remove.missing Whether to remove cells not expressing any measured genes, and genes not
#'   expressed in any cells (if take.gene.union = T, removes only genes not expressed in any
#'   dataset) (default TRUE).
#'
#'  Select informative genes:
#' @param var.thresh Variance threshold. Main threshold used to identify variable genes. Genes with
#'   expression variance greater than threshold (relative to mean) are selected.
#'   (higher threshold -> fewer selected genes). Accepts single value or vector with separate
#'   var.thresh for each dataset. (default 0.1)
#' @param alpha.thresh Alpha threshold. Controls upper bound for expected mean gene expression
#'   (lower threshold -> higher upper bound). (default 0.99)
#' @param num.genes Number of genes to find for each dataset. Optimises the value of var.thresh
#'   for each dataset to get this number of genes. Accepts single value or vector with same length
#'   as number of datasets (optional, default=NULL).
#' @param tol Tolerance to use for optimization if num.genes values passed in (default 0.0001).
#' @param datasets.use List of datasets to include for discovery of highly variable genes. 
#'   (default 1:length(object@raw.data))
#' @param combine How to combine variable genes across experiments. Either "union" or "intersect".
#'   (default "union")
#' @param keep.unique Keep genes that occur (i.e., there is a corresponding column in raw.data) only
#'    in one dataset (default FALSE)
#' @param capitalize Capitalize gene names to match homologous genes (ie. across species)
#'   (default FALSE)
#' @param do.plot Display log plot of gene variance vs. gene expression for each dataset.
#'   Selected genes are plotted in green. (default FALSE)
#' @param cex.use Point size for plot.
#'
#'  Scale genes by root-mean-square across cells:
#' @param remove.missing Whether to remove cells from scale.data with no gene expression
#'   (default TRUE)
#'
#' Remove cells/genes with no expression across any genes/cells:
#' @param use.cols Treat each column as a cell (default TRUE)
#'
#' Perform iNMF on scaled datasets:
#' @param k Inner dimension of factorization (number of factors). Run suggestK to determine
#'   appropriate value; a general rule of thumb is that a higher k will be needed for datasets with
#'   more sub-structure.
#' @param lambda Regularization parameter. Larger values penalize dataset-specific effects more
#'   strongly (ie. alignment should increase as lambda increases). Run suggestLambda to determine
#'   most appropriate value for balancing dataset alignment and agreement (default 5.0).
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh.
#'   (default 1e-4)
#' @param max.iters Maximum number of block coordinate descent iterations to perform (default 100).
#' @param nrep Number of restarts to perform (iNMF objective function is non-convex, so taking the
#'   best objective from multiple successive initializations is recommended). For easier
#'   reproducibility, this increments the random seed by 1 for each consecutive restart, so future
#'   factorizations of the same dataset can be run with one rep if necessary. (default 1)
#' @param H.init Initial values to use for H matrices. (default NULL)
#' @param W.init Initial values to use for W matrix (default NULL)
#' @param V.init Initial values to use for V matrices (default NULL)
#' @param rand.seed Random seed to allow reproducible results (default 1).
#' @param print.obj Print objective function values after convergence (default FALSE).
#'
#' Visually suggest appropiate k value:
#' @param k.test Set of factor numbers to test (default seq(5, 50, 5)).
#' @param lambda Lambda to use for all foctorizations (default 5).
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max.iters Maximum number of block coordinate descent iterations to perform
#' @param num.cores Number of cores to use for optimizing factorizations in parallel (default 1)
#' @param rand.seed Random seed for reproducibility (default 1).
#' @param gen.new Do not use optimizeNewK in factorizations. Results in slower factorizations.
#'   (default FALSE).
#' @param nrep Number restarts to perform at each k value tested (increase to produce
#'   smoother curve if results unclear) (default 1).
#' @param plot.log2 Plot log2 curve for reference on K-L plot (log2 is upper bound and con
#'   sometimes help in identifying "elbow" of plot). (default TRUE)
#' @param return.data Whether to return list of data matrices (raw) or dataframe (processed)
#'   instead of ggplot object (default FALSE).
#' @param return.raw If return.results TRUE, whether to return raw data (in format described below),
#'   or dataframe used to produce ggplot object. Raw data is list of matrices of K-L divergences
#'   (length(k.test) by n_cells). Length of list corresponds to nrep. (default FALSE)
#'
#' Visually suggest appropriate lambda value:
#' @param k Number of factors to use in test factorizations. See optimizeALS documentation.
#' @param lambda.test Vector of lambda values to test. If not given, use default set spanning
#'   0.25 to 60
#' @param rand.seed Random seed for reproducibility (default 1).
#' @param num.cores Number of cores to use for optimizing factorizations in parallel (default 1).
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max.iters Maximum number of block coordinate descent iterations to perform
#' @param knn_k Number of nearest neighbors for within-dataset knn in quantileAlignSNF (default 20).
#' @param k2 Horizon parameter for quantileAlignSNF (default 500).
#' @param ref_dataset Reference dataset for quantileAlignSNF (defaults to larger dataset).
#' @param resolution Resolution for quantileAlignSNF (default 1).
#' @param gen.new Do not use optimizeNewLambda in factorizations. Recommended to set TRUE
#'   when looking at only a small range of lambdas (ie. 1:7) (default FALSE)
#' @param nrep Number restarts to perform at each lambda value tested (increase to produce
#'   smoother curve if results unclear) (default 1).
#' @param return.data Whether to return list of data matrices (raw) or dataframe (processed)
#'   instead of ggplot object (default FALSE).
#' @param return.raw If return.results TRUE, whether to return raw data (in format described below),
#'   or dataframe used to produce ggplot object. Raw data is matrix of alignment values for each
#'   lambda value tested (each column represents a different rep for nrep).(default FALSE)
#'
#' Quantile align (normalize) factor loadings:
#' @param knn_k Number of nearest neighbors for within-dataset knn graph (default 20).
#' @param k2 Horizon parameter for shared nearest factor graph. Distances to all but the k2 nearest
#'   neighbors are set to 0 (cuts down on memory usage for very large graphs). (default 500)
#' @param prune.thresh Minimum allowed edge weight. Any edges below this are removed (given weight
#'  0) (default 0.2)
#' @param ref_dataset Name of dataset to use as a "reference" for normalization. By default,
#'   the dataset with the largest number of cells is used.
#' @param min_cells Minimum number of cells to consider a cluster shared across datasets (default 2)
#' @param quantiles Number of quantiles to use for quantile normalization (default 50).
#' @param nstart Number of times to perform Louvain community detection with different random
#'   starts (default 10).
#' @param resolution Controls the number of communities detected. Higher resolution -> more
#'   communities. (default 1)
#' @param dims.use Indices of factors to use for shared nearest factor determination (default
#'   1:ncol(H[[1]])).
#' @param dist.use Distance metric to use in calculating nearest neighbors (default "CR").
#' @param center Centers the data when scaling factors (useful for less sparse modalities like
#'   methylation data). (default FALSE)
#' @param small.clust.thresh Extracts small clusters loading highly on single factor with fewer
#'   cells than this before regular alignment (default 0 -- no small cluster extraction).
#' @param id.number Number to use for identifying edge file (when running in parallel)
#'   (generates random value by default).
#' @param print.mod Print modularity output from clustering algorithm (default FALSE).
#' @param print.align.summary Print summary of clusters which did not align normally (default FALSE).
#'
#' @return liger object with H, H.norm, W, and V slots sets.
#' @return Matrix of alignment plot results if indicated or ggplot object. Plots alignment vs. lambda to console.
#' @return Matrix of K-L plot results if indicated or ggplot object. Plots K-L divergence vs. k to console.
#'
#' @importFrom liger createLiger 
#' @importFrom BiocGenerics normalize
#' @importFrom liger selectGenes
#' @importFrom liger scaleNotCenter
#' @importFrom liger removeMissingObs
#' @importFrom liger optimizeALS
#' @importFrom liger suggestK
#' @importFrom liger suggestLambda
#' @importFrom liger quantileAlignSNF
#'
#' @export

liger_reduce_dims <- function(sce,
							make.sparse = T, 
							take.gene.union = F, 
                            remove.missing = T,
							var.thresh = 0.1, 
							alpha.thresh = 0.99, 
							num.genes = NULL,
							tol = 0.0001, 
							datasets.use = 1:length(object@raw.data), 
							combine = "union",
							keep.unique = F, 
							capitalize = F, 
							do.plot = F, 
							cex.use = 0.3
							remove.missing = T,
							use.cols = T,
							k = 20,
							lambda = 5.0,
							thresh = 1e-4,
							max.iters = 100,
							nrep = 1,
							H.init = NULL,
							W.init = NULL,
							V.init = NULL,
							rand.seed = 1,
							print.obj = FALSE,
							k.test = seq(5, 50, 5), 
							num.cores = 1, 
							gen.new = F, 
							plot.log2 = T,
							return.data = F, 
							return.raw = F,
							lambda.test = NULL, 
							knn_k = 20, 
							k2 = 500, 
							ref_dataset = NULL,
                            resolution = 1, 
							dims.use = 1:ncol(H[[1]]),
							dist.use = FALSE,
							prune.thresh = 0.2,
							min_cells = 2,
							quantiles = 50,
							nstart = 10,
							center = FALSE,
							small.clust.thresh = 0,
							id.number = NULL,
							print.mod = FALSE,
							print.align.summary = FALSE,
                            ...) {

  # Extract sparse matrix from sce
  mat <- sce@assays$data
  
  # Make a Liger object. Pass in the sparse matrix.  
  ligerex <- createLiger(raw.data = mat, make.sparse = make.sparse, take.gene.union = take.gene.union, 
                         remove.missing = remove.missing)

  ### preprocessing steps
  
  # Normalize the data to control for different numbers of UMIs per cell
  ligerex <- BiocGenerics::normalize(ligerex)

  # Select variable (informative) genes 
  ligerex <- liger::selectGenes(ligerex, var.thresh = var.thresh, alpha.thresh = alpha.thresh, num.genes = num.genes,
                         tol = tol, datasets.use = datasets.use, combine = combine,
                         keep.unique = keep.unique, capitalize = capitalize, do.plot = do.plot, cex.use = cex.use)

  # Scale the data by root-mean-square across cells
  ligerex = liger::scaleNotCenter(ligerex, remove.missing = remove.missing)

  # Remove cells/genes with no expression across any genes/cells
  ligerex <- liger::removeMissingObs(ligerex, use.cols = use.cols)
  
  ### Factorization
  
  # Perform iNMF on scaled datasets
  ligerex <- liger::optimizeALS(ligerex, k = k, lambda = lambda, thresh = thresh,
						 max.iters = max.iters, nrep = nrep, H.init = H.init, W.init = W.init,
						 V.init = V.init, rand.seed = rand.seed, print.obj = print.obj)
  
  # Visually suggest appropiate k value
  ligerex <- liger::suggestK(ligerex, k.test = k.test, lambda = lambda, thresh = thresh, max.iters = max.iters,
                         num.cores = num.cores, rand.seed = rand.seed, gen.new = gen.new, nrep = nrep, plot.log2 = plot.log2,
                         return.data = return.data, return.raw = return.raw)
  
  # Visually suggest appropriate lambda value
  ligerex <- liger::suggestLambda(ligerex, k = k, lambda.test = lambda.test, rand.seed = rand.seed, num.cores = num.cores,
                         thresh = thresh, max.iters = max.iters, knn_k = knn_k, k2 = k2, ref_dataset = ref_dataset,
                         resolution = resolution, gen.new = gen.new, nrep = nrep, return.data = return.data, return.raw = return.raw)
  
  ### Quantile Alignment/Normalization
  
  # Quantile align (normalize) factor loadings
  ligerex <- liger::quantileAlignSNF(ligerex, knn_k = knn_k, k2 = k2, prune.thresh = prune.thresh, ref_dataset = ref_dataset, min_cells = min_cells,
						 quantiles = quantiles, nstart = nstart, resolution = resolution, dims.use = dims.use, dist.use = dist.use, center = center,
						 small.clust.thresh = small.clust.thresh, id.number = id.number, print.mod = print.mod, print.align.summary = print.align.summary)  

  return(ligerex)

}