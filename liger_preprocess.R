################################################################################
#' Preprocessing steps for Liger dimensionality reduction  
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
#' @return liger object.
#'
#' @importFrom liger createLiger 
#' @importFrom liger normalize
#' @importFrom liger selectGenes
#' @importFrom liger scaleNotCenter
#' @importFrom liger removeMissingObs
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
                            ...) {

  # Extract sparse matrix from sce
  mat <- sce@assays$data$counts
  
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
  
  return(ligerex)

}
