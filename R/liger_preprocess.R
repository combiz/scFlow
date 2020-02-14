################################################################################
#' Preprocessing steps for Liger dimensionality reduction
#'
#'
#' Split merged object into multiple sce objects and extract sparse matrices:
#' @param sce merged SingleCellExperiment objects
#'
#' Make a Liger object:
#' @param raw.data List of expression matrices (gene by cell)
#'   each from an sce object.
#' @param take.gene.union Whether to fill out raw.data matrices with union
#'   of genes across all datasets (filling in 0 for missing data)
#'   (requires make.sparse=T) (default FALSE).
#' @param remove.missing Whether to remove cells not expressing any
#'   measured genes, and genes not
#'   expressed in any cells (if take.gene.union = T, removes only genes
#'   not expressed in any dataset) (default TRUE).
#'
#'  Select informative genes:
#' @param num.genes Number of genes to find for each dataset.
#'   Set to 3000 as default.
#' @param combine How to combine variable genes across experiments.
#'   Either "union" or "intersect".
#'   (default "union")
#' @param keep.unique Keep genes that occur (i.e., there is a
#'   corresponding column in raw.data) only in one dataset (default FALSE)
#' @param capitalize Capitalize gene names to match homologous genes
#'   (ie. across species)
#'   (default FALSE)
#' @param do.plot Display log plot of gene variance vs. gene expression
#'   for each dataset.
#'   Selected genes are plotted in green. (default FALSE)
#' @param cex.use Point size for plot.
#'
#'  Scale genes by root-mean-square across cells:
#' @param remove.missing Whether to remove cells from scale.data
#'   with no gene expression (default TRUE)
#'
#' Remove cells/genes with no expression across any genes/cells:
#' @param use.cols Treat each column as a cell (default TRUE)
#'
#' @return liger preprocessed object.
#'
#' @importFrom liger createLiger
#' @importFrom liger normalize
#' @importFrom liger selectGenes
#' @importFrom liger scaleNotCenter
#' @importFrom liger removeMissingObs
#'
#' @export

liger_preprocess <- function(sce,
                             unique_id_var = "manifest",
                             take_gene_union = F,
                             remove.missing = T,
                             num_genes = 3000,
                             combine = "union",
                             keep_unique = F,
                             capitalize = F,
                             do_plot = F,
                             cex_use = 0.3,
                             use_cols = T,
                             ...) {

  fargs <- as.list(environment())
  fargs <- fargs[fargs = c(
    "unique_id_var",
    "take_gene_union",
    "remove.missing",
    "num_genes",
    "combine",
    "keep_unique",
    "capitalize",
    "do_plot",
    "cex_use",
    "use_cols"
  )]

  do.call(.check_sce_for_liger, c(sce = sce, fargs))

  # Split merged sce object into multiple objects and extract sparse matrices
  dataset_list <- list()
  mat_list <- list()
  manifests <- unique(sce@colData[, unique_id_var])
  for (mnft in manifests) {
    dataset_name <- paste0("dataset_", mnft)
    dataset_list[[dataset_name]] <-
      sce[, sce@colData$manifest == mnft]
    mat_list[[dataset_name]] <-
      sce@assays$data$counts[, sce@colData$manifest == mnft]
  }

  # Make a Liger object. Pass in the sparse matrix.
  ligerex <- createLiger(
    raw.data = mat_list, take.gene.union = take_gene_union,
    remove.missing = remove.missing
  )

  ligerex@parameters$liger_params$liger_preprocess <- fargs


  ### preprocessing steps

  # Normalize the data to control for different numbers of UMIs per cell
  ligerex <- liger::normalize(ligerex)

  # Select variable (informative) genes

  ligerex <- liger::selectGenes(ligerex,
    num.genes = num_genes, combine = combine, keep.unique = keep_unique,
    capitalize = capitalize, do.plot = do_plot, cex.use = cex_use
  )

  # Scale the data by root-mean-square across cells
  ligerex <- liger::scaleNotCenter(ligerex, remove.missing = remove.missing)

  # Remove cells/genes with no expression across any genes/cells
  ligerex <- liger::removeMissingObs(ligerex, use.cols = use_cols)


  return(ligerex)
}

################################################################################
#' Helper function to check SCE before running liger
#'
#' @param sce SingleCellExperiment
#' @param ... args
#'
#' @return return 1 if completed
#' @family helper functions
#' @importFrom SummarizedExperiment colData
#' @importFrom assertthat is.scalar assert_that
#' @keywords internal
.check_sce_for_liger <- function(sce, ...) {

  fargs <- list(...)

  assertthat::is.scalar(fargs$unique_id_var)

  assertthat::assert_that(
    fargs$unique_id_var %in% names(SummarizedExperiment::colData(sce)))

  min_cells_per_id <- 5
  assertthat::assert_that(
    min(table(droplevels(sce[[unique_id_var]]))) >= min_cells_per_id,
    msg = sprintf("Need at least %s cells per id.", min_cells_per_id))

  return(1)

}

