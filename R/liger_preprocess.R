################################################################################
#' Preprocessing steps for Liger dimensionality reduction
#'
#' Split merged object into multiple sce objects and extract sparse matrices:
#' @param sce SingleCellExperiment object or merged objects
#' @param k Inner dimension of factorization (number of factors).
#' @param unique_id_var the colData variable identifying unique samples.
#' Default is "manifest".
#'
#' Make a Liger object:
#' @param take_gene_union Whether to fill out raw.data matrices with union
#'   of genes across all datasets (filling in 0 for missing data)
#'   (requires make.sparse=T) (default FALSE).
#' @param remove.missing Whether to remove cells not expressing any
#'   measured genes, and genes not
#'   expressed in any cells (if take.gene.union = T, removes only genes
#'   not expressed in any dataset) (default TRUE).
#'
#'  Select informative genes:
#' @param num_genes Number of genes to find for each dataset.
#'   Set to 3000 as default.
#' @param combine How to combine variable genes across experiments.
#'   Either "union" or "intersect".
#'   (default "union")
#' @param capitalize Capitalize gene names to match homologous genes
#'   (ie. across species)
#'   (default FALSE)
#'  Scale genes by root-mean-square across cells:
#'
#' Remove cells/genes with no expression across any genes/cells:
#' @param use_cols Treat each column as a cell (default TRUE)
#' @param num_cores Number of cores used on user's machine to run function.
#' Default is 1.
#' @param ... Additional arguments.
#'
#' @return liger preprocessed object.
#'
#' @family Data integration
#'
#' @importFrom rliger createLiger normalize selectGenes
#' @importFrom rliger scaleNotCenter removeMissingObs
#' @importFrom parallel mclapply
#' @importFrom cli cli_alert
#'
#' @export

liger_preprocess <- function(sce,
                             k,
                             unique_id_var = "manifest",
                             take_gene_union = F,
                             remove.missing = T,
                             num_genes = 2000,
                             combine = "union",
                             capitalize = F,
                             use_cols = T,
                             num_cores = 1,
                             ...) {
  fargs <- as.list(environment())
  fargs <- fargs[fargs = c(
    "k",
    "unique_id_var",
    "take_gene_union",
    "remove.missing",
    "num_genes",
    "combine",
    "capitalize",
    "use_cols",
    "num_cores"
  )]

  do.call(.check_sce_for_liger, c(sce = sce, fargs))
  # Split merged sce object into multiple objects and extract sparse matrices

  cli::cli_alert("Extracting sparse matrices")

  id <- as.character(unique(sce[[unique_id_var]]))

  mat_list <- list()

  for(i in id){
    tmp <- subset(sce, ,get(unique_id_var) == i)
    mat_list[[i]] <- tmp@assays@data$counts
    }

  names(mat_list) <- paste0("dataset_", id)

  # Make a Liger object. Pass in the sparse matrix.
  cli::cli_alert("Creating LIGER object")
  ligerex <- rliger::createLiger(
    raw.data = mat_list, take.gene.union = take_gene_union,
    remove.missing = remove.missing
  )

  ligerex@parameters$liger_params$liger_preprocess <- fargs
  ### preprocessing steps
  # Normalize the data to control for different numbers of UMIs per cell
  cli::cli_alert("Normalizing data")
  ligerex <- rliger::normalize(ligerex)

  # Select variable (informative) genes
  cli::cli_alert("Selecting variable (informative) genes")
  ligerex <- rliger::selectGenes(ligerex,
    num.genes = num_genes,
    combine = combine,
    capitalize = capitalize
  )

  # Scale the data by root-mean-square across cells
  cli::cli_alert("Scaling data")
  ligerex <- rliger::scaleNotCenter(ligerex, remove.missing = remove.missing)
  # Remove cells/genes with no expression across any genes/cells
  cli::cli_alert("Removing non-expressed genes")
  ligerex <- rliger::removeMissingObs(ligerex, use.cols = use_cols)
  ########## Selecting and storing variable genes for each dataset ##########
  cli::cli_alert("Selecting and storing variable genes for each dataset")

  var.genes_per_dataset <- parallel::mclapply(
    mat_list,
    mc.cores = num_cores,
    function(mat) {
      single_mat <- mat
      single_ligerex <- rliger::createLiger(
        raw.data = list(single_dataset = single_mat),
        remove.missing = remove.missing
      )
      single_ligerex <- rliger::normalize(single_ligerex)
      single_ligerex <- rliger::selectGenes(single_ligerex,
        num.genes = num_genes,
        combine = combine,
        capitalize = capitalize
      )
      single_ligerex_var_genes <- single_ligerex@var.genes
    }
  )

  ligerex@agg.data$var.genes_per_dataset <- var.genes_per_dataset

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
    fargs$unique_id_var %in% names(SummarizedExperiment::colData(sce))
  )

  min_cells_per_id <- 5
  assertthat::assert_that(
    min(table(droplevels(as.factor(sce[[fargs$unique_id_var]])))) >= min_cells_per_id,
    msg = sprintf("Need at least %s cells per id.", min_cells_per_id)
  )

  return(1)
}
