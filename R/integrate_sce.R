################################################################################
#' Integrate datasets within a SingleCellExperiment Object
#'
#' Produces a reducedDim slot representing common factors across datasets
#'
#' @param sce a SingleCellExperiment object or merged sce objects
#' @param method the integration method to use. Set to liger as default.
#' See `liger_preprocess` and `liger_reduce_dims` functions for parameter options.
#' @param ... Additional arguments. Such as k, inner dimension of factorization
#' (number of factors).
#'
#' @return sce SingleCellExperiment object annotated with reducedDims
#'
#' @importFrom cli cli_h1 cli_h2 cli_h3 cli_alert_success
#'
#' @family Data integration
#'
#' @export

integrate_sce <- function(sce,
                          method = "Liger",
                          ...) {
  fargs <- list(...)
  integration_methods <- c("Liger")
  assertthat::assert_that(
    method %in% c(integration_methods),
    msg = sprintf(
      "Available integration methods are: %s",
      paste0(integration_methods, collapse = ",")
    )
  )
  cli::cli_h1("Integrating Datasets")
  # Reduce dimensions with Liger
  if (method == "Liger") {
    cli::cli_h2(
      "Running Linked Inference of Genomic Experimental Relationships (LIGER)"
    )
    # Preprocess with Liger
    cli::cli_h3("Pre-processing SingleCellExperiment for LIGER")
    ligerex <- do.call(liger_preprocess, c(list(sce = sce), fargs))
    sce@metadata$liger_params$liger_preprocess <-
      ligerex@parameters$liger_params$liger_preprocess
    sce@metadata$liger_var.genes <-
      ligerex@var.genes
    sce@metadata$dataset_integration$var.genes_per_dataset <-
      ligerex@agg.data$var.genes_per_dataset
    # Reduce dimensions with Liger
    cli::cli_h3("Computing integrated factors with LIGER")
    ligerex <- do.call(liger_reduce_dims, c(list(ligerex = ligerex), fargs))
    sce@metadata$liger_params$liger_reduce_dims <-
      ligerex@parameters$liger_params$liger_reduce_dims

    # Order ligerex@H.norm barcodes to match the orders of original sce

    idx <- match(sce$barcode, rownames(ligerex@H.norm))
    ligerex@H.norm <- ligerex@H.norm[idx, ]

    SingleCellExperiment::reducedDim(sce, "Liger") <- ligerex@H.norm
    cli::cli_alert_success(
      "Successfully computed integrated factors with LIGER"
    )
  }
  return(sce)
}
