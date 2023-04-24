################################################################################
#' Find singlets in a SingleCellExperiment
#'
#' Runs a doublet/multiplet identification algorithm on the SingleCellExperiment
#' and returns a SingleCellExperiment annotated with is_singlet and, depending
#' on the algorithm, additional annotations / dimensionality reductions.
#'
#' @param sce a SingleCellExperiment object
#' @param singlet_find_method the method to use for identifying singlets
#' @param ... additional parameters passed to singlet finding algorithm.
#'
#' @return sce a SingleCellExperiment object annotated for singlets
#'
#' @family annotation functions
#' @importFrom SummarizedExperiment rowData colData
#' @export

find_singlets <- function(sce,
                          singlet_find_method,
                          ...) {
  fargs <- list(...)

  cli::cli_h1("Finding Singlets")

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  if (sce@metadata$scflow_steps$singlets_annotated == 1) {
    warning(cli::cli_alert_danger(
      "find_singlets was previously run on this data. Re-running."
    ))
  }

  singlet_find_methods_l <- list()
  singlet_find_methods_l[["doubletfinder"]] <- "run_doubletfinder"

  if (!singlet_find_method %in% names(singlet_find_methods_l)) {
    valid_methods <- paste(names(singlet_find_methods_l), collapse = ", ")
    stop(cli::cli_alert_danger(c(
      "{'{singlet_find_method}'} is not a valid singlet finding method. ",
      "Try: {'{valid_methods}'}."
    )))
  }

  if (singlet_find_method == "doubletfinder") {
    fargs <- list(
      pK = NULL,
      pca_dims = 10,
      vars_to_regress_out = "nCount_RNA",
      var_features = 2000,
      doublet_rate = 0,
      dpk = 8, # estimated doublets per thousand cells
      num_cores = max(1, future::availableCores() - 2)
    )
    inargs <- list(...)
    fargs[names(inargs)] <- inargs # override defaults if provided
    sce <- run_doubletfinder(
      sce,
      pK = fargs$pK,
      pca_dims = fargs$pca_dims,
      vars_to_regress_out = fargs$vars_to_regress_out,
      var_features = fargs$var_features,
      doublet_rate = fargs$doublet_rate,
      dpk = fargs$dpk,
      num_cores = fargs$num_cores
    )
  }
  ## this approach results in a bug with ScaleData
  # sce <- do.call(
  #  singlet_find_methods_l[[singlet_find_method]],
  #  c(list(sce = sce), fargs)
  # )

  sce@metadata$scflow_steps$singlets_annotated <- 1
  sce@metadata$scflow_steps$singlets_method <- singlet_find_method

  return(sce)
}
