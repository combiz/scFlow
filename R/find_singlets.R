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
#' @import cli Matrix SummarizedExperiment dplyr SingleCellExperiment
#' @export

find_singlets <- function(sce,
                          singlet_find_method,
                          ...) {

  args <- list(...)

  cat(cli::rule("Finding Singlets", line = 2), "\r\n")

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  if(sce@metadata$scflow_steps$singlets_annotated == 1) {
    warning(cli::cli_alert_danger(
      "find_singlets was previously run on this data. Re-running."))
  }

  singlet_find_methods_l <- list()
  singlet_find_methods_l[["doubletfinder"]] <- "run_doubletfinder"

  if (!singlet_find_method %in% names(singlet_find_methods_l)) {
    valid_methods <- paste(names(singlet_find_methods_l), collapse = ", ")
    stop(cli::cli_alert_danger(c(
      "{'{singlet_find_method}'} is not a valid singlet finding method. ",
      "Try: {'{valid_methods}'}."))
    )
  }

  sce <- do.call(
    singlet_find_methods_l[[singlet_find_method]],
    list(sce, args)
  )

  sce@metadata$scflow_steps$singlets_annotated <- 1
  sce@metadata$scflow_steps$singlets_method <- singlet_find_method

  return(sce)

}
