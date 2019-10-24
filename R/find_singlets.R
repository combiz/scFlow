################################################################################
#' Find singlets in a SingleCellExperiment
#'
#' Runs a doublet/multiplet identification algorithm on the SingleCellExperiment
#' and returns a SingleCellExperiment annotated with is_singlet and, depending
#' on the algorithm, additional annotations / dimensionality reductions.
#'
#' @param sce a SingleCellExperiment object
#' @param singlet_find_method the method to use for identifying singlets
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

  cat(cli::rule("Finding Singlets", line = 1), "\r\n")

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
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

  #find_singlet_fn <- get(singlet_find_methods_l[[singlet_find_method]])

  sce <- do.call(
    singlet_find_methods_l[[singlet_find_method]],
    list(sce, args)
  )
  #sce <- find_singlet_fn(sce)

  return(sce)

}


