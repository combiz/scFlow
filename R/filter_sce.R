################################################################################
#' Filter SingleCellExperiment according to QC metrics

#' Quality control metrics previously calculated are used to subset the
#' SingleCellExperiment for both cells and genes passing QC criteria previously
#' specified in the `annotate_sce_genes` and `annotate_sce_cells` functions
#'
#' @param sce a SingleCellExperiment object
#' @param filter_genes if set `FALSE`, genes will not be filtered
#' @param filter_cells if set `FALSE`, cells will not be filtered
#' @param keep_mito set `FALSE` to remove mitochondrial genes
#' @param keep_ribo set `FALSE` to remove ribosomal genes
#'
#' @return sce a SingleCellExperiment object filtered for QC passed cells and
#'   genes, with QC metrics annotations removed
#'
#' @family annotation functions
#' @import cli Matrix SummarizedExperiment dplyr SingleCellExperiment
#' @export
filter_sce <- function(sce,
                       filter_genes = FALSE,
                       filter_cells = TRUE,
                       keep_mito = FALSE,
                       keep_ribo = TRUE) {

  cat(cli::rule("Filtering SingleCellExperiment", line = 1), "\r\n")

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }


  if (filter_genes == TRUE) {

    if (keep_mito == FALSE) {
      n_mito <- sum(rowData(sce)$qc_metric_is_mito)
      sce <- sce[
        SummarizedExperiment::rowData(sce)$qc_metric_is_mito == FALSE, ]
      cli::cli_alert_success(
        "{.emph {n_mito}} mitochondrial genes were dropped. \r\n")
    }

    if (keep_ribo == FALSE) {
      n_ribo <- sum(rowData(sce)$qc_metric_is_ribo)
      sce <- sce[
        SummarizedExperiment::rowData(sce)$qc_metric_is_ribo == FALSE, ]
      cli::cli_alert_success(
        "{.emph {n_ribo}} ribosomal genes were dropped. \r\n")
    }

    n_genes <- length(rowData(sce)$qc_metric_is_expressive)
    sce <- sce[SummarizedExperiment::rowData(sce)$qc_metric_is_ribo == FALSE, ]
    cli::cli_alert_success(
      "{.emph {n_non_expressive}} non-expressive genes were dropped. \r\n")

  }

  if (filter_cells == TRUE) {

    n_qc_failed <- sum(!sce$passed_qc)
    sce <- sce[, sce$passed_qc == TRUE]
    cli::cli_alert_success(
      "{.emph {n_qc_failed}} cells failed QC and were dropped. \r\n"
    )

  } else {

    if (filter_genes == FALSE) {
      stop(cli::cli_alert_danger("Nothing to filter."))
    }

  }

  cli::cli_alert_success(c(
    "SingleCellExperiment was filtered successfully to ",
    "{.emph {dim(sce)[[1]]}} genes and ",
    "{.emph {dim(sce)[[2]]}} cells. \r\n")
  )

  return(sce)
}
