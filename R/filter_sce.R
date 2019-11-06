################################################################################
#' Filter SingleCellExperiment according to QC metrics

#' Quality control metrics previously calculated are used to subset the
#' SingleCellExperiment for both cells and genes passing QC criteria previously
#' specified in the `annotate_sce_genes` and `annotate_sce_cells` functions
#'
#' @param sce a SingleCellExperiment object
#' @param filter_genes if set `FALSE`, genes will not be filtered
#' @param filter_cells if set `FALSE`, cells will not be filtered
#'
#' @return sce a SingleCellExperiment object filtered for QC passed cells and
#'   genes, with QC metrics annotations removed
#'
#' @family annotation functions
#' @import cli Matrix dplyr SingleCellExperiment
#' @importFrom SummarizedExperiment rowData colData
#' @export
filter_sce <- function(sce,
                       filter_genes = FALSE,
                       filter_cells = TRUE) {

  cat(cli::rule("Filtering SingleCellExperiment", line = 2), "\r\n")

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  cli::cli_text(c(
    "Pre-filtered SingleCellExperiment contains ",
    "{.emph {dim(sce)[[1]]}} genes and ",
    "{.emph {dim(sce)[[2]]}} cells.")
  )

  if (filter_genes == TRUE) {

    n_mito <- sum(SummarizedExperiment::rowData(sce)$qc_metric_is_mito)
    if (drop_mito == TRUE) {
      sce <-
        sce[SummarizedExperiment::rowData(sce)$qc_metric_is_mito == FALSE, ]
      cli::cli_alert_success(
        "{.emph {n_mito}} mitochondrial genes were dropped. \r\n")
    } else {
      cli::cli_alert_info(
        "{.emph {n_mito}} mitochondrial genes were retained. \r\n")
    }

    n_ribo <- sum(SummarizedExperiment::rowData(sce)$qc_metric_is_ribo)
    if (drop_ribo == TRUE) {
      sce <-
        sce[SummarizedExperiment::rowData(sce)$qc_metric_is_ribo == FALSE, ]
      cli::cli_alert_success(
        "{.emph {n_ribo}} ribosomal genes were dropped. \r\n")
    } else {
      cli::cli_alert_info(
        "{.emph {n_ribo}} ribosomal genes were retained. \r\n")
    }

    n_unmapped <-
      sum(!SummarizedExperiment::rowData(sce)$qc_metric_ensembl_mapped)
    if (drop_unmapped == TRUE) {
      sce <-
        sce[SummarizedExperiment::rowData(sce)$qc_metric_ensembl_mapped == TRUE, ]
      cli::cli_alert_success(
        "{.emph {n_unmapped}} unmapped genes were dropped. \r\n")
    } else {
      cli::cli_alert_info(
        "{.emph {n_unmapped}} unmapped genes were retained. \r\n")
    }

    n_non_expressive <- sum(
      !(SummarizedExperiment::rowData(sce)$qc_metric_is_expressive)
    )
    sce <-
      sce[SummarizedExperiment::rowData(sce)$qc_metric_is_expressive == TRUE, ]

    sce@metadata$scflow_steps$genes_filtered <- 1

    cli::cli_alert_success(
      "{.emph {n_non_expressive}} non-expressive genes were dropped. \r\n")

  }

  if (filter_cells == TRUE) {

    n_qc_failed <- sum(!sce$qc_metric_passed)

    cli::cli_alert_success(
      "{.emph {n_qc_failed}} cells failed QC and were dropped: - \r\n"
    )

    cli::cli_div(theme = list(
      span.pass = list(color = "green"),
      span.fail = list(color = "red")
    ))

    cli::cli_dl(c(
      "Library size QC" =
        paste0("{.pass PASS: {sum(sce$qc_metric_min_library_size, na.rm=T)}}/",
          "{.fail FAIL: {sum(!sce$qc_metric_min_library_size, na.rm=T)}}"),
      "Number of genes QC" =
        paste0("{.pass PASS: {sum(sce$qc_metric_min_features, na.rm=T)}}/",
               "{.fail FAIL: {sum(!sce$qc_metric_min_features, na.rm=T)}}"),
      "Mitochondria counts proportion QC" =
        paste0("{.pass PASS: {sum(sce$qc_metric_pc_mito_ok, na.rm=T)}}/",
               "{.fail FAIL: {sum(!sce$qc_metric_pc_mito_ok, na.rm=T)}}"),
      "Ribosomal counts proportion QC" =
        paste0("{.pass PASS: {sum(sce$qc_metric_pc_ribo_ok, na.rm=T)}}/",
               "{.fail FAIL: {sum(!sce$qc_metric_pc_ribo_ok, na.rm=T)}}"),
      "{.strong Overall cell QC}" =
        paste0("{.pass PASS: {sum(sce$qc_metric_passed, na.rm=T)}}/",
               "{.fail FAIL: {sum(!sce$qc_metric_passed, na.rm=T)}}")
      )
    )

    # subset
    sce <- sce[, sce$qc_metric_passed == TRUE]

    sce@metadata$scflow_steps$cells_filtered <- 1

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
