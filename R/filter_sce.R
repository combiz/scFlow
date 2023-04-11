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
#' @importFrom SummarizedExperiment rowData colData
#' @export
filter_sce <- function(sce,
                       filter_genes = TRUE,
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
    n_mito_keep <- sum(!SummarizedExperiment::rowData(sce)$qc_metric_mito_keep)
    cli::cli_alert_success(
        "{.emph {n_mito_keep}/{n_mito}} mitochondrial genes were dropped. \r\n")

    n_ribo <- sum(SummarizedExperiment::rowData(sce)$qc_metric_is_ribo)
    n_ribo_keep <- sum(!SummarizedExperiment::rowData(sce)$qc_metric_ribo_keep)
    cli::cli_alert_success(
      "{.emph {n_ribo_keep}/{n_ribo}} ribosomal genes were dropped. \r\n")

    n_unmapped <- sum(!SummarizedExperiment::rowData(sce)$qc_metric_ensembl_mapped)
    n_unmapped_keep <- sum(!SummarizedExperiment::rowData(sce)$qc_metric_mapped_keep)
    cli::cli_alert_success(
      "{.emph {n_unmapped_keep}/{n_unmapped}} unmapped genes were dropped. \r\n")

    n_non_expressive <- sum(
      !(SummarizedExperiment::rowData(sce)$qc_metric_is_expressive)
    )

    cli::cli_alert_success(
      "{.emph {n_non_expressive}} non-expressive genes were dropped. \r\n")

    n_genes_passed <- sum(SummarizedExperiment::rowData(sce)$qc_metric_gene_passed)
    n_genes_total <- dim(sce)[[1]]

    cli::cli_alert_success(
      "{.emph {n_genes_passed}/{n_genes_total}} genes passed QC. \r\n")

    sce <- sce[SummarizedExperiment::rowData(sce)$qc_metric_gene_passed, ]
    sce@metadata$scflow_steps$genes_filtered <- 1

    #update total_counts and total_features_by_counts
    sce$total_counts <- Matrix::colSums(SingleCellExperiment::counts(sce))
    sce$total_features_by_counts <- Matrix::colSums(
      SingleCellExperiment::counts(sce) > 0
    )

    sce$qc_metric_min_library_size <-
      sce$total_counts >= sce@metadata$qc_params$min_library_size
    sce$qc_metric_min_features <-
      sce$total_features_by_counts >= sce@metadata$qc_params$min_features
    sce$qc_metric_passed <- and_list_fn(
      sce$qc_metric_passed,
      sce$qc_metric_min_library_size,
      sce$qc_metric_min_features
    )

  }

  if (filter_cells == TRUE) {

    n_qc_failed <- sum(!sce$qc_metric_passed)

    cli::cli_alert_success(
      "{.emph {n_qc_failed}} cells failed QC and were dropped: - \r\n"
    )

    if(sce@metadata$scflow_steps$singlets_annotated){
      n_multiplets <- sum(!sce$is_singlet)
      cli::cli_alert_info("Singlet annotations found!")
      cli::cli_alert_success(
        "{.emph {n_multiplets}} multiplets were dropped. \r\n")
    }

    if(sce@metadata$scflow_steps$emptydrops_annotated){
      n_empty_drops <- sum(sce$is_empty_drop)
      cli::cli_alert_info("EmptyDrops annotations found!")
      cli::cli_alert_success(
        "{.emph {n_empty_drops}} empty drops were dropped. \r\n")
    }

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

    sce@metadata$qc_summary$qc_cells_n_cells_passed <- dim(sce)[[2]]


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

and_list_fn <- function(...) Reduce("&", list(...))
