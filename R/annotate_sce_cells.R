################################################################################
#' Add basic cell-wise annotations for a SingleCellExperiment
#'
#' Calculates the following QC metrics:
#' * qc_metric_total_counts - was the ensembl_gene_id found in biomaRt
#' * qc_metric_total_features_by_counts - is the gene mitochondrial
#'
#' @param sce a SingleCellExperiment object
#' @param min_library_size the minimum number of counts per cell
#' @param min_features the minimum number of features per cell (i.e. the minimum
#'   number of genes with >0 counts)
#' @param max_mito the maximum proportion of counts mapping to
#'   mitochondrial genes (0 - 1)
#' @param max_ribo the maximum proportion of counts mapping to
#'   ribosomal genes (0 - 1)
#' @param min_counts the minimum number of counts per cell in min_cells
#' @param min_cells the minimum number of cells with min_counts
#'
#' @return sce a SingleCellExperiment object annotated with cell QC metrics
#'
#' @family annotation functions
#' @import cli Matrix SummarizedExperiment dplyr SingleCellExperiment
#' @export
annotate_sce_cells <- function(sce,
                               min_library_size = 300,
                               min_features = 100,
                               max_mito = 0.10,
                               max_ribo = 1.00,
                               min_counts = 2,
                               min_cells = 2) {
  cat(cli::rule("Annotating SingleCellExperiment cells", line = 1), "\r\n")

  if (typeof(sce) != "S4") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  before_coldata_colnames <- colnames(SummarizedExperiment::colData(sce))
  before_rowdata_colnames <- colnames(SummarizedExperiment::rowData(sce))

  sce$total_counts <- Matrix::colSums(SingleCellExperiment::counts(sce))
  sce$total_features_by_counts <- Matrix::colSums(
    SingleCellExperiment::counts(sce) > 0
  )

  sce$qc_metric_min_library_size <- sce$total_counts >= min_library_size

  sce$qc_metric_min_features <- sce$total_features_by_counts >= min_features

  sce$pc_mito <- Matrix::colSums(
    SingleCellExperiment::counts(
      sce[SummarizedExperiment::rowData(sce)$qc_metric_is_mito == 1, ]
    )
  ) / sce$total_counts

  sce$qc_metric_pc_mito_ok <- sce$pc_mito <= max_mito

  sce$pc_ribo <- Matrix::colSums(
    SingleCellExperiment::counts(
      sce[SummarizedExperiment::rowData(sce)$qc_metric_is_ribo == 1, ]
    )
  ) / sce$total_counts

  sce$qc_metric_pc_ribo_ok <- sce$pc_ribo <= max_ribo

  and_list_fn <- function(...) Reduce("&", list(...))
  sce$passed_qc <- and_list_fn(
    sce$qc_metric_min_library_size,
    sce$qc_metric_min_features,
    sce$qc_metric_pc_mito_ok,
    sce$qc_metric_pc_ribo_ok
  )

  # fast qc for expressive genes with counts from only qc passed cells
  mat <- SingleCellExperiment::counts(sce)
  mat <- mat >= min_counts
  qc_metric_n_cells_expressing <- Matrix::rowSums(mat[, sce$passed_qc])
  SummarizedExperiment::rowData(sce)$qc_metric_n_cells_expressing <-
    qc_metric_n_cells_expressing
  qc_metric_is_expressive <- qc_metric_n_cells_expressing >= min_cells
  SummarizedExperiment::rowData(sce)$qc_metric_is_expressive <-
    qc_metric_is_expressive

  cli::cli_alert_success(
    "SingleCellExperiment cells were successfully annotated with: \r\n"
  )

  cli::cli_ul(setdiff(
    colnames(SummarizedExperiment::colData(sce)),
    before_coldata_colnames
  ))

  cli::cli_alert_success(
    "SingleCellExperiment genes were successfully annotated with: \r\n"
  )

  cli::cli_ul(setdiff(
    colnames(SummarizedExperiment::rowData(sce)),
    before_rowdata_colnames
  ))

  return(sce)
}
