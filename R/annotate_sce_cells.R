################################################################################
#' Add basic cell-wise annotations for a SingleCellExperiment
#'
#' Calculates the following QC metrics:
#' * qc_metric_total_counts - was the ensembl_gene_id found in biomaRt
#' * qc_metric_total_features_by_counts - is the gene mitochondrial
#'
#' @param sce a SingleCellExperiment object
#' @param ... See \code{\link{annotate_sce}}
#'
#' min_library_size the minimum number of counts per cell
#' min_features the minimum number of features per cell (i.e. the minimum
#'   number of genes with >0 counts)
#' max_mito the maximum proportion of counts mapping to
#'   mitochondrial genes (0 - 1)
#' max_ribo the maximum proportion of counts mapping to
#'   ribosomal genes (0 - 1)
#' min_counts the minimum number of counts per cell in min_cells
#' min_cells the minimum number of cells with min_counts
#'
#' @return sce a SingleCellExperiment object annotated with cell QC metrics
#'
#' @family annotation functions
#' @import cli Matrix SummarizedExperiment dplyr SingleCellExperiment
#' @export
annotate_sce_cells <- function(sce, ...) {

  cat(cli::rule("Annotating SingleCellExperiment cells", line = 1), "\r\n")

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  args <- list(...)

  sce$total_counts <- Matrix::colSums(SingleCellExperiment::counts(sce))
  sce$total_features_by_counts <- Matrix::colSums(
    SingleCellExperiment::counts(sce) > 0
  )

  sce$qc_metric_min_library_size <- sce$total_counts >= args$min_library_size

  sce$qc_metric_min_features <- sce$total_features_by_counts >= args$min_features

  sce$pc_mito <- Matrix::colSums(
    SingleCellExperiment::counts(
      sce[SummarizedExperiment::rowData(sce)$qc_metric_is_mito == 1, ]
    )
  ) / sce$total_counts

  sce$qc_metric_pc_mito_ok <- sce$pc_mito <= args$max_mito

  sce$pc_ribo <- Matrix::colSums(
    SingleCellExperiment::counts(
      sce[SummarizedExperiment::rowData(sce)$qc_metric_is_ribo == 1, ]
    )
  ) / sce$total_counts

  sce$qc_metric_pc_ribo_ok <- sce$pc_ribo <= args$max_ribo

  and_list_fn <- function(...) Reduce("&", list(...))
  sce$passed_qc <- and_list_fn(
    sce$qc_metric_min_library_size,
    sce$qc_metric_min_features,
    sce$qc_metric_pc_mito_ok,
    sce$qc_metric_pc_ribo_ok
  )

  # fast qc for expressive genes with counts from only qc passed cells
  mat <- SingleCellExperiment::counts(sce)
  mat <- mat >= args$min_counts
  qc_metric_n_cells_expressing <- Matrix::rowSums(mat[, sce$passed_qc])
  SummarizedExperiment::rowData(sce)$qc_metric_n_cells_expressing <-
    qc_metric_n_cells_expressing
  qc_metric_is_expressive <- qc_metric_n_cells_expressing >= args$min_cells
  SummarizedExperiment::rowData(sce)$qc_metric_is_expressive <-
    qc_metric_is_expressive

  return(sce)

}
