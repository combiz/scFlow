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
#' min_ribo
#' max_ribo the maximum proportion of counts mapping to
#'   ribosomal genes (0 - 1)
#' min_counts the minimum number of counts per cell in min_cells
#' min_cells the minimum number of cells with min_counts
#'
#' @return sce a SingleCellExperiment object annotated with cell QC metrics
#'
#' @family annotation functions
#' @import cli Matrix dplyr SingleCellExperiment
#' @importFrom SummarizedExperiment rowData colData
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

  sce$qc_metric_min_features <-
    sce$total_features_by_counts >= args$min_features

  if (length(sce$pc_mito) == 0) {
    sce$pc_mito <- Matrix::colSums(
      SingleCellExperiment::counts(
        sce[SummarizedExperiment::rowData(sce)$qc_metric_is_mito == 1, ]
      )
    ) / sce$total_counts
  } else {
    print("Retaining previously calculated percentage mitochondrial counts.")
  }

  sce$qc_metric_pc_mito_ok <- sce$pc_mito <= args$max_mito

  if (length(sce$pc_ribo) == 0) {
    sce$pc_ribo <- Matrix::colSums(
      SingleCellExperiment::counts(
        sce[SummarizedExperiment::rowData(sce)$qc_metric_is_ribo == 1, ]
      )
    ) / sce$total_counts
  } else {
    print("Retaining previously calculated percentage ribosomal counts.")
  }

  sce$qc_metric_pc_ribo_ok <-
    (sce$pc_ribo <= args$max_ribo) & (sce$pc_ribo >= args$min_ribo)

  and_list_fn <- function(...) Reduce("&", list(...))
  sce$qc_metric_passed <- and_list_fn(
    sce$qc_metric_min_library_size,
    sce$qc_metric_min_features,
    sce$qc_metric_pc_mito_ok,
    sce$qc_metric_pc_ribo_ok
  )

  if (sce@metadata$scflow_steps$emptydrops_annotated == 1) {
    sce$qc_metric_passed <- sce$qc_metric_passed & !sce$is_empty_drop
  }

  # fast qc for expressive genes with counts from only qc passed cells
  mat <- SingleCellExperiment::counts(sce)
  mat <- mat >= args$min_counts
  qc_metric_n_cells_expressing <- Matrix::rowSums(mat[, sce$qc_metric_passed])
  SummarizedExperiment::rowData(sce)$qc_metric_n_cells_expressing <-
    qc_metric_n_cells_expressing
  qc_metric_is_expressive <- qc_metric_n_cells_expressing >= args$min_cells
  SummarizedExperiment::rowData(sce)$qc_metric_is_expressive <-
    qc_metric_is_expressive

  # final row qc flag
  SummarizedExperiment::rowData(sce)$qc_metric_gene_passed <- and_list_fn(
    SummarizedExperiment::rowData(sce)$qc_metric_is_expressive,
    SummarizedExperiment::rowData(sce)$qc_metric_mapped_keep,
    SummarizedExperiment::rowData(sce)$qc_metric_mito_keep,
    SummarizedExperiment::rowData(sce)$qc_metric_ribo_keep
  )

  sce@metadata$scflow_steps$cells_annotated <- 1

  return(sce)

}
