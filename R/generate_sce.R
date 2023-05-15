################################################################################
#' Generate a SingleCellExperiment object from a matrix and metadata
#'
#' Generate a SingleCellExperiment with a sparse matrix and a one-row dataframe
#' containing metadata for the sample.
#'
#' @param mat a sparse feature-barcode matrix
#' @param metadata a single row dataframe with sample metadata
#'
#' @return sce a SingleCellExperiment object annotated with sample metadata
#'
#' @family annotation functions
#' @importFrom cli cli_h1 cli_alert_danger cli_alert_success cli_text
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom stringr str_replace
#' @export
generate_sce <- function(mat, metadata) {

  cli::cli_h1("Generating SingleCellExperiment")

  if (typeof(mat) != "S4") {
    stop(cli::cli_alert_danger("A sparse Matrix::dgTMatrix is required."))
  }

  if (!is.data.frame(metadata)) {
    stop(cli::cli_alert_danger("Metadata should be a data.frame"))
  } else {
    if (dim(metadata)[[1]] != 1) {
      stop(cli::cli_alert_danger("The metadata data.frame should be dim 1xn"))
    }
  }

  if (is.null(colnames(mat))) {
    stop(cli::cli_alert_danger("Matrix coldata is missing."))
  }

  if (is.null(rownames(mat))) {
    stop(cli::cli_alert_danger("Matrix rowdata is missing."))
  }

  barcode <- paste(metadata[["manifest"]],
                   colnames(mat),
                   sep = "_") # unique barcode

  colnames(mat) <- barcode

  col_data <- cbind(barcode, metadata, row.names = NULL)

  rownames(mat) <- stringr::str_replace(
    rownames(mat),
    pattern = "\\..*",
    replacement = ""
  )

  mat <- mat[!duplicated(rownames(mat)), ]

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = mat),
    colData = data.frame(col_data),
    rowData = data.frame(ensembl_gene_id = rownames(mat))
  )

  SummarizedExperiment::rowData(sce)$ensembl_gene_id <- as.character(
    SummarizedExperiment::rowData(sce)$ensembl_gene_id
  )
  sce <- sce[order(rownames(sce)), ]

  cli::cli_alert_success("SingleCellExperiment object generated")
  cli::cli_text(c(
    dim(sce)[[1]], " gene x ", dim(sce)[[2]],
    " cells, annotated with ",
    dim(SummarizedExperiment::colData(sce))[[2]], " metadata variables (incl. barcode).")
  )

  sce@metadata$metadata <- metadata
  sce@metadata$scflow_steps$cells_annotated <- 0
  sce@metadata$scflow_steps$genes_annotated <- 0
  sce@metadata$scflow_steps$singlets_annotated <- 0
  sce@metadata$scflow_steps$singlets_method <- ""
  sce@metadata$scflow_steps$cells_filtered <- 0
  sce@metadata$scflow_steps$genes_filtered <- 0
  sce@metadata$scflow_steps$emptydrops_annotated <- 0


  return(sce)

}
