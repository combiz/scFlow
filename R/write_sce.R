################################################################################
#' Write a SingleCellExperiment to a folder
#'
#' Write a \linkS4class{SingleCellExperiment} into a folder: -
#' * Feature-barcode-matrix files: -
#'   - 'barcodes.tsv.gz'
#'   - 'features.tsv.gz'
#'   - 'matrix.mtx.gz'
#' * Metadata: -
#'   - sce_rowdata.tsv
#'   - sce_coldata.tsv
#'   - scecoldata_classes.tsv
#' * reducedDim(sce, x) (if present): -
#'   - ReducedDim_x.tsv
#'
#' @param sce a SingleCellExperiment Object
#' @param folder_path path to save the SingleCellExperiment
#'
#' @family import and export functions
#'
#' @import cli Matrix SummarizedExperiment dplyr SingleCellExperiment
#' @importFrom R.utils gzip
#'
#' @export
write_sce <- function(sce, folder_path) {

  cat(cli::rule("Writing SingleCellExperiment", line = 1), "\r\n")

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger(
      "Expected SingleCellExperiment, received {{class(sce)}}."))
  }

  if (!dir.exists(file.path(folder_path))) {
    cli::cli_text("Creating dir: {.path {folder_path}}")
    dir.create(file.path(folder_path))
  }

  # write metadata to tsv files
  cli::cli_text("Writing: {.path sce-rowdata.tsv}")
  write.table(
    as.data.frame(rowData(sce)),
    file = file.path(folder_path, "sce-rowdata.tsv"),
    quote = FALSE,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )

  cli::cli_text("Writing: {.path 'scecoldata.tsv'}")
  write.table(
    as.data.frame(colData(sce)),
    file = file.path(folder_path, "sce-coldata.tsv"),
    quote = FALSE,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )

  cli::cli_text("Writing: {.path 'scecoldata_classes.tsv'}")
  col_classes <- sapply(colData(sce), class)
  write.table(
    col_classes,
    file = file.path(folder_path, "scecoldata_classes.tsv"),
    sep = "\t",
    col.names = FALSE,
    row.names = T
  )

  # write reducedDim data to tsv files
  for (rd in names(reducedDims(sce))){
    rd_filename <- file.path(paste0("ReducedDim_", rd, ".tsv"))
    cli::cli_text("Writing: {.path {rd_filename}}")
    write.table(as.data.frame(SingleCellExperiment::reducedDim(sce, rd)),
                file = file.path(folder_path, rd_filename),
                quote = FALSE, sep = "\t",
                col.names = TRUE, row.names = FALSE
    )
  }

  # save the matrix
  write_feature_barcode_matrix(
    counts(sce),
    folder_path
  )

  cli::cli_alert_success(
    "Saved SingleCellExperiment to {.file {folder_path}}"
  )

}
