################################################################################
#' Write a feature-barcode matrix
#'
#' Write a feature-barcode matrix into a folder with .gz compressed files.  This
#' folder contains 'barcodes.tsv.gz', 'features.tsv.gz',and 'matrix.mtx.gz'.
#'
#' @param mat a matrix
#' @param folder_path path to save the feature barcode matrix files
#' @param overwrite_files set to `TRUE` to overwrite files if they already
#' exist.
#'
#' @family import and export functions
#'
#' @importFrom cli cli_alert_danger cli_text cli_alert_success cli_h1 cli_h2
#' @importFrom utils write.table
#' @importFrom Matrix writeMM
#' @importFrom R.utils gzip
#'
#' @export
write_sparse_matrix <- function(mat,
                                folder_path,
                                overwrite_files = TRUE) {

  cli::cli_h2("Writing feature-barcode matrix")

  mat <- as(mat, "dgTMatrix")

  if (class(mat) != "dgTMatrix") {
    cli::cli_alert_danger("Expected dgTMatrix, received {{class(mat)}}.")
    stop("Expected dgTMatrix.")
  }

  if (!dir.exists(file.path(folder_path))) {
    cli::cli_text("Creating dir: {.path {folder_path}} \r\n")
    dir.create(file.path(folder_path))
  }

  paths_l <- list()
  paths_l[["barcodes_path"]] <-  file.path(folder_path, "barcodes.tsv.gz")
  paths_l[["features_path"]] <-  file.path(folder_path, "features.tsv.gz")
  paths_l[["matrix_path"]] <-  file.path(folder_path, "matrix.mtx")

  # check files dont exist or throw error
  if (overwrite_files == FALSE) {
    if (!any(purrr::map_lgl(paths_l, file.exists))) {
      stop(cli::cli_alert_danger(c(
        "Files already exist.",
        "Pass overwrite_files=TRUE to overwrite."
      )))
    }
  }

  cli::cli_text("Writing: {.path {paths_l$barcodes_path}}")
  utils::write.table(
    as.character(colnames(mat)),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    file = gzfile(paths_l$barcodes_path))

  cli::cli_text("Writing: {.path {paths_l$features_path}}")
  utils::write.table(
    as.character(rownames(mat)),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    file = gzfile(paths_l$features_path))

  cli::cli_text("Writing: {.path {paths_l$matrix_path}}")
  Matrix::writeMM(
    obj = as(mat, "dgTMatrix"),
    file = paths_l$matrix_path)

  cli::cli_text("Compressing: {.path {paths_l$matrix_path}}")
  R.utils::gzip(
    paths_l$matrix_path,
    overwrite = overwrite_files)

  cli::cli_alert_success("Saved sparse matrix to {.file {folder_path}}")

}
