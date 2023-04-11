################################################################################
#' Read a feature-barcode matrix
#'
#' Reads a feature-barcode matrix folder into a sparse matrix.  This folder is
#' output by cellranger and must contain a 'barcodes.tsv.gz', 'features.tsv.gz',
#' and 'matrix.mtx.gz' file.
#'
#' @param folder_path path to folder containing feature-barcode matrix files
#'
#' @return mat a sparse matrix with rownames and colnames
#'
#'
#' @family import and export functions
#'
#' @importFrom Matrix Matrix
#' @importFrom cli cli_h1 cli_h2
#'
#' @export
read_sparse_matrix <- function(folder_path) {

  cli::cli_h2("Reading feature-barcode matrix")

  paths_l <- list()
  paths_l[["barcodes_path"]] <-  file.path(folder_path, "barcodes.tsv.gz")
  paths_l[["features_path"]] <-  file.path(folder_path, "features.tsv.gz")
  paths_l[["matrix_path"]] <-  file.path(folder_path, "matrix.mtx.gz")

  # check all files exist or throw error
  if (!all(purrr::map_lgl(paths_l, file.exists))) {
    stop(cli::cli_alert_danger("Files not found.  \
                               Specify a valid folder path."))
  }

  cli::cli_text("Reading: {.path {paths_l$barcodes_path}}")
  col_data <- read.csv(gzfile(paths_l$barcodes_path),
                       sep = "",
                       header = FALSE,
                       stringsAsFactors = FALSE)

  cli::cli_text("Reading: {.path {paths_l$features_path}}")
  row_data <- read.delim(gzfile(paths_l$features_path),
                         header = FALSE,
                         stringsAsFactors = FALSE)

  cli::cli_text("Reading: {.path {paths_l$matrix_path}}")
  mat <- Matrix::readMM(gzfile(paths_l$matrix_path))

  cli::cli_alert_success(sprintf("Imported sparse matrix: %s cols x %s rows.",
      dim(col_data)[[1]],
      dim(row_data)[[1]])
  )

  rownames(mat) <- as.character(row_data[, 1]) # ensembl ids
  colnames(mat) <- as.character(col_data[, 1]) # barcodes

  return(mat)

}
