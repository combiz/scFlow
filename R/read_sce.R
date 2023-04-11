################################################################################
#' Read a SingleCellExperiment from a folder
#'
#' Read a \linkS4class{SingleCellExperiment} from a folder: -
#' * Feature-barcode-matrix files: -
#'   - 'barcodes.tsv.gz'
#'   - 'features.tsv.gz'
#'   - 'matrix.mtx.gz'
#' * Metadata: -
#'   - 'sce_rowdata.tsv'
#'   - 'sce_coldata.tsv'
#'   - 'scecoldata_classes.tsv'
#' * reducedDim(sce, x) (if present): -
#'   - 'ReducedDim_x.tsv'
#'
#' @param folder_path path to save the SingleCellExperiment
#' @param read_metadata Enable import of metadata if previously saved.
#'
#' @return sce a SingleCellExperiment object
#'
#' @family import and export functions
#'
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim
#' @importFrom cli rule cli_alert_danger col_green cli_alert_info cli_alert_success cli_text cli_h1
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom R.utils gzip
#' @importFrom tools file_path_sans_ext
#' @importFrom qs qread
#' @importFrom future availableCores
#'
#' @export
read_sce <- function(folder_path, read_metadata = FALSE) {

  cli::cli_h1("Reading SingleCellExperiment")

  ## Check core files are present
  paths_l <- list()
  paths_l[["all_coldata"]] <-  file.path(folder_path, "sce-coldata.tsv")
  paths_l[["all_rowdata"]] <-  file.path(folder_path, "sce-rowdata.tsv")
  paths_l[["col_classes"]] <-  file.path(folder_path, "scecoldata_classes.tsv")
  paths_l[["barcodes_path"]] <-  file.path(folder_path, "barcodes.tsv.gz")
  paths_l[["features_path"]] <-  file.path(folder_path, "features.tsv.gz")
  paths_l[["matrix_path"]] <-  file.path(folder_path, "matrix.mtx.gz")
  if(read_metadata) {
    paths_l[["metadata_path"]] <-  file.path(folder_path, "metadata.qs")
    assertthat::assert_that(file.exists(paths_l[["metadata_path"]]),
                            msg = "Metadata was not found.")
  }

  # check all files exist or throw error
  if (!(all(purrr::map_lgl(paths_l, file.exists)))) {
    stop("Files not found.  Specify a valid folder path.")
  }

  all_counts <- read_sparse_matrix(folder_path)

  cli::cli_h2("Reading colData/rowData")
  cli::cli_text("Reading: {.path {paths_l$all_rowdata}}")
  all_rowdata <- read.delim(
    file = paths_l$all_rowdata
  )

  cli::cli_text("Reading: {.path {paths_l$col_classes}}")
  col_classes <- read.delim(
    paths_l$col_classes,
    header = FALSE
  )
  cc <- as.character(col_classes[, 2])
  names(cc) <- as.character(col_classes[, 1])

  cli::cli_text("Reading: {.path {paths_l$all_coldata}}")
  all_coldata <- read.delim(
    file = paths_l$all_coldata,
    colClasses = cc
  )

  cli::cli_text("Generating SingleCellExperiment")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = all_counts),
    rowData = all_rowdata,
    colData = all_coldata
  )

  rm(all_coldata, all_counts, all_rowdata)
  rd_files_l <- list.files(
    folder_path)[startsWith(list.files(folder_path), prefix = "ReducedDim_")]

  # if there are reduced dimension data, import them
  if (length(rd_files_l) > 0) {
    for (rd_file in rd_files_l) {
      rdname <- tools::file_path_sans_ext(gsub("ReducedDim_", "", rd_file))
      cli::cli_text("Reading: {.path {rd_file}}")
      SingleCellExperiment::reducedDim(sce, rdname) <- as.matrix(read.delim(
        file = file.path(folder_path, rd_file)), nrow = ncol(sce))
    }
  }

  if(read_metadata){
    cli::cli_h2("Appending metadata")
    metadata_path <- paths_l[["metadata_path"]]
    cli::cli_text("Reading: {.path {metadata_path}}")
    metadata <- qs::qread(
      metadata_path,
      nthreads = max(1, future::availableCores()-1))
    sce@metadata <- metadata
  }

  cli::cli_alert_success("Imported SingleCellExperiment.")


  return(sce)

}
