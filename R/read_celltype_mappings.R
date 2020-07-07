################################################################################
#' Read Celltype Mappings
#'
#' @param file_path path to read the tsv
#'
#' @family import and export functions
#'
#' @importFrom cli cli_alert_danger cli_text cli_alert_success
#' @importFrom utils read.table
#'
#' @export
read_celltype_mappings <- function(file_path) {

  cli::cli_h1("Reading Celltype Mappings")

  # check all files exist or throw error
  if (!file.exists(file_path)) {
    stop("File not found.  Specify a valid path.")
  }


  cli::cli_text("Reading: {.path {file_path}}")
  celltype_mappings <- read.delim(
    file = file_path, stringsAsFactors = FALSE
  )

  cli::cli_alert_success(
    "Imported celltype mappings from {.file {file_path}}")

  return(celltype_mappings)

}
