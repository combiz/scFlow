################################################################################
#' Write Celltype Mappings
#'
#' @param sce a singlecellexperiment
#' @param folder_path path to save the tsv
#' @param overwrite_files set to `TRUE` to overwrite files if they already
#' exist.
#'
#' @family import and export functions
#'
#' @importFrom cli cli_alert_danger cli_text cli_alert_success
#' @importFrom utils write.table
#' @importFrom Matrix writeMM
#' @importFrom R.utils gzip
#'
#' @export
write_celltype_mappings <- function(sce,
                                    folder_path,
                                    overwrite_files = TRUE) {

  cli::cli_h1("Writing Celltype Mappings")

  assertthat::assert_that(class(sce@metadata$mappings) == "data.frame")

  if (!dir.exists(file.path(folder_path))) {
    cli::cli_text("Creating dir: {.path {folder_path}} \r\n")
    dir.create(file.path(folder_path))
  }

  paths_l <- list()
  paths_l[["mappings_fp"]] <-  file.path(folder_path, "celltype_mappings.tsv")

  # check files dont exist or throw error
  if (overwrite_files == FALSE) {
    if (!any(purrr::map_lgl(paths_l, file.exists))) {
      stop(cli::cli_alert_danger(c(
        "Files already exist.",
        "Pass overwrite_files=TRUE to overwrite."
      )))
    }
  }

  cli::cli_text("Writing: {.path {paths_l$mappings_fp}}")
  utils::write.table(
    sce@metadata$mappings,
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    file = paths_l$mappings_fp)

  cli::cli_alert_success(
    "Saved celltype mappings to {.file {paths_l$mappings_fp}}")

}
