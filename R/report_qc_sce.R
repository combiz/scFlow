################################################################################
#' Generate plots and a QC report for a SingleCellExperiment
#'
#' @param sce a SingleCellExperiment object
#' @param report_folder_path folder path to save the report
#'
#' @return sce a annotated SingleCellExperiment object
#'
#' @family annotation functions
#' @import cli Matrix SummarizedExperiment dplyr SingleCellExperiment purrr
#' @import ggplot2
#' @importFrom rmarkdown render
#' @export
#'
report_qc_sce <- function(sce,
                          report_folder_path = getwd()) {

  if(!class(sce) == "SingleCellExperiment"){
    stop("expecting singlecellexperiment")
  }

  cat(cli::rule(
    "Generating QC report for SingleCellExperiment", line = 2),
    "\r\n")

  metadata_tmp_path <- file.path(tempdir(), "metadata.rds")

  cli::cli_text("Writing data for QC report...")
  saveRDS(sce@metadata, metadata_tmp_path)

  cli::cli_text("Generating QC report...")
  rmarkdown::render(
    # for dev use file.path(getwd(),
    # "inst/rmarkdown/templates/quality-control/skeleton/skeleton.Rmd")
    system.file(
      "rmarkdown/templates/quality-control/skeleton/skeleton.Rmd",
      package = "scflow"),
    params = list(
      metadata_path = metadata_tmp_path
    ),
    output_dir = report_folder_path,
    output_file = "qc_report_scflow"
  )

  cli::cli_alert_success(
    "QC report succesfully generated. \r\n")

  return(sce)
}

