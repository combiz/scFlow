################################################################################
#' Generate plots and a QC report for a SingleCellExperiment
#'
#' @param sce a SingleCellExperiment object
#' @param report_folder_path folder path to save the report
#' @param report_file filename for report (without an extension)
#'
#' @return sce a annotated SingleCellExperiment object
#'
#' @family annotation functions
#' @import cli Matrix dplyr SingleCellExperiment purrr
#' @import ggplot2
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom rmarkdown render
#' @importFrom tools file_path_sans_ext
#' @export
#'
report_merged_sce <- function(sce,
                              report_folder_path = getwd(),
                              report_file = "merged_report_scflow") {

  if(!class(sce) == "SingleCellExperiment"){
    stop("expecting singlecellexperiment")
  }

  report_file <- tools::file_path_sans_ext(report_file)

  cat(cli::rule(
    "Generating merged QC report for SingleCellExperiment", line = 2),
    "\r\n")

  metadata_tmp_path <- file.path(tempdir(), "metadata.rds")

  cli::cli_text("Writing temp files for report...")
  saveRDS(
    sce@metadata,
    metadata_tmp_path
  )

  krd <- file.path(tempdir(), "krdqc")
  intd <- file.path(tempdir(), "idqc")
  dir.create(krd, showWarnings = FALSE)
  dir.create(intd, showWarnings = FALSE)

  cli::cli_text("Generating merged QC report...")
  rmarkdown::render(
    # for dev use file.path(getwd(),
    # "inst/rmarkdown/templates/quality-control/skeleton/skeleton.Rmd")
    system.file(
      "rmarkdown/templates/merged-quality-control/skeleton/skeleton.Rmd",
      package = "scflow"),
    params = list(
      metadata_path = metadata_tmp_path
    ),
    output_dir = report_folder_path,
    output_file = report_file,
    knit_root_dir = krd,
    intermediates_dir = intd,
    quiet = TRUE
  )

  cli::cli_text(c(
    "{cli::col_green(symbol$tick)} merged QC report succesfully generated: ",
    "{.file {file.path(report_folder_path, 'merged_qc_report_scflow.html')}}")
  )

  return(sce)

}

