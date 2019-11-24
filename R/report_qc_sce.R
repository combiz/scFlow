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
report_qc_sce <- function(sce,
                          report_folder_path = getwd(),
                          report_file = "qc_report_scflow") {

  if(!class(sce) == "SingleCellExperiment"){
    stop("expecting singlecellexperiment")
  }

  report_file <- tools::file_path_sans_ext(report_file)

  cat(cli::rule(
    "Generating QC report for SingleCellExperiment", line = 2),
    "\r\n")

  metadata_tmp_path <- file.path(tempdir(), "metadata.rds")

  sce@metadata$qc_plots <- lapply(
    sce@metadata$qc_plots,
    function(p) {
      p$plot_env$sce <- NULL
      return(p)
    }
  )

  cli::cli_text("Writing temp files for QC report...")
  saveRDS(
    sce@metadata[!(names(sce@metadata) %in% "qc_plot_data")],
    metadata_tmp_path
  )

  #biblio_tmp_path <- file.path(tempdir(), "scflow_references.bib")
  #bib2df::df2bib(sce@metadata$citations, biblio_tmp_path)

  krd <- file.path(tempdir(), "krdqc")
  intd <- file.path(tempdir(), "idqc")
  dir.create(krd)
  dir.create(intd)

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
    output_file = report_file,
    knit_root_dir = krd,
    intermediates_dir = intd,
    quiet = TRUE
  )

  cli::cli_text(c(
    "{cli::col_green(symbol$tick)} QC report succesfully generated: ",
    "{.file {file.path(report_folder_path, 'qc_report_scflow.html')}}")
  )

}

