################################################################################
#' Generate plots and a QC report for a imppacted pathway analysis
#'
#' @param res Pahtway enrichment result
#' @param report_folder_path folder path to save the report
#' @param report_file filename for report (without an extension)
#'
#' @family annotation functions
#' @import cli Matrix dplyr purrr
#' @import ggplot2
#' @importFrom rmarkdown render
#' @importFrom tools file_path_sans_ext
#' @export
#'
report_impacted_pathway <- function(res,
                          report_folder_path = getwd(),
                          report_file = "impacted_pathway_report_scflow") {

  report_file <- tools::file_path_sans_ext(report_file)

  cat(cli::rule(
    "Generating report for impacted pathway", line = 2),
    "\r\n")

  metadata_tmp_path <- file.path(tempdir(), "metadata.rds")

  res$plot <- lapply(
    res$plot,
    function(p) {
      p$plot_env$dt <- NULL
      return(p)
    }
  )

  cli::cli_text("Writing temp files for report...")
  saveRDS(
    res,
    metadata_tmp_path
  )

  krd <- file.path(tempdir(), "krdqc")
  intd <- file.path(tempdir(), "idqc")
  dir.create(krd, showWarnings = FALSE)
  dir.create(intd, showWarnings = FALSE)

  cli::cli_text("Generating QC report...")
  rmarkdown::render(
    # for dev use file.path(getwd(),
    # "inst/rmarkdown/templates/impacted-pathway/skeleton/skeleton.Rmd")
    system.file(
      "rmarkdown/templates/impacted-pathway/skeleton/skeleton.Rmd",
      package = "scFlow"),
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
    "{cli::col_green(symbol$tick)} Report succesfully generated: ",
    "{.file {file.path(report_folder_path, 'impacted_pathway_report_scflow.html')}}")
  )

}

