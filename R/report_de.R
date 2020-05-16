################################################################################
#' Generate a report for differential expression analysis
#'
#' @param res Differential expression result table from perform_de() function.
#' @param report_folder_path folder path to save the report.
#' @param report_file filename for report (without an extension).
#'
#' @family differential gene expression
#'
#' @import ggplot2
#' @importFrom cli cli_h2 cli_text col_green
#' @importFrom rmarkdown render
#' @importFrom tools file_path_sans_ext
#'
#' @export
#'
report_de <- function(res,
                      report_folder_path = getwd(),
                      report_file = "de_report_scflow") {
  report_file <- tools::file_path_sans_ext(report_file)

  cli::cli_h2("Generating report for differential expression analysis")

  metadata_tmp_path <- file.path(tempdir(), "metadata.qs")

  cli::cli_text("Writing temp files for report...")
  qs::qsave(
    res,
    metadata_tmp_path
  )

  krd <- file.path(tempdir(), "krdqc")
  intd <- file.path(tempdir(), "idqc")
  dir.create(krd, showWarnings = FALSE)
  dir.create(intd, showWarnings = FALSE)

  cli::cli_text("Generating differential expression analysis report...")
  rmarkdown::render(
    system.file(
      "rmarkdown/templates/differential-expression/skeleton/skeleton.Rmd",
      package = "scFlow"
    ),
    params = list(
      metadata_path = metadata_tmp_path
    ),
    output_dir = report_folder_path,
    output_file = report_file,
    knit_root_dir = krd,
    intermediates_dir = intd,
    quiet = TRUE
  )

  report_file_name <- paste(report_file, ".html", sep = "")

  cli::cli_text(c(
    "{cli::col_green(symbol$tick)} Report succesfully generated: ",
    "{.file {file.path(report_folder_path, report_file_name)}}"
  ))
}
