################################################################################
#' Generate a report for impacted pathway analysis
#'
#' @param res Pathway enrichment result from find_impacted_pathway() function.
#' @param report_folder_path folder path to save the report.
#' @param report_file filename for report (without an extension).
#'
#' @family Impacted pathway analysis
#'
#' @import ggplot2
#' @importFrom cli cli_h2 cli_text col_green
#' @importFrom rmarkdown render
#' @importFrom tools file_path_sans_ext
#'
#' @export
#'
report_impacted_pathway <- function(res,
                                    report_folder_path = getwd(),
                                    report_file = "ipa_report_scflow") {
  report_file <- tools::file_path_sans_ext(report_file)

  cli::cli_h2("Generating report for impacted pathway")

  metadata_tmp_path <- file.path(report_folder_path, "metadata.rds")

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

  cli::cli_text("Generating impacted pathway analysis report...")
  rmarkdown::render(
    system.file(
      "rmarkdown/templates/impacted-pathway/skeleton/skeleton.Rmd",
      package = "scFlow"
    ),
    params = list(
      metadata_path = metadata_tmp_path
    ),
    output_dir = report_folder_path,
    output_file = report_file,
    knit_root_dir = report_folder_path,
    intermediates_dir = report_folder_path,
    quiet = TRUE
  )

  report_file_name <- paste(report_file, ".html", sep = "")

  cli::cli_text(c(
    "{cli::col_green(symbol$tick)} Report succesfully generated: ",
    "{.file {file.path(report_folder_path, report_file_name)}}"
  ))

  unlink(metadata_tmp_path)


}
