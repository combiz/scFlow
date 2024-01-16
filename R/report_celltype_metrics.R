################################################################################
#' Generate A Report of Cell-Type Metrics for a SingleCellExperiment
#'
#' @param sce a SingleCellExperiment object
#' @param report_folder_path folder path to save the report
#' @param report_file filename for report (without an extension)
#'
#' @return sce a annotated SingleCellExperiment object
#'
#' @family annotation functions
#' @import ggplot2
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom rmarkdown render
#' @importFrom tools file_path_sans_ext
#' @importFrom formattable formattable icontext
#' @importFrom assertthat assert_that
#' @importFrom ggridges geom_density_ridges_gradient
#' @importFrom cli cli_h1 col_green symbol
#' @importFrom qs qsave
#' @export
#'
report_celltype_metrics <- function(sce,
                                    report_folder_path = getwd(),
                                    report_file = "celltype_metrics_scflow") {

  assertthat::assert_that(class(sce) == "SingleCellExperiment")
  assertthat::assert_that(
    class(sce@metadata$scflow_steps) == "list",
    msg = "Before producing the report, run annotate_celltype_metrics()")
  assertthat::assert_that(
    sce@metadata$scflow_steps$celltype_plots_annotated == TRUE,
    msg = "Before producing the report, run annotate_celltype_metrics()")

  report_file <- tools::file_path_sans_ext(report_file)

  cli::cli_h1("Generating Report for Cell-Type Metrics")

  metadata_tmp_path <- file.path(report_folder_path, "metadata.qs")

  cli::cli_text("Writing temp files for report...")
  qs::qsave(
    sce@metadata,
    metadata_tmp_path
  )


  cli::cli_text("Generating report...")
  rmarkdown::render(
    # for dev use file.path(getwd(),
    # "inst/rmarkdown/templates/quality-control/skeleton/skeleton.Rmd")
    system.file(
      "rmarkdown/templates/celltype-metrics/skeleton/skeleton.Rmd",
      package = "scFlow"),
    params = list(
      metadata_path = metadata_tmp_path
    ),
    output_dir = report_folder_path,
    output_file = report_file,
    knit_root_dir = report_folder_path,
    intermediates_dir = report_folder_path,
    quiet = TRUE
  )

  file.remove(metadata_tmp_path)

  cli::cli_text(c(
    "{cli::col_green(cli::symbol$tick)} Report succesfully generated: ",
    "{.file {file.path(report_folder_path, paste0(report_file, '.html'))}}")
  )

  unlink(metadata_tmp_path)

  return(sce)

}
