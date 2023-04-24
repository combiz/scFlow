################################################################################
#' Generate a report for dataset integration, dims reduction, and clustering
#'
#' @param sce a SingleCellExperiment object
#' @param report_folder_path folder path to save the report
#' @param report_file filename for report (without an extension)
#' @param categorical_covariates list of categorical variables
#' @return sce SingleCellExperiment object annotated with reducedDims
#'
#' @family Data integration
#' @import ggplot2
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom rmarkdown render
#' @importFrom tools file_path_sans_ext
#' @importFrom formattable formattable icontext
#' @export
#'
report_integrated_sce <- function(sce,
                                  report_folder_path = getwd(),
                                  report_file = "integrate_scFlow",
                                  categorical_covariates = list()) {
  if (!class(sce) == "SingleCellExperiment") {
    stop("expecting singlecellexperiment")
  }
  report_file <- tools::file_path_sans_ext(report_file)
  cat(cli::rule(
    "Generating Report for Dataset Integration, Dimension Reduction, and Clustering",
    line = 2),
    "\r\n")
  metadata_tmp_path <- file.path(report_folder_path, "metadata.qs")
  cli::cli_text("Writing temp files for report...")
  qs::qsave(
    sce@metadata,
    metadata_tmp_path
  )


  cli::cli_text("Generating Dataset Integration, Dimension Reduction, and Clustering report...")
  rmarkdown::render(
    system.file(
      "rmarkdown/templates/integrate/skeleton/skeleton.Rmd",
      package = "scFlow"),
    params = list(
      metadata_path = metadata_tmp_path,
      categorical_covariates = categorical_covariates
    ),
    output_dir = report_folder_path,
    output_file = report_file,
    knit_root_dir = report_folder_path,
    intermediates_dir = report_folder_path,
    quiet = TRUE
  )
  cli::cli_text(c(
    "Report successfully generated: ",
    "{.file {file.path(report_folder_path, 'integrate_scFlow.html')}}")
  )

  unlink(metadata_tmp_path)

  return(sce)
}
