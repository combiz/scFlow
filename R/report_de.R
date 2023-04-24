################################################################################
#' Generate a report for differential expression analysis
#'
#' @param res Differential expression result table from perform_de() function.
#' @param logFC_threshold fold change up/down threshold.
#' @param padj_cutoff the adjusted pvalue cutoff threshold.
#' @param n_label number of genes to be labelled in the volcano plot.
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
                      logFC_threshold = 0.25,
                      padj_cutoff = 0.05,
                      n_label = 5,
                      report_folder_path = getwd(),
                      report_file = "de_report_scflow") {
  report_file <- tools::file_path_sans_ext(report_file)

  cli::cli_h2("Generating report for differential expression analysis")

  p <- volcano_plot(
    dt = res,
    logFC_threshold = logFC_threshold,
    padj_cutoff = padj_cutoff,
    n_label = n_label
  )

  DGEs <- c(
    res %>%
      filter(padj <= padj_cutoff, logFC >= logFC_threshold) %>%
      pull(gene) %>%
      length(),
    res %>%
      filter(padj <= padj_cutoff, logFC <= -logFC_threshold) %>%
      pull(gene) %>%
      length()
  )

  names(DGEs) <- c("Up", "Down")

  attr(res, "de_result") <- DGEs

  attr(res, "plot") <- p

  attr(res, "report_params") <- setNames(
    c(logFC_threshold, padj_cutoff),
    c("logFC_threshold", "padj_cutoff")
  )


  metadata_tmp_path <- file.path(report_folder_path, "metadata.qs")

  cli::cli_text("Writing temp files for report...")
  qs::qsave(
    res,
    metadata_tmp_path
  )

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
