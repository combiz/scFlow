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
    "Generating QC plots and report for SingleCellExperiment", line = 2),
    "\r\n")

  # run all functions starting with .qc_plot_ on sce
  cli::cli_text("Generating QC plots and appending to metadata...")
  all_scflow_fns <- ls(getNamespace("scflow"), all.names=TRUE)
  qc_plot_fns <- all_scflow_fns[startsWith(all_scflow_fns, ".qc_plot_")]
  for (fn in qc_plot_fns) { sce <- get(fn)(sce) }

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
    "QC plots and report succesfully generated. \r\n")

  return(sce)
}

#' x axis barcode rank, y axis total counts
#' @keywords internal
.qc_plot_count_depth_distribution <- function(sce) {

  dt <- dplyr::as_tibble(data.frame(total_counts = sce$total_counts)) %>%
    dplyr::filter(total_counts > 10) %>%
    dplyr::arrange(-total_counts) %>%
    dplyr::mutate(barcode_rank = as.integer(rownames(.)))

  p <- ggplot2::ggplot(dt) +
    geom_point(aes(x = barcode_rank, y = total_counts))+
    scale_y_continuous(trans = 'log10') +
    scale_x_continuous(
      limits=c(0, max(dt$barcode_rank)),
      breaks = seq(0, max(dt$barcode_rank), by = 10000))+
    geom_hline(
      yintercept = sce@metadata$qc_params$min_library_size,
      linetype = "solid",
      color = "red")+
    labs(x = "Barcode rank", y = "Count depth")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(angle = 90),
          axis.title=element_text(size=16),
          legend.text=element_text(size=10),
          plot.title = element_text(size = 18, hjust = 0.5))

  sce@metadata$qc_plots$count_depth_distribution <- p
  sce@metadata$qc_plot_data$count_depth_distribution <- dt

  return(sce)

}

#' x axis count depth, y axis number of genes
#' @keywords internal
.qc_plot_features_vs_count_depth <- function(sce) {

  dt <- dplyr::as_tibble(data.frame(
    total_features = sce$total_features_by_counts,
    total_counts = sce$total_counts,
    pc_mito = sce$pc_mito)) %>%
    dplyr::filter(total_counts > 0)

  p <- ggplot2::ggplot(dt) +
    geom_point(aes(x = total_counts, y = total_features, colour = pc_mito), size = .01)+
    geom_hline(
      yintercept = sce@metadata$qc_params$min_features,
      linetype = "solid",
      color = "red") +
    geom_vline(
      xintercept = sce@metadata$qc_params$min_library_size,
      linetype = "solid",
      color = "red") +
    scale_colour_gradientn(limits = c(0, 0.10),
                           colours = c("darkblue", "red"),
                           na.value = "red",
                           name = "Relative mito") +
    labs(x = "Count depth", y = "Number of genes")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title=element_text(size=16),
          legend.text=element_text(size=10),
          plot.title = element_text(size = 18, hjust = 0.5))

  sce@metadata$qc_plots$number_genes_vs_count_depth <- p
  sce@metadata$qc_plot_data$number_genes_vs_count_depth <- dt

  return(sce)

}


#' x axis count depth, y axis number of genes
#' @keywords internal
.qc_plot_count_depth_histogram <- function(sce) {

  dt <- dplyr::as_tibble(data.frame(
    total_counts = sce$total_counts)) %>%
    filter(total_counts <= quantile(total_counts, c(.99)) ) %>%
    filter(total_counts > 10 )

  p <- ggplot2::ggplot(dt) +
    geom_histogram(aes(x = total_counts), bins = 100) +
    geom_vline(
      xintercept = sce@metadata$qc_params$min_library_size,
      linetype = "solid",
      color = "red") +
    labs(x = "Count depth", y = "Frequency")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title=element_text(size=16),
          legend.text=element_text(size=10),
          plot.title = element_text(size = 18, hjust = 0.5))

  sce@metadata$qc_plots$count_depth_histogram <- p
  sce@metadata$qc_plot_data$count_depth_histogram <- dt

  return(sce)

}


#' x axis count depth, y axis number of genes
#' @keywords internal
.qc_plot_number_genes_histogram <- function(sce) {

  dt <- dplyr::as_tibble(data.frame(
    total_features = sce$total_features_by_counts)) %>%
    filter(total_features <= quantile(total_features, c(.99))) %>%
    filter(total_features >= 10 )

  p <- ggplot2::ggplot(dt) +
    geom_histogram(aes(x = total_features), bins = 100) +
    geom_vline(
      xintercept = sce@metadata$qc_params$min_features,
      linetype = "solid",
      color = "red") +
    labs(x = "Number of genes", y = "Frequency")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title=element_text(size=16),
          legend.text=element_text(size=10),
          plot.title = element_text(size = 18, hjust = 0.5))

  sce@metadata$qc_plots$number_genes_histogram <- p
  sce@metadata$qc_plot_data$number_genes_histogram <- dt

  return(sce)

}



#' x axis count depth, y axis number of genes
#' @keywords internal
.qc_plot_mito_fraction_histogram <- function(sce) {

  dt <- dplyr::as_tibble(data.frame(
    pc_mito = sce$pc_mito)) %>%
    filter(pc_mito > 0 )

  p <- ggplot2::ggplot(dt) +
    geom_histogram(aes(x = pc_mito), bins = 100) +
    geom_vline(
      xintercept = sce@metadata$qc_params$max_mito,
      linetype = "solid",
      color = "red") +
    labs(x = "Fraction mitochondrial counts", y = "Frequency")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title=element_text(size=16),
          legend.text=element_text(size=10),
          plot.title = element_text(size = 18, hjust = 0.5))

  sce@metadata$qc_plots$mito_fraction_histogram <- p
  sce@metadata$qc_plot_data$mito_fraction_histogram <- dt

  return(sce)

}


#' x axis count depth, y axis number of genes
#' @keywords internal
.qc_plot_ribo_fraction_histogram <- function(sce) {

  dt <- dplyr::as_tibble(data.frame(
    pc_ribo = sce$pc_ribo)) %>%
    filter(pc_ribo > 0 )

  p <- ggplot2::ggplot(dt) +
    geom_histogram(aes(x = pc_ribo), bins = 100) +
    geom_vline(
      xintercept = sce@metadata$qc_params$max_ribo,
      linetype = "solid",
      color = "red") +
    labs(x = "Fraction ribosomal counts", y = "Frequency")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title=element_text(size=16),
          legend.text=element_text(size=10),
          plot.title = element_text(size = 18, hjust = 0.5))

  sce@metadata$qc_plots$ribo_fraction_histogram <- p
  sce@metadata$qc_plot_data$ribo_fraction_histogram <- dt

  return(sce)

}

#' temp
#' @keywords internal
.generate_rmd_report <- function(sce) {


}

