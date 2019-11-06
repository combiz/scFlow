################################################################################
#' Annotate a SingleCellExperiment with gene names and QC metrics
#'
#' Adds biomaRt annotations (e.g. gene, gene_biotype) and QC metric annotations.
#'
#' @section Quality control options and thresholds:
#'
#' In addition to calculating QC metrics and annotating gene information, this
#' function adds boolean (TRUE/FALSE) indicators of which cells/genes met the QC
#' criteria.  This enables QC reports, plots, and various QC-related tables to
#' be saved before filtering with the  [filter_sce()] function.
#'
#' @section Annotations:
#'
#' With the default settings, the SingleCellExperiment object is annotated with:
#'
#' Cell-level annotations
#' * total_counts - sum of counts across all genes
#' * total_features_by_counts - total number of unique genes with expression >0
#' * qc_metric_min_library_size - did the cell have at least min_library_size
#' counts
#' * qc_metric_min_features - did the cell have counts >0 in at least
#' min_features number of cells?
#' * pc_mito - percentage of counts mapping to mitochondrial genes in this cell
#' * qc_metric_pc_mito_ok was pc_mito <= the max_mito cutoff?
#' * pc_ribo - percentage of counts mapping to ribosomal genes in this cell
#' * qc_metric_pc_ribo_ok was pc_ribo <= the max_ribo cutoff?
#' * qc_metric_passed - did the cell pass all of the cell QC tests
#'
#' Gene-level annotations
#' * gene - official gene name
#' * gene_biotype - protein_coding, lncRNA, pseudogene, etc.
#' * qc_metric_ensembl_mapped - was the ensembl_gene_id found in biomaRt
#' * qc_metric_is_mito - is the gene mitochondrial
#' * qc_metric_is_ribo - is the gene ribosomal
#' * qc_metric_n_cells_expressing - number of cells with at least min_counts
#' * qc_metric_is_expressive - did at least min_cells have min_counts?
#'
#' @param sce a SingleCellExperiment object
#' @param min_library_size the minimum number of counts per cell
#' @param min_features the minimum number of features per cell (i.e. the minimum
#'   number of genes with >0 counts)
#' @param max_mito the maximum proportion of counts mapping to
#'   mitochondrial genes (0 - 1)
#' @param min_ribo the minimum proportion of counts mapping to
#'   ribosomal genes (0 - 1)
#' @param max_ribo the maximum proportion of counts mapping to
#'   ribosomal genes (0 - 1)
#' @param min_counts the minimum number of counts per cell in min_cells
#' @param min_cells the minimum number of cells with min_counts
#' @param drop_unmapped set `TRUE` to remove unmapped ensembl_gene_id
#' @param drop_mito set `TRUE` to remove mitochondrial genes
#' @param drop_ribo set `TRUE` to remove ribosomal genes
#' @param ensembl_mapping_file a local tsv file with ensembl_gene_id and
#'   additional columns for mapping ensembl_gene_id to gene info.  If
#'   not provided, the biomaRt db is queried (slower).
#' @param annotate_genes optionally skip gene annotation with FALSE
#' @param annotate_cells optionally skip cell annotation with FALSE
#'
#' @return sce a annotated SingleCellExperiment object
#'
#' @family annotation functions
#' @import cli Matrix dplyr SingleCellExperiment purrr
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#' @export
annotate_sce <- function(sce,
                         min_library_size = 300,
                         min_features = 100,
                         max_mito = 0.10,
                         min_ribo = 0.00,
                         max_ribo = 1.00,
                         min_counts = 2,
                         min_cells = 2,
                         drop_unmapped = TRUE,
                         drop_mito = TRUE,
                         drop_ribo = FALSE,
                         annotate_genes = TRUE,
                         annotate_cells = TRUE,
                         ensembl_mapping_file = NULL) {

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  cat(cli::rule("Annotating SingleCellExperiment", line = 2), "\r\n")

  before_coldata_colnames <- colnames(sce@colData)
  before_rowdata_colnames <- colnames(SummarizedExperiment::rowData(sce))

  # add the qc parameters to the metadata
  qc_params <- setdiff(names(formals(annotate_sce)),
                       c("sce", "ensembl_mapping_file",
                         "annotate_genes", "annotate_cells")) #not these args

  qc_params_l <- purrr::map(qc_params, ~ get(.))
  qc_params_l <- purrr::set_names(qc_params_l, qc_params)
  sce@metadata[["qc_params"]] <- qc_params_l

  if (annotate_genes) {
    sce <- annotate_sce_genes(
      sce,
      drop_unmapped = drop_unmapped,
      drop_mito = drop_mito,
      drop_ribo = drop_ribo,
      ensembl_mapping_file = ensembl_mapping_file)
  }
  if (annotate_cells) {
    sce <- annotate_sce_cells(
      sce,
      min_library_size = min_library_size,
      min_features = min_features,
      max_mito = max_mito,
      min_ribo = min_ribo,
      max_ribo = max_ribo,
      min_counts = min_counts,
      min_cells = min_cells
    )
  } else {
    if (!annotate_genes) {
      stop(cli::cli_alert_danger("Nothing to do. Specify gene/cell/both."))
    }
  }

  if (annotate_cells) {
    cli::cli_alert_success(
      "SingleCellExperiment cells were successfully annotated with: \r\n"
    )

    cli::cli_ul(setdiff(
      colnames(sce@colData),
      before_coldata_colnames
    ))
  }

  if (annotate_genes) {
    cli::cli_alert_success(
      "SingleCellExperiment genes were successfully annotated with: \r\n"
    )

    cli::cli_ul(setdiff(
      colnames(SummarizedExperiment::rowData(sce)),
      before_rowdata_colnames
    ))
  }

  cat(cli::rule(
    "Generating QC plots for SingleCellExperiment", line = 2),
    "\r\n")

  # generate plots - run all functions starting with .qc_plot_ on sce
  cli::cli_text("Generating QC plots and appending to metadata.")
  all_scflow_fns <- ls(getNamespace("scflow"), all.names=TRUE)
  qc_plot_fns <- all_scflow_fns[startsWith(all_scflow_fns, ".qc_plot_")]
  for (fn in qc_plot_fns) { sce <- get(fn)(sce) }

  # generate qc summary table
  cli::cli_text("Generating QC summary table and appending to metadata.")
  sce <- .qc_append_summary_table(sce)

  return(sce)

}

#' generate table of qc results
#' @keywords internal
.qc_append_summary_table <- function(sce) {

  # qc params to data frame
  qc_params_df <- map_df(sce@metadata$qc_params, ~ .) %>%
    dplyr::rename_all(~ paste0("qc_params_",.))

  # gene qc results to data frame
  genes_qc <- list()
  genes_qc$n_genes <- dim(sce)[[1]]
  genes_qc$n_mito <-
    sum(SummarizedExperiment::rowData(sce)$qc_metric_is_mito, na.rm=T)
  genes_qc$n_ribo <-
    sum(SummarizedExperiment::rowData(sce)$qc_metric_is_ribo, na.rm=T)
  genes_qc$n_unmapped <-
    sum(!SummarizedExperiment::rowData(sce)$qc_metric_ensembl_mapped, na.rm=T)

  genes_qc$n_expressive <- sum(
    SummarizedExperiment::rowData(sce)$qc_metric_is_expressive, na.rm=T)

  genes_qc$n_nonexpressive <- sum(
    !SummarizedExperiment::rowData(sce)$qc_metric_is_expressive, na.rm=T)

  genes_qc$n_unmapped_dropped <- sum(
    !SummarizedExperiment::rowData(sce)$qc_metric_mapped_keep, na.rm=T)

  genes_qc$n_mito_dropped <- sum(
    !SummarizedExperiment::rowData(sce)$qc_metric_mito_keep, na.rm=T)

  genes_qc$n_ribo_dropped <- sum(
    !SummarizedExperiment::rowData(sce)$qc_metric_ribo_keep, na.rm=T)

  genes_qc$n_genes_passed <-
    sum(SummarizedExperiment::rowData(sce)$qc_metric_gene_passed, na.rm=T)
  genes_qc$n_genes_failed <-
    sum(!SummarizedExperiment::rowData(sce)$qc_metric_gene_passed, na.rm=T)

  genes_qc_df <- purrr::map_df(genes_qc, ~ .) %>%
    dplyr::rename_all(~ paste0("qc_genes_",.))

  # cell qc results to data frame

  cells_qc <- list()
  cells_qc$n_cells <- dim(sce)[[2]]
  cells_qc$n_library_size_passed <-
    sum(sce$qc_metric_min_library_size, na.rm=T)
  cells_qc$n_library_size_failed <-
    sum(!sce$qc_metric_min_library_size, na.rm=T)
  cells_qc$n_min_features_passed <-
    sum(sce$qc_metric_min_features, na.rm=T)
  cells_qc$n_min_features_failed <-
    sum(!sce$qc_metric_min_features, na.rm=T)
  cells_qc$n_mito_fraction_passed <-
    sum(sce$qc_metric_pc_mito_ok, na.rm=T)
  cells_qc$n_mito_fraction_failed <-
    sum(!sce$qc_metric_pc_mito_ok, na.rm=T)
  cells_qc$n_ribo_fraction_passed <-
    sum(sce$qc_metric_pc_ribo_ok, na.rm=T)
  cells_qc$n_ribo_fraction_failed <-
    sum(!sce$qc_metric_pc_ribo_ok, na.rm=T)
  cells_qc$n_cells_passed <-
    sum(sce$qc_metric_passed, na.rm=T)
  cells_qc$n_cells_failed <-
    sum(!sce$qc_metric_passed, na.rm=T)

  cells_qc_df <- purrr::map_df(cells_qc, ~ .) %>%
    dplyr::rename_all(~ paste0("qc_cells_",.))

  qc_summary <- cbind(
    #sce@metadata$metadata,
    qc_params_df,
    genes_qc_df,
    cells_qc_df
  )

  rownames(qc_summary) <- NULL
  sce@metadata$qc_summary <- qc_summary

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
    filter(total_counts <= stats::quantile(total_counts, c(.99)) ) %>%
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
      xintercept = sce@metadata$qc_params$min_ribo,
      linetype = "solid",
      color = "red") +
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
