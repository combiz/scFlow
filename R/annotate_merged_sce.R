################################################################################
#' Annotate a post-merge SingleCellExperiment with plots
#'
#' @param sce a SingleCellExperiment object
#' @param plot_vars the colData variable(s) to generate plots for
#' @param unique_id_var the colData variable identifying unique samples
#' @param facet_vars the colData variable(s) to facet/subset by
#' @param outlier_vars the colData variable(s) to apply adaptive thresholding
#' @param outlier_mads the minimum median absolute deviations for thresholding
#'
#' @return sce a annotated SingleCellExperiment object
#'
#' @family annotation functions
#' @import cli Matrix dplyr SingleCellExperiment purrr
#' @import ggplot2
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom rmarkdown render
#' @importFrom purrr map_lgl
#' @importFrom tools file_path_sans_ext
#' @importFrom tidyr pivot_longer
#' @export
annotate_merged_sce <- function(sce,
                                plot_vars = c("total_features_by_counts",
                                            "total_counts", "pc_mito",
                                            "pc_ribo"),
                                unique_id_var = "manifest",
                                facet_vars = NULL,
                                outlier_vars = c("total_features_by_counts",
                                                 "total_counts"),
                                outlier_mads = 3) {

  cat(cli::rule("Annotating Merged SingleCellExperiment", line = 2), "\r\n")
  cat(cli::rule("Appending Merge Summary Plots", line = 1), "\r\n")

  plot_vars <- unique(
    c("total_features_by_counts", "total_counts", "pc_mito", "pc_ribo",
      plot_vars)
  )

  # generate cell metadata plots
  merged_plots_l <- list()
  merged_plots_data_l <- list()
  cli::cli_alert_success(
    "SingleCellExperiment successfully annotated with merge summary plots: \r\n"
  )

  for (pv in plot_vars) {
    # generate plot without faceting
    merged_plots_l[[pv]][[pv]] <- .generate_merge_summary_plot(
      sce,
      plot_var = pv,
      unique_id_var = unique_id_var,
      facet_var = NULL,
      plot_points = TRUE) # no facets

    # generate plot data table for export
    dt <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
      dplyr::select(unique(c(pv, unique_id_var))) %>%
      dplyr::group_by_at(unique_id_var) %>%
      dplyr::summarize(
        mean_avg = round(mean(!!rlang::sym(pv), na.rm = TRUE), digits = 3),
        stdev_mean = round(sd(!!rlang::sym(pv), na.rm = TRUE), digits = 3),
        median_avg = round(median(!!rlang::sym(pv), na.rm = TRUE), digits = 3),
        mad = round(mad(!!rlang::sym(pv), na.rm = TRUE), digits = 3)
        ) %>%
      mutate(z = scale(mean_avg)[, 1])

    merged_plots_data_l[[pv]][[pv]] <- dt

    cli::cli_ul(sprintf("sce@metadata$merged_plots$%s$%s", pv, pv))

    for (fv in facet_vars) {
      # generate plot with faceting
      plot_name <- paste(pv, fv, sep = "_vs_")
      merged_plots_l[[pv]][[plot_name]] <- .generate_merge_summary_plot(
        sce,
        plot_var = pv,
        unique_id_var = unique_id_var,
        facet_var = fv,
        plot_points = TRUE)

      # generate plot data table for export
      dt <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
        dplyr::select(unique(c(pv, unique_id_var, fv))) %>%
        dplyr::group_by_at(fv) %>%
        dplyr::summarize(
          mean_avg = round(
            mean(!!rlang::sym(pv), na.rm = TRUE), digits = 3),
          stdev_mean = round(
            sd(!!rlang::sym(pv), na.rm = TRUE), digits = 3),
          median_avg = round(
            median(!!rlang::sym(pv), na.rm = TRUE), digits = 3),
          mad = round(
            mad(!!rlang::sym(pv), na.rm = TRUE), digits = 3)
          ) %>%
        mutate(z = scale(mean_avg)[, 1])

      merged_plots_data_l[[pv]][[plot_name]] <- dt

      cli::cli_ul(sprintf("sce@metadata$merged_plots$%s$%s", pv, plot_name))
    }
  }

  sce@metadata$merged_plots <- merged_plots_l
  sce@metadata$merged_plots_data <- merged_plots_data_l

  cat(cli::rule(
    "Generating Pseudobulking Matrices and Plots", line = 1), "\r\n")
  # generate whole sample pseudobulk and plots
  x <- Sys.time()
  message(sprintf("Performing Pseudobulking of Cells by %s", unique_id_var))
  pbsce <- .generate_pbsce_by_id(
    sce,
    unique_id_var = unique_id_var
    )
  time_taken <- difftime(Sys.time(), x, units = "mins")
  cli::cli_alert_success(sprintf(
    "Pseudobulking completed (%.2f minutes taken). \r\n",
    time_taken)
  )

  cli::cli_alert_success(
    "SingleCellExperiment successfully annotated with pseudobulk plots: \r\n"
  )
  # generate reduced dimension pseudobulk sample plots
  sce@metadata$pseudobulk_rd_plots <- list()
  for (rd_method in SingleCellExperiment::reducedDimNames(pbsce)) {
    p <- plot_reduced_dim(
      pbsce,
      feature_dim = unique_id_var,
      reduced_dim = rd_method,
      size = 6,
      alpha = 1, label_clusters = TRUE
      )
    p$plot_env$sce <- NULL
    sce@metadata$pseudobulk_rd_plots[[rd_method]] <- p
    cli::cli_ul(sprintf("sce@metadata$pseudobulk_rd_plots$%s", rd_method))
  }

  # dendrogram and heatmap of expressed genes
  pbsce <- .plot_heatmap_of_pbsce(pbsce, binarize = TRUE, trim_name = TRUE)
  sce@metadata$pseudobulk_plots <- pbsce@metadata$pseudobulk_plots
  sce@metadata$pseudobulk_data <- pbsce@metadata$pseudobulk_data
  cli::cli_ul("sce@metadata$pseudobulk_plots$combined_heatmap")

  # extend pseudobulk reducedDim coordinates to all cells and save to sce
  for (rd_method in SingleCellExperiment::reducedDimNames(pbsce)) {
    mat <- as.data.frame(SingleCellExperiment::reducedDim(pbsce, rd_method))
    mat$id <- as.character(rownames(mat))
    big_mat <- dplyr::left_join(
      data.frame(
        "id" = as.character(sce[[unique_id_var]]),
        stringsAsFactors = FALSE),
      mat,
      by = "id")
    rownames(big_mat) <- sce$barcode
    big_mat$id <- NULL
    big_mat <- as.matrix(big_mat)
    SingleCellExperiment::reducedDim(sce, paste0(rd_method, "_PB")) <- big_mat
  }

  pb_reddims <- SingleCellExperiment::reducedDimNames(
    sce)[purrr::map_lgl(SingleCellExperiment::reducedDimNames(sce),
                        ~ endsWith(., "_PB"))]
  gt <- cli::col_green(cli::symbol$tick)
  cli::cli_text(sprintf(c(
    gt,
    " Pseudobulk reducedDim matrices {.var (%s)} successfully appended ",
    "to SingleCellExperiment."), paste0(pb_reddims, collapse = ", ")))

  sce@metadata$scflow_steps$merged_annotated <- 1
  sce@metadata$merge_qc_params$plot_vars <- plot_vars
  sce@metadata$merge_qc_params$unique_id_var <- unique_id_var
  sce@metadata$merge_qc_params$facet_vars <- facet_vars

  cli::cli_alert_success("Done! \r\n")

  return(sce)

}

#' helper fn to generate a pseudobulk sce for post-merge qc
#'
#' @param sce a singlecellexperiment object
#' @param unique_id_var the colData variable identifying unique samples
#' @keywords internal
.generate_pbsce_by_id <- function(sce,
                                  unique_id_var = "manifest") {

  pbsce <- pseudobulk_sce(
    sce,
    sample_var = unique_id_var,
    celltype_var = unique_id_var,
    keep_vars = unique_id_var
  )

  colnames(pbsce) <- purrr::map_chr(
    colnames(pbsce), ~ strsplit(., "_")[[1]][[1]])

  as.data.frame(SummarizedExperiment::colData(pbsce))
  colnames(pbsce) <- pbsce[[unique_id_var]]

  pbsce <- reduce_dims_sce(
    pbsce,
    vars_to_regress_out = c("n_cells"),
    pca_dims = 5,
    n_neighbors = 2,
    reduction_methods = c("PCA", "UMAP")
  )

  return(pbsce)

}

#' helper fn to generate a heatmap and dendrogram for features of a pbsce
#'
#' @param pbsce a pseudobulked singlecellexperiment object
#' @param binarize if TRUE, counts are TRUE/FALSE if counts > 0
#' @param trim_name if TRUE, the pseudobulked name is trimmed to before _
#'
#' @import ggplot2
#' @importFrom SingleCellExperiment counts
#' @importFrom purrr map_chr
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom stats hclust
#' @importFrom ggdendro ggdendrogram
#' @importFrom ggpubr ggarrange
#' @keywords internal
.plot_heatmap_of_pbsce <- function(pbsce, binarize = TRUE, trim_name = TRUE) {

  if (binarize) {
    dt <- as.data.frame(SingleCellExperiment::counts(pbsce) > 0)
  } else {
    dt <- as.data.frame(SingleCellExperiment::counts(pbsce))
  }

  if (trim_name) {
    colnames(dt) <- purrr::map_chr(colnames(dt), ~ strsplit(., "_")[[1]][[1]])
  }

  clust <- stats::hclust(dist(t(dt)))

  dt_long <- dt %>%
    dplyr::mutate(ensembl_gene_id = rownames(.)) %>%
    tidyr::pivot_longer(-ensembl_gene_id)

  dt_long$name <- factor(x = dt_long$name,
                         levels = dt_long$name[clust$order],
                         ordered = TRUE)

  p1 <- ggdendro::ggdendrogram(clust, labels = FALSE, leaf_labels = FALSE) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "lines"))

  p2 <- ggplot(dt_long, aes(x = name, y = ensembl_gene_id)) +
    geom_tile(aes(fill = !value)) +
    scale_fill_viridis_d() +
    xlab("Sample") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 16),
          axis.title = element_text(size = 18),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_blank(),
          plot.margin = unit(c(-1, 1.5, 0, 1.5), "lines"))

  p <- ggpubr::ggarrange(p1, p2,
                         ncol = 1, nrow = 2,
                         heights = c(0.5, 2),
                         align = "none")

  pbsce@metadata$pseudobulk_data$dt_long <- dt_long
  pbsce@metadata$pseudobulk_data$dt <- dt
  pbsce@metadata$pseudobulk_data$clust <- clust

  pbsce@metadata$pseudobulk_plots$dendrogram <- p1
  pbsce@metadata$pseudobulk_plots$heatmap <- p2
  pbsce@metadata$pseudobulk_plots$combined_heatmap <- p

  return(pbsce)
}



#' helper fn to generate a merge summary plot for a SingleCellExperiment
#'
#' if plot_var is a factor, a bar plot is returned.
#' if plot_var is.numeric, a violin plot is returned.
#'
#' @param sce a singlecellexperiment object
#' @param plot_var the colData variable with values to plot
#' @param unique_id_var the colData variable identifying unique samples
#' @param facet_var the colData variable to facet/subset by
#' @param plot_points if TRUE plot individual values for violin plots
#' @keywords internal
.generate_merge_summary_plot <- function(sce,
                                         plot_var,
                                         unique_id_var,
                                         facet_var = NULL,
                                         plot_points = FALSE) {

  vars_to_keep <- unique(c(plot_var, unique_id_var, facet_var))
  dt <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    select(vars_to_keep)
  dt$plot_var <- dt[[plot_var]]
  if (!is.null(facet_var)) dt$facet_var <- as.factor(dt[[facet_var]])

  if (startsWith(plot_var, "pc_")) {
    scale_y <- scale_y_continuous(labels = scales::percent)
  } else {
    scale_y <- scale_y_continuous(breaks = scales::pretty_breaks())
  }

  if (class(dt$plot_var) == "factor") {

    p <- ggplot(dt, aes(x = plot_var)) +
      geom_bar() +
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      xlab(plot_var) +
      theme_bw() +
      theme(
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)
      )

    if (!is.null(facet_var)) {
      p <- p + facet_grid(~ get(facet_var))
    }

  }

  if (is.numeric(dt$plot_var)) {

    if (is.null(facet_var)) dt$facet_var <- dt[[unique_id_var]]

    p <- ggplot(dt, aes(x = facet_var, y = plot_var)) +
      geom_violin(fill = "grey40", trim = TRUE) +
      scale_y +
      ylab(plot_var) +
      xlab(facet_var) +
      theme_bw() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
            text = element_text(size = 16),
            axis.text = element_text(size = 16),
            axis.text.x = element_text(angle = 90, hjust = 0.5),
            axis.title = element_text(size = 18)
            )

      if (plot_points == TRUE) {
        p <- p + geom_jitter(size = .01, width = .2, alpha = .05)
      } else {
        p <- p + stat_summary(fun.y = median,
                              fun.ymin = function(x) max(0, mean(x) - sd(x)),
                              fun.ymax = function(x) mean(x) + sd(x),
                              geom = "crossbar", width = 0.1, fill = "white")
      }
  }

  p$plot_env$sce <- NULL
  return(p)
}