library(magrittr)
# 3.3.1 bad
devtools::install_version("ggplot2", version = "3.1.1", repos = "http://cran.us.r-project.org")

unique_id_var <- "individual"
celltype_var <- "cluster_celltype"
metric_var <- "pc_mito"
var <- "cluster_celltype"
facet_var <- "pc_mito"
group_by_var <- "diagnosis"

e <- new.env()
e$sce <- sce



sce <- do.call(
  .append_cell_metric_plots_sce_test,
  list(
    sce = sce,
    celltype_var = var,
    unique_id_var = unique_id_var,
    facet_var = facet_var,
    group_by_var = group_by_var,
    metric_var = metric_var
  ),
  envir = e,
  quote = TRUE
)

sce <- .append_cell_metric_plots_sce_test(
    sce = sce,
    celltype_var = var,
    unique_id_var = unique_id_var,
    facet_var = facet_var,
    group_by_var = group_by_var,
    metric_var = metric_var
)

weigh(x@metadata$celltype_annotations$metric_plots$pc_mito$cluster_celltype$metric_plot)

weigh(x@metadata$celltype_annotations$metric_plots$pc_mito$cluster_celltype$ridge_plot)


.append_cell_metric_plots_sce_test <- function(sce, celltype_var = var,
                                               unique_id_var = unique_id_var,
                                               facet_var = facet_var,
                                               group_by_var = group_by_var,
                                               metric_var = metric_var, ...) {

  # metric_var is pc_mito, total_counts, etc.
  #list2env(list(...), environment())

  print(environment())
  print(ls(environment()))
  dt <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::group_by(!!rlang::sym(celltype_var)) %>%
    dplyr::summarize(
      mean = mean(!!rlang::sym(metric_var)),
      sd = sd(!!rlang::sym(metric_var)),
      se = sd(!!rlang::sym(metric_var) / sqrt(dplyr::n())),
      median = median(!!rlang::sym(metric_var))
    ) %>%
    dplyr::mutate(ymin = mean - se,
                  ymax = mean + se)

  dt[[celltype_var]] <- forcats::fct_reorder(as.factor(dt[[celltype_var]]), -dt$mean)
  dt$metric_var <- metric_var # for export context

  p <- ggplot2::ggplot(
    dt,
    ggplot2::aes(
      x = .data[[celltype_var]], y = mean, ymin = ymin, ymax = ymax
    )
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_errorbar(#ggplot2::aes(ymin = ymin, ymax = ymax),
                           width = .2,
                           position = ggplot2::position_dodge(.9)
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::ylab(paste(metric_var)) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 16),
      axis.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 10)
    )

  dt2 <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::select(.data[[celltype_var]], .data[[metric_var]])

  # reorder based on the dt means
  dt2[[celltype_var]] <- factor(
    dt2[[celltype_var]],
    levels = levels(dt[[celltype_var]]) #from the previous dt
  )

  p2 <- ggplot2::ggplot(dt2,
                        ggplot2::aes(x = rlang::sym(metric_var),
                                     y = rlang::sym(celltype_var))) +
    ggridges::geom_density_ridges_gradient(
      scale = 3,
      rel_min_height = 0.01,
      colour = "white",
      fill = "grey30"
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 16, colour = "black"),
      axis.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 10),
      legend.position = "none"
    )

  p$plot_env <- rlang::new_environment()
  p2$plot_env <- rlang::new_environment()

  sce@metadata$celltype_annotations$
    metric_plots[[metric_var]][[celltype_var]] <- list()
  sce@metadata$celltype_annotations$
    metric_plots[[metric_var]][[celltype_var]]$metric_data <- dt
  sce@metadata$celltype_annotations$
    metric_plots[[metric_var]][[celltype_var]]$metric_plot <- p
  sce@metadata$celltype_annotations$
    metric_plots[[metric_var]][[celltype_var]]$ridge_data <- dt2
  sce@metadata$celltype_annotations$
    metric_plots[[metric_var]][[celltype_var]]$ridge_plot <- p2

  return(sce)
}

