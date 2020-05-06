annotate_celltype_plots <- function(sce,
                                    cluster_var = "clusters",
                                    celltype_var = "cluster_celltype",
                                    unique_id_var = "manifest",
                                    facet_vars = c("manifest", "group", "sex"),
                                    input_reduced_dim = "UMAP",
                                    metric_vars = c("pc_mito", "pc_ribo", "total_counts", "total_features_by_counts"),
                                    ...) {

  sce@metadata$celltype_annotations <- list()
  sce@metadata$celltype_annotations$reddim_plots <- list()
  sce@metadata$celltype_annotations$prop_plots <- list()
  sce@metadata$celltype_annotations$metric_plots <- list()

  sce@metadata$celltype_annotations$params <- c(as.list(environment()), list(...))
  sce@metadata$celltype_annotations$params$sce <- NULL

  # input validation checks
  var_present <- purrr::map_lgl(
    c(cluster_var, celltype_var, unique_id_var, facet_vars, metric_vars),
    ~ . %in% colnames(SummarizedExperiment::colData(sce))
  )
  assertthat::assert_that(
    all(var_present),
    msg = "Specified variables missing from colData."
  )

  var_classes <- purrr::map_lgl(
    c(cluster_var, celltype_var, unique_id_var, facet_vars),
    ~ class(SummarizedExperiment::colData(sce)[[.]]) %in% c("character", "factor")
  )
  assertthat::assert_that(
    all(var_classes),
    msg = "Specify categorical variables only."
  )

  var_cat_classes <- purrr::map_lgl(
    metric_vars,
    ~ is.numeric(SummarizedExperiment::colData(sce)[[.]])
  )
  assertthat::assert_that(
    all(var_cat_classes),
    msg = "Specify numeric metric variables only."
  )

  # generate the plots
  # cluster numbers
  sce@metadata$celltype_annotations$reddim_plots$cluster_var <- do.call(
    plot_reduced_dim,
    list(sce = sce,
         feature_dim = cluster_var,
         reduced_dim = input_reduced_dim,
         highlight_feature = NA,
         label_clusters = TRUE,
         ...)
  )

  # reddim cluster aliases
  sce@metadata$celltype_annotations$reddim_plots$celltype_var <- do.call(
    plot_reduced_dim,
    list(sce = sce,
         feature_dim = celltype_var,
         reduced_dim = input_reduced_dim,
         highlight_feature = NA,
         label_clusters = TRUE,
         ...)
  )

  # reddim custom coldata variables
  for (facet_var in facet_vars) {
    sce@metadata$celltype_annotations$reddim_plots[[facet_var]] <- do.call(
      plot_reduced_dim,
      list(sce = sce,
           feature_dim = facet_var,
           reduced_dim = input_reduced_dim,
           highlight_feature = NA,
           label_clusters = FALSE,
           ...)
    )
  }

  # proportion of cell types (rel/abs) by groups
  for (var in c(cluster_var, celltype_var)) {
    for (group_by_var in unique(c(facet_vars, unique_id_var, "All"))) {
      sce <- do.call(
        .append_celltype_prop_plots_sce,
        list(sce = sce,
             cluster_var = cluster_var,
             celltype_var = var,
             unique_id_var = unique_id_var,
             facet_var = facet_var,
             group_by_var = group_by_var)
      )
      sce@metadata$celltype_annotations$prop_plots[[group_by_var]][[var]] <-
        lapply(sce@metadata$celltype_annotations$prop_plots[[group_by_var]][[var]], .clean_ggplot_plot_env)
    }
  }


  for (var in c(cluster_var, celltype_var)) {
    for (metric_var in metric_vars) {
      sce <- do.call(
        .append_celltype_metric_plots_sce,
        list(sce = sce,
             celltype_var = var,
             unique_id_var = unique_id_var,
             facet_var = facet_var,
             group_by_var = group_by_var,
             metric_var = metric_var)
      )
      sce@metadata$celltype_annotations$metric_plots[[metric_var]][[var]] <-
        lapply(sce@metadata$celltype_annotations$metric_plots[[metric_var]][[var]], .clean_ggplot_plot_env)
    }
  }

  if (is.null(sce@metadata$scflow_steps) sce@metadata$scflow_steps <- list())
  sce@metadata$scflow_steps$celltype_report_plots_annotated <- TRUE

  return(sce)

}

.append_celltype_prop_plots_sce <- function(...) {

  list2env(list(...), environment())

  if(group_by_var == "All"){ sce$All <- "All" }

  dt <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::group_by(.data[[group_by_var]], .data[[celltype_var]]) %>%
    dplyr::tally() %>%
    dplyr::ungroup(.data[[group_by_var]], .data[[celltype_var]]) %>%
    dplyr::group_by(.data[[group_by_var]]) %>%
    dplyr::mutate(total_cells = sum(n), pc = n/total_cells)

  pal_values <- .get_d_palette("ggsci::nrc_npg", length(unique(dt[[celltype_var]])))
  scale_colours <- ggplot2::scale_fill_manual(values = pal_values)

  # absolute counts (n)
  p <- ggplot2::ggplot(dt,
                       ggplot2::aes(
                         x = .data[[group_by_var]], y = n,
                         fill = .data[[celltype_var]])) +
    ggplot2::geom_col() +
    scale_colours +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 16),
                   axis.title = ggplot2::element_text(size = 18),
                   legend.text = ggplot2::element_text(size = 10)
    )
  #p$plot_env$p <- NULL

  # relative counts (pc)
  p2 <- ggplot2::ggplot(dt,
                        ggplot2::aes(
                          x = .data[[group_by_var]], y = round(pc * 100, 2),
                          fill = .data[[celltype_var]])) +
    ggplot2::geom_col() +
    scale_colours +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::ylab("%") +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 16),
                   axis.title = ggplot2::element_text(size = 18),
                   legend.text = ggplot2::element_text(size = 10)
    )
  #p2$plot_env$p <- NULL
  #p2$plot_env$p2 <- NULL

  sce@metadata$celltype_annotations$prop_plots[[group_by_var]][[celltype_var]] <- list()
  sce@metadata$celltype_annotations$prop_plots[[group_by_var]][[celltype_var]]$prop_data <- dt
  sce@metadata$celltype_annotations$prop_plots[[group_by_var]][[celltype_var]]$absolute_cell_numbers <- p
  sce@metadata$celltype_annotations$prop_plots[[group_by_var]][[celltype_var]]$relative_cell_numbers <- p2

  if(group_by_var == "All"){ sce$All <- NULL }

  return(sce)

}

.clean_ggplot_plot_env <- function(p,
                                   drop_classes = c("SingleCellExperiment",
                                                    "ggplot")) {
  if ("ggplot" %in% class(p)) {
    env_classes <- lapply(p$plot_env, class)
    drop_idx <- purrr::map_lgl(env_classes, ~ any(drop_classes %in% .))
    drop_names <- names(drop_idx)
    for (drop_name in drop_names) {
      p$plot_env[[eval(drop_name)]] <- NULL
    }
    p$plot_env$... <- NULL
  }
  return(p)
}

.append_celltype_metric_plots_sce <- function(...) {

  # metric_var is pc_mito, total_counts, etc.
  list2env(list(...), environment())

  dt <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::group_by(.data[[celltype_var]]) %>%
    dplyr::summarize(
      mean = mean(.data[[metric_var]]),
      sd = sd(.data[[metric_var]]),
      se = sd(.data[[metric_var]]) / sqrt(dplyr::n()),
      median = median(.data[[metric_var]])
    )

  dt[[celltype_var]] <- forcats::fct_reorder(dt[[celltype_var]], -dt$mean)
  dt$metric_var <- metric_var # for export context

  p <- ggplot2::ggplot(dt,
                       ggplot2::aes(
                         x = .data[[celltype_var]], y = mean)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se),
                           width = .2,
                           position = ggplot2::position_dodge(.9)
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::ylab(paste(metric_var)) +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 16),
                   axis.title = ggplot2::element_text(size = 18),
                   legend.text = ggplot2::element_text(size = 10)
    )

  sce@metadata$celltype_annotations$metric_plots[[metric_var]][[celltype_var]] <- list()
  sce@metadata$celltype_annotations$metric_plots[[metric_var]][[celltype_var]]$metric_data <- dt
  sce@metadata$celltype_annotations$metric_plots[[metric_var]][[celltype_var]]$metric_plot <- p

  return(sce)

}


.get_d_palette <- function(pal_name = "ggsci::nrc_npg", n_colours = NULL) {
  palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")
  getPalette <- grDevices::colorRampPalette(palette_choice)
  if (n_colours <= length(palette_choice)) {
    pal_values <- palette_choice[1:n_colours]
  } else {
    pal_values <- getPalette(n_colours)
  }
  return(pal_values)
}

library(magrittr)
sce$age <- as.factor(sce$age)
sce$Cluster <- as.factor(sce$Cluster)
x <- annotate_celltype_plots(sce, cluster_var = "Cluster", unique_id_var = "individual", facet_var = c("age", "sex", "group"), metric_vars = c("pc_mito", "total_counts", "total_features_by_counts"))


