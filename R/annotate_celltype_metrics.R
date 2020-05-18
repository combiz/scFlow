################################################################################
#' Annotate a SingleCellExperiment With Cell-type Metrics
#'
#' @param sce a SingleCellExperiment object
#' @param cluster_var the colData variable with the cluster numbers
#' @param celltype_var the colData variable with the celltype alias
#' @param unique_id_var the colData variable with the unique sample name
#' @param facet_vars the colData variable(s) for grouped analyses
#' @param input_reduced_dim the reducedDim slot used for clustering
#' @param metric_vars the numeric colData variable(s) for metric comparisons
#'
#' @return sce a annotated SingleCellExperiment object
#'
#' @family annotation functions
#'
#' @importFrom purrr map_lgl
#' @importFrom assertthat assert_that
#' @importFrom cli cli_alert_danger cli_alert_success rule cli_text
#' @importFrom SummarizedExperiment colData
#' @export
annotate_celltype_metrics <- function(sce,
                                    cluster_var = "clusters",
                                    celltype_var = "cluster_celltype",
                                    unique_id_var = "manifest",
                                    facet_vars = c("manifest", "group", "sex"),
                                    input_reduced_dim = "UMAP",
                                    metric_vars = c("pc_mito",
                                                    "pc_ribo",
                                                    "total_counts",
                                                    "total_features_by_counts"),
                                    ...) {
  sce@metadata$celltype_annotations <- list()
  sce@metadata$celltype_annotations$reddim_plots <- list()
  sce@metadata$celltype_annotations$prop_plots <- list()
  sce@metadata$celltype_annotations$metric_plots <- list()

  sce@metadata$celltype_annotations$params <- c(
    as.list(environment()),
    list(...)
    )
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
    ~ class(
      SummarizedExperiment::colData(sce)[[.]]) %in% c("character", "factor")
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
    list(
      sce = sce,
      feature_dim = cluster_var,
      reduced_dim = input_reduced_dim,
      highlight_feature = NA,
      label_clusters = TRUE,
      ...
    )
  )

  # reddim cluster aliases
  sce@metadata$celltype_annotations$reddim_plots$celltype_var <- do.call(
    plot_reduced_dim,
    list(
      sce = sce,
      feature_dim = celltype_var,
      reduced_dim = input_reduced_dim,
      highlight_feature = NA,
      label_clusters = TRUE,
      ...
    )
  )

  # reddim custom coldata variables
  for (facet_var in union(facet_vars, unique_id_var)) {
    sce@metadata$celltype_annotations$reddim_plots[[facet_var]] <- do.call(
      plot_reduced_dim,
      list(
        sce = sce,
        feature_dim = facet_var,
        reduced_dim = input_reduced_dim,
        highlight_feature = NA,
        label_clusters = FALSE,
        ...
      )
    )
  }

  # proportion of cell types (rel/abs) by groups
  for (var in c(cluster_var, celltype_var)) {
    for (group_by_var in unique(c(facet_vars, unique_id_var, "all"))) {
      sce <- do.call(
        .append_cell_prop_plots_sce,
        list(
          sce = sce,
          cluster_var = cluster_var,
          celltype_var = var,
          unique_id_var = unique_id_var,
          facet_var = facet_var,
          group_by_var = group_by_var
        )
      )
      #sce@metadata$celltype_annotations$prop_plots[[group_by_var]][[var]] <-
      #  lapply(
      #    sce@metadata$celltype_annotations$prop_plots[[group_by_var]][[var]],
      #    .clean_ggplot_plot_env
      #    )
    }
  }


  for (var in c(cluster_var, celltype_var)) {
    for (metric_var in metric_vars) {
      sce <- do.call(
        .append_cell_metric_plots_sce,
        list(
          sce = sce,
          celltype_var = var,
          unique_id_var = unique_id_var,
          facet_var = facet_var,
          group_by_var = group_by_var,
          metric_var = metric_var
        )
      )
      sce@metadata$celltype_annotations$metric_plots[[metric_var]][[var]] <-
        lapply(
          sce@metadata$celltype_annotations$metric_plots[[metric_var]][[var]],
          .clean_ggplot_plot_env
          )
    }
  }

  if (is.null(sce@metadata$scflow_steps)) sce@metadata$scflow_steps <- list()
  sce@metadata$scflow_steps$celltype_plots_annotated <- TRUE

  return(sce)
}

################################################################################
#' Attach proportional celltype plots to sce
#'
#' @family helper
#'
#' @importFrom dplyr group_by tally ungroup group_by mutate
#' @importFrom ggplot2 scale_fill_manual ggplot aes geom_col scale_y_continuous
#' @importFrom ggplot2 coord_flip theme_bw theme xlab ylab
#' @importFrom ggplot2 element_blank element_text
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment colData
#'
#' @keywords internal
.append_cell_prop_plots_sce <- function(...) {
  list2env(list(...), environment())

  if (group_by_var == "all") {
    sce$all <- "all"
  }

  dt <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::group_by(.data[[group_by_var]], .data[[celltype_var]]) %>%
    dplyr::tally() %>%
    dplyr::ungroup(.data[[group_by_var]], .data[[celltype_var]]) %>%
    dplyr::group_by(.data[[group_by_var]]) %>%
    dplyr::mutate(total_cells = sum(n), pc = n / total_cells)

  pal_values <- .get_d_palette(
    "ggsci::nrc_npg", length(unique(dt[[celltype_var]]))
    )
  scale_colours <- ggplot2::scale_fill_manual(values = pal_values)

  # absolute counts (n)
  p <- ggplot2::ggplot(
    dt,
    ggplot2::aes(
      x = .data[[group_by_var]], y = n,
      fill = .data[[celltype_var]]
    )
  ) +
    ggplot2::geom_col() +
    scale_colours +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 16),
      axis.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_blank()
    )

  # relative counts (pc)
  p2 <- ggplot2::ggplot(
    dt,
    ggplot2::aes(
      x = .data[[group_by_var]], y = round(pc * 100, 2),
      fill = .data[[celltype_var]]
    )
  ) +
    ggplot2::geom_col() +
    scale_colours +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::ylab("%") +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 16),
      axis.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_blank()
    )

  p$plot_env <- rlang::new_environment()
  p2$plot_env <- rlang::new_environment()
  sce@metadata$celltype_annotations$
    prop_plots[[group_by_var]][[celltype_var]] <- list()
  sce@metadata$celltype_annotations$
    prop_plots[[group_by_var]][[celltype_var]]$prop_data <- dt
  sce@metadata$celltype_annotations$
    prop_plots[[group_by_var]][[celltype_var]]$absolute_cell_numbers <- p
  sce@metadata$celltype_annotations$
    prop_plots[[group_by_var]][[celltype_var]]$relative_cell_numbers <- p2

  if (group_by_var == "all") {
    sce$all <- NULL
  }

  return(sce)
}


################################################################################
#' Attach celltype associated metric plots to sce
#'
#' @family helper
#'
#' @importFrom dplyr group_by tally ungroup group_by mutate summarize
#' @importFrom forcats fct_reorder
#' @importFrom ggplot2 scale_fill_manual ggplot aes geom_col scale_y_continuous
#' @importFrom ggplot2 coord_flip theme_bw theme xlab ylab geom_bar
#' @importFrom ggplot2 geom_errorbar element_blank element_text
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment colData
#'
#' @keywords internal
.append_cell_metric_plots_sce <- function(...) {

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

  p <- ggplot2::ggplot(
    dt,
    ggplot2::aes(
      x = .data[[celltype_var]], y = mean
    )
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se),
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
                        ggplot2::aes(x = .data[[metric_var]],
                                     y = .data[[celltype_var]])) +
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

################################################################################
#' Retrieve a palette of n_colours discrete colours from paletteer
#'
#' Useful to remove large objects before writing to disk with qs or rds
#'
#' @family helper
#'
#' @importFrom paletteer paletteer_d
#' @importFrom grDevices colorRampPalette
#'
#' @keywords internal
.get_d_palette <- function(pal_name = "ggsci::nrc_npg", n_colours = NULL) {
  palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")
  get_palette <- grDevices::colorRampPalette(palette_choice)
  if (n_colours <= length(palette_choice)) {
    pal_values <- palette_choice[1:n_colours]
  } else {
    pal_values <- get_palette(n_colours)
  }
  return(pal_values)
}


################################################################################
#' Remove objects from a ggplot plot_env
#'
#' Useful to remove large objects before writing to disk with qs or rds
#'
#' @family helper
#'
#' @keywords internal
.clean_ggplot_plot_env <- function(p,
                                   drop_classes = c(
                                     "SingleCellExperiment",
                                     "ggplot"
                                   )) {
  if ("ggplot" %in% class(p)) {
    env_classes <- lapply(p$plot_env, class)
    drop_idx <- purrr::map_lgl(env_classes, ~ any(drop_classes %in% .))
    drop_names <- names(drop_idx[drop_idx])
    for (drop_name in drop_names) {
      p$plot_env[[eval(drop_name)]] <- NULL
    }
    p$plot_env$... <- NULL
  }
  return(p)
}

################################################################################
#' scale_y_continuous(labels=.fancy_scientific)
#'
#' @family helper
#'
#' @keywords internal
.fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # prevent 0 x xx
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2)
  l <- gsub("e\\+","e",l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}
