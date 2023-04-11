################################################################################
#' Plot Gene Expression Violin Plots Stratified by Group
#'
#' Plots the cells in reduced dimensionality space with log10 expression values
#' for a specified gene.
#'
#' Generates a 2d plot of cells with gene expression
#'
#' @param sce a SingleCellExperiment object
#' @param group_var The colData variable for x-axis groups
#' @param subset_var The colData variable to subset on
#' @param subset_group The specific subset_var group to subset
#' @param gene The gene of interest
#' @param var_order Optional re-ordering of subset_group factor levels
#' @param palette Optional custom palette
#' @param size point size
#' @param alpha point alpha
#' @param label_angle The angle of x-axis labels (e.g. 0, 45)
#'
#' @return p a ggplot object
#'
#' @family plotting functions
#' @import ggplot2
#' @importFrom dplyr select group_by summarize
#' @importFrom cli cli_alert cli_alert_warning
#' @importFrom assertthat assert_that
#' @importFrom paletteer paletteer_d
#' @importFrom Matrix t
#' @importFrom SingleCellExperiment reducedDim reducedDimNames counts
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom stats median
#'
#' @export
plot_violin <- function(sce,
                        group_var = "group",
                        subset_var = "cluster_celltype",
                        subset_group = "Oligo",
                        gene = "PLP1",
                        var_order = NULL,
                        palette = NULL,
                        alpha = .05,
                        size = .01,
                        label_angle = 0) {

  assertthat::assert_that(gene %in% SummarizedExperiment::rowData(sce)$gene)

  sce <- sce[SummarizedExperiment::rowData(sce)$gene == gene,]

  if (!is.null(subset_var)) {
    assertthat::assert_that(
      subset_var %in% names(SummarizedExperiment::colData(sce)),
      msg = "The subset_var is missing from colData(sce)")
    assertthat::assert_that(subset_group %in% sce[[subset_var]])
    cli::cli_alert("Subsetting {.var {subset_var}} == {.val {subset_group}}")
    sce <- sce[, sce[[subset_var]] == subset_group]
  }

  unique_groups <- unique(sce[[group_var]])
  n_groups <- length(unique_groups)
  assertthat::assert_that(n_groups <= 50, msg = "Too many groups!")

  if (!is.null(palette)) {
    assertthat::assert_that(
      length(palette) >= n_groups,
      msg = sprintf("%s colour(s) specified for %s groups.",
                    length(palette),
                    n_groups))
  }

  cli::cli_alert(c("Plotting gene {.strong {gene}} ",
                   "across {.val {n_groups}} groups: ",
                   "{.var {paste0(unique_groups, collapse = \", \")}}"))

  size_factors <- sce$total_counts/mean(sce$total_counts)
  expr <- as.numeric(
    log2(Matrix::t(
      Matrix::t(SingleCellExperiment::counts(sce))/size_factors) + 1))

  df <- data.frame(group = sce[[group_var]], expr = expr)

  if (class(df$group) != "factor") {
    cli::cli_alert_warning(
      "group_var {.var {group_var}} is not a factor! Coercing to factor..")
    df$group <- as.factor(df$group)
    }

  if (is.null(palette)) {
    if(n_groups <= 10) palette <- paletteer::paletteer_d("ggsci::default_aaas")
    if(n_groups > 10) palette <- paletteer::paletteer_d("ggsci::default_igv")
  }

  if (!is.null(var_order)) {
    assertthat::assert_that(all(var_order %in% unique_groups))
    df$group <- factor(df$group, levels = var_order)
  }

  # Violin plot + boxplot with median +/- stdev, point plot
  p <- ggplot(df, aes(x = group, y = expr)) +
    geom_violin(aes(fill = group), trim = TRUE) +
    scale_fill_manual(values = palette) +
    geom_jitter(size = size, width = .07, alpha = alpha) +
    stat_summary(fun.y = stats::median,
                 fun.ymin = function(x) max(0, mean(x) - sd(x)),
                 fun.ymax = function(x) mean(x) + sd(x),
                 geom = "crossbar", width = 0.06, fill = "white") +
    ggtitle(gene) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
          axis.text.x = element_text(size = 16, angle = label_angle))

  p$plot_env$sce <- NULL
  return(p)

}
