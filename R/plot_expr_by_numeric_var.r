################################################################################
#' Plot Gene Expression for Samples against Numerical Variable
#'
#' Equivalent of plot_violin with a numeric x-axis.
#'
#' Generates a scatter plot of log2(counts+1) expression against numeric var.
#' An additional plot of proportion of cells with non-zero counts per sample.
#' Fits a linear model with confidence intervals for each plot.
#'
#'
#' @param sce a SingleCellExperiment object
#' @param numeric_var The colData variable for x-axis groups. Default is p_tau
#' @param subset_var The colData variable to subset on
#' @param subset_group The specific subset_var group to subset
#' @param gene The gene of interest
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
#' @importFrom patchwork plot_layout
#' @importFrom scales pretty_breaks percent_format
#' @importFrom stats lm coef predict
#'
#' @export
plot_expr_by_numeric_var <- function(sce,
                        numeric_var = "p_tau",
                        subset_var = "cluster_celltype",
                        subset_group = "Oligo",
                        gene = "PLP1",
                        #palette = NULL,
                        alpha = .05,
                        size = .01,
                        label_angle = 0) {

  assertthat::assert_that(gene %in% SummarizedExperiment::rowData(sce)$gene)
  assertthat::assert_that(class(sce[[numeric_var]]) == "numeric")

  sce <- sce[SummarizedExperiment::rowData(sce)$gene == gene,]

  if (!is.null(subset_var)) {
    assertthat::assert_that(
      subset_var %in% names(SummarizedExperiment::colData(sce)),
      msg = "The subset_var is missing from colData(sce)")
    assertthat::assert_that(subset_group %in% sce[[subset_var]])
    cli::cli_alert("Subsetting {.var {subset_var}} == {.val {subset_group}}")
    sce <- sce[, sce[[subset_var]] == subset_group]
  }


  cli::cli_alert(c("Plotting gene {.strong {gene}} "))

  size_factors <- sce$total_counts/mean(sce$total_counts)
  expr <- as.numeric(
    log2(Matrix::t(
      Matrix::t(SingleCellExperiment::counts(sce))/size_factors) + 1))

  df <- data.frame(numeric_var = sce[[numeric_var]], expr = expr)

  #if (is.null(palette)) {
  #  if(n_groups <= 10) palette <- paletteer::paletteer_d("ggsci::default_aaas")
  #  if(n_groups > 10) palette <- paletteer::paletteer_d("ggsci::default_igv")
  #}

  #pred <- predict(lm(expr ~ numeric_var, df),  ## ALL cells
  pred <- predict(lm(expr ~ numeric_var, df[df$expr > 0, ]), # only expressive cells
                  se.fit = TRUE, interval = "confidence")
  limits <- as.data.frame(pred$fit)
  limits$numeric_var <- df[df$expr > 0,]$numeric_var

  # Main scatter plot
  p <- ggplot(df, aes(x = numeric_var, y = expr)) +
    #geom_point() +
    geom_jitter(size = size, width = .03, alpha = alpha) +
    stat_summary(data = df[df$expr > 0,], fun = median,
                 #fun.ymin = function(x) max(0, mean(x) - sd(x)),
                 #fun.ymax = function(x) mean(x) + sd(x),
                 geom = "point", fill = "white", size = 2) +
    #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_smooth(data = df[df$expr > 0,], method = "lm", colour = "#366092", se = TRUE, fill = "#366092", alpha = 0.3) +
    geom_line(data = limits, aes(x = numeric_var, y = lwr),
              linetype = 2, colour = "#366092") +
    geom_line(data = limits, aes(x = numeric_var, y = upr),
              linetype = 2, colour = "#366092") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    ggplot2::ylab(bquote(Log[2]*" (counts + 1)")) +
    ggplot2::xlab(numeric_var) +
    ggtitle(gene) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
          axis.text = element_text(size = 16, angle = label_angle),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(5.5, 5.5, 0, 5.5))

  # calculate zero counts proportion by x
  dt <- df %>%
    dplyr::group_by(numeric_var) %>%
    dplyr::mutate(is_zero = ifelse(expr == 0, "zero", "non_zero")) %>%
    dplyr::group_by(numeric_var, is_zero) %>%
    dplyr::tally() %>%
    tidyr::pivot_wider(names_from = c("is_zero"), values_from = "n") %>%
    dplyr::mutate(pc_zero = zero/(zero + non_zero),
                  pc_expressive = non_zero/(zero + non_zero))

  pred <- predict(lm(pc_expressive ~ numeric_var, dt),
                  se.fit = TRUE, interval = "confidence")
  limits <- as.data.frame(pred$fit)
  limits$numeric_var <- dt$numeric_var

  # lower plot of zero counts proportion
  zp <- ggplot(dt, aes(x = numeric_var, y = pc_expressive)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", colour = "#366092", se = TRUE, fill = "#366092", alpha = 0.3) +
    geom_line(data = limits, aes(x = numeric_var, y = lwr),
              linetype = 2, colour = "#366092") +
    geom_line(data = limits, aes(x = numeric_var, y = upr),
              linetype = 2, colour = "#366092") +
    xlab(numeric_var) +
    ylab("Non-zero (%)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::pretty_breaks(n = 4)) +
    theme_bw()+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(0, 5.5, 5.5, 5.5)
    )

  p <- p / zp + patchwork::plot_layout(heights = c(2.5, 1))

  p$plot_env$sce <- NULL
  return(p)

}
