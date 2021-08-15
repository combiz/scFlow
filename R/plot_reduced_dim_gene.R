################################################################################
#' Plot Cells Reduced Dimensions With Gene Expression
#'
#' Plots the cells in reduced dimensionality space with log10 expression values
#' for a specified gene.
#'
#' Generates a 2d plot of cells with gene expression
#'
#' @param sce a SingleCellExperiment object
#' @param reduced_dim the reducedDim slot for plotting
#' @param gene the gene of interest
#' @param size point size
#' @param alpha point alpha
#' @param palette off and on colours
#'
#' @return p a ggplot object
#'
#' @family plotting functions
#' @import ggplot2
#' @importFrom dplyr select group_by summarize
#' @importFrom assertthat assert_that
#' @importFrom paletteer paletteer_d
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom magrittr %>%
#' @importFrom grDevices colorRampPalette
#' @importFrom SummarizedExperiment rowData colData
#'
#' @export
plot_reduced_dim_gene <- function(sce,
                                  reduced_dim = "UMAP",
                                  gene = "PLP1",
                                  size = getOption(
                                    "scflow_reddimplot_pointsize",
                                    default = 0.1),
                                  alpha = getOption(
                                    "scflow_reddimplot_alpha",
                                    default = 0.2),
                                  palette = c("grey80", "#440154FF")) {

  assertthat::assert_that(
    reduced_dim %in% SingleCellExperiment::reducedDimNames(sce),
    msg = sprintf(
      "Invalid reducedDim (%s) specified.", reduced_dim)
  )

  assertthat::assert_that(
    all(gene %in% SummarizedExperiment::rowData(sce)$gene),
    msg = sprintf(
      "Gene %s is absent from the matrix.  Try another gene.", gene)
  )

  assertthat::assert_that(
    length(palette) == 2,
    msg = "Palette should contain exactly two colours."
  )

  cli::cli_alert(sprintf(
    "Plotting {.strong %s} with gene expression for {.var %s}",
    reduced_dim, gene)
  )

  sce <- sce[SummarizedExperiment::rowData(sce)$gene == gene,]
  # prepare data
  dt <- data.frame(
    "dim1" = SingleCellExperiment::reducedDim(sce, reduced_dim)[, 1],
    "dim2" = SingleCellExperiment::reducedDim(sce, reduced_dim)[, 2]
  )

  size_factors <- sce$total_counts/mean(sce$total_counts)
  dt$gene <- as.numeric(
    log2(Matrix::t(
      Matrix::t(SingleCellExperiment::counts(sce))/size_factors) + 1))
  #dt$gene <- log10(dt$gene) + 1


  p <- ggplot(data = dt) +
    geom_point(aes(x = dim1, y = dim2,
                   colour = gene),
               shape = 16, size = size, alpha = alpha) +
    scale_colour_gradientn(#limits = c(0, 1),
      colours = c(palette[1], palette[2]),
      na.value = palette[1])+
    ggtitle(gene)+
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      line = element_blank(),
      text = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 18, hjust = 0.5, colour = "black")
    )

  p$plot_env <- rlang::new_environment()
  return(p)

}
