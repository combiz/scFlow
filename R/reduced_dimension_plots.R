################################################################################
#' Plot Cells Reduced Dimensions With Feature Highlighting
#'
#' Generates a 2d plot of cells with feature highlighting
#'
#' @param sce a SingleCellExperiment object
#' @param feature_dim the colData variable of interest
#' @param reduced_dim the reducedDim slot for plotting
#' @param highlight_feature highlights a feature
#' @param label_clusters prints cluster labels
#' @param size point size
#' @param alpha point alpha
#'
#' @return p a ggplot object
#'
#' @family plotting functions
#' @import ggplot2
#' @importFrom dplyr select group_by summarize
#' @importFrom paletteer paletteer_d
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom magrittr %>%
#' @importFrom grDevices colorRampPalette
#' @importFrom SummarizedExperiment rowData colData
#'
#' @export
plot_reduced_dim <- function(sce,
                             feature_dim = "Cluster",
                             reduced_dim = "UMAP",
                             highlight_feature = NA,
                             label_clusters = FALSE,
                             size = 0.1,
                             alpha = 0.2) {

  assertthat::assert_that(
    feature_dim %in% names(SummarizedExperiment::colData(sce)),
    msg = sprintf(
      "The feature_dim %s is not present in colData(sce)", feature_dim)
  )

  assertthat::assert_that(
    reduced_dim %in% SingleCellExperiment::reducedDimNames(sce),
    msg = sprintf(
      "The reduced_dim %s is not present in reducedDim(sce)", reduced_dim)
  )

  # prepare data
  dt <- as.data.frame(
    cbind(colData(sce),
          "dim1" = SingleCellExperiment::reducedDim(sce, reduced_dim)[, 1],
          "dim2" = SingleCellExperiment::reducedDim(sce, reduced_dim)[, 2])
  )
  if (is.na(highlight_feature)) {
    dt$feature_dim <- dt[[feature_dim]]
  } else {
    dt$feature_dim <- (dt[[feature_dim]] == highlight_feature)
  }

  # colour palette
  colourCount <- length(unique(dt$feature_dim))
  palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")
  getPalette <- grDevices::colorRampPalette(palette_choice)
  if (colourCount <= length(palette_choice)) {
    if (colourCount == 2) {
      if(is.na(highlight_feature)){
        pal_values <- c(palette_choice[2], palette_choice[1]) #blue red
      } else {
        pal_values <- c("grey80", "red") #blue red
      }
    } else {
      pal_values <- palette_choice[1:colourCount]
    }
  } else {
    pal_values <- getPalette(colourCount)
  }

  if (is.double(dt$feature_dim) | is.integer(dt$feature_dim)){
    scale_colours <- ggplot2::scale_colour_gradientn(
      limits = c(0, mean(feature_dim)),
      colours=c("grey80", "red"),
      na.value="red"
    )

  } else {
    scale_colours <- ggplot2::scale_colour_manual(values = pal_values)
  }

  p <- ggplot2::ggplot(data = dt) +
    ggplot2::geom_point(
      ggplot2::aes(x = dt$dim1,
                   y = dt$dim2,
                   colour = dt$feature_dim),
      shape = 16, size = size, alpha = alpha) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 6,
                            alpha = 1))) +
    scale_colours +
    ggplot2::ylab(dt$dim1) +
    ggplot2::xlab(dt$dim2) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      line = ggplot2::element_blank(),
      text = ggplot2::element_blank(),
      title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(size = 18, hjust = 0.5)
    )

  # if, clusters labelled individually
  if (label_clusters == TRUE){
    cluster_centroids <- dt %>%
      dplyr::select(dim1, dim2, feature_dim) %>%
      dplyr::group_by(feature_dim) %>%
      dplyr::summarize(mean_dim1 = mean(dim1),
                       mean_dim2 = mean(dim2))

    p <- p +
      ggplot2::geom_text(
        data = cluster_centroids,
        ggplot2::aes(x = mean_dim1,
                     y = mean_dim2,
                     label = feature_dim)
      ) +
      ggplot2::theme(legend.position = "none")
  }

  p$plot_env$sce <- NULL
  return(p)

}

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
                                  size = 0.1,
                                  alpha = 0.2,
                                  palette = c("grey80", "red")) {

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

  # prepare data
  dt <- data.frame(
    "dim1" = SingleCellExperiment::reducedDim(sce, reduced_dim)[, 1],
    "dim2" = SingleCellExperiment::reducedDim(sce, reduced_dim)[, 2],
    "gene" = SingleCellExperiment::counts(sce)[which(
      SingleCellExperiment::rowData(sce)$gene == gene)[1], ]
  )

  dt$gene <- log10(dt$gene) + 1

  med_counts <- median(dt[dt$gene > 0,]$gene)
  med_total_counts <- median(sce$total_counts)

  p <- ggplot(data = dt) +
    geom_point(aes(x = dt$dim1, y = dt$dim2,
                   colour = dt$gene),
               shape = 16, size = size, alpha = alpha) +
    scale_colour_gradientn(#limits = c(0, 1),
                           colours = c(palette[1], palette[2]),
                           na.value = palette[1])+
    ggtitle(gene)+
    ylab(dt$dim1) +
    xlab(dt$dim2) +
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

  p$plot_env$sce <- NULL
  return(p)

}
