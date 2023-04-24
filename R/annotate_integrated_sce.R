################################################################################
#' Annotate integrated, reduced dimension,
#' and clustered SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object
#' @param categorical_covariates list of categorical variables
#' @param input_reduced_dim which reduced dim to use for integration and
#' clustering. Either "tSNE" or "UMAP" (default "tSNE").
#'
#' @return sce annotated SingleCellExperiment object
#'
#' @family annotation functions
#'
#' @import ggplot2
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom tools file_path_sans_ext
#' @importFrom formattable formattable icontext
#' @importFrom UpSetR upset fromList
#' @importFrom kBET kBET
#' @export

annotate_integrated_sce <- function(sce,
                                    categorical_covariates = list(),
                                    input_reduced_dim = "tSNE") {
  if (!class(sce) == "SingleCellExperiment") {
    stop("expecting singlecellexperiment")
  }
  cat(
    cli::rule(
      "Annotating Selected Variable Genes",
      line = 2
    ),
    "\r\n"
  )

  cli::cli_text("Generating Upset chart for selected variable genes...")

  sce <- .generate_upset_plot(sce)

  cat(
    cli::rule(
      "Annotating Batch Effect Correction by LIGER",
      line = 2
    ),
    "\r\n"
  )

  sce@metadata$dataset_integration$annotation$input_reduced_dim <-
    input_reduced_dim

  pca_reducedDim_plots <- list()
  liger_reducedDim_plots <- list()
  pca_kbet_plots <- list()
  liger_kbet_plots <- list()
  cli::cli_text("Generating tSNE/UMAP and kBET plots for each covariate...")

  for (variable in categorical_covariates) {
    cat(paste("* covariate:", variable, sep = " "), "\n")
    plot_pca <- plot_reduced_dim(sce,
      feature_dim = variable,
      reduced_dim = sprintf(
        "%s_PCA",
        input_reduced_dim
      )
    )

    plot_pca <- .grobify_ggplot(plot_pca)
    pca_reducedDim_plots[[variable]] <- plot_pca

    plot_liger <- plot_reduced_dim(sce,
      feature_dim = variable,
      reduced_dim = sprintf(
        "%s_Liger",
        input_reduced_dim
      )
    )

    plot_liger <- .grobify_ggplot(plot_liger)
    liger_reducedDim_plots[[variable]] <- plot_liger

    # sub-sampling data
    batch <- SummarizedExperiment::colData(sce)[[variable]]
    subset_id <- sample.int(
      n = length(batch), size = floor(0.1 * length(batch)),
      replace = FALSE
    )

    # calculate kBET

    data <- SingleCellExperiment::reducedDim(sce, "PCA")

    pca_kbet_plots[[variable]] <- .generate_kbet_plot(
      data = data,
      batch = batch,
      subset_id = subset_id,
      integration_method = "PCA"
    )


    data <- SingleCellExperiment::reducedDim(sce, "Liger")

    liger_kbet_plots[[variable]] <- .generate_kbet_plot(
      data = data,
      batch = batch,
      subset_id = subset_id,
      integration_method = "LIGER"
    )
  }

  sce@metadata$dataset_integration$
    batch_correction_plots$pca_reducedDim_plots <- pca_reducedDim_plots
  sce@metadata$dataset_integration$
    batch_correction_plots$liger_reducedDim_plots <- liger_reducedDim_plots
  sce@metadata$dataset_integration$
    batch_correction_plots$pca_kbet_plots <- pca_kbet_plots
  sce@metadata$dataset_integration$
    batch_correction_plots$liger_kbet_plots <- liger_kbet_plots

  cat(
    cli::rule(
      "Annotating Clustering",
      line = 2
    ),
    "\r\n"
  )

  cli::cli_text("Generating tSNE/UMAP plots for PCA and LIGER data...")
  cluster_pca <- plot_reduced_dim(sce,
    feature_dim = "clusters",
    reduced_dim = sprintf(
      "%s_PCA",
      input_reduced_dim
    ),
    label_clusters = TRUE
  )

  cluster_pca <- .grobify_ggplot(cluster_pca)
  sce@metadata$dataset_integration$clustering_plots$cluster_pca <- cluster_pca

  cluster_liger <- plot_reduced_dim(sce,
    feature_dim = "clusters",
    reduced_dim = sprintf(
      "%s_Liger",
      input_reduced_dim
    ),
    label_clusters = TRUE
  )

  cluster_liger <- .grobify_ggplot(cluster_liger)

  sce@metadata$dataset_integration$clustering_plots$cluster_liger <-
    cluster_liger

  return(sce)
}


#' @importFrom UpSetR upset fromList
#' @keywords internal
.generate_upset_plot <- function(sce) {
  upset_sets <- sce@metadata$dataset_integration$var.genes_per_dataset
  my_upset <- UpSetR::upset(
    UpSetR::fromList(upset_sets),
    nsets = length(upset_sets),
    sets.x.label = "Variable genes per dataset",
    text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1)
  )
  sce@metadata$dataset_integration$var.genes_plots$upset <- my_upset
  return(sce)
}

#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_bw scale_y_continuous
#' @importFrom kBET kBET
#' @keywords internal
.generate_kbet_plot <- function(data, batch, subset_id, integration_method) {
  batch_estimate <- kBET::kBET(data[subset_id, ], batch[subset_id],
    plot = FALSE, do.pca = FALSE
  )

  plot.data <- data.frame(
    class = rep(c("observed", "expected"),
      each = length(batch_estimate$stats$kBET.observed)
    ),
    data = c(
      batch_estimate$stats$kBET.observed,
      batch_estimate$stats$kBET.expected
    )
  )
  x_label <- sprintf(
    "P-value = %s",
    formatC(batch_estimate$summary$kBET.signif[1],
      format = "e",
      digits = 2
    )
  )
  kbet_plot <- ggplot2::ggplot(plot.data, ggplot2::aes(class, data)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(x = x_label, y = "Rejection rate") +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::theme(
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 13),
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
    ) +
    ggtitle(sprintf("kBET test results - %s", integration_method))

  kbet_plot <- .grobify_ggplot(kbet_plot)

  return(kbet_plot)
}
