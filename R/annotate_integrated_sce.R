################################################################################
#' Annotate integrated, reduced dimension,
#' and clustered SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object
#' @param categorical_covariates list of categorical variables
#'
#' @return sce annotated SingleCellExperiment object
#'
#' @family integration, dimension reduction, and clustering
#'
#' @import cli Matrix dplyr SingleCellExperiment purrr
#' @import ggplot2
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom tools file_path_sans_ext
#' @importFrom formattable formattable icontext
#' @export

annotate_integrated_sce <- function(sce,
                                    categorical_covariates = list()) {
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
  cli::cli_text("Generating Venn diagram for selected variable genes...")
  if (length(sce_merged@metadata$var.genes_per_dataset) < 11) {
    venn_sets <- sce@metadata$var.genes_per_dataset
    my_nv <- nVennR::plotVenn(venn_sets)
    my_nv <- nVennR::plotVenn(nVennObj = my_nv)
    sce@metadata$dataset_integration$var.genes_plots$venn <- my_nv
  } else {
    print("The number of datasets is too high for a Venn diagram.")
  }

  cli::cli_text("Generating Upset chart for selected variable genes...")

  upset_sets <- sce@metadata$var.genes_per_dataset
  my_upset <- UpSetR::upset(fromList(upset_sets),
    nsets = length(upset_sets),
    sets.x.label = "Variable genes per dataset",
    text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1)
  )
  sce@metadata$dataset_integration$var.genes_plots$upset <- my_upset

  cat(
    cli::rule(
      "Annotating Batch Effect Correction by LIGER",
      line = 2
    ),
    "\r\n"
  )

  pca_tsne_plots <- list()
  liger_tsne_plots <- list()
  pca_kbet_plots <- list()
  liger_kbet_plots <- list()
  cli::cli_text("Generating tSNE and kBET plots for each covariate...")

  for (variable in categorical_covariates) {
    cat(paste("• covariate:", variable, sep = " "), "\n")
    plot_pca <- plot_reduced_dim(sce,
      feature_dim = variable,
      reduced_dim = "tSNE_PCA", size = 1, alpha = 0.3
    )
    tsne_pca <- plot_pca + ggtitle("tSNE using PCA data") +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.8))
    pca_tsne_plots[[variable]] <- tsne_pca

    plot_liger <- plot_reduced_dim(sce,
      feature_dim = variable,
      reduced_dim = "tSNE_Liger", size = 0.75, alpha = 0.3
    )
    tsne_liger <- plot_liger + ggtitle("tSNE using LIGER data") +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.8))
    liger_tsne_plots[[variable]] <- tsne_liger

    data <- SingleCellExperiment::reducedDim(sce, "PCA")
    batch <- SummarizedExperiment::colData(sce)[[variable]]
    subset_id <- sample.int(
      n = length(batch), size = floor(0.1 * length(batch)),
      replace = FALSE
    )
    batch_estimate_pca <- kBET::kBET(data[subset_id, ], batch[subset_id],
      plot = FALSE
    )
    plot.data <- data.frame(
      class = rep(c("observed", "expected"),
        each = length(batch_estimate_pca$stats$kBET.observed)
      ),
      data = c(
        batch_estimate_pca$stats$kBET.observed,
        batch_estimate_pca$stats$kBET.expected
      )
    )
    x_label <- sprintf(
      "P-value = %s",
      formatC(batch_estimate_pca$summary$kBET.signif[1], format = "e",
              digits = 2)
    )
    kbet_pca <- ggplot(plot.data, aes(class, data)) +
      geom_boxplot() +
      labs(x = x_label, y = "Rejection rate") +
      theme_bw() +
      scale_y_continuous(limits = c(0, 1)) +
      theme(
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(), panel.background = element_blank()
      ) +
      ggtitle("kBET test results - PCA") +
      theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
    pca_kbet_plots[[variable]] <- kbet_pca

    data <- SingleCellExperiment::reducedDim(sce, "Liger")
    batch <- SummarizedExperiment::colData(sce)[[variable]]
    subset_id <- sample.int(
      n = length(batch), size = floor(0.1 * length(batch)),
      replace = FALSE
    )
    batch_estimate_liger <- kBET::kBET(data[subset_id, ], batch[subset_id],
      plot = FALSE
    )
    plot.data <- data.frame(
      class = rep(c("observed", "expected"),
        each = length(batch_estimate_liger$stats$kBET.observed)
      ),
      data = c(
        batch_estimate_liger$stats$kBET.observed,
        batch_estimate_liger$stats$kBET.expected
      )
    )
    x_label <- sprintf(
      "P-value = %s",
      formatC(batch_estimate_liger$summary$kBET.signif[1],
        format = "e",
        digits = 2
      )
    )
    kbet_liger <- ggplot(plot.data, aes(class, data)) +
      geom_boxplot() +
      labs(x = x_label, y = "Rejection rate") +
      theme_bw() +
      scale_y_continuous(limits = c(0, 1)) +
      theme(
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(), panel.background = element_blank()
      ) +
      ggtitle("kBET test results - LIGER") +
      theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
    liger_kbet_plots[[variable]] <- kbet_liger
  }

  sce@metadata$dataset_integration$batch_correction_plots$pca_tsne_plots <-
    pca_tsne_plots
  sce@metadata$dataset_integration$batch_correction_plots$liger_tsne_plots <-
    liger_tsne_plots
  sce@metadata$dataset_integration$batch_correction_plots$pca_kbet_plots <-
    pca_kbet_plots
  sce@metadata$dataset_integration$batch_correction_plots$liger_kbet_plots <-
    liger_kbet_plots

  cat(
    cli::rule(
      "Annotating Clustering",
      line = 2
    ),
    "\r\n"
  )

  cli::cli_text("Generating UMAP plots for PCA and LIGER data...")
  plot_pca <- plot_reduced_dim(sce,
    feature_dim = "clusters",
    reduced_dim = "UMAP_PCA", size = 0.75, alpha = 0.3
  )
  umap_pca <- plot_pca + ggtitle("UMAP using PCA data (colored by cluster)") +
    theme(plot.title = element_text(size = 15, face = "bold"))

  sce@metadata$dataset_integration$clustering_plots$umap_pca <- umap_pca

  plot_liger <- plot_reduced_dim(sce,
    feature_dim = "clusters",
    reduced_dim = "UMAP_Liger", size = 0.75, alpha = 0.3
  )
  umap_liger <- plot_liger +
    ggtitle("UMAP using LIGER data (colored by cluster)") +
    theme(plot.title = element_text(size = 15, face = "bold"))

  sce@metadata$dataset_integration$clustering_plots$umap_liger <- umap_liger
  return(sce)
}