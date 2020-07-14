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
#' @family integration, dimension reduction, and clustering
#'
#' @import cli Matrix dplyr SingleCellExperiment purrr
#' @import ggplot2
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom tools file_path_sans_ext
#' @importFrom formattable formattable icontext
#' @importFrom nVennR plotVenn
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
  cli::cli_text("Generating Venn diagram for selected variable genes...")
  if (length(sce@metadata$dataset_integration$var.genes_per_dataset) < 11) {
    venn_sets <- sce@metadata$dataset_integration$var.genes_per_dataset
    my_nv <- nVennR::plotVenn(venn_sets)
    my_nv <- nVennR::plotVenn(nVennObj = my_nv)
    sce@metadata$dataset_integration$var.genes_plots$venn <- my_nv
  } else {
    print("The number of datasets is too high for a Venn diagram.")
  }

  cli::cli_text("Generating Upset chart for selected variable genes...")

  upset_sets <- sce@metadata$dataset_integration$var.genes_per_dataset
  my_upset <- UpSetR::upset(
    UpSetR::fromList(upset_sets),
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

  sce@metadata$dataset_integration$annotation$input_reduced_dim <-
    input_reduced_dim

  pca_reducedDim_plots <- list()
  liger_reducedDim_plots <- list()
  pca_kbet_plots <- list()
  liger_kbet_plots <- list()
  cli::cli_text("Generating tSNE/UMAP and kBET plots for each covariate...")

  for (variable in categorical_covariates) {
    cat(paste("• covariate:", variable, sep = " "), "\n")
    plot_pca <- plot_reduced_dim(sce,
                                 feature_dim = variable,
                                 reduced_dim = sprintf("%s_PCA",
                                                       input_reduced_dim),
                                 size = 1, alpha = 0.3)
    # reducedDim_pca <- plot_pca +
    # ggtitle(sprintf("%s using PCA data", input_reduced_dim)) +
    # theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.8))
    # pca_reducedDim_plots[[variable]] <- reducedDim_pca
    pca_reducedDim_plots[[variable]] <- plot_pca

    plot_liger <- plot_reduced_dim(sce,
                                   feature_dim = variable,
                                   reduced_dim = sprintf("%s_Liger",
                                                         input_reduced_dim),
                                   size = 0.75, alpha = 0.3)
    # reducedDim_liger <- plot_liger +
    # ggtitle(sprintf("%s using LIGER data", input_reduced_dim)) +
    # theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.8))
    # liger_reducedDim_plots[[variable]] <- reducedDim_liger
    liger_reducedDim_plots[[variable]] <- plot_liger

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

  sce@metadata$dataset_integration$batch_correction_plots$pca_reducedDim_plots <-
    pca_reducedDim_plots
  sce@metadata$dataset_integration$batch_correction_plots$liger_reducedDim_plots <-
    liger_reducedDim_plots
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

  cli::cli_text("Generating tSNE/UMAP plots for PCA and LIGER data...")
  sce <- cluster_sce(sce, k = 50, reduction_method = sprintf("%s_PCA",
                                                             input_reduced_dim))
  pca_clusters <- sce$clusters

  plot_pca <- plot_reduced_dim(sce,
                               feature_dim = "clusters",
                               reduced_dim = sprintf("%s_PCA",
                                                     input_reduced_dim),
                               size = 0.75, alpha = 0.3, label_clusters = TRUE
  )
  # cluster_pca <- plot_pca +
  # ggtitle(sprintf("%s using PCA data (colored by cluster)",
  # input_reduced_dim)) +
  # theme(plot.title = element_text(size = 15, face = "bold"))

  cluster_pca <- plot_pca
  sce@metadata$dataset_integration$clustering_plots$cluster_pca <- cluster_pca

  sce <- cluster_sce(sce, k = 50, reduction_method = sprintf("%s_Liger",
                                                             input_reduced_dim))

  liger_clusters <- sce$clusters

  plot_liger <- plot_reduced_dim(sce,
                                 feature_dim = "clusters",
                                 reduced_dim = sprintf("%s_Liger",
                                                       input_reduced_dim),
                                 size = 0.75, alpha = 0.3,
                                 label_clusters = TRUE
  )
  # cluster_liger <- plot_liger +
  # ggtitle(sprintf("%s using LIGER data (colored by cluster)",
  # input_reduced_dim)) +
  # theme(plot.title = element_text(size = 15, face = "bold"))
  cluster_liger <- plot_liger

  sce@metadata$dataset_integration$clustering_plots$cluster_liger <-
    cluster_liger

  sankey_data <- as.data.frame.matrix(table(pca_clusters, liger_clusters)) %>%
    rownames_to_column %>%
    tidyr::gather(key = "key", value = "value", -rowname) %>%
    dplyr::filter(value > 0)
  colnames(sankey_data) <- c("pca_cluster", "liger_cluster", "overlap")
  sankey_data$liger_cluster <- paste(sankey_data$liger_cluster, " ", sep = "")
  nodes <- data.frame(name = c(as.character(sankey_data$pca_cluster),
                             as.character(sankey_data$liger_cluster)) %>%
                        unique())
  sankey_data$pca_cluster_id <-
    match(sankey_data$pca_cluster, nodes$name) - 1
  sankey_data$liger_cluster_id <-
    match(sankey_data$liger_cluster, nodes$name) - 1
  my_colour_scale <-
    'd3.scaleOrdinal() .range(["red","cornflowerblue","green","blue","darkorange","grey"])'
  sankey_plot <- networkD3::sankeyNetwork(Links = sankey_data, Nodes = nodes,
                                          Source = "pca_cluster_id",
                                          Target = "liger_cluster_id",
                                          Value = "overlap", NodeID = "name",
                                          sinksRight = TRUE,
                                          colourScale = my_colour_scale,
                                          nodeWidth = 80, fontSize = 13,
                                          nodePadding = 20)

  sce@metadata$dataset_integration$clustering_plots$sankey_plot <- sankey_plot

  pca_proportional_barplots <- list()
  liger_proportional_barplots <- list()

  for (variable in categorical_covariates) {

    pca_barplot_data <-
      data.frame("clusters" = pca_clusters, "group" = sce@colData[[variable]])
    clusters <- pca_barplot_data$clusters
    group <- pca_barplot_data$group

    pca_barplot <- ggplot(pca_barplot_data, aes(x = clusters, fill = group)) +
      geom_bar(position = "fill") + theme_bw() +
      theme(
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
      theme(legend.text = element_text(size = 13))

    pca_proportional_barplots[[variable]] <- pca_barplot

    liger_barplot_data <-
      data.frame("clusters" = liger_clusters,
                 "group" = sce@colData[[variable]])
    clusters <- liger_barplot_data$clusters
    group <- liger_barplot_data$group

    liger_barplot <-
      ggplot(liger_barplot_data, aes(x = clusters, fill = group)) +
      geom_bar(position = "fill") + theme_bw() +
      theme(
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
      theme(legend.text = element_text(size = 13))

    liger_proportional_barplots[[variable]] <- liger_barplot
  }
  sce@metadata$dataset_integration$clustering_plots$pca_proportional_barplots <-
    pca_proportional_barplots
  sce@metadata$dataset_integration$clustering_plots$liger_proportional_barplots <-
    liger_proportional_barplots

  return(sce)
}
