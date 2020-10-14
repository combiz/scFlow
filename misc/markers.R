library(scFlow)
library(monocle3)
sce <- read_sce("~/Documents/junk/MS_Custom_Mapped_SCE/")
idx <- as.numeric(caret::createDataPartition(sce$individual, p = .05, list = FALSE)) # 15% subset
mini_sce <- sce[, idx]
sce <- mini_sce


################################################################################
#' Model Changes in Celltype Frequencies with Dirichlet Multinomial Regression
#'
#' @param sce A SingleCellExperiment object
#' @param by_vars The colData variable(s) to group cells by
#' @param fraction_expressing Top marker genes should be expressed in a minimum
#' of this fraction of cells (default 0.10)
#' @param top_n The top_n genes to use for plotting / subset table generation
#' @param n_cores The number of cores to use
#'
#' @return results_l A list of results
#'
#' @family Further analyses
#'
#' @importFrom cli cli_h1 cli_alert
#' @importFrom monocle3 top_markers plot_genes_by_group
#' @importFrom dplyr filter group_by top_n pull
#'
#' @export
find_marker_genes <- function(sce,
                              by_vars = c("cluster_celltype", "clusters"),
                              fraction_expressing = 0.10,
                              top_n = 5,
                              max_point_size = 3,
                              n_cores = future::availableCores()
                              ) {

  assertthat::assert_that(class(sce) == "SingleCellExperiment")
  assertthat::assert_that(
    all(by_vars %in% colnames(SummarizedExperiment::colData(sce)))
  )

  cli::cli_h1("Finding Top Marker Genes")

  cds <- scFlow:::.sce_to_cds(sce)

  results <- list()

  for(by_var in by_vars) {

    cli::cli_alert("Finding markers for {.emph {by_var}}")

    marker_test_res  <- monocle3::top_markers(
      cds,
      group_cells_by = by_var,
      cores = n_cores
      )

    top_specific_markers <- marker_test_res %>%
      dplyr::filter(fraction_expressing >= fraction_expressing) %>%
      dplyr::group_by(cell_group) %>%
      dplyr::top_n(top_n, pseudo_R2)

    top_specific_marker_ids <-
      unique(top_specific_markers %>% dplyr::pull(gene_id))

    marker_plot <- monocle3::plot_genes_by_group(
      cds,
      top_specific_marker_ids,
      group_cells_by = by_var,
      ordering_type = "cluster_row_col",
      max.size = max_point_size
      )

    results[[by_var]] <- list(
      marker_test_res = marker_test_res,
      top_specific_markers = top_specific_markers,
      marker_plot = marker_plot
    )
  }

  return(results)

}
