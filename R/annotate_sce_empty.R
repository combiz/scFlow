################################################################################
#' Annotate a SingleCellExperiment with EmptyDrops predictions
#'
#' Calculates the following QC metric:
#' * is_empty_drop - is the barcode an empty drop
#'
#' @param sce a SingleCellExperiment object
#' @param lower see [dropletUtils::emptyDrops()]
#' @param retain see [dropletUtils::emptyDrops()]
#' @param alpha minimum FDR for [dropletUtils::emptyDrops()]
#' @param ... additional arguments for uniformity testing
#'
#' @return sce a SingleCellExperiment object annotated with emptyDrops metrics
#'
#' @family annotation functions
#' @importFrom cli cli_alert_danger rule cli_alert_success cli_text
#' @importFrom uniftest kolmogorov.unif.test
#' @importFrom DropletUtils emptyDrops
#'
#' @export
annotate_sce_empty <- function(sce,
                               lower = 100,
                               retain = NULL,
                               alpha_cutoff = 0.001,
                               ...) {

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  if (sce@metadata$scflow_steps$cells_filtered == 1) {
    stop(cli::cli_alert_danger(
      "Ambient RNA modelling must be performed on an unfiltered SCE."))
  }

  cat(cli::rule("Annotating Empty Droplets", line = 2), "\r\n")
  cat(cli::rule(
    "Running Simulations", line = 1), "\r\n")

  # add citations and print
  #sce <- .append_citation_sce(sce, key = c("emptydrops"))

  # test ambient to evaluate performance with histogram
  genes_idx <- rowSums(SingleCellExperiment::counts(sce)) > 0 # drop non-exp
  ed_out <- DropletUtils::emptyDrops(
    SingleCellExperiment::counts(sce[genes_idx, ]),
    lower = lower,
    test.ambient = TRUE
  )

  # for uniformity testing
  pval_dt <- data.frame(
    "pval" = ed_out$PValue[ed_out$Total <= lower & ed_out$Total > 0])

  unif_test <- .perform_ks_uniformity_test(pval_dt$pval)

  # plot histogram
  p <- .plot_emptydrops_distribution(pval_dt)
  sce@metadata$emptydrops_plots$emptydrops_hist <- p

  cat(cli::rule(
    "Running Simulations to Identify Empty Droplets", line = 1), "\r\n")
  # main fn - actual droplet calling
  ed_out <- DropletUtils::emptyDrops(
    SingleCellExperiment::counts(sce[genes_idx, ]),
    lower = lower,
    test.ambient = FALSE
  )

  # process ed_out for storage and plotting
  dt <- as.data.frame(ed_out)
  dt$is_empty_drop <- TRUE
  dt$is_empty_drop[which(dt$FDR <= alpha_cutoff)] <- FALSE
  sce$is_empty_drop <- dt$is_empty_drop
  dt <- tidyr::drop_na(dt)
  dt <- dt[dt$Total >= lower, ]

  # plot -log likelihood by counts
  p <- .plot_emptydrops_loglike_vs_counts(dt)
  sce@metadata$emptydrops_plots$log_likelihood_by_total_counts <- p

  # save plot data
  sce@metadata$emptydrops_plots_data$log_likelihood_by_total_counts <- dt

  # save params
  emptydrops_params <- list()
  emptydrops_params$lower <- lower
  emptydrops_params$retain <- retain
  emptydrops_params$alpha_cutoff <- alpha_cutoff
  emptydrops_params$emptydrops_found <- sum(sce$is_empty_drop)

  # save calculated params
  emptydrops_params$uniformity_method <- unif_test$method[[1]]
  emptydrops_params$uniformity_statistic <- unif_test$statistic[[1]]
  emptydrops_params$uniformity_pval <- unif_test$p.value[[1]]

  sce@metadata$emptydrops_params <- emptydrops_params

  emptydrops_params <- purrr::map_df(
    unlist(emptydrops_params), ~ .) %>%
    dplyr::rename_all(~ paste0("emptydrops_",.))

  sce@metadata$emptydrops_params <- emptydrops_params

  cli::cli_text(c("Identified {.val {sum(sce$is_empty_drop)}} empty droplets ",
                  "from {.val {length(sce$is_empty_drop)}} barcodes."))
  cli::cli_alert_success("EmptyDrops algorithm ran successfully.")
  cat(cli::rule(line = 1), "\r\n")

  sce@metadata$scflow_steps$emptydrops_annotated <- 1
  #sce$qc_metric_passed <- sce$qc_metric_passed & !sce$is_empty_drop

  return(sce)

}

#' helper fn - Kolmogorov-Smirnov Uniformity Test
#' @importFrom uniftest kolmogorov.unif.test
#' @importFrom cli cli_alert_warning rule
#' @export
#' @keywords internal
.perform_ks_uniformity_test <- function(x, nrepl = 1000) {
  cat(cli::rule("Evaluating Distribution Uniformity", line = 1), "\r\n")
  unif_test <- uniftest::kolmogorov.unif.test(x, nrepl = nrepl, k = 0)

  if(unif_test$p.value <= 0.05) {
    cli::cli_alert_warning(c(
      "Non-uniform distribution ",
      "{.val (Kolmogorov-Smirnov P-Value {unif_test$p.value})}: ",
      "consider decreasing {.var lower} parameter."))
  }
  return(unif_test)
}

#' helper fn - histogram plot to confirm uniform distribution
#' @export
#' @keywords internal
.plot_emptydrops_distribution <- function(dt) {

  brx <- pretty(range(dt$pval),
                n = nclass.Sturges(dt$pval), min.n = 1)

  p <- ggplot(data = dt) +
    geom_histogram(aes(x = pval), breaks = brx, fill = "grey80", colour = "black") +
    xlab("P-Value") +
    ylab("Frequency") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 16),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(hjust = 0.5, vjust = 0.5))

  return(p)

}



#' helper fn - histogram plot to confirm uniform distribution
#' @export
#' @keywords internal
.plot_emptydrops_loglike_vs_counts <- function(dt) {

  p <- ggplot2::ggplot(dt) +
    geom_point(
      aes(x = Total/1000, y = -LogProb/1000, colour = is_empty_drop),
      size = .01, alpha = 0.1) +
    scale_colour_manual(values = c("darkblue", "red")) +
    labs(x = bquote("Total counts (x"*10^3*")"),
         y = bquote("-Log likelihood (x"*10^3*")"))+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title=element_text(size=16),
          legend.text=element_text(size=10),
          plot.title = element_text(size = 18, hjust = 0.5))

  return(p)

}



