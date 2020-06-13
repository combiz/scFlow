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
#' @import cli
#' @importFrom uniftest kolmogorov.unif.test
#' @importFrom Matrix rowSums colSums
#' @importFrom DropletUtils emptyDrops
#'
#' @export
find_cells <- function(sce,
                       lower = 100,
                       retain = NULL,
                       alpha_cutoff = 0.001,
                       niters = 10000,
                       ...) {
  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  if (sce@metadata$scflow_steps$cells_filtered == 1 |
    sce@metadata$scflow_steps$cells_annotated == 1) {
    stop(cli::cli_alert_danger(
      "Ambient RNA modelling must be performed before annotation/filtering."
    ))
  }

  cli::cli_h1("Annotating Empty Droplets")
  sce <- .append_citation_sce(sce, key = c("emptydrops"))
  cli::cli_h2("Running Simulations")
  cli::cli_text(
    "Generating ambient model with barcodes with fewer than ",
    "{.val {lower}} counts."
  )

  # add citations and print


  # test ambient to evaluate performance with histogram
  cli::cli_text(
    "Performing {.val {prettyNum(niters, big.mark = ", ")}} ",
    "Monte Carlo iterations..."
  )
  genes_idx <- Matrix::rowSums(SingleCellExperiment::counts(sce)) > 0
  ed_out <- DropletUtils::emptyDrops(
    SingleCellExperiment::counts(sce[genes_idx, ]),
    lower = lower,
    retain = retain,
    test.ambient = TRUE,
    niters = niters
  )

  # for uniformity testing
  pval_dt <- data.frame(
    "pval" = ed_out$PValue[ed_out$Total <= lower & ed_out$Total > 0]
  )

  unif_test <- .perform_ks_uniformity_test(pval_dt$pval)

  # plot histogram
  p <- .plot_emptydrops_distribution(pval_dt)
  sce@metadata$emptydrops_plots$emptydrops_hist <- p

  cli::cli_h2(
    "Running Simulations to Identify Empty Droplets"
  )
  cli::cli_text(
    "Performing {.val {prettyNum(niters, big.mark = ", ")}} ",
    "Monte Carlo iterations..."
  )
  # main fn - actual droplet calling
  ed_out <- DropletUtils::emptyDrops(
    SingleCellExperiment::counts(sce[genes_idx, ]),
    lower = lower,
    retain = retain,
    test.ambient = FALSE,
    niters = niters
  )

  # p values limited by niters?
  lim_dt <- as.data.frame(table(
    Significant = ed_out$FDR <= alpha_cutoff,
    Limited = ed_out$Limited
  ))

  n_limited <-
    lim_dt[lim_dt$Significant == FALSE & lim_dt$Limited == TRUE, ]$Freq
  if (n_limited > 0) {
    cli::cli_alert_warning(c(
      "For {.val {n_limited}} non-significant barcodes ",
      "the p-value was limited by the number of iterations ",
      "(iters: {.val {niters}}).  Consider increasing {.var niters}."
    ))
  }

  # process ed_out for storage and plotting
  dt <- as.data.frame(ed_out)
  dt$is_empty_drop <- TRUE
  dt$is_empty_drop[which(dt$FDR <= alpha_cutoff)] <- FALSE
  sce$is_empty_drop <- dt$is_empty_drop
  dt <- tidyr::drop_na(dt)
  dt <- dt[dt$Total >= lower, ]

  # plot -log likelihood by counts
  p <- .plot_emptydrops_loglike(dt)
  sce@metadata$emptydrops_plots$log_likelihood_by_total_counts <- p

  # save plot data
  sce@metadata$emptydrops_plots_data$log_likelihood_by_total_counts <- dt

  # save params
  emptydrops_params <- list()
  emptydrops_params$lower <- lower
  emptydrops_params$retain <- retain
  emptydrops_params$alpha_cutoff <- alpha_cutoff
  emptydrops_params$niters <- niters

  # save calculated params
  emptydrops_params$emptydrops_found <- sum(sce$is_empty_drop)
  emptydrops_params$cells_found <- sum(!sce$is_empty_drop)
  emptydrops_params$n_limited <- n_limited
  emptydrops_params$uniformity_method <- unif_test$method[[1]]
  emptydrops_params$uniformity_statistic <- unif_test$statistic[[1]]
  emptydrops_params$uniformity_pval <- unif_test$p.value[[1]]

  sce@metadata$emptydrops_params <- emptydrops_params

  cli::cli_text(c(
    "Identified {.val {sum(!sce$is_empty_drop)}} cells and ",
    "{.val {sum(sce$is_empty_drop)}} empty droplets ",
    "from {.val {length(sce$is_empty_drop)}} barcodes."
  ))
  cli::cli_alert_success("EmptyDrops algorithm ran successfully.")
  cat(cli::rule(line = 1), "\r\n")

  sce@metadata$scflow_steps$emptydrops_annotated <- 1

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

  if (unif_test$p.value <= 0.05) {
    cli::cli_alert_warning(c(
      "Non-uniform distribution ",
      "{.val (Kolmogorov-Smirnov P-Value {unif_test$p.value})}: ",
      "consider decreasing {.var lower} parameter."
    ))
  } else {
    cli::cli_alert_success(c(
      "Distribution looks uniform ",
      "{.val (Kolmogorov-Smirnov P-Value {unif_test$p.value})}: ",
      "emptyDrops should be performant."
    ))
  }
  return(unif_test)
}

#' helper fn - histogram plot to confirm uniform distribution
#' @export
#' @keywords internal
.plot_emptydrops_distribution <- function(dt) {
  brx <- pretty(range(dt$pval),
    n = nclass.Sturges(dt$pval), min.n = 1
  )

  p <- ggplot(data = dt) +
    geom_histogram(aes(x = pval),
                   breaks = brx, fill = "grey80", colour = "black") +
    xlab("P-Value") +
    ylab("Frequency") +
    theme_bw() +
    theme(
      legend.position = "none",
      text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(hjust = 0.5, vjust = 0.5)
    )

  p <- .grobify_ggplot(p)
  return(p)
}



#' helper fn - histogram plot to confirm uniform distribution
#' @export
#' @keywords internal
.plot_emptydrops_loglike <- function(dt) {
  p <- ggplot2::ggplot(dt) +
    geom_point(
      aes(x = Total / 1000, y = -LogProb / 1000, colour = is_empty_drop),
      size = .01, alpha = 0.1
    ) +
    scale_colour_manual(values = c("darkblue", "red")) +
    labs(
      x = bquote("Total counts (x" * 10^3 * ")"),
      y = bquote("-Log likelihood (x" * 10^3 * ")")
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 18, hjust = 0.5)
    )

  p <- .grobify_ggplot(p)
  return(p)
}
