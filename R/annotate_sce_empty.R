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

  cat(cli::rule("Annotating Empty Droplets", line = 2), "\r\n")
  cat(cli::rule("Running Simulations", line = 1), "\r\n")

  # add citations and print
  sce <- .append_citation_sce(sce, key = c("emptydrops"))

  ed_out <- DropletUtils::emptyDrops(
    SingleCellExperiment::counts(sce),
    lower = lower,
    test.ambient = TRUE
  )

  dt <- data.frame(
    "pval" = ed_out$PValue[ed_out$Total <= lower & ed_out$Total > 0])

  # a quantitative approach to confirming uniformity of distribution
  cat(cli::rule("Evaluating Distribution Uniformity", line = 1), "\r\n")
  unif_test <- kolmogorov.unif.test(dt$pval, nrepl = 2000, k = 0)

  if(unif_test$p.value <= 0.05) {
    cli::cli_alert_warning(c(
    "Non-uniform distribution ",
    "{.val (Kolmogorov-Smirnov P-Value {unif_test$p.value})}: ",
    "consider decreasing {.var lower} parameter."))
  }

  # plot histogram
  p <- .plot_emptydrops_distribution(dt)
  sce@metadata$emptydrops_plots$emptydrops_hist <- p

  # annotate
  sce$is_empty_drop <- FALSE # using recycling
  sce$is_empty_drop[which(ed_out$FDR <= alpha_cutoff)] <- TRUE

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
  sce$qc_metric_passed <- sce$qc_metric_passed & !sce$is_empty_drop

  return(sce)

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
