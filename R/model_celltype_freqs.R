################################################################################
#' Model Changes in Celltype Frequencies with Dirichlet Multinomial Regression
#'
#' @param sce a SingleCellExperiment object
#' @param unique_id_var the unique identifier variable for each sample
#' @param celltype_var the colData variable specifying celltype or subtype
#' @param dependent_var the name of the colData variable for contrasts
#' @param ref_class the class of dependent_var used as reference
#' @param var_order Optional re-ordering of subset_group factor levels. Default
#' NULL.
#' @param confounding_vars Additional confounding variables from colData.
#' @param ... Additional arguments
#'
#' @return results_l a list of results
#'
#' @family Celltype annotation
#'
#' @importFrom cli cli_h1 cli_alert
#' @importFrom DirichletReg DR_data DirichReg
#' @importFrom stats as.formula fisher.test p.adjust relevel
#'
#' @export
model_celltype_freqs <- function(sce,
                                 unique_id_var = "manifest",
                                 celltype_var = "cluster_celltype",
                                 dependent_var = "group",
                                 ref_class = "Control",
                                 var_order = NULL,
                                 confounding_vars = NULL,
                                 ...) {

  if (is.null(var_order)){ var_order <- levels(as.factor(sce[[dependent_var]]))}
  fargs <- c(as.list(environment()), list(...))

  cli::cli_h1("Modelling Cell-type Frequencies")

  covariates <- do.call(.retrieve_covariates, fargs)

  cli::cli_alert("Covariates retrieved for {.val {nrow(covariates)}} samples")

  mat <- do.call(.tally_cells, fargs)

  cli::cli_alert(
    "Cell frequencies calculated across {.val {dim(mat)[[2]]}} cell-types"
  )

  prop_mat <- prop.table(mat, margin = 1)
  df <- as.data.frame(prop_mat)
  df <- df[order(rownames(df)), ]
  df$counts <- DirichletReg::DR_data(df)
  df[[unique_id_var]] <- rownames(df)
  df <- cbind(df, covariates)

  cli::cli_h2("Fitting Dirichlet Model")

  if (!is.null(confounding_vars)){

    model_formula <- stats::as.formula(sprintf("counts ~ %s + %s",
                                               dependent_var,
                                               paste(confounding_vars, collapse = " + ")))
    cli::cli_alert(
      "Fitting model: {.var {scFlow:::.formula_to_char(model_formula)}}"
    )

  } else {

  model_formula <- stats::as.formula(sprintf("counts ~ %s", dependent_var))
  cli::cli_alert(
    "Fitting model: {.var {scFlow:::.formula_to_char(model_formula)}}"
  )

  }
  fit <- do.call(
    DirichletReg::DirichReg,
    #list(formula = model_formula, data = df, model = "alternative")
    list(formula = model_formula, data = df)
  )

  cli::cli_alert("Post-processing model")
  pvals <- do.call(
    .process_dirichlet_fit,
    c(
      list(fit = fit),
      fargs
    )
  )


  cli::cli_alert_success("Dirichlet Model fit successfully.")

  results <- list()
  results$counts_mat <- mat
  results$prop_counts_mat <- prop_mat
  results$DR_data_df <- df
  df$counts <- NULL
  results$dirichlet_fit <- fit
  results$dirichlet_pvals <- pvals

  pvals <- pvals %>%
    dplyr::filter(get(dependent_var) %in%
                    setdiff(unique(df[[dependent_var]]), ref_class))

  results$dirichlet_plot_table <- do.call(
    .prepare_dirichlet_plot_table,
    c(
      list(
        df = as.data.frame(df),
        celltypes = unique(sce[[celltype_var]]),
        pvals = pvals
      ),
      fargs
    )
  )

  # combined plots with faceting
  results$dirichlet_plot <- do.call(
    .plot_dirichlet_results,
    c(list(
      df = results$dirichlet_plot_table,
      n_groups = length(unique(df[[dependent_var]])),
      facet_plots = TRUE
    ), fargs)
  )
  #results$dirichlet_plot$plot_env$sce <- NULL

  # individual celltype plots for dirichlet
  results$dirichlet_plots_by_celltype <- list()
  for (celltype in unique(results$dirichlet_plot_table[[celltype_var]])) {
    results$dirichlet_plots_by_celltype[[celltype]] <- do.call(
      .plot_dirichlet_results,
      c(list(
        df = results$dirichlet_plot_table[
          results$dirichlet_plot_table[[celltype_var]] == celltype,],
        n_groups = length(unique(df[[dependent_var]])),
        facet_plots = TRUE
      ), fargs)
    )
  }

  results$counts_df <- do.call(
    .prepare_fisher_counts_table,
    c(
      list(
        df = cbind(as.data.frame(mat), covariates),
        celltypes = unique(sce[[celltype_var]]),
        pvals = pvals
      ),
      fargs
    )
  )

  results$fisher_df <- do.call(
    .model_fisher_celltype,
    c(
      list(
        counts_df = results$counts_df
      ),
      fargs
    )
  )

  results$unique_id_plot_table <- do.call(
    .prepare_unique_id_var_plot_table,
    c(
      list(
        df = results$DR_data_df,
        celltypes = unique(sce[[celltype_var]])
      ),
      fargs
    )
  )

  results$unique_id_plots_by_celltype <- list()
  for (celltype in unique(results$unique_id_plot_table[[celltype_var]])) {
    results$unique_id_plots_by_celltype[[celltype]] <- do.call(
      .plot_unique_id_var,
      c(list(
        df = results$unique_id_plot_table[
          results$unique_id_plot_table[[celltype_var]] == celltype,],
        n_groups = length(unique(df[[dependent_var]])),
        facet_plots = FALSE
      ), fargs)
    )
  }

  results$fargs <- fargs
  results$fargs$sce <- NULL

  return(results)
}


################################################################################
#' Prepare data
#'
#' Arrange proportional cell numbers by celltype, group, individual
#'
#' @family helper
#'
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom assertthat assert_that
#'
#' @keywords internal
.prepare_unique_id_var_plot_table <- function(df,
                                              celltypes,
                                              celltype_var,
                                              dependent_var,
                                              unique_id_var,
                                              var_order = NULL,
                                              ...) {

  x <- df
  x$counts <- NULL

  if (!is.null(var_order)) {
    assertthat::assert_that(all(var_order %in% x[[dependent_var]]))
    x[[dependent_var]] <- factor(
      x[[dependent_var]],
      levels = var_order
    )
  }

  x <- x[order(x[[dependent_var]]), ]
  x[[unique_id_var]] <- factor(x[[unique_id_var]], levels = x[[unique_id_var]])
  x <- x %>% tidyr::pivot_longer(
    # cols = setdiff(colnames(x),
    #                c(dependent_var, unique_id_var)),
    cols = all_of(celltypes),
    names_to = celltype_var,
    values_to = "cells_pc"
  )

  x <- as.data.frame(x)

  return(x)
}


################################################################################
#' Plot relative celltype proportions by sample
#'
#' @family helper
#'
#' @importFrom paletteer paletteer_d
#'
#' @keywords internal
.plot_unique_id_var <- function(df,
                                n_groups,
                                celltype_var,
                                dependent_var,
                                unique_id_var,
                                facet_plots = FALSE,
                                ...) {
  fargs <- list(...)

  if (is.null(fargs$palette)) {
    if (n_groups <= 10) palette <- paletteer::paletteer_d("ggsci::default_aaas")
    if (n_groups > 10) palette <- paletteer::paletteer_d("ggsci::default_igv")
  } else {
    palette <- fargs$palette
  }

  p <- ggplot(df, aes(x = .data[[unique_id_var]], y = cells_pc)) +
    geom_col(aes(fill = .data[[dependent_var]]), colour = "black") +
    ylab("Relative Proportion") +
    xlab(NULL) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
      axis.title = element_text(size = 18),
      axis.text.y = element_text(size = 16, colour = "black"),
      axis.text.x = element_text(size = 16, colour = "black", angle = 45, hjust = 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.position = "top",
      legend.justification = "left",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_text(size = 16),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black")
    )


  if(facet_plots == TRUE) {
    p <- p +
      facet_wrap(~ .data[[celltype_var]], ncol = 3)
  }

  p <- ggplot2::ggplotGrob(p)
  p <- ggpubr::as_ggplot(p)
  return(p)
}

################################################################################
#' Prepare data for `.model_fisher_celltype`
#'
#' Sums absolute cell numbers by celltypes and groups
#'
#' @family helper
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate_if select group_by summarise_each left_join mutate
#' @importFrom assertthat assert_that
#'
#' @keywords internal
.prepare_fisher_counts_table <- function(df,
                                         pvals,
                                         celltypes,
                                         celltype_var,
                                         dependent_var,
                                         var_order = NULL,
                                         ...) {
  x <- tidyr::pivot_longer(df,
                           cols = all_of(celltypes),
                           names_to = celltype_var) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::select(!!dependent_var, !!celltype_var, value) %>%
    dplyr::group_by(!!(as.name(dependent_var)),
                    !!(as.name(celltype_var))) %>%
    dplyr::summarise(
      tibble::tibble(
        dplyr::across(where(is.numeric), ~sum(.x), .names = "sum")
      )) %>%
    dplyr::left_join(pvals, by = c(celltype_var, dependent_var))

  if (!is.null(var_order)) {
    assertthat::assert_that(all(var_order %in% x[[dependent_var]]))
    x[[dependent_var]] <- factor(
      x[[dependent_var]],
      levels = var_order
    )
  }

  x <- as.data.frame(x)

  return(x)
}

################################################################################
#' Prepare data for `.plot_dirichlet_results`
#'
#' @family helper
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate_if select group_by summarise_each left_join mutate n
#' @importFrom assertthat assert_that
#'
#' @keywords internal
.prepare_dirichlet_plot_table <- function(df,
                                          pvals,
                                          celltypes,
                                          celltype_var,
                                          dependent_var,
                                          var_order = NULL,
                                          ...) {
  x <- tidyr::pivot_longer(df,
                           cols = all_of(celltypes),
                           names_to = celltype_var) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::select(!!dependent_var, !!celltype_var, value) %>%
    dplyr::group_by(!!(as.name(dependent_var)),
                    !!(as.name(celltype_var))) %>%
    dplyr::summarise(
      tibble::tibble(
        dplyr::across(where(is.numeric), ~mean(.x), .names = "mean"),
        dplyr::across(where(is.numeric), ~sd(.x), .names = "sd"),
        dplyr::across(where(is.numeric), ~ sd(.x) / sqrt(dplyr::n()), .names = "se")
      )) %>%
    dplyr::left_join(pvals, by = c(celltype_var, dependent_var)) %>%
    dplyr::mutate(label = case_when(
      is.na(label) ~ "",
      !is.na(label) ~ label
    ))

  if (!is.null(var_order)) {
    assertthat::assert_that(all(var_order %in% x[[dependent_var]]))
    x[[dependent_var]] <- factor(
      x[[dependent_var]],
      levels = var_order
    )
  }

  x <- as.data.frame(x)

  return(x)
}

################################################################################
#' Plot relative celltype proportions with Dirichlet model significance
#'
#' @family helper
#'
#' @importFrom paletteer paletteer_d
#'
#' @keywords internal
.plot_dirichlet_results <- function(df,
                                    n_groups,
                                    celltype_var,
                                    dependent_var,
                                    facet_plots = FALSE,
                                    ...) {
  fargs <- list(...)

  if (is.null(fargs$palette)) {
    if (n_groups <= 10) palette <- paletteer::paletteer_d("ggsci::default_aaas")
    if (n_groups > 10) palette <- paletteer::paletteer_d("ggsci::default_igv")
  } else {
    palette <- fargs$palette
  }

  p <- ggplot(df, aes(x = .data[[dependent_var]], y = mean)) +
    geom_col(aes(fill = .data[[dependent_var]]), colour = "black") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
      width = .2,
      position = position_dodge(.9)
    ) +
    geom_text(aes(y = (mean + se) * 1.05, label = label), size = 6) +
    ylab("Relative Proportion") +
    xlab(NULL) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
      axis.title = element_text(size = 18),
      axis.text.y = element_text(size = 16, colour = "black"),
      axis.text.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.position = "top",
      legend.justification = "left",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_text(size = 16),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black")
    )

  if(facet_plots == TRUE) {
    p <- p +
      facet_grid(~ .data[[celltype_var]], scales = "free_y", switch = "x")
  }

  p <- ggplot2::ggplotGrob(p)
  p <- ggpubr::as_ggplot(p)
  #p$plot_env$sce <- NULL
  #p$plot_env$... <- NULL
  #p$plot_env$fargs$sce <- NULL

  return(p)
}

################################################################################
#' Tally the cells in each celltype for each unique sample
#'
#' @param unique_id_var the unique identifier variable for each sample
#' @param celltype_var the colData variable specifying celltype or subtype
#'
#' @return mat a matrix sample x celltypes with absolute cell numbers
#'
#' @family helper
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom dplyr select group_by count mutate_if
#' @importFrom tidyr pivot_wider
#'
#' @keywords internal
.tally_cells <- function(sce,
                         unique_id_var = "manifest",
                         celltype_var = "cluster_celltype",
                         ...) {
  mat <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::select(!!unique_id_var, !!celltype_var) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::group_by(!!(as.name(unique_id_var))) %>%
    dplyr::count(!!(as.name(celltype_var))) %>%
    tidyr::pivot_wider(names_from = !!(as.name(celltype_var)), values_from = "n") %>%
    as.data.frame()

  mat <- mat[order(mat[[unique_id_var]]), ]
  rownames(mat) <- mat[[unique_id_var]]
  mat[[unique_id_var]] <- NULL
  mat[is.na(mat)] <- 0
  mat <- as.matrix(mat)

  return(mat)
}


################################################################################
#' Retrieve the covariates from a SCE colData
#'
#' @param sce a SingleCellExperiment object
#' @param unique_id_var the unique identifier variable for each sample
#' @param dependent_var the name of the colData variable for contrasts
#' @param ref_class the class of dependent_var used as reference
#'
#' @return mat a dataframe with covariates ordered by unique_id_var
#'
#' @family helper
#'
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom dplyr select
#'
#' @keywords internal
.retrieve_covariates <- function(sce,
                                 unique_id_var,
                                 dependent_var,
                                 confounding_vars,
                                 ref_class,
                                 ...) {
  covariates <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::select(!!unique_id_var, !!dependent_var, !!confounding_vars) %>%
    unique()

  rownames(covariates) <- covariates[[unique_id_var]]
  covariates <- covariates[order(covariates[[unique_id_var]]), ]
  covariates[[unique_id_var]] <- NULL
  covariates[[dependent_var]] <- stats::relevel(
    as.factor(covariates[[dependent_var]]),
    ref = ref_class
  )

  return(covariates)
}

################################################################################
#' Process the fit model to generate an accessible summary data.frame
#'
#' @param fit a model generated by `DirichletReg::DirichReg`
#' @param dependent_var the name of the colData variable for contrasts
#' @param celltype_var the colData variable specifying celltype or subtype
#'
#' @return pvals a long format data frame with accessible pvalues
#'
#' @family helper
#'
#' @importFrom tibble rownames_to_column tibble
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr case_when
#'
#' @keywords internal
.process_dirichlet_fit <- function(fit, dependent_var, celltype_var, ...) {
  u <- summary(fit)
  pvals <- u$coef.mat[grep("Intercept", rownames(u$coef.mat), invert = TRUE), 4]
  v <- names(pvals)
  pvals <- matrix(pvals, ncol = length(u$varnames))
  rownames(pvals) <- gsub(as.name(dependent_var), "", v[1:nrow(pvals)])
  colnames(pvals) <- u$varnames

  pvals <- as.data.frame(t(pvals)) %>%
    tibble::rownames_to_column(var = celltype_var) %>%
    tidyr::pivot_longer(
      cols = rownames(pvals),
      names_to = dependent_var,
      values_to = "pval"
    ) %>%
    dplyr::mutate(
      padj = p.adjust(pval, method = "bonferroni"),
      label = dplyr::case_when(
        padj <= 0.001 ~ "***",
        padj <= 0.01 ~ "**",
        padj <= 0.05 ~ "*",
        padj > 0.05 ~ "",
        is.na(padj) ~ ""
      )
    ) %>%
    as.data.frame()

  return(pvals)
}

################################################################################
#' Run a fisher exact test on absolute cell numbers for each celltype
#'
#' @param counts_df absolute cell number counts stratified
#' @param dependent_var the name of the colData variable for contrasts
#' @param celltype_var the colData variable specifying celltype or subtype
#' @param ref_class the class of dependent_var used as reference
#'
#' @return df a data frame of cell numbers with fisher p values adjusted
#'
#' @family helper
#'
#' @importFrom  dplyr mutate case_when
#'
#' @keywords internal
.model_fisher_celltype <- function(counts_df,
                                   dependent_var,
                                   celltype_var,
                                   ref_class,
                                   ...) {
  df <- NULL

  for (celltype in unique(counts_df[[celltype_var]])) {
    ct <- counts_df[counts_df[[celltype_var]] == celltype, ]
    not_ct <- counts_df[counts_df[[celltype_var]] != celltype, ]
    for (contrast_class in setdiff(ct[[dependent_var]], ref_class)) {
      contab <- matrix(0, nrow = 2, ncol = 2)
      contab[1, 1] <- ct[ct[[dependent_var]] == ref_class, ]$sum
      contab[1, 2] <- sum(not_ct[not_ct[[dependent_var]] == ref_class, ]$sum)
      contab[2, 1] <- ct[ct[[dependent_var]] == contrast_class, ]$sum
      contab[2, 2] <- sum(not_ct[not_ct[[dependent_var]] == contrast_class, ]$sum)
      rownames(contab) <- c(ref_class, contrast_class)
      colnames(contab) <- c("cells", "other_cells")
      res <- stats::fisher.test(contab)
      contab_df <- as.data.frame(contab)
      contab_df[[dependent_var]] <- rownames(contab_df)
      contab_df[[celltype_var]] <- celltype
      contab_df$pval <- NA
      contab_df[contab_df[[dependent_var]] == contrast_class, ]$pval <- res$p.value
      df <- unique(rbind(df, contab_df))
    }
  }
  df$padj <- stats::p.adjust(df$pval, method = "bonferroni")
  df <- df %>%
    dplyr::mutate(
      label = dplyr::case_when(
        padj <= 0.001 ~ "***",
        padj <= 0.01 ~ "**",
        padj <= 0.05 ~ "*",
        padj > 0.05 ~ "",
        is.na(padj) ~ ""
      )
    )
  rownames(df) <- NULL
  df <- as.data.frame(df)
  return(df)
}
