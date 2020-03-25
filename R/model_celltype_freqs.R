################################################################################
#' Model Changes in Celltype Frequencies with Dirichlet Multinomial Regression
#'
#' @param sce a SingleCellExperiment object
#' @param unique_id_var the unique identifier variable for each sample
#' @param celltype_var the colData variable specifying celltype or subtype
#' @param dependent_var the name of the colData variable for contrasts
#' @param ref_class the class of dependent_var used as reference
#'
#' @return results_l a list of results
#'
#' @family Further analyses
#'
#' @importFrom cli cli_h1 cli_alert
#' @importFrom DirichletReg DR_data DirichReg
#'
#' @export
model_celltype_freqs <- function(sce,
                                 unique_id_var = "manifest",
                                 celltype_var = "cluster_celltype",
                                 dependent_var = "group",
                                 ref_class = "Control",
                                 var_order = NULL,
                                 ...) {

  fargs <- c(as.list(environment()), list(...))

  cli::cli_h1("Modelling Cell-type Frequencies")

  covariates <- do.call(.retrieve_covariates, fargs)

  cli::cli_alert("Covariates retrieved for {.val {nrow(covariates)}} samples")

  mat <- do.call(.tally_cells, fargs)

  cli::cli_alert(
    "Cell frequencies calculated across {.val {dim(mat)[[2]]}} cell-types")

  prop_mat <- prop.table(mat, margin = 1)
  df <- as.data.frame(prop_mat)
  df$counts <- DirichletReg::DR_data(df)
  df[[unique_id_var]] <- as.factor(rownames(df))
  df <- cbind(df, covariates)

  cli::cli_h2("Fitting Dirichlet Model")
  model_formula <- as.formula(sprintf("counts ~ %s", dependent_var))
  cli::cli_alert(
    "Fitting model: {.var {scFlow:::.formula_to_char(model_formula)}}")
  fit <- do.call(
    DirichletReg::DirichReg,
    list(formula = model_formula, data = df)
    )

  cli::cli_alert("Post-processing model")
  pvals <- do.call(.process_dirichlet_fit,
                   c(list(fit = fit),
                     fargs))

  cli::cli_alert_success("Dirichlet Model fit successfully.")

  results <- list()
  results$counts_mat <- mat
  results$prop_counts_mat <- prop_mat
  df$counts <- NULL
  results$df <- df
  results$fit <- fit
  results$pvals <- pvals

  results$plot_table <- do.call(.prepare_dirichlet_plot_table,
                                c(list(
                                  df = df,
                                  celltypes = unique(sce[[celltype_var]])
                                ), fargs))

  results$plot <- do.call(
    .plot_dirichlet_results,
    c(list(
      df = results$plot_table,
      n_groups = length(unique(df[[dependent_var]]))
      ), fargs)
  )

  return(results)

}

.prepare_dirichlet_plot_table <- function(df,
                                          celltypes,
                                          celltype_var,
                                          dependent_var,
                                          var_order = NULL,
                                          ...) {

  x <- tidyr::pivot_longer(results$df, cols = all_of(celltypes), names_to = celltype_var) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::select(!!dependent_var, !!celltype_var, value) %>%
    dplyr::group_by(!!(as.name(dependent_var)), !!(as.name(celltype_var))) %>%
    dplyr::summarise_each(funs(mean, sd, se = sd(.) / sqrt(n()))) %>%
    dplyr::left_join(pvals$pvals, by = c(celltype_var, dependent_var)) %>%
    dplyr::mutate(label = case_when(
      is.na(label) ~ "",
      !is.na(label) ~ label))

  if (!is.null(var_order)) {
    assertthat::assert_that(all(var_order %in% x[[dependent_var]]))
    x[[dependent_var]] <- factor(
      x[[dependent_var]], levels = var_order)
  }
  return(x)

}


.plot_dirichlet_results <- function(df, n_groups, individual_plots = FALSE, ...) {

  fargs <- list(...)

  if (is.null(fargs$palette)) {
    if(n_groups <= 10) palette <- paletteer::paletteer_d("ggsci::default_aaas")
    if(n_groups > 10) palette <- paletteer::paletteer_d("ggsci::default_igv")
  } else {
    palette <- fargs$palette
  }

  p <- ggplot(df, aes(x = group, y = mean)) +
    geom_col(aes(fill = .data[[dependent_var]]), colour = "black") +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                  position=position_dodge(.9)) +
    geom_text(aes(y = (mean+se) * 1.05, label = label), size = 5) +
    ylab("Relative Proportion") +
    xlab(NULL) +
    facet_grid(~ .data[[celltype_var]], scales = "free_y", switch = "x") +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
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
          panel.border = element_rect(colour = "black"))

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
#' @importFrom dplyr select group_by count
#' @importFrom tidyr pivot_wider
#'
#' @keywords internal
.tally_cells <- function(sce,
                         unique_id_var = "manifest",
                         celltype_var = "cluster_celltype",
                         ...) {

  mat <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::select(!!unique_id_var, !!celltype_var) %>%
    dplyr::group_by(!!(as.name(unique_id_var))) %>%
    dplyr::count(!!(as.name(celltype_var))) %>%
    tidyr::pivot_wider(names_from = !!(as.name(celltype_var)), values_from = "n") %>%
    as.data.frame()

  mat <- mat[order(mat[[unique_id_var]]),]
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
                                 ref_class,
                                 ...
                                 ) {

  covariates <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::select(!!unique_id_var, !!dependent_var) %>%
    unique()

  rownames(covariates) <- covariates[[unique_id_var]]
  covariates <- covariates[order(covariates[[unique_id_var]]),]
  covariates[[unique_id_var]] <- NULL
  covariates[[dependent_var]] <- relevel(
    covariates[[dependent_var]], ref = ref_class
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
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
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
      label = case_when(
        pval <= 0.001 ~ "***",
        pval <= 0.01 ~ "**",
        pval <= 0.05 ~ "*",
        pval > 0.05 ~ "")) %>%
    as.data.frame()

  return(pvals)

}
