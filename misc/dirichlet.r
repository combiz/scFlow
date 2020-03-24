
library(dplyr)

unique_id_var <- "manifest"
celltype_var <- "cluster_celltype"
dependent_var <- "group"
ref_class <- "Control"
confounding_vars <- c("sex", "age")


.tally_cells <- function(sce,
                         unique_id_var = "manifest",
                         celltype_var = "cluster_celltype") {

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

.retrieve_covariates <- function(sce, unique_id_var, dependent_var, confounding_vars, ref_class) {

  # retrieve covariates
  covariates <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::select(!!unique_id_var, !!dependent_var, !!confounding_vars) %>%
    unique()

  rownames(covariates) <- covariates[[unique_id_var]]
  covariates <- covariates[order(covariates[[unique_id_var]]),]
  covariates[[unique_id_var]] <- NULL
  # define the reference class
  covariates[[dependent_var]] <- relevel(
    covariates[[dependent_var]], ref = ref_class
    )

  return(covariates)

}

.process_dirichlet_fit <- function(fit, dependent_var) {

  u <- summary(fit)
  pvals <- u$coef.mat[grep("Intercept", rownames(u$coef.mat), invert = TRUE), 4]
  v <- names(pvals)
  pvals <- matrix(pvals, ncol = length(u$varnames))
  rownames(pvals) <- gsub(as.name(dependent_var), "", v[1:nrow(pvals)])
  colnames(pvals) <- u$varnames
  fit$pvals <- pvals
  return(fit)

}

covariates <- .retrieve_covariates(sce,
                                   unique_id_var = unique_id_var,
                                   dependent_var = dependent_var,
                                   confounding_vars = confounding_vars,
                                   ref_class = ref_class
                                   )

mat <- .tally_cells(sce,
                    unique_id_var = unique_id_var,
                    celltype_var = celltype_var)

prop_mat <- prop.table(mat, margin = 1) # relative proportions
df <- as.data.frame(prop_mat)
df$counts <- DirichletReg::DR_data(df)
df[[unique_id_var]] <- as.factor(rownames(df))
df <- cbind(df, covariates)

model_formula <- as.formula(sprintf("counts ~ %s", dependent_var))
fit <- do.call(
  DirichletReg::DirichReg,
  list(formula = model_formula, data = df)
  )

fit <- .process_dirichlet_fit(fit, dependent_var)

df$counts <- NULL



x <- tidyr::pivot_longer(df, cols = intersect(colnames(df), colnames(mat)))

library(ggplot2)
ggplot(x, aes(x = group, y = value)) +
         geom_col(aes(fill = name)) +
  facet_grid(~ name, scales = "free_y")


fit$pvals
p.adjust(fit$pvals[1,], method = "bonferroni") < 0.05

dirichlet_regression = function(counts, covariates, formula) {

  # Dirichlet multinomial regression to detect changes in cell frequencies
  # formula is not quoted, example: counts ~ condition
  # counts is a [samples x cell types] matrix
  # covariates holds additional data to use in the regression
  #
  # Example:
  # counts = do.call(cbind, tapply(seur@data.info$orig.ident, seur@ident, table))
  # covariates = data.frame(condition=gsub('[12].*', '', rownames(counts)))
  # res = dirichlet_regression(counts, covariates, counts ~ condition)

  # Calculate regression
  counts = as.data.frame(counts)
  counts$counts = DR_data(counts)
  data = cbind(counts, covariates)
  fit = DirichReg(counts ~ condition, data)

  # Get p-values
  u = summary(fit)
  pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
  v = names(pvals)
  pvals = matrix(pvals, ncol=length(u$varnames))
  rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
  colnames(pvals) = u$varnames
  fit$pvals = pvals

  fit
}
