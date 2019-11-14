perform_de <- function(sce,
                       de_method = "MASTZLM",
                       min_counts = 1,
                       min_cells_pc = 0.10,
                       rescale_numerics = TRUE,
                       dependent_var = "group",
                       ref_class = "Control",
                       confounding_vars = c("individual", "cngeneson", "sex", "age", "PMI", "RIN", "seqdate", "pc_mito"),
                       random_effects_var = NULL,
                       fc_threshold = 1.1,
                       pval_cutoff = 0.05
                       ) {

  sce <- .preprocess_sce_for_de(
    sce = sce,
    min_counts = min_counts,
    min_cells_pc = min_cells_pc,
    rescale_numerics = rescale_numerics,
    dependent_var = dependent_var,
    ref_class = ref_class
  )

  if(de_method == "MASTZLM") {
    de_fn <- .perform_de_with_mast
    de_results <- do.call(
      de_fn,
      list(
        sce = sce,
        fc_threshold = fc_threshold,
        pval_cutoff
      )
    )
  }
  dependent_var_ref_class
}




################################################################################
#' Preprocess the SCE for DE analysis
#'
#' Subset genes, set factor reference level, and optionally pseudobulk
#'
#' @param ...
#'
#' @return sce a SingleCellExperiment object
#'
#' @family differential gene expression
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom scater normalize
#'
#' @keywords internal
.preprocess_sce_for_de <- function(...) {

  args <- list(...)
  sce <- args$sce

  sce <- .filter_sce_genes_for_de(sce,
                                  min_counts = args$min_counts,
                                  min_cells_pc = args$min_cells_pc
  )

  sce <- scater::normalize(sce)

  if (args$rescale_numerics == TRUE) {
    # recale numerics avoids issues with glmer
    nums <- unlist(lapply(SummarizedExperiment::colData(sce), is.numeric))
    SummarizedExperiment::colData(sce)[, nums] <- scale(
      as.data.frame(SummarizedExperiment::colData(sce)[, nums])
    )
  }

  # from MAST tutorial
  cdr2 <- colSums(SingleCellExperiment::counts(sce) > 0)
  SummarizedExperiment::colData(sce)$cngeneson <- as.numeric(scale(cdr2))

  # define the reference class
  sce[[args$dependent_var]] <- relevel(
    sce[[args$dependent_var]],
    ref = args$ref_class
  )

  return(sce)

}


################################################################################
#' Perform DE with MAST
#'
#'
#' @param sce a SingleCellExperiment object
#'
#' @return sce a SingleCellExperiment object
#'
#' @family differential gene expression
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom scater normalize
#'
#' @keywords internal
.perform_de_with_mast <- function(...) {

  message("##### Starting MAST Zero-Inflated Linear Mixed Model DE")



  sca <- as(sce, "SingleCellAssay")

  # setup model variables
  de_vars_l <- list()
  de_vars_l[["model"]] <- list()
  # de model variables
  de_vars_l$model$dependent_var <- "group"
  de_vars_l$model$confounding_vars <- c(
    "sex",
    "age", "pc_mito", "PMI", "total_features_by_counts", "seqdate"
  )

  de_vars_l$model$random_effects_vars <- "individual"
  # within dependent group reference variable
  de_vars_l$dep_var_ref_class <- "Control"

  # add prefixes for estimable vars
  de_vars_l$model_prefixed <- map(
    de_vars_l$model, function(vars) {
      .add_prefix_to_model_vars(vars, "sce$")
    }
  )

  sce$group <- relevel(sce$group, ref = de_vars_l$dep_var_ref_class)

  model_formula <- as.formula("~ group + individual + cngeneson + sex + age + PMI + RIN + seqdate + pc_mito")  #for bayesglm

  # fit model
  message("Fitting model\n")
  #options(warn=-1) #temporary silencing
  x <- Sys.time()
  zlmCond <- zlm(
    formula = model_formula,
    sca = sca,
    method = "bayesglm", # note: glmer requires a random effects var
    ebayes = FALSE,
    parallel = TRUE
  )

  gc()

  ## test each contrast separately
  # obtain the group names for non-controls
  dependent_var_names <- unique(sce[[de_vars_l$model$dependent_var]])

  dependent_var_names <- droplevels(
    dependent_var_names[dependent_var_names != de_vars_l$dep_var_ref_class]
  )

  contrasts <- map_chr(dependent_var_names, ~ paste0(de_vars_l$model$dependent_var, .))

  # contrasts <- .add_prefix_to_model_vars(contrasts, "sce$")

  results_l <- list()

  for (ctrast in contrasts) {
    message(sprintf("Running DE for %s\n", ctrast))

    # only test the condition coefficient.
    summaryCond <- summary(zlmCond, doLRT = ctrast)

    # create data table of results
    summaryDt <- summaryCond$datatable

    # primerid is numbered for H, so retrieve rowname (ensembl)
    row_num2name <- seq_along(rownames(counts(sca)))
    names(row_num2name) <- rownames(counts(sca))

    hurdle_pvals <- summaryDt %>%
      dplyr::filter(contrast == ctrast & component == "H") %>%
      dplyr::select(primerid, `Pr(>Chisq)`)

    hurdle_pvals$primerid <- map_chr(
      hurdle_pvals$primerid,
      function(x) names(row_num2name[as.numeric(x)])
    )

    logFCcoefs <- summaryDt %>%
      dplyr::filter(contrast == ctrast & component == "logFC") %>%
      dplyr::select(primerid, coef, ci.hi, ci.lo)

    # merge
    fcHurdle <- merge(hurdle_pvals, # hurdle P values
                      logFCcoefs,
                      by = "primerid"
    ) # logFC coefficients

    fcHurdle <- fcHurdle %>% dplyr::rename(
      ensembl_gene_id = primerid,
      pval = `Pr(>Chisq)`,
      logFC = coef
    )

    # append gene names
    ensembl_res <- EnsemblGeneVersionVectorToGeneNameVector(
      fcHurdle$ensembl_gene_id,
      mappings_file = ensembl_mappings
    ) %>%
      dplyr::rename(gene = external_gene_name)

    fcHurdle <- merge(fcHurdle, ensembl_res, by = "ensembl_gene_id")

    fcHurdle$padj <- p.adjust(fcHurdle$pval, "fdr", n = dim(sce)[[1]])

    fcHurdle$contrast <- ctrast

    model_formula_string <- as.character(Reduce(paste, deparse(model_formula)))
    fcHurdle$model <- gsub(" ", "", model_formula_string, fixed = TRUE) # no whitespace

    results_l[[ctrast]] <- fcHurdle
  }

  message("Merging results")
  results <- Reduce(rbind, results_l)

  results <- results %>%
    dplyr::filter(padj <= pval_cutoff) %>%
    dplyr::filter(gene_biotype == "protein_coding") %>%
    dplyr::filter(abs(logFC) > log2(fc_threshold)) %>%
    dplyr::arrange(-logFC)

  message("Done!  Returning results")

  return(results)
}


################################################################################
#' Subset a SingleCellExperiment for genes meeting expressivity criteria
#'
#' This function is applied after the subsetting of a SCE (e.g. celltype).
#' It's useful to subset before differential expression analysis.
#'
#' @param sce a SingleCellExperiment object
#'
#' @return sce a SingleCellExperiment object
#'
#' @family differential gene expression
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom cli cli_text
#'
#' @keywords internal
.filter_sce_genes_for_de <- function(sce,
                                     min_counts = 1L,
                                     min_cells_pc = 0.1) {

  n_genes_before <- dim(sce)[[1]]

  min_cells <- floor(dim(sce)[[2]] * min_cells_pc)

  keep_genes <- apply(
    SingleCellExperiment::counts(sce), 1,
    function(x) {length(x[x >= min_counts]) > min_cells}
  )

  sce <- sce[keep_genes, ]

  n_genes_after <- dim(sce)[[1]]

  cli::cli_text(c(
    "Selected {n_genes_after} from {n_genes_before} genes ",
    "with >{min_counts} count(s) in >{min_cells_pc * 100}% of cells."
  ))

  return(sce)

}
