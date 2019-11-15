################################################################################
#' Perform Differential Gene Expression on a SingleCellExperiment
#'
#' @param sce a SingleCellExperiment object
#' @param min_counts minimum number of counts
#' @param min_cells_pc percentage of cells with min_counts for gene selection
#' @param rescale_numerics rescaling numerics may improve model
#' @param dependent_var the name of the colData variable for contrasts
#' @param ref_class the class of dependent_var used as reference
#' @param confounding_vars the independent variables of the model
#' @param random_effects_var variable(s) to model as random effects
#' @param fc_threshold fold change up/down threshold
#' @param pval_cutoff the adjusted pvalue cutoff threshold
#'
#' @return results_l a list of DE table results
#'
#' @family differential gene expression
#'
#' @importFrom cli cli_text
#'
#' @export
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

  pp_sce <- .preprocess_sce_for_de(
    sce = sce,
    min_counts = min_counts,
    min_cells_pc = min_cells_pc,
    rescale_numerics = rescale_numerics,
    dependent_var = dependent_var,
    ref_class = ref_class
  )

  print(sce$cngeneson)

  if(de_method == "MASTZLM") {
    de_fn <- .perform_de_with_mast
    de_results <- do.call(
      de_fn,
      list(
        sce = pp_sce,
        fc_threshold = fc_threshold,
        pval_cutoff = pval_cutoff,
        dependent_var = dependent_var,
        ref_class = ref_class,
        confounding_vars = confounding_vars,
        random_effects_var = random_effects_var,
        fc_threshold = fc_threshold,
        pval_cutoff = pval_cutoff
      )
    )
  }

  return(de_results)
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
#' @importFrom cli cli_text
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
  sce$cngeneson <- as.numeric(scale(cdr2))

  cli::cli_text(c(
    "Setting '{.emph {args$ref_class}}' ",
    "as the reference class of '{.strong {args$dependent_var}}'.")
  )

  # define the reference class
  sce[[args$dependent_var]] <- relevel(
    sce[[args$dependent_var]],
    ref = args$ref_class
  )

  return(sce)

}


################################################################################
#' Perform DE with MAST#'
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
#' @importFrom methods as
#' @importFrom magrittr %>%
#' @import MAST
#'
#' @keywords internal
.perform_de_with_mast <- function(...) {

  args <- list(...)

  message("##### Starting MAST Zero-Inflated Linear Mixed Model DE")

  sca <- as(sce, "SingleCellAssay")

  message("Generating model formula")
  model_formula <- do.call(.generate_model_from_vars, args)

  # fit model
  message("Fitting model\n")
  print(sce$cngeneson)
  #options(warn=-1) #temporary silencing
  x <- Sys.time()
  zlmCond <- MAST::zlm(
    formula = model_formula,
    sca = sca,
    method = "bayesglm", # note: glmer requires a random effects var
    ebayes = FALSE,
    parallel = TRUE
  )
  message(Sys.time() - x)

  gc()

  ## test each contrast separately
  # obtain the group names for non-controls
  dependent_var_names <- unique(sce[[args$dependent_var]])

  dependent_var_names <- droplevels(
    dependent_var_names[dependent_var_names != args$ref_class]
  )

  contrasts <- purrr::map_chr(
    dependent_var_names,
    ~ paste0(args$dependent_var, .)
    )

  results_l <- list()

  for (ctrast in contrasts) {
    message(sprintf("Running DE for %s\n", ctrast))

    # only test the condition coefficient.
    summaryCond <- suppressWarnings(summary(zlmCond, doLRT = ctrast))

    # create data table of results
    summaryDt <- summaryCond$datatable

    # primerid is numbered for H, so retrieve rowname (ensembl)
    row_num2name <- seq_along(rownames(SingleCellExperiment::counts(sca)))
    names(row_num2name) <- rownames(SingleCellExperiment::counts(sca))

    hurdle_pvals <- summaryDt %>%
      dplyr::filter(contrast == ctrast & component == "H") %>%
      dplyr::select(primerid, `Pr(>Chisq)`)

    hurdle_pvals$primerid <- purrr::map_chr(
      hurdle_pvals$primerid,
      ~ names(row_num2name[as.numeric(.)])
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
    ensembl_res <- map_ensembl_gene_id(fcHurdle$ensembl_gene_id) %>%
      dplyr::rename(gene = external_gene_name)

    fcHurdle <- merge(fcHurdle, ensembl_res, by = "ensembl_gene_id")

    fcHurdle$padj <- p.adjust(fcHurdle$pval, "fdr", n = dim(sce)[[1]])

    fcHurdle$contrast <- ctrast
    fcHurdle$reference <- args$ref_class

    model_formula_string <- as.character(Reduce(paste, deparse(model_formula)))
    fcHurdle$model <- gsub(" ", "", model_formula_string, fixed = TRUE) # no whitespace

    results <- fcHurdle %>%
      dplyr::filter(padj <= args$pval_cutoff) %>%
      dplyr::filter(gene_biotype == "protein_coding") %>%
      dplyr::filter(abs(logFC) > log2(args$fc_threshold)) %>%
      dplyr::arrange(-logFC)

    element_name <- paste(args$ref_class, ctrast, sep = "_vs_")
    results_l[[element_name]] <- results
  }

  merged_name <- paste0(args$ref_class, "_Merged_Results")
  message("Merging results")
  results_l[[merged_name]] <- Reduce(rbind, results_l)

  message("Done!  Returning results")

  return(results_l)
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



################################################################################
#' Generate a model formula from dependent, confounding, and random vars
#'
#' @param dependent_var the dependent variable
#' @param confounding_vars the confounding variables
#' @param re_vars the random effect variables (optional)
#'
#' @return model_formula a model formula
#'
#' @family differential gene expression
#'
#' @importFrom purrr map_chr
#'
#' @keywords internal
.generate_model_from_vars <- function(...) {

  args <- list(...)

  if (!is.null(args$random_effects_var)) {
    message("Generating random effects model")
    args$random_effects_var <- purrr::map_chr(
      args$random_effects_var,
      function(var) sprintf("(1|%s)", var)
    )
    args$random_effects_var <- paste(args$random_effects_var, collapse = " + ")

    model_formula <- as.formula(
      sprintf(
        "~ %s + %s + %s",
        args$dependent_var,
        args$random_effects_var,
        paste(args$confounding_vars, collapse = " + ")
      )
    )
  } else {
    message("Generating mixed effects model")
    model_formula <- as.formula(
      sprintf(
        "~ %s + %s",
        args$dependent_var,
        paste(args$confounding_vars, collapse = " + ")
      )
    )
  }

  return(model_formula)

}
