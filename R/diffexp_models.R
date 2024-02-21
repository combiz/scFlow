################################################################################
#' Perform Differential Gene Expression on a SingleCellExperiment
#'
#' @param sce a SingleCellExperiment object
#' @param de_method The differential gene expression method.
#' @param mast_method If `de_method` is "MASTZLM" then mast_method should be
#' provided. Possible values are "glm", "glmer", "bayesglm". Default is "glm".
#' For "glmer" and "random_effects_var" should be provided.
#' @param min_counts minimum number of counts
#' @param min_cells_pc percentage of cells with min_counts for gene selection
#' @param rescale_numerics rescaling numerics may improve model
#' @param dependent_var the name of the colData variable for contrasts
#' @param ref_class the class of dependent_var used as reference
#' @param confounding_vars the independent variables of the model
#' @param random_effects_var variable(s) to model as random effects
#' @param interaction_vars two or more variables to model as interacting
#' @param unique_id_var the colData variable identifying unique samples
#' @param species human or mouse
#' @param parallel enable parallel processing
#' @param ... advanced options
#'
#' @return results_l a list of DE table results
#'
#' @family differential gene expression
#'
#' @importFrom cli cli_text rule cli_h1 cli_h2 cli_h3
#'
#'
#' @export
perform_de <- function(sce,
                       de_method = "MASTZLM",
                       mast_method = "glm",
                       min_counts = 1,
                       min_cells_pc = 0.10,
                       rescale_numerics = TRUE,
                       dependent_var = "group",
                       ref_class = "Control",
                       confounding_vars = c("individual",
                                            "cngeneson",
                                            "sex",
                                            "age",
                                            "PMI",
                                            "RIN",
                                            "seqdate",
                                            "pc_mito"),
                       random_effects_var = NULL,
                       interaction_vars = NULL,
                       unique_id_var = "individual",
                       species = getOption(
                         "scflow_species",
                         default = "human"),
                       parallel = TRUE,
                       ...) {
  fargs <- c(as.list(environment()), list(...))

  cli::cli_h1("Starting Differential Gene Expression")

  cli::cli_h2("Pre-processing SingleCellExperiment")

  # preprocess
  fargs$sce <- do.call(.preprocess_sce_for_de, fargs)

  assertthat::assert_that(
    dim(fargs$sce)[1] >= 100,
    msg = "Less than 500 genes passed the min_counts and min_cells_pc threshold;
    Consider decreasing both values!"
  )

  # select method
  if (de_method == "MASTZLM") {
    cli::cli_h2("MAST Zero-inflated Linear Model")
    de_fn <- .perform_de_with_mast
  } else {
    stop("Invalid method specified.")
  }

  # perform differential expression
  de_results <- do.call(de_fn, fargs)

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
#' @rawNamespace import(scater, except = "normalize")
#' @importFrom SingleCellExperiment counts tpm
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom Matrix colSums
#' @importFrom cli cli_text
#' @importFrom sctransform vst
#' @importFrom limma normalizeQuantiles
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom stats relevel
#'
#' @keywords internal
.preprocess_sce_for_de <- function(...) {
  fargs <- list(
    quantile_norm = FALSE,
    sctransform = FALSE,
    rescale_numerics = TRUE,
    pseudobulk = FALSE,
    subset_var = NULL, # for variance calculation, subset on
    subset_class = NULL # for variance calculation, subset class
  )
  inargs <- list(...)
  fargs[names(inargs)] <- inargs
  sce <- fargs$sce

  # run this before subset to avoid non-conformable array error
  if (fargs$sctransform) {
    cli::cli_alert("Normalizing with sctransform")
    #sce$log10_total_counts <- log10(sce$total_counts)
    vst_out <- sctransform::vst(
      umi = as(SingleCellExperiment::counts(sce), "CsparseMatrix"),
      cell_attr = as.data.frame(SummarizedExperiment::colData(sce)),
      latent_var = c("log_umi"),#c("log10_total_counts"),
      batch_var = NULL,
      latent_var_nonreg = NULL,
      n_genes = 2000,
      n_cells = NULL,
      return_gene_attr = TRUE,
      return_cell_attr = TRUE,
      method = "poisson",
      do_regularize = TRUE,
      residual_type = "pearson",
      show_progress = TRUE
    )
    mat <- vst_out$y
    ridx <- rownames(sce) %in% rownames(mat)
    sce <- sce[ridx, ]
    SingleCellExperiment::normcounts(sce) <- mat
  }

  if (fargs$quantile_norm) {
    cli::cli_alert("Quantile normalizing")
    #mat <- as(SingleCellExperiment::counts(sce), "dgCMatrix")
    new_mat <- matrix(0, nrow = dim(mat)[1], ncol = dim(mat)[2])
    #new_mat <- Matrix::Matrix(0, nrow = dim(mat)[1], ncol = dim(mat)[2], sparse = TRUE)
    for(unique_id in unique(sce[[fargs$unique_id_var]])) {
      print(sprintf("normalizing sample %s", unique_id))
      idx <- sce[[fargs$unique_id_var]] == unique_id
      submat <- as.matrix(mat[, idx])
      qn_submat <- preprocessCore::normalize.quantiles(submat)
      new_mat[, idx] <- qn_submat
    }
    mat <- new_mat
    SingleCellExperiment::normcounts(sce) <- mat
  }

  sce <- .filter_sce_genes_for_de(sce,
    min_counts = fargs$min_counts,
    min_cells_pc = fargs$min_cells_pc
  )


  if (fargs$rescale_numerics == TRUE) {
    cli::cli_alert("Rescaling numeric variables")
    # recale numerics avoids issues with glmer
    nums <- unlist(lapply(SummarizedExperiment::colData(sce), is.numeric))
    SummarizedExperiment::colData(sce)[, nums] <- scale(
      as.data.frame(SummarizedExperiment::colData(sce)[, nums])
    )
  }

  # from MAST tutorial
  cdr2 <- Matrix::colSums(SingleCellExperiment::counts(sce) > 0)
  sce$cngeneson <- as.numeric(scale(cdr2))

  #calculate genes with high inter-individual variance
  variable_genes <- .get_variance_explained(
    sce,
    variable = fargs$unique_id_var,
    subset_var = fargs$subset_var,
    subset_class = fargs$subset_class)
  sce@metadata$variable_genes <- variable_genes


  if (fargs$pseudobulk) {
    cli::cli_alert("Pseudobulking")
    sce <- pseudobulk_sce(
      sce,
      keep_vars = unique(c(
        fargs$dependent_var,
        fargs$confounding_vars,
        fargs$random_effects_var,
        fargs$interaction_vars)),
      sample_var = fargs$unique_id_var
      )
    sce@metadata$variable_genes <- variable_genes

    if(fargs$quantile_norm) {
      cli::cli_alert("Quantile normalizing merged")
      mat <- SingleCellExperiment::normcounts(sce)
      mat <- preprocessCore::normalize.quantiles(mat)
      SingleCellExperiment::normcounts(sce) <- mat
    }
  }

  cli::cli_text(c(
    "Setting '{.emph {fargs$ref_class}}' ",
    "as the reference class of '{.strong {fargs$dependent_var}}'."
  ))

  # define the reference class
  if (!is.numeric(sce[[fargs$dependent_var]])) {
    sce[[fargs$dependent_var]] <- stats::relevel(
      as.factor(sce[[fargs$dependent_var]]),
      ref = fargs$ref_class
    )
  }

  return(sce)
}


################################################################################
#' Perform DE with MAST
#'
#' @param sce a SingleCellExperiment object
#'
#' @return sce a SingleCellExperiment object
#'
#' @family differential gene expression
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom methods as
#' @importFrom magrittr %>%
#' @importFrom MAST zlm summary SceToSingleCellAssay
#' @importFrom cli cli_text cli_alert_info cli_alert_danger cli_alert_success
#'
#' @keywords internal
.perform_de_with_mast <- function(...) {
  fargs <- list(
    mast_method = "bayesglm",
    ebayes = FALSE,
    force_run = FALSE,
    nAGQ = 0,
    parallel = TRUE,
    append_info = NULL # attached to results info column
  )
  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  sce <- fargs$sce
  SingleCellExperiment::tpm(sce) <- scater::calculateTPM(sce, exprs_values = "counts")
  SingleCellExperiment::logcounts(sce) <- log2(SingleCellExperiment::tpm(sce) + 1)

  suppressMessages(sca <- MAST::SceToSingleCellAssay(sce))

  message("Generating model formula")

  # generate formula with and without prefix
  prefixes <- c("", "sce$")
  mod_formulae <- purrr::map(
    prefixes,
    ~ do.call(
      .generate_model_from_vars,
      list(
        sce = sce,
        dependent_var = fargs$dependent_var,
        confounding_vars = fargs$confounding_vars,
        random_effects_var = fargs$random_effects_var,
        interaction_vars = fargs$interaction_vars,
        prefix = .
      )
    )
  )

  # main, without prefixes
  model_formula <- mod_formulae[[1]]

  cli::cli_text(c(
    "Model formula: ",
    .formula_to_char(model_formula)
  ))

  # mod_check <- do.call(.check_model, list(model_formula = mod_formulae[[2]]))
  if (is.null(fargs$random_effects_var)) {
    is_full_rank <- .check_model(mod_formulae[[2]])
    if (!fargs$force_run) {
      assertthat::assert_that(
        is_full_rank,
        msg = "A full rank model specification is required."
      )
    } else {
      cli::cli_alert_info("Forcing run, ignoring model full rank.")
    }
  }

  # fit model
  message("Fitting model\n")

  if (fargs$mast_method == "glmer") {
    fit_args_D <- list(nAGQ = fargs$nAGQ)
  } else {
    fit_args_D <- list()
  }

  # options(warn=-1) #temporary silencing
  x <- Sys.time()
  zlmCond <- MAST::zlm(
    formula = model_formula,
    sca = sca,
    exprs_value = 'logcounts',
    method = fargs$mast_method, # note: glmer requires a random effects var
    ebayes = fargs$ebayes,
    parallel = fargs$parallel,
    fitArgsD = fit_args_D
  )
  message(Sys.time() - x)
  zlmCond@hookOut <- NULL

  if(is.numeric(sce[[fargs$dependent_var]])) {
    contrasts <- fargs$dependent_var
  } else {
    ## test each contrast separately
    # obtain the group names for non-controls
    dependent_var_names <- unique(sce[[fargs$dependent_var]])

    dependent_var_names <- droplevels(
      dependent_var_names[dependent_var_names != fargs$ref_class]
    )

    contrasts <- purrr::map_chr(
      dependent_var_names,
      ~ paste0(fargs$dependent_var, .)
    )
  }

  results_l <- list()

  for (ctrast in contrasts) {
    message(sprintf("Running DE for %s\n", ctrast))
    gc()

    # only test the condition coefficient.
    # suppress warnings is a temporary fix for missing namespace in mast
    summaryCond <- suppressWarnings(MAST::summary(zlmCond, doLRT = ctrast))

    # create data table of results
    summaryDt <- summaryCond$datatable

    hurdle_pvals <- summaryDt %>%
      dplyr::filter(contrast == ctrast & component == "H") %>%
      dplyr::select(primerid, `Pr(>Chisq)`)


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

    assertthat::assert_that(
      !all(is.na(fcHurdle$pval)),
      msg = "Singularities prevented p-value calculations. Revise model."
      )

    # append gene names
    ensembl_res <- map_ensembl_gene_id(
      fcHurdle$ensembl_gene_id,
      mappings = c("external_gene_name", "gene_biotype"),
      ensembl_mapping_file = fargs$ensembl_mapping_file,
      species = fargs$species
    ) %>%
      dplyr::rename(gene = external_gene_name)

    fcHurdle <- merge(fcHurdle, ensembl_res, by = "ensembl_gene_id")

    fcHurdle$padj <- stats::p.adjust(
      fcHurdle$pval, "fdr", n = dim(sca)[[1]])

    fcHurdle$contrast <- ctrast
    fcHurdle$reference <- fargs$ref_class
    fcHurdle$FCRO <- order(abs(2^fcHurdle$logFC))

    model_formula_string <- .formula_to_char(model_formula)
    fcHurdle$model <- gsub(" ", "", model_formula_string,
                           fixed = TRUE) # no whitespace

    if(!is.null(fargs$append_info)){ fcHurdle$info <- fargs$append_info }

    results <- fcHurdle %>%
      dplyr::arrange(padj)

    if (!is.null(sce@metadata$variable_genes)) {
    results <- dplyr::left_join(
      results,
      sce@metadata$variable_genes,
      by = "ensembl_gene_id")
    }

    element_name <- paste(ctrast, fargs$ref_class, sep = "_vs_")

    de_params <- list(
      celltype = unique(fargs$sce$cluster_celltype),
      de_method = fargs$de_method,
      mast_method = fargs$mast_method,
      pseudobulk = ifelse(isTRUE(fargs$sce@metadata$scflow_steps$pseudobulk), "Yes", "No"),
      min_counts = fargs$min_counts,
      min_cells_pc = fargs$min_cells_pc,
      rescale_numerics = fargs$rescale_numerics,
      dependent_var = fargs$dependent_var,
      ref_class = fargs$ref_class,
      confounding_vars = fargs$confounding_vars,
      random_effects_var = fargs$random_effects_var,
      interaction_vars = fargs$interaction_vars,
      cells_per_group = table(
        as.data.frame(SingleCellExperiment::colData(fargs$sce))[[fargs$dependent_var]]),
      n_genes = dim(sce)[[1]],
      model = gsub(" ", "", model_formula_string, fixed = TRUE),
      model_full_rank = ifelse(is.null(fargs$random_effects_var), is_full_rank, NA),
      contrast_name = element_name
    )

    if (fargs$de_method != "MASTZLM") {
      de_params$mast_method <- NULL
    }

    de_params <- unlist(de_params)

    attr(results, "de_params") <- de_params

    results_l[[element_name]] <- results
  }

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
    function(x) {
      length(x[x >= min_counts]) > min_cells
    }
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
#' @param sce the data
#' @param dependent_var the dependent variable
#' @param confounding_vars the confounding variables
#' @param random_effects_var the random effect variables (optional)
#' @param interaction_vars variables which interact (optional)
#' @param prefix a variable prefix (e.g. "sce$") (optional)
#'
#' @return model_formula a model formula
#'
#' @family differential gene expression
#'
#' @importFrom purrr map_chr
#' @importFrom stats as.formula p.adjust
#'
#' @keywords internal
.generate_model_from_vars <- function(sce,
                                      dependent_var,
                                      confounding_vars,
                                      random_effects_var = NULL,
                                      interaction_vars = NULL,
                                      prefix = NULL) {
  if (!is.null(prefix)) {
    dependent_var <- purrr::map_chr(
      dependent_var, ~ paste0(prefix, .)
    )
    if (!is.null(confounding_vars)) {
      confounding_vars <- purrr::map_chr(
        confounding_vars, ~ paste0(prefix, .)
      )
    }
    if (!is.null(random_effects_var)) {
      random_effects_var <- purrr::map_chr(
        random_effects_var, ~ paste0(prefix, .)
      )
    }
    if (!is.null(interaction_vars)) {
      interaction_vars <- purrr::map_chr(
        interaction_vars, ~ paste0(prefix, .)
      )
    }
  }

  if (!is.null(random_effects_var)) {
    random_effects_var <- purrr::map_chr(
      random_effects_var, ~ sprintf("(1|%s)", .)
    )
    random_effects_var <- paste(random_effects_var, collapse = " + ")

    model_formula <- stats::as.formula(
      sprintf(
        "~ %s + %s%s%s%s%s",
        dependent_var,
        random_effects_var,
        ifelse(!is.null(confounding_vars), " + ", ""),
        paste(confounding_vars, collapse = " + "),
        ifelse(!is.null(interaction_vars), " + ", ""),
        ifelse(
          !is.null(interaction_vars),
          sprintf("( %s )", paste(interaction_vars, collapse = ":")
        ), "")
      )
    )
  } else {
    model_formula <- stats::as.formula(
      sprintf(
        "~ %s%s%s%s%s",
        dependent_var,
        ifelse(!is.null(confounding_vars), " + ", ""),
        paste(confounding_vars, collapse = " + "),
        ifelse(!is.null(interaction_vars), " + ", ""),
        ifelse(
          is.null(interaction_vars), "",
          sprintf("( %s )", paste(interaction_vars, collapse = ":"))
        )
      )
    )
  }

  return(model_formula)
}

#' check model is full rank
#' @importFrom limma is.fullrank nonEstimable
#' @importFrom cli cli_text cli_ul cli_alert_danger cli_alert_success
#' @importFrom stats model.matrix
#' @importFrom purrr walk
#' @keywords internal
.check_model <- function(model_formula) {
  mod0 <- stats::model.matrix(model_formula)
  if (limma::is.fullrank(mod0)) {
    cli::cli_alert_success("Model formula is full rank")
    return(TRUE)
  } else {
    cli::cli_alert_danger(
      "The model formula is not full rank."
    )
    cli::cli_text("The following coefficients are non-estimable: ")
    purrr::walk(limma::nonEstimable(mod0), ~ cli::cli_ul(.))
  }
  return(FALSE)
}

#' model formula to char/string
#' @keywords internal
.formula_to_char <- function(model_formula) {
  as.character(Reduce(paste, deparse(model_formula)))
}

#' Get the variance explained by a variable for all genes
#' @rawNamespace import(SingleCellExperiment, except = "cpm")
#' @importFrom assertthat assert_that
#' @importFrom cli cli_h2 cli_alert
#' @importFrom SingleCellExperiment normcounts logcounts
#' @importFrom dplyr mutate rename dense_rank desc
#' @importFrom tidyr drop_na
#' @keywords internal
.get_variance_explained <- function(sce,
                                    variable = "individual",
                                    subset_var = NULL,
                                    subset_class = NULL) {

  assertthat::assert_that(length(variable) == 1)
  assertthat::assert_that(class(sce) == "SingleCellExperiment")
  start_time <- Sys.time()
  cli::cli_h2("Calculating variance explained by {.var {variable}}")
  cli::cli_alert("Normalizing counts")
  sce_in <- sce
  if (!is.null(subset_var)) {
    assertthat::assert_that(
      all(subset_var %in% colnames(SummarizedExperiment::colData(sce)))
      )
    assertthat::assert_that(
      subset_class %in% unique(sce[[subset_var]])
    )
    cli::cli_alert(c("Subsetting cells where {.var {subset_var}} ",
                     "is {.var {subset_class}}"))
    sce <- sce[, sce[[subset_var]] == subset_class]
  }
  SingleCellExperiment::normcounts(sce) <-
    scater::normalizeCounts(sce, log = FALSE)
  SingleCellExperiment::logcounts(sce) <-
    log2(SingleCellExperiment::normcounts(sce) + 1)
  cli::cli_alert(c(
    "Calculating variance explained by {.var {variable}} ",
    "for {.val {dim(sce)[[1]]}} genes"
  ))
  #remove unused levels - necessary as each cell type passed to DEG separately
  sce[[variable]] <- droplevels(sce[[variable]])
  vemat <- scater::getVarianceExplained(sce, variables = variable)
  vedf <- as.data.frame.matrix(vemat) %>%
    dplyr::mutate(
      ensembl_gene_id = rownames(.)
    )
  vedf <- tidyr::drop_na(vedf) %>%
    dplyr::rename(variance = 1) %>%
    dplyr::mutate(variance_rank = dplyr::dense_rank(dplyr::desc(variance)))

  end_time <- Sys.time()
  cli::cli_alert(c("Time taken: ", end_time - start_time))

  return(vedf)
}
