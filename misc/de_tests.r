library(parallel)
options(mc.cores = 12)
library(scFlow)
library(magrittr)
sce <- read_sce("~/Documents/junk/enriched/final_sce/")
sce$seqdate <- as.factor(sce$seqdate)
sce_all <- sce
sce_subset <- sce[, sce$cluster_celltype == "Micro" & sce$brain_region == "SSC"] #EC

colnames(SummarizedExperiment::colData(sce))
plot_violin(sce_subset, group_var = "individual", subset_group = "Micro", gene = "DPYD", label_angle = 90)

?plot_violin

# working
fargs <-  list(
  sce = sce_subset,
  de_method = "MASTZLM",
  min_counts = 1,
  min_cells_pc = 0.10,
  rescale_numerics = TRUE,
  dependent_var = "group",
  ref_class = "Control",
  confounding_vars = c(#"individual",
    "cngeneson",
    "sex",
    "age",
    "PMI",
    "RIN",
    "seqdate",
    "pc_mito"),
  random_effects_var = NULL,
  fc_threshold = 1.1,
  pval_cutoff = 0.05,
  ensembl_mapping_file = "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv"
)

# working experimenting
fargs <-  list(
  sce = sce_subset,
  de_method = "MASTZLM",
  min_counts = 1,
  min_cells_pc = 0.10,
  rescale_numerics = TRUE,
  dependent_var = "group",
  ref_class = "Control",
  confounding_vars = c(#"individual",
    "cngeneson",
    "sex",
    "age",
    "PMI",
    "RIN",
    "seqdate",
    "pc_mito"),
  random_effects_var = NULL,
  fc_threshold = 1.05,
  pval_cutoff = 0.1,
  ensembl_mapping_file = "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv"
)


# working experimenting 2
fargs <-  list(
  sce = sce_subset,
  de_method = "MASTZLM",
  min_counts = 1,
  min_cells_pc = 0.10,
  rescale_numerics = TRUE,
  dependent_var = "group",
  ref_class = "Control",
  confounding_vars = c(#"individual",
    "cngeneson",
    "individual",
    #"sex",
    #"age",
    #"PMI",
    #"RIN",
    #"seqdate",
    "pc_mito"),
  random_effects_var = NULL,
  fc_threshold = 1.05,
  pval_cutoff = 0.1,
  ensembl_mapping_file = "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv"
)

# glmer
fargs <-  list(
  sce = sce_subset,
  de_method = "MASTZLM",
  ebayes = FALSE,
  mast_method = "glmer",
  min_counts = 1,
  min_cells_pc = 0.10,
  rescale_numerics = TRUE,
  dependent_var = "group",
  ref_class = "Control",
  confounding_vars = c(
    "cngeneson",
    "sex",
    "age",
    "PMI",
    "RIN",
    "seqdate",
    "pc_mito"
    ),
  random_effects_var = "individual",
  fc_threshold = 1.05,
  pval_cutoff = 0.1,
  ensembl_mapping_file = "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv"
)
de_results <- do.call(perform_de, fargs)


#sce_pp <- do.call(scFlow:::.preprocess_sce_for_de, fargs)
#sce <- sce_pp





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
perform_de_test <- function(sce,
                       de_method = "MASTZLM",
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
                       fc_threshold = 1.1,
                       pval_cutoff = 0.05,
                       ...) {
  fargs <- c(as.list(environment()), list(...))

  cli::cli_h1("Starting Differential Gene Expression")

  cli::cli_h2("Pre-processing SingleCellExperiment")

  # preprocess
  fargs$sce <- do.call(scFlow:::.preprocess_sce_for_de, fargs)

  # select method
  if (de_method == "MASTZLM") {
    cli::cli_h2("MAST Zero-inflated Linear Model")
    de_fn <- .perform_de_with_mast_test
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
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom scater librarySizeFactors normalize logNormCounts
#' @importFrom Matrix colSums
#' @importFrom cli cli_text
#'
#' @keywords internal
.preprocess_sce_for_de_test <- function(...) {
  fargs <- list(...)
  sce <- fargs$sce

  sce <- scFlow:::.filter_sce_genes_for_de(sce,
                                  min_counts = fargs$min_counts,
                                  min_cells_pc = fargs$min_cells_pc
  )

  #sce@int_colData$size_factor <- scater::librarySizeFactors(sce)
  sce <- scater::logNormCounts(sce)

  if (fargs$rescale_numerics == TRUE) {
    # recale numerics avoids issues with glmer
    nums <- unlist(lapply(SummarizedExperiment::colData(sce), is.numeric))
    SummarizedExperiment::colData(sce)[, nums] <- scale(
      as.data.frame(SummarizedExperiment::colData(sce)[, nums])
    )
  }

  # from MAST tutorial
  cdr2 <- Matrix::colSums(SingleCellExperiment::counts(sce) > 0)
  sce$cngeneson <- as.numeric(scale(cdr2))

  cli::cli_text(c(
    "Setting '{.emph {fargs$ref_class}}' ",
    "as the reference class of '{.strong {fargs$dependent_var}}'."
  ))

  # define the reference class
  sce[[fargs$dependent_var]] <- relevel(
    sce[[fargs$dependent_var]],
    ref = fargs$ref_class
  )

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
#' @importFrom scater normalize
#' @importFrom methods as
#' @importFrom magrittr %>%
#' @importFrom MAST zlm
#' @importFrom cli cli_text cli_alert_info cli_alert_danger cli_alert_success
#'
#' @keywords internal
.perform_de_with_mast_test <- function(...) {
  fargs <- list(
    mast_method = "bayesglm",
    ebayes = FALSE,
    force_run = FALSE
  )
  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  sce <- fargs$sce
  sca <- as(sce, "SingleCellAssay")

  message("Generating model formula")

  # generate formula with and without prefix
  prefixes <- c("", "sce$")
  mod_formulae <- purrr::map(
    prefixes,
    ~ do.call(
      scFlow:::.generate_model_from_vars,
      list(
        sce = sce,
        dependent_var = fargs$dependent_var,
        confounding_vars = fargs$confounding_vars,
        random_effects_var = fargs$random_effects_var,
        prefix = .
      )
    )
  )

  # main, without prefixes
  model_formula <- mod_formulae[[1]]

  cli::cli_text(c(
    "Model formula: ",
    scFlow:::.formula_to_char(model_formula)
  ))

  # mod_check <- do.call(.check_model, list(model_formula = mod_formulae[[2]]))
  if (is.null(fargs$random_effects_var)) {
    is_full_rank <- scFlow:::.check_model(mod_formulae[[2]])
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

  # options(warn=-1) #temporary silencing
  x <- Sys.time()
  zlmCond <- MAST::zlm(
    formula = model_formula,
    sca = sca,
    method = fargs$mast_method, # note: glmer requires a random effects var
    ebayes = fargs$ebayes,
    parallel = TRUE
  )
  message(Sys.time() - x)

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

  results_l <- list()

  for (ctrast in contrasts) {
    print(ctrast)
    message(sprintf("Running DE for %s\n", ctrast))
    gc()
    print("!!!!!!!did gc")
    # only test the condition coefficient.
    summaryCond <- MAST::summary(zlmCond, doLRT = ctrast)
    print("!!!!!!!did summaryCond")
    # create data table of results
    print(summaryCond)
    dim(summaryCond)
    #summaryDt <- summaryCond$datatable
    summaryDt <- summaryCond[["datatable"]]
    print("!!!!!!!did summaryCond")

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
    print("!!!!!!!did hurdles")
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
    print("!!!!!!!did stuff")
    assertthat::assert_that(
      !all(is.na(fcHurdle$pval)),
           msg = "Singularities prevented pval calculations. Revise model.")

    # append gene names
    ensembl_res <- map_ensembl_gene_id(
      fcHurdle$ensembl_gene_id,
      mappings = c("external_gene_name", "gene_biotype"),
      ensembl_mapping_file = fargs$ensembl_mapping_file
    ) %>%
      dplyr::rename(gene = external_gene_name)

    print("!!!!!!!finished ensem")
    fcHurdle <- merge(fcHurdle, ensembl_res, by = "ensembl_gene_id")

    fcHurdle$padj <- p.adjust(
      fcHurdle$pval, "fdr", n = dim(sca)[[1]])
    table(fcHurdle$padj <= 0.05)

    fcHurdle$contrast <- ctrast
    fcHurdle$reference <- fargs$ref_class

    print("!!!!!!!did fargs stuff")
    model_formula_string <- scFlow:::.formula_to_char(model_formula)
    fcHurdle$model <- gsub(" ", "", model_formula_string,
                           fixed = TRUE) # no whitespace


    p <- scFlow:::.volcano_plot(
      dt = fcHurdle,
      fc_threshold = fargs$fc_threshold,
      pval_cutoff = fargs$pval_cutoff,
      n_label = fargs$n_label
    )
    print("!!!!!!!did volcano")
    results <- fcHurdle %>%
      dplyr::filter(padj <= fargs$pval_cutoff) %>%
      dplyr::filter(gene_biotype == "protein_coding") %>%
      dplyr::filter(abs(logFC) > log2(fargs$fc_threshold)) %>%
      dplyr::arrange(-logFC)

    DGEs <- c(sum(results$logFC > 0), sum(results$logFC < 0))
    names(DGEs) <- c("Up", "Down")

    element_name <- paste(fargs$ref_class, ctrast, sep = "_vs_")

    de_params <- list(
      celltype = unique(fargs$sce$cluster_celltype),
      de_method = fargs$de_method,
      mast_method = fargs$mast_method,
      pseudobulk = fargs$sce@metadata$scflow_steps$pseudobulk,
      min_counts = fargs$min_counts,
      min_cells_pc = fargs$min_cells_pc,
      rescale_numerics = fargs$rescale_numerics,
      dependent_var = fargs$dependent_var,
      ref_class = fargs$ref_class,
      confounding_vars = fargs$confounding_vars,
      random_effects_var = fargs$random_effects_var,
      fc_threshold = fargs$fc_threshold,
      pval_cutoff = fargs$pval_cutoff,
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

    attr(results, "de_parameters") <- de_params

    attr(results, "de_result") <- DGEs

    attr(results, "plot") <- p

    results_l[[element_name]] <- results
  }

  message("Done!  Returning results")

  return(results_l)
}
