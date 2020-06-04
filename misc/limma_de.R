fargs <- list(
  min_counts = 1,
  min_cells_pc = 0.10,
  rescale_numerics = TRUE,
  dependent_var = "sex",
  ref_class = "male",
  confounding_vars = c(#"individual",
    #"cngeneson",
    #"sex"#,
    "age",
    #"PMI",
    #"RIN",
    #"seqdate",
    "pc_mito"),
  #),
  random_effects_var = NULL,
  fc_threshold = 1,
  pval_cutoff = 1,
  ensembl_mapping_file = "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv",
  n_label = 12,
  gene_biotype = "protein_coding",
  pseudobulk = TRUE,
  sctransform = TRUE,
  quantile_norm = TRUE
)

sce[[fargs$dependent_var]] <- relevel(
  sce[[fargs$dependent_var]],
  ref = fargs$ref_class
)

sce <- sce_subset
sce <- do.call(scFlow:::.preprocess_sce_for_de, c(list(sce = sce), fargs))
sce
fargs$sce <- sce
res <- do.call(.perform_de_with_limma, fargs)
sce
sce <- sce_subset

#sce <- sce_subset
################################################################################
#' Perform DE with Limma
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
#' @importFrom MAST zlm summary
#' @importFrom cli cli_text cli_alert_info cli_alert_danger cli_alert_success
#'
#' @keywords internal
.perform_de_with_limma <- function(...) {
  fargs <- list(
    force_run = FALSE,
    n_label = 10
  )
  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  sce <- fargs$sce

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

  # with prefixes
  model_formula <- mod_formulae[[2]]

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

  model_mat <- stats::model.matrix(model_formula)
  # fit model
  message("Fitting model\n")

  # options(warn=-1) #temporary silencing
  x <- Sys.time()
  mat <- as(SingleCellExperiment::counts(sce), "dgCMatrix")

  if (fargs$sctransform) {
    cli::cli_alert("Normalizing with sctransform")
    sce$log10_total_counts <- log10(sce$total_counts)
    vst_out <- sctransform::vst(
      umi = mat,
      cell_attr = as.data.frame(SummarizedExperiment::colData(sce)),
      latent_var = c("log10_total_counts"),
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
    mat <- limma::normalizeQuantiles(mat)
    SingleCellExperiment::normcounts(sce) <- mat
  }
  if (fargs$pseudobulk) {
    sce <- pseudobulk_sce(sce, keep_vars = unique(c(fargs$dependent_var, fargs$confounding_vars, fargs$random_effects_var)))
  }
  fit <- limma::lmFit(mat, model_mat)#, #block = sce$individual,
  #correlation = dupcor$consensus)

  fit <- limma::eBayes(fit)

  message(Sys.time() - x)

  ## test each contrast separately
  # obtain the group names for non-controls
  dependent_var_names <- unique(sce[[fargs$dependent_var]])

  dependent_var_names <- droplevels(
    dependent_var_names[dependent_var_names != fargs$ref_class]
  )

  #contrasts <- setdiff(colnames(fit$design), "(Intercept)")

  contrasts <- purrr::map_chr(
    dependent_var_names,
    ~ paste0("sce$", fargs$dependent_var, .)
  )
  print(contrasts)
  print("colnames:")
  print(colnames(fit$design))

  results_l <- list()

  for (ctrast in contrasts) {
    message(sprintf("Running DE for %s\n", ctrast))
    gc()
    print(ctrast)
    tt <- limma::topTable(fit, coef = ctrast,
                          adjust = "fdr", number = 10000000)
    tt$ensembl_gene_id <- rownames(tt)

    # append gene names
    ensembl_res <- map_ensembl_gene_id(
      tt$ensembl_gene_id,
      mappings = c("external_gene_name", "gene_biotype"),
      ensembl_mapping_file = fargs$ensembl_mapping_file
    ) %>%
      dplyr::rename(gene = external_gene_name)

    tt <- dplyr::left_join(tt, ensembl_res, by = "ensembl_gene_id")
    tt <- tt %>%
      dplyr::rename(padj = adj.P.Val, pval = P.Value) %>%
      dplyr::filter(gene_biotype %in% fargs$gene_biotype)

    tt$contrast <- ctrast
    tt$reference <- fargs$ref_class
    tt$FCRO <- order(abs(2^tt$logFC))

    model_formula_string <- scFlow:::.formula_to_char(model_formula)
    tt$model <- gsub(" ", "", model_formula_string,
                     fixed = TRUE) # no whitespace

    p <- scFlow:::.volcano_plot(
      dt = tt,
      fc_threshold = fargs$fc_threshold,
      pval_cutoff = fargs$pval_cutoff,
      n_label = fargs$n_label
    )

    results <- tt %>%
      dplyr::filter(padj <= fargs$pval_cutoff) %>%
      dplyr::filter(gene_biotype %in% fargs$gene_biotype) %>%
      dplyr::filter(abs(logFC) >= log2(fargs$fc_threshold)) %>%
      dplyr::arrange(FCRO)

    DGEs <- c(sum(results$logFC > 0), sum(results$logFC < 0))
    names(DGEs) <- c("Up", "Down")

    element_name <- paste(fargs$ref_class, ctrast, sep = "_vs_")
    element_name <- gsub("\\$", "", element_name)

    de_params <- list(
      celltype = unique(fargs$sce$cluster_celltype),
      de_method = fargs$de_method,
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

    #if (fargs$de_method != "limma") {
    #  de_params$mast_method <- NULL
    #}

    de_params <- unlist(de_params)

    attr(results, "de_parameters") <- de_params

    attr(results, "de_result") <- DGEs

    attr(results, "plot") <- p

    results_l[[element_name]] <- results
  }

  message("Done!  Returning results")

  return(results_l)
}
