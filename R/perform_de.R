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
#' @importFrom cli cli_text rule
#'
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
                       pval_cutoff = 0.05,
                       ...
) {

  fargs <- c(as.list(environment()), list(...))

  cat(cli::rule(
    "Starting Differential Gene Expression",
    line = 2))
  cat("\r\n")

  cat(cli::rule(
    "Pre-processing SingleCellExperiment",
    line = 1))
  cat("\r\n")

  # preprocess
  fargs$sce <- do.call(.preprocess_sce_for_de, fargs)

  # select method
  if(de_method == "MASTZLM") {
    cat(cli::rule(
      "MAST Zero-inflated Linear Model",
      line = 1))
    cat("\r\n")
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
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom scater librarySizeFactors normalize
#' @importFrom cli cli_text
#'
#' @keywords internal
.preprocess_sce_for_de <- function(...) {

  fargs <- list(...)
  sce <- fargs$sce

  sce <- .filter_sce_genes_for_de(sce,
                                  min_counts = fargs$min_counts,
                                  min_cells_pc = fargs$min_cells_pc
  )

  sce@int_colData$size_factor <- scater::librarySizeFactors(sce)
  sce <- scater::normalize(sce)

  if (fargs$rescale_numerics == TRUE) {
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
    "Setting '{.emph {fargs$ref_class}}' ",
    "as the reference class of '{.strong {fargs$dependent_var}}'.")
  )

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
.perform_de_with_mast <- function(...) {

  fargs <- list(
    mast_method = "bayesglm",
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
      .generate_model_from_vars,
      list(
        sce = sce,
        dependent_var = fargs$dependent_var,
        confounding_vars = fargs$confounding_vars,
        random_effects_var = fargs$random_effects_var,
        prefix = .)
    ))

  # main, without prefixes
  model_formula <- mod_formulae[[1]]

  cli::cli_text(c(
    "Model formula: ",
    .formula_to_char(model_formula))
  )

  #mod_check <- do.call(.check_model, list(model_formula = mod_formulae[[2]]))
  if (is.null(fargs$random_effects_var)) {
    is_full_rank <- .check_model(mod_formulae[[2]])
    if(!fargs$force_run){
      assertthat::assert_that(
        is_full_rank,
        msg = "A full rank model specification is required.")
    } else {
      cli::cli_alert_info("Forcing run, ignoring model full rank.")
    }
  }

  # fit model
  message("Fitting model\n")
  #options(warn=-1) #temporary silencing
  x <- Sys.time()
  zlmCond <- MAST::zlm(
    formula = model_formula,
    sca = sca,
    method = fargs$mast_method, # note: glmer requires a random effects var
    ebayes = FALSE,
    parallel = TRUE
  )
  message(Sys.time() - x)

  gc()

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
    ensembl_res <- map_ensembl_gene_id(
      fcHurdle$ensembl_gene_id,
      mappings = c("external_gene_name", "gene_biotype"),
      ensembl_mapping_file = fargs$ensembl_mapping_file) %>%
      dplyr::rename(gene = external_gene_name)

    fcHurdle <- merge(fcHurdle, ensembl_res, by = "ensembl_gene_id")

    fcHurdle$padj <- p.adjust(fcHurdle$pval, "fdr", n = dim(sca)[[1]])

    fcHurdle$contrast <- ctrast
    fcHurdle$reference <- fargs$ref_class

    model_formula_string <- .formula_to_char(model_formula)
    fcHurdle$model <- gsub(" ", "", model_formula_string, fixed = TRUE) # no whitespace

    p <- .volcano_plot(dt = fcHurdle,
                       fc_threshold = fargs$fc_threshold,
                       pval_cutoff = fargs$pval_cutoff,
                       n_label = fargs$n_label)

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
      min_counts = fargs$min_counts,
      min_cells_pc = fargs$min_cells_pc,
      rescale_numerics = fargs$rescale_numerics,
      dependent_var = fargs$dependent_var,
      ref_class = fargs$ref_class,
      confounding_vars = fargs$confounding_vars,
      random_effects_var = fargs$random_effects_var,
      fc_threshold = fargs$fc_threshold,
      pval_cutoff = fargs$pval_cutoff,
      cells_per_group = table(SingleCellExperiment::colData(fargs$sce)[fargs$dependent_var]),
      n_genes = dim(sce)[[1]],
      model = gsub(" ", "", model_formula_string, fixed = TRUE),
      model_full_rank = is_full_rank,
      contrast_name = element_name
    )

    if (fargs$de_method != "MASTZLM"){
      de_params$mast_method <- NULL
    }

    de_params <- unlist(de_params)

    attr(results, "de_parameters") <- de_params

    attr(results, "de_result") <- DGEs

    attr(results, "plot") <- p

    results_l[[element_name]] <- results

  }

  merged_name <- paste0(fargs$ref_class, "_Merged_Results")
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
#' @param sce the data
#' @param dependent_var the dependent variable
#' @param confounding_vars the confounding variables
#' @param random_effects_var the random effect variables (optional)
#' @param prefix a variable prefix (e.g. "sce$") (optional)
#'
#' @return model_formula a model formula
#'
#' @family differential gene expression
#'
#' @importFrom purrr map_chr
#'
#' @keywords internal
.generate_model_from_vars <- function(sce,
                                      dependent_var,
                                      confounding_vars,
                                      random_effects_var = NULL,
                                      prefix = NULL) {

  if(!is.null(prefix)) {
    dependent_var <- purrr::map_chr(
      dependent_var, ~ paste0(prefix, .))
    if(!is.null(confounding_vars)) {
      confounding_vars <- purrr::map_chr(
        confounding_vars, ~ paste0(prefix, .))
    }
    if(!is.null(random_effects_var)){
      random_effects_var <- purrr::map_chr(
        random_effects_var, ~ paste0(prefix, .))
    }
  }

  # allows a model without confounding variables
  if (!is.null(confounding_vars)) {
    plus_or_blank <- " + "
  } else {
    plus_or_blank <- ""
  }

  if (!is.null(random_effects_var)) {
    random_effects_var <- purrr::map_chr(
      random_effects_var, ~ sprintf("(1|%s)", .)
    )
    random_effects_var <- paste(random_effects_var, collapse = " + ")

    model_formula <- as.formula(
      sprintf(
        "~ %s + %s%s %s",
        dependent_var,
        random_effects_var,
        plus_or_blank,
        paste(confounding_vars, collapse = " + ")
      )
    )
  } else {

    model_formula <- as.formula(
      sprintf(
        "~ %s%s%s",
        dependent_var,
        plus_or_blank,
        paste(confounding_vars, collapse = " + ")
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
  if(limma::is.fullrank(mod0)) {
    cli::cli_alert_success("Model formula is full rank")
    return(TRUE)
  } else {
    cli::cli_alert_danger(
      "The model formula is not full rank.")
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


#' volcano plot
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr %>% filter top_n
#'
#' @keywords internal

.volcano_plot <- function(dt = fcHurdle,
                          fc_threshold = 2,
                          pval_cutoff = 0.001,
                          n_label = 10){
  dt <- dt[!is.na(dt$padj), ]
  dt <- dt[!is.nan(dt$logFC), ]
  dt <- dt[order(dt$padj, decreasing = FALSE), ]
  dt$de <- "Not sig"
  dt$de[dt$padj <= pval_cutoff & dt$logFC > log2(fc_threshold)] <- "Up"
  dt$de[dt$padj <= pval_cutoff & dt$logFC < -log2(fc_threshold)] <- "Down"
  dt$de <- factor(dt$de, levels = c("Up", "Down", "Not sig"))
  top_up <- dt %>%
    dplyr::filter(de == "Up" & gene_biotype == "protein_coding") %>%
    dplyr::top_n(n_label, wt = -padj)
  top_down <- dt %>%
    dplyr::filter(de == "Down" & gene_biotype == "protein_coding") %>%
    dplyr::top_n(n_label, wt = -padj)
  top_de <- rbind(top_up, top_down)
  dt$label <- NA
  dt$label[dt$gene %in% top_de$gene] <- "Yes"


  ggplot2::ggplot(dt) +
    geom_point(aes(x=logFC, y=-log10(padj), fill = de, colour=de), show.legend = T, alpha = 0.5) +
    ggrepel::geom_text_repel(data = dt,
                             aes(logFC, y = -log10(padj), label = ifelse(label == "Yes", gene, "")), max.iter = 1000, size = 3, na.rm = TRUE)+
    #ggtitle("High-inflammation vs control") +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    #scale_y_continuous(breaks = c(0, 10, 20, 50, 100, 200, 300)) +
    #scale_x_continuous(limits = c(-30, 30))+
    geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)), linetype = 2, size = 0.2, alpha = 0.5)+
    geom_hline(yintercept = pval_cutoff, linetype = 2, size = 0.2, alpha = 0.5)+
    scale_colour_manual(name = NULL, aesthetics = c("colour", "fill"), values = c("red", "blue", "grey"), label = c("Up", "Down", "Not sig"))+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    theme(axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA),
          legend.text = element_text(size = 12, colour = "black"),
          legend.key = element_blank(),
          plot.margin = margin(c(1,1,1,1), unit = "cm"))+
    coord_cartesian(xlim = c(-6,6), ylim = c(0, 150), expand = TRUE, clip = "on")
}
