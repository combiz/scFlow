fargs <- list(
  sce = sce,
  min_counts = 1,
  min_cells_pc = 0.10,
  rescale_numerics = TRUE,
  dependent_var = "diagnosis",
  ref_class = "control",
  confounding_vars = c("individual",
                       "diagnosis",
    "cngeneson",
    "sex",#,
    "age",
    #"PMI",
    "RIN",
    #"seqdate",
    "pc_mito"
    #),
  ),
  random_effects_var = NULL,
  fc_threshold = 1,
  pval_cutoff = 1,
  ensembl_mapping_file = "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv",
  n_label = 12,
  gene_biotype = "protein_coding",
  pseudobulk = TRUE,
  #pseudobulk = FALSE,
  sctransform = TRUE,
  quantile_norm = TRUE,
  unique_id_var = "individual"
)
sce_all <- sce
#sce_subset <- sce_all[, sce_all$cluster_celltype == "Micro"]
sce_subset <- sce_all[, sce_all$cluster_celltype == "Oligo"]
sce <- sce_subset
#sce[[fargs$dependent_var]] <- relevel(
#  sce[[fargs$dependent_var]],
#  ref = fargs$ref_class
#)

# ms
sce$sex <- relevel(
  sce$sex,
  ref = "M"
)

sce$group <- relevel(
  sce$group,
  ref = "Low"
)

fargs$sce <- sce
sce_pp <- do.call(scFlow:::.preprocess_sce_for_de, fargs)

sce <- sce_pp
fargs$sce <- sce_pp
res <- do.call(.perform_de_with_limma, fargs)
sce
sce <- sce_subset

idx <- as.numeric(caret::createDataPartition(sce$individual, p = .10, list = FALSE)) # 15% subset
mini_sce <- sce[, idx]
sce <- mini_sce

dt <- results_l$male_vs_diagnosiscontrol
dt <- results_l$control_vs_diagnosiscase
dt <- results_l$Control_vs_sexM
dt <- results_l$Control_vs_sexF
dt <- results_l$Control_vs_diagnosisMS
dt <- results_l$Control_vs_groupHigh
dt <- results_l$Control_vs_groupLow
dt$old_padj <- dt$padj
dt$padj <- dt$pval
scFlow:::.volcano_plot(res$male_vs_scesexfemale, fc_threshold = 1.05, pval_cutoff = 0.05, n_label = 12)
dt$gene <- dt$gene.x
scFlow:::.volcano_plot(dt, fc_threshold = 1.05, pval_cutoff = 0.5, n_label = 12)
dt$gene <- dt$gene.x
dt <- results_l$control_vs_diagnosiscase
dt <- results_l$Control_vs_sexF
dt <- results_l$Control_vs_scesexF
dt <- results_l$sce.diagnosisMS
dt <- results_l$sce.groupHigh
de_genes <- dt %>%
  dplyr::filter(padj <= 0.05 & 2^abs(logFC) >= 1.05)

de_genes <- dt %>%
  dplyr::filter(padj <= 0.2 & 2^abs(logFC) >= 1.05)

de_genes %>% dplyr::group_by(logFC > 0) %>% dplyr::tally()

chr_mappings <- map_ensembl_gene_id(dt[1:50,]$ensembl_gene_id, mappings = c("chromosome_name"))
#de_genes <- dplyr::left_join(de_genes, chr_mappings)
table(de_genes$chromosome_name[1:50])
head(de_genes)
dt <- results_l$male_vs_sexfemale
dt <- results_l$sce.sexF
dt_l <- rbind(dt, dt_l)

write.table(dt, "~/Documents/junk/ms_tmp/MS_Oligo_control_vs_High_dream.tsv", row.names = F, sep = "\t")

## y chr plots
y_dt <- dt_l[dt_l$gene %in% c("UTY", "PRKY", "USP9Y", "NLGN4Y"),]
ggplot(y_dt, aes(x = model, y = padj)) + geom_col() + coord_flip() + facet_grid(. ~ gene)

### dream #############

library(variancePartition) # NEEDED TO WORK!
#model_formula <- as.formula("~ sex + (1 | individual) + cngeneson + age + pc_mito + pc_ribo + (1 | diagnosis)") # nope
#model_formula <- as.formula("~ sex + (1 | individual) + cngeneson + pc_mito + pc_ribo + (1 | diagnosis)")
model_formula <- as.formula("~ sex + (1 | individual) + cngeneson + pc_mito + pc_ribo") # best so far
model_formula <- as.formula("~ diagnosis + (1 | individual) + cngeneson + pc_mito + pc_ribo")
model_formula <- as.formula("~ group + (1 | individual) + cngeneson + pc_mito + pc_ribo")
model_formula <- as.formula("~ diagnosis") # pb
model_formula <- as.formula("~ sex") # pb
model_formula <- as.formula("~ sex + pc_mito + cngeneson") # pb
model_formula <- as.formula("~ diagnosis + pc_mito + sex") # pb

sce$cngeneson
sce <- sce_pp

#sce_c2 <- sce[,sce$clusters == "2"]
#sce_c10 <- sce[,sce$clusters == "10"]
#sce <- sce_c10
#sce <- sce_c2

human <- biomaRt::useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
entrez_table = biomaRt::getBM(attributes=c("ensembl_gene_id","chromosome_name","hgnc_symbol"), mart=human)

plot_reduced_dim(sce, feature_dim = "individual", reduced_dim = "UMAP_Liger", size = 4)
library(SingleCellExperiment)
sce <- integrate_sce(sce, unique_id_var = "individual", k = 20)
sce <- reduce_dims_sce(sce, input_reduced_dim = "Liger", reduction_methods = "UMAP", k = 30, res = 0.01)


## pb tst
mm <- model.matrix(~ 0 + sce_subset$individual)
colnames(mm) <- levels(sce_subset$individual)
pb_matrix <- SingleCellExperiment::counts(sce_subset) %*% mm
pb_matrix[1,]

x <- sce_subset[, sce_subset$individual == "PDC016"]
x <- SingleCellExperiment::counts(sce_subset[,sce_subset$individual == "PDC016"])
x <- SingleCellExperiment::counts(sce)
x <- rowSums(x)
x[1]
##

mat <- as(SingleCellExperiment::counts(sce), "dgCMatrix")
vobjDream <- variancePartition::voomWithDreamWeights( mat, model_formula, as.data.frame(SummarizedExperiment::colData(sce)) )
fit <- variancePartition::dream( vobjDream, model_formula, as.data.frame(SummarizedExperiment::colData(sce)) )
if (!fargs$random_effects_var) {
  fit <- limma::eBayes(fit)
}

plot_reduced_dim(sce, feature_dim = "individual", reduced_dim = "PCA", size = 4)
pca_mat <- SingleCellExperiment::reducedDim(sce, "PCA")
plot(sce$cngeneson, pca_mat[, 1])

# model colinearity plots
#bad model
C <- variancePartition::canCorPairs( ~ sex + cngeneson + individual + age + pc_mito + pc_ribo + diagnosis + Braak_tangle_stage + PM_delay + RIN, as.data.frame(SummarizedExperiment::colData(sce)))
#simplified model without colinear vars
C <- variancePartition::canCorPairs( ~ sex + cngeneson + individual + pc_mito + pc_ribo, as.data.frame(SummarizedExperiment::colData(sce)))
C <- variancePartition::canCorPairs( ~ group + cngeneson + individual + pc_mito + pc_ribo, as.data.frame(SummarizedExperiment::colData(sce)))
# Plot correlation matrix
variancePartition::plotCorrMatrix( C )

# variance
varPart <- variancePartition::fitExtractVarPartModel( vobjDream, model_formula, as.data.frame(SummarizedExperiment::colData(sce)), showWarnings=FALSE )
variancePartition::plotVarPart( sortCols(varPart), label.angle=60 )
plotPercentBars( varPart[1:15,] )
var_genes <- map_ensembl_gene_id(rownames(varPart), ensembl_mapping_file = fargs$ensembl_mapping_file)
varpart_ensembl <- rownames(varPart)
rownames(varPart) <- var_genes$external_gene_name
prob_varPart <- rbind(varPart[1:15,], varPart[rownames(varPart) %in% c("LINGO1", "RASGEF1B", "SLC26A3","INO80D", "TMEM163"),])
plotPercentBars( prob_varPart )
varPart

######################

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
  model_mat <- model.matrix(~ sce$sex + sce$cngeneson + sce$pc_mito)
  model_mat <- model.matrix(~ sce$diagnosis + sce$pc_mito + sce$sex)
  model_mat <- stats::model.matrix(model_formula)
  # fit model
  message("Fitting model\n")

  # options(warn=-1) #temporary silencing
  x <- Sys.time()

  if (fargs$sctransform | fargs$quantile_norm) {
    mat <- SingleCellExperiment::normcounts(sce)
  } else {
    mat <- as(SingleCellExperiment::counts(sce), "dgCMatrix")
  }

  #dupcor <- limma::duplicateCorrelation(
  #  mat[sample(1:dim(mat)[1], 400, replace = TRUE),],
  #  model_mat,
  #  block = sce[[unique_id_var]]
  #  )

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
    #~ paste0("sce$", fargs$dependent_var, .)
    ~ paste0(fargs$dependent_var, .)
  )
  print(contrasts)
  print("colnames:")
  print(colnames(fit$design))

  contrasts <- c("sexfemale")
  #contrasts <- c("diagnosiscase")
  contrasts <- c("diagnosisMS")
  contrasts <- c("sexF")
  contrasts <- c("sce$sexF")
  contrasts <- c("groupHigh", "groupLow")
  contrasts <- c("groupHigh")
  contrasts <- "diagnosiscase"

    colnames(fit$coefficients)

colnames(fit$coefficients)

  results_l <- list()

  for (ctrast in contrasts) {
    message(sprintf("Running DE for %s\n", ctrast))
    gc()
    print(ctrast)
    tt <- limma::topTable(fit, coef = ctrast,
                          adjust = "fdr", number = 10000000, sort.by = "p")
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
