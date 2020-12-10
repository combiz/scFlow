# Attaching some column metadata.
library(edgeR)

model_mat
#d <- edgeR::DGEList(counts = mat, group = sce[[fargs$dependent_var]])

d <- scran::convertTo(sce, type = "edgeR")

#d <- edgeR::DGEList(counts = mat, group = sce$sex, samples = as.data.frame(SummarizedExperiment::colData(sce)), norm.factors = 2^rnorm())
y
d
#mm <- model.matrix(~ sce$sex + sce$pc_mito + sce$cngeneson + sce$diagnosis)
mm <- model.matrix(~ sce$diagnosis + sce$pc_mito + sce$cngeneson + sce$sex)
mm <- model.matrix(~ sce$group + sce$pc_mito + sce$cngeneson + sce$sex)
colnames(mm) <- make.names(colnames(mm))

d <- calcNormFactors(d)

d <- estimateGLMCommonDisp(d,mm)
d <- estimateGLMTrendedDisp(d,mm)
d <- estimateGLMTagwiseDisp(d,mm)

d <- estimateDisp(d, mm, trend="none")
plotBCV(d, cex=1)
fit <- glmFit(d, mm, robust=TRUE, abundance.trend=FALSE)
plotQLDisp(fit, cex = 1)
results <- glmQLFTest(fit, coef = ncol(mm))
#fit <- glmQLFTest(fit, coef = ncol(mm))
topTags(results)
colnames(results$coefficients)
colnames(fit$coefficients)
#contrasts <- c("sce.diagnosisMS")
contrasts <- c("sce.groupHigh")

results_l <- list()
for (contrast in contrasts) {
  message(sprintf("Performing DE for contrast: %s", contrast))
  con <- makeContrasts(contrasts = contrast, levels = colnames(mm))
  lrt <- glmLRT(fit,contrast=con)
  #lrt <- glmQLFit(fit,contrast=con)
  res <- as.data.frame(topTags(lrt,nrow(d)))
  res$ensembl_gene_id <- rownames(res)
  res$padj <- p.adjust(res$PValue, method = "fdr", n = length(res$PValue))
  results_l[[contrast]] <- res
  results_l[[contrast]]$contrast <- contrast
  # append gene names
  # append gene names
  ensembl_res <- map_ensembl_gene_id(
    res$ensembl_gene_id,
    mappings = c("external_gene_name", "gene_biotype"),
    ensembl_mapping_file = fargs$ensembl_mapping_file
  ) %>%
    dplyr::rename(gene = external_gene_name)
  results_l[[contrast]] <- merge(results_l[[contrast]],
                             ensembl_res, by = 'ensembl_gene_id')
}
