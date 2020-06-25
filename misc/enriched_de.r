options(mc.cores = max(1, future::availableCores() - 2))

param = BiocParallel::SnowParam(max(1, future::availableCores() - 2), "SOCK", progressbar = TRUE)
BiocParallel::register(param)

library(parallel)
#options(mc.cores = 10)
  library(scFlow)
library(magrittr)
#sce <- read_sce("~/Documents/junk/enriched/final_sce/")
sce <- read_sce("~/Documents/junk/enriched_impute/")


##########

class(sce$p.Tau)
#sce$seqdate <- as.factor(sce$seqdate)
#sce_all <- sce
#sce_subset <- sce_all[, sce_all$cluster_celltype == "Micro" & sce_all$brain_region == "EC"] #EC
#sce <- sce_subset
#sce_subset <- sce_all[, sce_all$cluster_celltype == "Astro" & sce_all$brain_region == "SSC"] #EC
#sce_ssc <- sce[, sce$brain_region == "SSC"]
#sce_ec <- sce[, sce$brain_region == "EC"]


# glmer
fargs <-  list(
  sce = sce_subset,
  de_method = "MASTZLM",
  ebayes = FALSE,
  mast_method = "glmer",
  min_counts = 1,
  min_cells_pc = 0.10,
  rescale_numerics = TRUE,
  #dependent_var = "diagnosis_region",
  dependent_var = "diagnosis",
  ref_class = "control",
  #ref_class = "control_EC",
  confounding_vars = c(
    "cngeneson",
    "sex",
    "age",
    #"PM_delay",
    #"RIN",
    #"seqdate",
    "pc_mito",
    "pc_ribo"
  ),
  random_effects_var = "individual",
  fc_threshold = 1,
  #pval_cutoff = 0.1,
  pval_cutoff = 1.0,
  ensembl_mapping_file = "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv"
)
de_results <- do.call(perform_de, fargs)

#qs::qsave(de_results$decontrol_vs_diagnosiscase, "~/Documents/junk/enriched_de/Micro_SSC_decontrol_vs_diagnosiscase.qs")
#write.table(de_results$control_vs_diagnosiscase, "~/Documents/junk/enriched_de/Micro_SSC__control_vs_diagnosiscase.tsv", row.names = F, sep = "\t")

micro_ec <- qs::qread("~/Documents/junk/enriched_de/Micro_EC_decontrol_vs_diagnosiscase.qs")
micro_ssc <- qs::qread("~/Documents/junk/enriched_de/Micro_SSC_decontrol_vs_diagnosiscase.qs")
astro_ec <- qs::qread("~/Documents/junk/enriched_de/astro_ec_control_vs_diagnosiscase.qs")
astro_ssc <- qs::qread("~/Documents/junk/enriched_de/Astro_SSC_decontrol_vs_diagnosiscase.qs")

micro_ec <- read.delim("~/Documents/junk/enriched_de/Micro_EC__control_vs_diagnosiscase.tsv")
micro_ssc <- read.delim("~/Documents/junk/enriched_de/Micro_SSC__control_vs_diagnosiscase.tsv")
astro_ssc <- read.delim("~/Documents/junk/enriched_de/Astro_SSC_control_vs_diagnosiscase.tsv")
astro_ec <- read.delim("~/Documents/junk/enriched_de/astro_ec_control_vs_diagnosiscase.tsv")
dt <- micro_ssc

scFlow:::.volcano_plot(dt, fc_threshold = 1.05, pval_cutoff = 0.05, n_label = 12)
de_genes <- dt %>%
  dplyr::filter(padj <= 0.05 & 2^abs(logFC) >= 1.05)

de_genes %>% dplyr::group_by(logFC > 0) %>% dplyr::tally()

plot_violin(
  sce_ec,
  group_var = "diagnosis",
  subset_var = "cluster_celltype",
  subset_group = "Astro",
  gene = "INO80D",
  var_order = c("control", "case")
  )


write.table(de_genes, "de_genes.tsv", row.names = F, sep = "\t")

ipa_res <- find_impacted_pathways(de_genes)

de

# bayesglm
fargs <-  list(
  sce = sce_subset,
  de_method = "MASTZLM",
  ebayes = FALSE,
  mast_method = "bayesglm",
  min_counts = 1,
  min_cells_pc = 0.10,
  rescale_numerics = TRUE,
  #dependent_var = "diagnosis_region",
  dependent_var = "sex",
  ref_class = "male",
  #ref_class = "control_SSC",
  confounding_vars = c(
    "individual",
    "cngeneson",
    #"sex",
    "age",
    "pc_mito",
    "pc_ribo"
  ),
  #random_effects_var = "individual",
  random_effects_var = NULL,
  #fc_threshold = 1.05,
  fc_threshold = 1,
  pval_cutoff = 1,
  #pval_cutoff = 0.1,
  ensembl_mapping_file = "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv",
  force_run = TRUE
)
#de_results <- do.call(perform_de_test, fargs)
#sce_pb <- pseudobulk_sce(sce_subset, keep_vars = c("individual", "age", "sex", "diagnosis"))
de_results <- do.call(perform_de, fargs)

#write.table(de_results$control_vs_diagnosiscase, "~/Documents/junk/enriched_de/Micro_SSC_BAYESGLM_DIAGNOSISONLY_control_vs_diagnosiscase.tsv", row.names = F, sep = "\t")

#dt <- read.delim("~/Documents/junk/enriched_de/Astro_SSC_BAYESGLM_control_vs_diagnosiscase.tsv")

dt <- de_results$control_vs_diagnosiscase
dt <- de_results$male_vs_sexfemale

scFlow:::.volcano_plot(dt, fc_threshold = 2^0.25, pval_cutoff = 0.05, n_label = 12)
de_genes <- dt %>%
  dplyr::filter(padj <= 0.05 & 2^abs(logFC) >= 2^0.25)

chr_mappings <- map_ensembl_gene_id(de_genes$ensembl_gene_id, mappings = c("chromosome_name"))
de_genes <- dplyr::left_join(de_genes, chr_mappings)
table(de_genes$chromosome_name)

de_genes %>% dplyr::group_by(logFC > 0) %>% dplyr::tally()

write.table(de_genes, "de_genes.tsv", row.names = F, sep = "\t")



########### DCESEQ2

expt_contrasts <- list()
expt_contrasts[["numerators"]] <- "male"
expt_contrasts[["denominators"]] <- "female"

pb_mat <- as.matrix(SingleCellExperiment::normcounts(pb))
pb_de <- pb
SingleCellExperiment::counts(pb_de) <- pb_mat
dds <- DESeq2::DESeqDataSet(pb_de, design = ~ sex)
dds <- DESeq2::DESeq(dds)

contrasts <- DESeq2::resultsNames(dds)[startsWith(DESeq2::resultsNames(dds), "sex")]
res_l <- list()
for (contrast in contrasts) {
  res_l[[contrast]] <- DESeq2::lfcShrink(dds, coef = contrast, type = "apeglm")
  res_l[[contrast]]$contrast <- contrast
  res_l[[contrast]]$ensembl_gene_id <- rownames(res_l[[contrast]])

  # append gene names
  ensembl_res <- map_ensembl_gene_id(res_l[[contrast]]$ensembl_gene_id) %>%
    dplyr::rename(gene = external_gene_name)
  res_l[[contrast]] <- dplyr::left_join(data.frame(res_l[[contrast]]),
                             ensembl_res, by = "ensembl_gene_id")
}
res <- Reduce(rbind, res_l)
dt <- res
dt$logFC <- dt$log2FoldChange


scFlow:::.volcano_plot(dt, fc_threshold = 2^0.25, pval_cutoff = 0.05, n_label = 12)
de_genes <- dt %>%
  #dplyr::filter(padj <= 0.05 & 2^abs(logFC) >= 2^0.25)
  dplyr::filter(padj <= 0.2 & 2^abs(logFC) >= 2^0.25)

chr_mappings <- map_ensembl_gene_id(de_genes$ensembl_gene_id, mappings = c("chromosome_name"))
de_genes <- dplyr::left_join(de_genes, chr_mappings)
table(de_genes$chromosome_name)

de_genes %>% dplyr::group_by(logFC > 0) %>% dplyr::tally()


#### edgeR

expt_contrasts <- list()
expt_contrasts[["numerators"]] <- "male"
expt_contrasts[["denominators"]] <- "female"


pb_mat <- as.matrix(SingleCellExperiment::counts(pb))
pb_de <- pb
SingleCellExperiment::counts(pb_de) <- pb_mat
d <- edgeR::DGEList(counts=pb_mat, group=pb$sex)
d <- calcNormFactors(d)

contrasts <- DESeq2::resultsNames(dds)[startsWith(DESeq2::resultsNames(dds), "sex")]
res_l <- list()
for (contrast in contrasts) {
  res_l[[contrast]] <- DESeq2::lfcShrink(dds, coef = contrast, type = "apeglm")
  res_l[[contrast]]$contrast <- contrast
  res_l[[contrast]]$ensembl_gene_id <- rownames(res_l[[contrast]])

  # append gene names
  ensembl_res <- map_ensembl_gene_id(res_l[[contrast]]$ensembl_gene_id) %>%
    dplyr::rename(gene = external_gene_name)
  res_l[[contrast]] <- dplyr::left_join(data.frame(res_l[[contrast]]),
                                        ensembl_res, by = "ensembl_gene_id")
}
res <- Reduce(rbind, res_l)
dt <- res
dt$logFC <- dt$log2FoldChange


scFlow:::.volcano_plot(dt, fc_threshold = 2^0.25, pval_cutoff = 0.05, n_label = 12)
de_genes <- dt %>%
  #dplyr::filter(padj <= 0.05 & 2^abs(logFC) >= 2^0.25)
  dplyr::filter(padj <= 0.2 & 2^abs(logFC) >= 2^0.25)

chr_mappings <- map_ensembl_gene_id(de_genes$ensembl_gene_id, mappings = c("chromosome_name"))
de_genes <- dplyr::left_join(de_genes, chr_mappings)
table(de_genes$chromosome_name)

de_genes %>% dplyr::group_by(logFC > 0) %>% dplyr::tally()

