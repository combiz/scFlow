options(mc.cores = 10)
library(parallel)
options(mc.cores = 10)
library(scFlow)
library(magrittr)
#seu <- readRDS("~/Documents/Amy_Seurat/GliaEnrichmentExpt_20200113_AS.Rds")
#sce <- Seurat::as.SingleCellExperiment(seu)
#write_sce(sce, "~/Documents/Amy_Seurat")
#sce <- read_sce("~/Documents/junk/enriched_with_ptau")
sce <- read_sce("~/Documents/Amy_Seurat")


sce_all <- sce

#sce_all$cluster_celltype <- "Other"
#sce_all$cluster_celltype[sce_all$seurat_clusters %in% c("1", "6", "7", "9")] <- "Micro"
#sce_all$cluster_celltype[sce_all$seurat_clusters %in% c("0", "2", "3", "4", "5", "10")] <- "Astro"
#sce <- sce_all

## ~ diagnosis
#for (region in c("EC", "SSC", "ALL")) {
for (celltype in c("Micro", "Astro")) {

  for (region in c("ALL")) {

#      if(region != "ALL") {
#        sce_subset <- sce_all[, sce_all$cluster_celltype == celltype & sce_all$brain_region == region]
#        confounders <- c(
#          "cngeneson",
#          #"sex",
#          #"age",
#          "pc_mito"
#          #"pc_ribo"
#        )
#      } else {
#        sce_subset <- sce_all[, sce_all$cluster_celltype == celltype]
#        confounders <- c(
#          "cngeneson",
#          #"sex",
#          #"age",
#          "pc_mito"
#          #"pc_ribo"
#          #"brain_region"
#        )
#      }
#
#      de_results <- perform_de(
#        sce = sce_subset,
#        de_method = "MASTZLM",
#        ebayes = FALSE,
#        mast_method = "glmer",
#        min_counts = 1,
#        min_cells_pc = 0.10,
#        rescale_numerics = TRUE,
#        dependent_var = "diagnosis",
#        ref_class = "Control",
#        confounding_vars = confounders,
#        random_effects_var = "individual",
#        fc_threshold = 1,
#        pval_cutoff = 1,
#        ensembl_mapping_file = "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv",
#        unique_id_var = "individual",
#        subset_var = "diagnosis",
#        subset_class = "Control"
#      )
#
#      #for (contrast in names(de_results)) {
#      contrast <- "diagnosis"
#      file_name <- sprintf("%s_%s_%s_glmer.tsv", celltype, region, contrast)
#      write.table(
#        #de_results[[contrast]],
#        de_results[[1]],
#        file.path("~/Documents/junk/enriched_de_seu", file_name),
#        col.names = TRUE, row.names = FALSE, sep = "\t"
#      )
#      #}
    }


# ~ p_tau
#for (region in c("EC", "SSC", "ALL")) {
  for (region in c("ALL")) {
      if(region != "ALL") {
        sce_subset <- sce_all[, sce_all$cluster_celltype == celltype & sce_all$brain_region == region]
        confounders <- c(
          "cngeneson",
          #"sex",
          #"age",
          "pc_mito"
          #"pc_ribo"
        )
      } else {
        sce_subset <- sce_all[, sce_all$cluster_celltype == celltype]
        confounders <- c(
          "cngeneson",
          #"sex",
          #"age",
          "pc_mito"
          #"pc_ribo"
          #"brain_region"
        )
      }

      de_results <- perform_de(
        sce = sce_subset,
        de_method = "MASTZLM",
        ebayes = FALSE,
        mast_method = "glmer",
        min_counts = 1,
        min_cells_pc = 0.10,
        rescale_numerics = TRUE,
        dependent_var = "p_tau",
        ref_class = NULL,
        confounding_vars = confounders,
        random_effects_var = "individual",
        fc_threshold = 1,
        pval_cutoff = 1,
        ensembl_mapping_file = "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv",
        unique_id_var = "individual",
        subset_var = "diagnosis",
        subset_class = "Control"
      )

      #for (contrast in names(de_results)) {
      contrast <- "p_tau"
      file_name <- sprintf("%s_%s_%s_glmer.tsv", celltype, region, contrast)
      write.table(
        #de_results[[contrast]],
        de_results[[1]],
        file.path("~/Documents/junk/enriched_de_seu", file_name),
        col.names = TRUE, row.names = FALSE, sep = "\t"
        )
      #}
    }

}



##############################



#################


#library(ggplot2)
#pc_cutoff <- 0.02
#rank_cutoff <- dim(vedf)[[1]] * pc_cutoff
#vedf$drop <- vedf$rank <= rank_cutoff
#ggplot(vedf, aes(x = rank, y = individual)) +
#  geom_point(aes(colour = drop)) +
#  ggrepel::geom_text_repel(data = vedf[1:25,], aes(label = ifelse(drop == TRUE, as.character(.data[["gene"]]), "")),
#    max.iter = 1000, size = 3, na.rm = TRUE
#  ) +
#  theme_bw() +
#  theme(
#    axis.text = element_text(size = 16),
#    axis.title = element_text(size = 18)
#  )

## compare
#de_results <- read.delim("~/Documents/junk/enriched_de/ALL_MICRO_GLMER_PTAU.tsv")
#de_results <- de_results[de_results$padj <= 0.05,]

#hvgs <- vedf[vedf$drop == TRUE, ]
#dim(hvgs)

#sum(hvgs$gene %in% de_results$gene)
