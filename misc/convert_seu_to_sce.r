library("biomaRt")

# import seurat data and convert to sce

# save the mappings from biomart
human <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
attrib_hum <- listAttributes(human)
entrez_df <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "ensembl_gene_id_version"), mart = human)

rn_genes <- data.frame(external_gene_name = rownames(sce), stringsAsFactors = FALSE)

final_df <- dplyr::left_join(rn_genes, dplyr::distinct(entrez_df, external_gene_name, .keep_all = T))
all(rn_genes$external_gene_name == final_df$external_gene_name) # check all aligned

SummarizedExperiment::rowData(sce)$gene <- final_df$external_gene_name
SummarizedExperiment::rowData(sce)$ensembl_gene_id <- final_df$ensembl_gene_id

# drop version from gene name
SummarizedExperiment::rowData(sce)$gene <-
  unlist(lapply(stringr::str_split(SummarizedExperiment::rowData(sce)$gene, "[.]"), "[[",1))

# subset unique only
unique_idx <- !duplicated(SummarizedExperiment::rowData(sce)$gene)
sce <- sce[unique_idx,]

# subset not na
sce <- sce[!is.na(SummarizedExperiment::rowData(sce)$ensembl_gene_id),]

rownames(sce) <- SummarizedExperiment::rowData(sce)$ensembl_gene_id

# standard aliases
sce$individual <- sce$Case
sce$diagnosis <- as.factor(sce$Disease)
sce$pc_mito <- sce$percent.mito
sce$brain_region <- sce$Region

# combine new metadata
  new_metadata <- read.delim("~/Documents/junk/Enriched_Histology.csv", sep = ",")

  imputedata <- new_metadata
  imputedata$individual <- gsub("/", "_", imputedata$individual, fixed = TRUE) # for seu only

  old_df <- as.data.frame(SummarizedExperiment::colData(sce))
  old_df$individual <- as.character(old_df$individual)
  imputedata$individual <- as.character(imputedata$individual)
  new_df <- dplyr::left_join(old_df, imputedata, by = c("individual", "brain_region"))
  #table(new_df$braak, new_df$Braak_tangle_stage)

  # check everything there
  x <- paste0(sce$individual, sce$brain_region)
  y <- paste0(imputedata$individual, imputedata$brain_region)
  unique(x[!(x %in% y)])

  # check
  dim(new_df)

  for (var in c("braak", "amyloid_beta", "p_tau", "gfap", "iba1", "hla", "gpnmb")) {
    sce[[var]] <- new_df[[var]]
  }


# save
write_sce(sce, "~/Documents/Amy_Seurat")

