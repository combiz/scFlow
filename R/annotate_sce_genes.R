################################################################################
#' Add basic gene-wise annotations for a SingleCellExperiment
#'
#' Adds biomaRt annotations (e.g. gene, gene_biotype) and calculates the
#' following QC metrics:
#' * qc_metric_ensembl_mapped - was the ensembl_gene_id found in biomaRt
#' * qc_metric_is_mito - is the gene mitochondrial
#' * qc_metric_is_ribosomal - is the gene ribosomal
#'
#' @param sce a SingleCellExperiment object with ensembl_gene_id rowData
#' @param ensembl_mapping_file a local tsv file with ensembl_gene_id and
#'   additional columns for mapping ensembl_gene_id to gene info.  If
#'   not provided, the biomaRt db is queried (slower).
#'
#' @return sce a SingleCellExperiment object annotated with gene data
#'
#' @family annotation functions
#' @import cli Matrix dplyr
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment rowData colData
#'
#' @export
annotate_sce_genes <- function(sce,
                               ensembl_mapping_file = NULL) {

  cat(cli::rule("Annotating SingleCellExperiment genes", line = 1), "\r\n")

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  if (!("ensembl_gene_id" %in%
        colnames(SummarizedExperiment::rowData(sce)))) {
    stop(cli::cli_alert_danger("The rowData is missing ensembl_gene_id."))
  }

  # annotate rowdata with biomart data
  mapped_ensembl_ids <- map_ensembl_gene_id(
    SummarizedExperiment::rowData(sce)$ensembl_gene_id,
    mappings_filepath = ensembl_mapping_file
  )

  mapped_ensembl_ids$ensembl_gene_id <- factor(
    mapped_ensembl_ids$ensembl_gene_id,
    levels = sort(union(
      levels(mapped_ensembl_ids$ensembl_gene_id),
      levels(SummarizedExperiment::rowData(sce)$ensembl_gene_id)))
  )

  SummarizedExperiment::rowData(sce) <- dplyr::full_join(
    data.frame(SummarizedExperiment::rowData(sce)),
    mapped_ensembl_ids,
    by = "ensembl_gene_id",
    all = TRUE
    ) %>% dplyr::rename(gene = external_gene_name)

  if(!all(
    SummarizedExperiment::rowData(sce)$ensembl_gene_id == rownames(sce))){
    stop(cli::cli_alert_danger("Fatal error: Misaligned new_rowdata."))
  }

  # annotate rowdata with qc_metric_ensembl_mapped
  SummarizedExperiment::rowData(sce)$qc_metric_ensembl_mapped <-
    !(is.na(SummarizedExperiment::rowData(sce)$gene)) + 0

  # annotate rowdata with qc_metric_mitochondrial
  qc_metric_is_mito <- startsWith(
    as.character(SummarizedExperiment::rowData(sce)$gene), "MT-") + 0
  qc_metric_is_mito[is.na(qc_metric_is_mito)] <- 0 # for non-mapped genes
  SummarizedExperiment::rowData(sce)$qc_metric_is_mito <- qc_metric_is_mito

  # annotate rowdata with qc_metric_is_ribo
  qc_metric_is_ribo <- startsWith(
    as.character(SummarizedExperiment::rowData(sce)$gene), "RPS") + 0
  qc_metric_is_ribo[is.na(qc_metric_is_ribo)] <- 0
  SummarizedExperiment::rowData(sce)$qc_metric_is_ribo <-
    qc_metric_is_ribo

  return(sce)

}

