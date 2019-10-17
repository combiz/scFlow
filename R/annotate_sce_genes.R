################################################################################
#' Add basic gene-wise annotations for a SingleCellExperiment
#'
#' Adds biomaRt annotations (e.g. gene, gene_biotype) and calculates the
#' following QC metrics:
#' * qc_metric_ensembl_mapped - was the ensembl_gene_id found in biomaRt
#' * qc_metric_is_mito - is the gene mitochondrial
#'
#' @param sce a SingleCellExperiment object with ensembl_gene_id rowData
#'
#' @return sce a SingleCellExperiment object annotated with gene data
#'
#' @family annotation functions
#' @import cli Matrix SingleCellExperiment dplyr
#' @export
annotate_sce_genes <- function(sce,
                               ensembl_mapping_file = NULL) {

  if (typeof(sce) != "S4") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  if (!("ensembl_gene_id" %in% colnames(rowData(sce)))) {
    stop(cli::cli_alert_danger("The rowData is missing ensembl_gene_id."))
  }

  before_rowdata_colnames <- colnames(rowData(sce))

  # annotate rowdata with biomart data
  mapped_ensembl_ids <- map_ensembl_gene_id(
    rowData(sce)$ensembl_gene_id,
    mappings_filepath = ensembl_mapping_file
  )

  mapped_ensembl_ids$ensembl_gene_id <- factor(
    mapped_ensembl_ids$ensembl_gene_id,
    levels = sort(union(
      levels(mapped_ensembl_ids$ensembl_gene_id),
      levels(rowData(sce)$ensembl_gene_id)))
  )

  rowData(sce) <- dplyr::full_join(
    data.frame(rowData(sce)),
    mapped_ensembl_ids,
    by = "ensembl_gene_id",
    all = TRUE
  ) %>% dplyr::rename(gene = external_gene_name)

  # annotate rowdata with qc_metric_ensembl_mapped
  rowData(sce)$qc_metric_ensembl_mapped <- is.na(rowData(sce)$gene) + 0

  # annotate rowdata with qc_metric_mitochondrial
  qc_metric_is_mito <- startsWith(
    as.character(rowData(sce)$gene), "MT-") + 0
  qc_metric_is_mito[is.na(qc_metric_is_mito)] <- 0 # for non-mapped genes
  rowData(sce)$qc_metric_is_mito <- qc_metric_is_mito

  cli::cli_alert_success("SingleCellExperiment genes were successfully \
                         annotated with: ")
  cat(paste(setdiff(colnames(rowData(sce)),
                    before_rowdata_colnames), collapse = ", "))

  return(sce)

}

