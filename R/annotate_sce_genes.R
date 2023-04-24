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
#' @param drop_unmapped set `TRUE` to remove unmapped ensembl_gene_id
#' @param drop_mito set `TRUE` to remove mitochondrial genes
#' @param drop_ribo set `TRUE` to remove ribosomal genes
#' @param ensembl_mapping_file a local tsv file with ensembl_gene_id and
#'   additional columns for mapping ensembl_gene_id to gene info.  If
#'   not provided, the biomaRt db is queried (slower).
#' @param species The biological species of the sample.#'
#'
#' @return sce a SingleCellExperiment object annotated with gene data
#'
#' @family annotation functions
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom dplyr left_join rename
#'
#' @export
annotate_sce_genes <- function(sce,
                               drop_unmapped = TRUE,
                               drop_mito = TRUE,
                               drop_ribo = FALSE,
                               ensembl_mapping_file = NULL,
                               species = getOption(
                                 "scflow_species",
                                 default = "human")) {

  cat(cli::rule("Annotating SingleCellExperiment genes", line = 2), "\r\n")

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  if (!("ensembl_gene_id" %in%
        colnames(SummarizedExperiment::rowData(sce)))) {
    stop(cli::cli_alert_danger("The rowData is missing ensembl_gene_id."))
  }

  # annotate rowdata with biomart data
  mapping_results <- map_ensembl_gene_id(
    SummarizedExperiment::rowData(sce)$ensembl_gene_id,
    ensembl_mapping_file = ensembl_mapping_file,
    species = species
  )

  mapping_results$ensembl_gene_id <- as.character(
    mapping_results$ensembl_gene_id
  )

  SummarizedExperiment::rowData(sce) <- dplyr::left_join(
    data.frame(SummarizedExperiment::rowData(sce)),
    mapping_results,
    by = "ensembl_gene_id",
    ) %>% dplyr::rename(gene = external_gene_name)

  if (!all(
    SummarizedExperiment::rowData(sce)$ensembl_gene_id == rownames(sce))) {
    stop(cli::cli_alert_danger("Fatal error: Misaligned new_rowdata."))
  }

  # annotate rowdata with qc_metric_ensembl_mapped
  SummarizedExperiment::rowData(sce)$qc_metric_ensembl_mapped <-
    !(is.na(SummarizedExperiment::rowData(sce)$gene)) + 0

  # annotate rowdata with qc_metric_mitochondrial
  qc_metric_is_mito <- grepl(
    "^mt|^MT", as.character(SummarizedExperiment::rowData(sce)$gene)
    )
  #qc_metric_is_mito <- startsWith(
  #  as.character(SummarizedExperiment::rowData(sce)$gene), mito_prefix) + 0
  qc_metric_is_mito[is.na(qc_metric_is_mito)] <- 0 # for non-mapped genes
  SummarizedExperiment::rowData(sce)$qc_metric_is_mito <- qc_metric_is_mito

  # annotate rowdata with qc_metric_is_ribo
  qc_metric_is_ribo <- grepl(
    "^RPS|^Rps|^RPL|^Rpl", as.character(SummarizedExperiment::rowData(sce)$gene)
  )
  #qc_metric_is_ribo <- startsWith(
  #  as.character(SummarizedExperiment::rowData(sce)$gene), "RPS") + 0
  qc_metric_is_ribo[is.na(qc_metric_is_ribo)] <- 0
  SummarizedExperiment::rowData(sce)$qc_metric_is_ribo <-
    qc_metric_is_ribo

  # boolean logic to obtain keep flags
  SummarizedExperiment::rowData(sce)$qc_metric_mapped_keep <-
    (SummarizedExperiment::rowData(sce)$qc_metric_ensembl_mapped |
       !drop_unmapped)
  SummarizedExperiment::rowData(sce)$qc_metric_mito_keep <-
    !(qc_metric_is_mito & drop_mito)
  SummarizedExperiment::rowData(sce)$qc_metric_ribo_keep <-
    !(qc_metric_is_ribo & drop_ribo)

  sce@metadata$scflow_steps$genes_annotated <- 1

  return(sce)

}
