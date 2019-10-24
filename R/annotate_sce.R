################################################################################
#' Add basic annotations and qc metrics for a SingleCellExperiment
#'
#' Adds biomaRt annotations (e.g. gene, gene_biotype) and QC metric annotations.
#'
#' Cell-level annotations
#' * total_counts - sum of counts across all genes
#' * total_features_by_counts - total number of unique genes with expression >0
#' * qc_metric_min_library_size - did the cell have at least min_library_size
#' counts
#' * qc_metric_min_features - did the cell have counts >0 in at least
#' min_features number of cells?
#' * pc_mito - percentage of counts mapping to mitochondrial genes in this cell
#' * qc_metric_pc_mito_ok was pc_mito <= the max_mito cutoff?
#' * pc_ribo - percentage of counts mapping to ribosomal genes in this cell
#' * qc_metric_pc_ribo_ok was pc_ribo <= the max_ribo cutoff?
#' * passed_qc - did the cell pass all of the cell QC tests
#'
#' Gene-level annotations
#' * gene - official gene name
#' * gene_biotype - protein_coding, lncRNA, pseudogene, etc.
#' * qc_metric_ensembl_mapped - was the ensembl_gene_id found in biomaRt
#' * qc_metric_is_mito - is the gene mitochondrial
#' * qc_metric_is_ribo - is the gene ribosomal
#' * qc_metric_n_cells_expressing - number of cells with at least min_counts
#' * qc_metric_is_expressive - did at least min_cells have min_counts?
#'
#' @param sce a SingleCellExperiment object
#' @param min_library_size the minimum number of counts per cell
#' @param min_features the minimum number of features per cell (i.e. the minimum
#'   number of genes with >0 counts)
#' @param max_mito the maximum proportion of counts mapping to
#'   mitochondrial genes (0 - 1)
#' @param max_ribo the maximum proportion of counts mapping to
#'   ribosomal genes (0 - 1)
#' @param min_counts the minimum number of counts per cell in min_cells
#' @param min_cells the minimum number of cells with min_counts
#' @param ensembl_mapping_file a local tsv file with ensembl_gene_id and
#'   additional columns for mapping ensembl_gene_id to gene info.  If
#'   not provided, the biomaRt db is queried (slower).
#'
#' @return sce a annotated SingleCellExperiment object
#'
#' @family annotation functions
#' @import cli Matrix SummarizedExperiment dplyr SingleCellExperiment purrr
#' @export
annotate_sce <- function(sce,
                         min_library_size = 300,
                         min_features = 100,
                         max_mito = 0.10,
                         max_ribo = 1.00,
                         min_counts = 2,
                         min_cells = 2,
                         annotate_genes = TRUE,
                         annotate_cells = TRUE,
                         ensembl_mapping_file = NULL) {

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  before_coldata_colnames <- colnames(SummarizedExperiment::colData(sce))
  before_rowdata_colnames <- colnames(SummarizedExperiment::rowData(sce))

  # add the qc parameters to the metadata
  qc_params <- setdiff(names(formals(annotate_sce)),
                       c("sce", "ensembl_mapping_file",
                         "annotate_genes", "annotate_cells")) #not these args

  metadata(sce)[["qc_params"]] <- purrr::map(qc_params, ~ get(.)) %>%
    purrr::set_names(qc_params)

  if (annotate_genes) {
    sce <- annotate_sce_genes(sce, ensembl_mapping_file)
  }
  if (annotate_cells) {
    sce <- annotate_sce_cells(sce,
                              min_library_size = min_library_size,
                              min_features = min_features,
                              max_mito = max_mito,
                              max_ribo = max_ribo,
                              min_counts = min_counts,
                              min_cells = min_cells
    )
  } else {
    if(!annotate_genes){
      stop(cli::cli_alert_danger("Nothing to do. Specify gene/cell/both."))
    }
  }

  if(annotate_cells){
    cli::cli_alert_success(
      "SingleCellExperiment cells were successfully annotated with: \r\n"
    )

    cli::cli_ul(setdiff(
      colnames(SummarizedExperiment::colData(sce)),
      before_coldata_colnames
    ))
  }

  if(annotate_genes) {
    cli::cli_alert_success(
      "SingleCellExperiment genes were successfully annotated with: \r\n"
    )

    cli::cli_ul(setdiff(
      colnames(SummarizedExperiment::rowData(sce)),
      before_rowdata_colnames
    ))
  }

  return(sce)

}

