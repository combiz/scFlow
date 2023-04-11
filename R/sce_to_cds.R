################################################################################
#' Convert a SingleCellExperiment object o a CellDataSet (Monocle)
#'
#' @param sce a SingleCellExperiment object
#'
#' @return cds a CellDataSet object
#'
#' @family internal helper
#'
#' @importFrom SingleCellExperiment counts reducedDims
#' @importFrom SummarizedExperiment rowData
#' @importFrom monocle3 new_cell_data_set
#' @keywords internal
.sce_to_cds <- function(sce) {
  phenoData <- data.frame(
    SummarizedExperiment::colData(sce),
    stringsAsFactors = FALSE
  )

  featureData <- data.frame(
    gene_short_name = SummarizedExperiment::rowData(sce)$gene,
    ensembl_gene_id = SummarizedExperiment::rowData(sce)$ensembl_gene_id,
    stringsAsFactors = FALSE
  )

  rownames(featureData) <- rownames(sce)

  cds <- monocle3::new_cell_data_set(
    expression_data = SingleCellExperiment::counts(sce),
    cell_metadata = phenoData,
    gene_metadata = featureData
  )

  SingleCellExperiment::reducedDims(cds) <-
    SingleCellExperiment::reducedDims(sce)

  return(cds)
}
