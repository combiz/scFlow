################################################################################
#' Add basic cell-wise annotations for a SingleCellExperiment
#'
#' Calculates the following QC metrics:
#' * qc_metric_total_counts - was the ensembl_gene_id found in biomaRt
#' * qc_metric_total_features_by_counts - is the gene mitochondrial
#'
#' @param sce a SingleCellExperiment object
#'
#' @return sce a SingleCellExperiment object annotated with cell data
#'
#' @family annotation functions
#' @import cli Matrix SingleCellExperiment dplyr
#' @export
annotate_sce_cells <- function(sce, ls_threshold = 300, tf_threshold = 100, mt_threshold = .10) {

  if(is.null(colData(sce)$total_counts)){
    sce$total_counts <- colSums(counts(sce))
  }
  if(is.null(colData(sce)$total_features_by_counts)){
    sce$total_features_by_counts <- colSums(counts(sce) > 0)
  }
  # doublet and triplet filter
  #upper_cutoff <- mean((sce$total_counts)) + 2*sd((sce$total_counts))
  #keep.cells.singlets <- sce$total_counts <= upper_cutoff
  #plot_and_tabulate_singlet_filter(sce, upper_cutoff, QCPLOTDIR)

  # library size filter (total counts)
  keep.cells.librarysize <- sce$total_counts >= ls_threshold
  # cells with low numbers of expressed genes filter
  keep.cells.expressive <- sce$total_features_by_counts > tf_threshold
  # annotate mito and filter
  MT.genes <- startsWith(rowData(sce)$gene, "MT-")+0
  MT.genes[is.na(MT.genes)] <- 0 # for deprecated ensembl_ids bug
  sce$pc_mito <- colSums(counts(sce[MT.genes == 1,])) / sce$total_counts
  keep.cells.not.excess.MT.genes <- sce$pc_mito < mt_threshold

  # logical AND
  ## define a function
  AND1 <- function (...)  Reduce("&", list(...))
  sce$keep <- AND1(keep.cells.librarysize,
                   keep.cells.expressive,
                   #keep.cells.singlets,
                   keep.cells.not.excess.MT.genes)

  return(sce)

}
