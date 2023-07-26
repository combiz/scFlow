################################################################################
#' Generate a Pseudo-bulk SingleCellExperiment from a SingleCellExperiment
#'
#' Generate a pseudo-bulk of a SingleCellExperiment by summing counts of all
#' cells from the same `sample_var` and `celltype_var`.  The SCE is then
#' re-annotated and library sizes calculated based on cell numbers.
#'
#' @param sce a SingleCellExperiment object
#' @param sample_var the variable name identifying unique samples
#' @param celltype_var the variable name identifying cell types
#' @param assay_name the assay name for pseudobulking (e.g. "counts")
#' @param keep_vars metadata variables to keep from the SingleCellExperiment
#'
#' @return pb_sce a pseudobulk SingleCellExperiment object
#'
#' @family differential gene expression
#'
#' @rawNamespace import(scater, except = "normalize")
#' @importFrom cli cli_text rule
#' @importFrom SummarizedExperiment colData rowData assays
#' @importFrom dplyr select distinct group_by tally select left_join
#' @importFrom magrittr %>% set_colnames
#' @importFrom stringr str_split
#' @importFrom purrr map_int map_df
#' @importFrom future availableCores
#' @importFrom assertthat are_equal
#' @importFrom tidyselect all_of
#' @importFrom SingleCellExperiment SingleCellExperiment counts
#' @importFrom edgeR cpm
#'
#' @export
pseudobulk_sce <- function(sce,
                           sample_var = "individual",
                           celltype_var = "cluster_celltype",
                           assay_name = "counts",
                           keep_vars = c("individual", "group", "sex",
                                         "age", "PMI", "RIN", "seqdate")
                           ) {

  # coldata variables to be merged to pseudobulk "cells" later
  keep_vars <- unique(c(sample_var, keep_vars))

  sce_cd <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::select(tidyselect::all_of(keep_vars)) %>%
    dplyr::distinct()

  # generate individual/cluster factor on which to pseudobulk
  sce$pseudobulk_id <- as.factor(paste(
    sce[[sample_var]],
    sce[[celltype_var]],
    sep = "||"
  ))

  # fast matrix multiplication method
  mm <- model.matrix(~ 0 + sce$pseudobulk_id)
  pb_matrix <- SingleCellExperiment::counts(sce) %*% mm
  pb_matrix <- Matrix::Matrix(pb_matrix, sparse = TRUE)
  colnames(pb_matrix) <- levels(sce$pseudobulk_id)
  pb_cpm_matrix <- edgeR::cpm(pb_matrix, log = F)

  if("normcounts" %in% names(SummarizedExperiment::assays(sce))) {
    pb_normcounts_matrix <- SingleCellExperiment::normcounts(sce) %*% mm
    colnames(pb_normcounts_matrix) <- levels(sce$pseudobulk_id)
  }

  # create lookup table for cellnumbers by individual/cluster
  n_lookup <- data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::group_by(pseudobulk_id) %>%
    dplyr::tally() %>%
    dplyr::group_by(pseudobulk_id) %>%
    dplyr::select(n, pseudobulk_id)

  lup <- n_lookup$n
  names(lup) <- n_lookup$pseudobulk_id

  # rebuild basic cell identifiers
  pb_cd <- data.frame(
    Reduce(
      rbind, strsplit(colnames(pb_matrix), "||", fixed = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  pb_cd <- as.data.frame(unclass(pb_cd)) # chr to factor trick

  pb_cd <- magrittr::set_colnames(
    pb_cd,
    c(sample_var, celltype_var)
    )

  pb_cd$pseudobulk_id <- colnames(pb_matrix)
  # if sample_var equals celltype_var
  pb_cd <- pb_cd[, !duplicated(colnames(pb_cd))]

  # append cell numbers to each pseudobulk sample
  pb_cd$n_cells <- purrr::map_int(pb_cd$pseudobulk_id, ~ lup[.])
  pb_cd[[sample_var]] <- as.factor(pb_cd[[sample_var]])
  pb_cd <- pb_cd[order(colnames(pb_matrix)), ]

  n_samples <- nrow(pb_cd)

  pb_cd[[sample_var]] <- as.character(pb_cd[[sample_var]])
  sce_cd[[sample_var]] <- as.character(sce_cd[[sample_var]])

  # append the rest of the sample information
  pb_cd <- dplyr::left_join(
    pb_cd, sce_cd,
    by = sample_var)

  # ensure the order of the coldata matches the matrix
  pb_cd <- pb_cd[match(colnames(pb_matrix), pb_cd$pseudobulk_id),]

  assertthat::are_equal(
    nrow(pb_cd), n_samples,
    msg = c("keep_vars variable(s) specified that are not ",
            "common across all cells of sample_var."))

  # basic rowdata
  pb_rd <- data.frame(SummarizedExperiment::rowData(sce)) %>%
    dplyr::select(ensembl_gene_id, gene)

  assays_l <- list(counts = pb_matrix,
                   cpm = pb_cpm_matrix)
  if("normcounts" %in% names(SummarizedExperiment::assays(sce))) {
    assays_l[["normcounts"]] <- pb_normcounts_matrix
  }

  pb_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_l,
    colData = pb_cd,
    rowData = pb_rd
  )

  # size factors are set to the number of cells
  pb_sce <-scater::addPerCellQC(pb_sce)
  SingleCellExperiment::sizeFactors(pb_sce) <- pb_sce$n_cells

  pb_sce@metadata$scflow_steps <- list()
  pb_sce@metadata$scflow_steps$pseudobulk <- TRUE

  return(pb_sce)
}
