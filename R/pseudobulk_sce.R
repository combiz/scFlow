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
#' @importFrom cli cli_text rule
#' @importFrom SummarizedExperiment colData rowData assays
#' @importFrom dplyr select distinct group_by tally select
#' @importFrom magrittr %>% set_colnames
#' @importFrom stringr str_split
#' @importFrom purrr map_int map_df
#' @importFrom assertthat assert_that
#' @importFrom scater librarySizeFactors normalize calculateQCMetrics
#' @importFrom SingleCellExperiment SingleCellExperiment
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
    dplyr::select(keep_vars) %>%
    dplyr::distinct()

  # generate individual/cluster factor on which to pseudobulk
  sce$pseudobulk_id <- paste(
    sce[[sample_var]],
    sce[[celltype_var]],
    sep = "_"
  )

  pb_matrix <- scater::sumCountsAcrossCells(
    sce,
    sce$pseudobulk_id,
    exprs_values = assay_name
  )

  rownames(pb_matrix) <- SummarizedExperiment::rowData(sce)$ensembl_gene_id

  # create lookup table for cellnumbers by individual/cluster
  n_lookup <- data.frame(SummarizedExperiment::colData(sce)) %>%
    dplyr::group_by(pseudobulk_id) %>%
    dplyr::tally() %>%
    dplyr::group_by(pseudobulk_id) %>%
    dplyr::select(n, pseudobulk_id)

  lup <- n_lookup$n
  names(lup) <- n_lookup$pseudobulk_id

  # rebuild basic cell identifiers
  pb_cd <- purrr::map_df(stringr::str_split(colnames(pb_matrix), "_"),
                  ~ magrittr::set_colnames(
                    data.frame(.[[1]], .[[2]]),
                    c(sample_var, celltype_var))
                  )
  pb_cd$pseudobulk_id <- colnames(pb_matrix)
  pb_cd <- as.data.frame(unclass(pb_cd)) # chr to factor trick

  # append cell numbers to each pseudobulk sample
  pb_cd$n_cells <- purrr::map_int(pb_cd$pseudobulk_id, ~ lup[.])
  pb_cd[[sample_var]] <- as.factor(pb_cd[[sample_var]])

  n_samples <- nrow(pb_cd)

  # append the rest of the sample information
  pb_cd <- merge(pb_cd, sce_cd,
                 by = sample_var,
                 all.x = TRUE, all.y = FALSE)

  assertthat::assert_that(
    nrow(pb_cd) == n_samples,
    msg = c("keep_vars variable(s) specified that are not ",
            "common across all cells of sample_var."))

  # basic rowdata
  pb_rd <- data.frame(SummarizedExperiment::rowData(sce)) %>%
    dplyr::select(ensembl_gene_id, gene)

  pb_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = pb_matrix),
    colData = pb_cd,
    rowData = pb_rd
  )

  names(SummarizedExperiment::assays(pb_sce)) <- assay_name

  # size factors are set to the number of cells
  pb_sce@int_colData$size_factor <- scater::librarySizeFactors(pb_sce)
  pb_sce <- scater::normalize(pb_sce, return_log = FALSE)
  pb_sce <- scater::calculateQCMetrics(pb_sce)

  return(pb_sce)
}