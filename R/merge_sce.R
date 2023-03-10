################################################################################
#' Merge Multiple SingleCellExperiment Objects
#'
#' @param sce_l a list of SingleCellExperiment objects, or, folder_paths to
#' SingleCellExperiment objects to be read in with \code{\link{read_sce}}
#' @param ensembl_mapping_file path to the mappings tsv file
#' @param species human or mouse
#'
#' @return sce a annotated SingleCellExperiment object
#'
#' @family annotation functions
#' @rawNamespace import(scater, except = "normalize")
#' @importFrom cli cli_alert_danger rule cli_text cli_alert_success
#' @importFrom purrr map map_lgl map_chr
#' @importFrom SingleCellExperiment counts
#' @importFrom Matrix colSums rowSums
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom dplyr left_join union rename
#' @importFrom english words
#' @importFrom tools toTitleCase
#' @export
#'
merge_sce <- function(sce_l, ensembl_mapping_file = NULL,
                      species = getOption(
                        "scflow_species",
                        default = "human")) {

  cat(cli::rule("Merging SingleCellExperiments", line = 2), "\r\n")

  if (!class(sce_l) %in% c("list", "character")) {
    stop(cli::cli_alert_danger(
      "Expected a list or vector, received a {.emph {class(sce_l)}}.")
    )
  }
  if (all((purrr::map_lgl(sce_l, is.character)))) {
    if (all(purrr::map_lgl(sce_l, dir.exists))) {
      sce_l <- purrr::map(sce_l, read_sce)
    } else {
      stop(cli::cli_alert_danger(
        "Path(s) not found.")
      )
    }
  } else if (!all(purrr::map_lgl(sce_l,
                                 ~ class(.) == "SingleCellExperiment"))) {
    stop(cli::cli_alert_danger(
      "Provide a list of SingleCellExperiments or folder paths.")
    )
  }

  cat(cli::rule("Starting merge", line = 1), "\r\n")
  cli::cli_text("Merging {.var {length(sce_l)}} SingleCellExperiments")
  ensembl_gene_id_l <- purrr::map(sce_l,
                                  ~ as.character(rowData(.)$ensembl_gene_id))
  n_sce_genes <- paste0(purrr::map_int(ensembl_gene_id_l, length), collapse = ", ")
  cli::cli_text("Genes in each SingleCellExperiment: {.var {n_sce_genes}}")
  union_ensembl_gene_id <- as.character(Reduce(dplyr::union, ensembl_gene_id_l))
  cli::cli_text(c("Union of genes across all SingleCellExperiments: ",
                  "{.var {length(union_ensembl_gene_id)}}"))
  cli::cli_text("Equalizing to {.var {length(union_ensembl_gene_id)}} genes.")
  sce_expanded_l <-
    purrr::map(sce_l,
               ~ .expand_sce_rows(., union_ensembl_gene_id))

  sce <- Reduce(cbind, sce_expanded_l)

  # drop previous cell annotations
  keep_idx <- !(startsWith(
    colnames(SummarizedExperiment::colData(sce)),
    "qc_metric"))

  new_rd <- map_ensembl_gene_id(
    names(sce),
    ensembl_mapping_file = ensembl_mapping_file) %>%
    dplyr::rename(gene = external_gene_name)

  SummarizedExperiment::rowData(sce) <- dplyr::left_join(
    as.data.frame(SummarizedExperiment::rowData(sce)),
    new_rd, by = "ensembl_gene_id"
  )

  SummarizedExperiment::colData(sce) <-
    SummarizedExperiment::colData(sce)[, keep_idx]

  # recalculate cell metrics
  sce$total_counts <- Matrix::colSums(SingleCellExperiment::counts(sce))

  sce$total_features_by_counts <- Matrix::colSums(
    SingleCellExperiment::counts(sce) > 0
  )

  # recalculate gene metrics
  SummarizedExperiment::rowData(sce)$total_counts <-
    Matrix::rowSums(SingleCellExperiment::counts(sce))

  SummarizedExperiment::rowData(sce)$n_cells_by_counts <-
    Matrix::rowSums(SingleCellExperiment::counts(sce) > 0)

  n_merged <- tools::toTitleCase(english::words(length(sce_l)))
  cli::cli_alert_success("{.val {n_merged}} SingleCellExperiment were merged.")

  sce@metadata$scflow_steps <- list()
  sce@metadata$scflow_steps$merged <- 1

  return(sce)

}

#' Helper function to expand a SingleCellExperiment to include specified genes
#'
#' The function identifies the missing ensembl_gene_id's (if any), and
#' generates new zero'd matrix rows, rbinds them to the matrix, then sorts
#' the matrix by ensembl_gene_id.  rowData (except ensembl_gene_id) and metadata
#' are dropped to facilitate subsequent cbind-ing of SingleCellExperiments.
#'
#' @param sce a SingleCellExperiment
#' @param union_ensembl_gene_id the full list of ensembl_gene_id to be included
#' in the new expanded SingleCellExperiment
#'
#' @importFrom SummarizedExperiment rowData colData
#' @return sce a annotated SingleCellExperiment object with any
#' @keywords internal
.expand_sce_rows <- function(sce, union_ensembl_gene_id) {

  missing_ids <- union_ensembl_gene_id[
    !union_ensembl_gene_id %in% rownames(counts(sce))]

  if (length(missing_ids) > 0) {
    zeros_mat <- matrix(
      data = rep(0, length(missing_ids) * dim(counts(sce))[[2]]),
      nrow = length(missing_ids),
      dimnames = list(missing_ids)
    )
    mat <- rbind(SingleCellExperiment::counts(sce), zeros_mat)
    mat <- mat[order(rownames(mat)), ]
    mat <- as(mat, "TsparseMatrix")

    newsce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = mat),
      colData = data.frame(sce@colData),
      rowData = data.frame(ensembl_gene_id = rownames(mat))
    )

    SingleCellExperiment::reducedDims(newsce) <-
      SingleCellExperiment::reducedDims(sce)

  } else {
    newsce <- sce
    SingleCellExperiment::counts(sce) <-
      as(SingleCellExperiment::counts(sce), "TsparseMatrix")
  }

  SummarizedExperiment::rowData(newsce)$ensembl_gene_id <-
    as.character(SummarizedExperiment::rowData(newsce)$ensembl_gene_id)

  newsce@metadata <- list() # drop for merge

  return(newsce)

}
