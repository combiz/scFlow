################################################################################
#' Merge Multiple SingleCellExperiment Objects
#'
#' @param sce_l a list of SingleCellExperiment objects, or, folder_paths to
#' SingleCellExperiment objects to be read in with \code{\link{read_sce}}
#' @param ensembl_mapping_file path to the mappings tsv file
#'
#' @return sce a annotated SingleCellExperiment object
#'
#' @family annotation functions
#' @import cli Matrix dplyr SingleCellExperiment purrr
#' @importFrom SummarizedExperiment rowData colData
#' @export
#'
merge_sce <- function(sce_l, ensembl_mapping_file = NULL) {

  cat(cli::rule("Merging SingleCellExperiments", line = 2), "\r\n")

  if (!class(sce_l) %in% c("list", "character")) {
    stop(cli::cli_alert_danger(
      "Expected a list or vector, received a {.emph {class(sce_l)}}.")
    )
  }
  if (all((purrr::map_lgl(sce_l, is.character)))) {
    if (all(purrr::map_lgl(sce_l, dir.exists))) {
      sce_l <- map(sce_l, read_sce)
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
  n_sce_genes <- paste0(map_int(ensembl_gene_id_l, length), collapse = ", ")
  cli::cli_text("Genes in each SingleCellExperiment: {.var {n_sce_genes}}")
  union_ensembl_gene_id <- as.character(Reduce(dplyr::union, ensembl_gene_id_l))
  cli::cli_text(c("Union of genes across all SingleCellExperiments: ",
                  "{.var {length(union_ensembl_gene_id)}}"))
  cli::cli_text("Equalizing to {.var {length(union_ensembl_gene_id)}} genes.")
  sce_expanded_l <-
    purrr::map(sce_l,
               ~ .expand_sce_rows(., union_ensembl_gene_id))

  sce <- Reduce(cbind, sce_expanded_l)

  sce <- annotate_sce(
    sce,
    annotate_cells = FALSE,
    ensembl_mapping_file = ensembl_mapping_file)

  n_merged <- tools::toTitleCase(english::words(length(sce_l)))
  cli::cli_alert_success("{.val {n_merged}} SingleCellExperiment were merged.")

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
#' @import purrr SingleCellExperiment
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
    mat <- rbind(counts(sce), zeros_mat)
    mat <- mat[order(rownames(mat)), ]

    newsce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = mat),
      colData = data.frame(sce@colData),
      rowData = data.frame(ensembl_gene_id = rownames(mat))
    )

    SingleCellExperiment::reducedDims(newsce) <-
      SingleCellExperiment::reducedDims(sce)

  } else {
    newsce <- sce
  }

  SummarizedExperiment::rowData(newsce)$ensembl_gene_id <-
    as.character(SummarizedExperiment::rowData(newsce)$ensembl_gene_id)

  metadata(newsce) <- list() # drop for merge

  return(newsce)

}
