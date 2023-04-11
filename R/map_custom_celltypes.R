################################################################################
#' Map Custom Celltype Annotations from a DataFrame
#'
#' @param sce a SingleCellExperiment
#' @param mappings a dataframe with a "cluster" column and additional columns with
#' celltype annotation data
#' @param cols specifies the subset of columns to annotate with.
#' defaults to NULL or all columns
#' @param clusters_colname the name of the colData column with cluster number
#'
#' @return sce a SingleCellExperiment object with custom annotations
#' @family clustering and dimensionality reduction
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom assertthat assert_that
#' @importFrom cli cli_text cli_alert_success
#' @importFrom purrr map_chr
#' @export
map_custom_celltypes <- function(sce,
                                 mappings,
                                 cols = NULL,
                                 clusters_colname = "clusters") {

  assertthat::assert_that(
    all(clusters_colname %in% names(SummarizedExperiment::colData(sce))))

  assertthat::assert_that(all(
    clusters_colname %in% colnames(mappings))
  )

  if (!is.null(cols)) {
    assertthat::assert_that(all(
      cols %in% colnames(mappings)), msg = "Invalid cols specified.")
  } else {
    cols <- names(mappings)
  }

  cols <- cols[cols != clusters_colname]

  for (var in cols) {
    mappings_lookup <- mappings[[var]]
    names(mappings_lookup) <- mappings[[clusters_colname]]

    sce[[var]] <- purrr::map_chr(
      as.character(sce[[clusters_colname]]), ~ mappings_lookup[.])

    #sce[[var]] <- as.factor(as.character(sce[[var]]))

    cli::cli_text("Appended {.var {var}} to SingleCellExperiment")
  }

  return(sce)
}
