################################################################################
#' Retrieve the metadata for a sample from a samplesheet
#'
#' Retrieves a single row from a samplesheet table where the provided
#' `unique_id` parameter is matched in the `id_colname` column of the
#' samplesheet.
#'
#' @param unique_id unique value for this sample in the `id_colname` column
#' @param id_colname the column name containing the `unique_id`
#' @param samplesheet_path path to a tsv file containing a samplesheet
#'
#' @return metadata the metadata for the specified sample
#'
#'
#' @family annotation functions
#' @import cli
#' @export
retrieve_sample_metadata <- function(unique_id,
                                     id_colname,
                                     samplesheet_path) {

  if (!file.exists(samplesheet_path)) {
    stop(cli::cli_alert_danger("Samplesheet not found.\
                                Check that the path is correct."))
  } else {
    samplesheet <- read.delim(samplesheet_path)
  }

  if (!(id_colname %in% colnames(samplesheet))) {
    stop(cli::cli_alert_danger(sprintf(
    "Column '%s' is not present in the samplesheet.", id_colname))
    )
  }

  if (!(unique_id %in% samplesheet[[id_colname]])) {
    stop(cli::cli_alert_danger(sprintf(
      "unique_id '%s' is not present in column %s of the samplesheet.",
      unique_id,
      id_colname))
    )
  }

  if (sum(unique_id %in% samplesheet[[id_colname]]) > 1) {
    stop(cli::cli_alert_danger("unique_id is not unique."))
  }

  metadata <- samplesheet[samplesheet[[id_colname]] == unique_id, ]

  if (dim(metadata)[[1]] != 1) {
    stop(cli::cli_alert_danger("Something went wrong, more than one row."))
  } else {
    cli::cli_alert_success(sprintf("%s metadata variables loaded for %s",
                           dim(metadata)[[2]],
                           unique_id))
    cli::cli_text(paste(colnames(metadata), collapse = ", "))
  }

  return(metadata)

}
