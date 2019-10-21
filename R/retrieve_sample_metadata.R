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
#' @importFrom tools toTitleCase
#' @importFrom english words
#' @export
retrieve_sample_metadata <- function(unique_id,
                                     id_colname,
                                     samplesheet_path,
                                     ...) {

  args_l <- list(colClasses = NA)
  inargs <- list(...)
  args_l[names(inargs)] <- inargs

  cat(cli::rule("Retrieving sample metadata", line = 1), "\r\n")

  if (!file.exists(samplesheet_path)) {
    stop(cli::cli_alert_danger("Samplesheet not found.\
                                Check that the path is correct."))
  } else {
    cat("Reading", cli::col_green(c(samplesheet_path, " \r\n")))
    samplesheet <- read.delim(
      samplesheet_path,
      colClasses = args_l$colClasses)
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
    stop(cli::cli_alert_danger("Something went wrong, more than one match."))
  } else {
    cli::cli_alert_success(sprintf("%s metadata variables loaded for %s='%s'",
                           tools::toTitleCase(
                             english::words(dim(metadata)[[2]])),
                           id_colname,
                           unique_id)
    )
    #cli::cli_text(paste(colnames(metadata), collapse = ", "))
    for (var in names(metadata)){
      cli::cli_text(c("{.strong {var}}: {.var {metadata[[var]][[1]]}} ",
        "{.emph ({class(metadata[[var]][[1]])}})")
      )
    }
  }

  return(metadata)

}
