################################################################################
#' Read the metadata for a sample from a samplesheet
#'
#' Retrieves a single row from a samplesheet table where the provided
#' `unique_key` parameter is matched in the `key_colname` column of the
#' samplesheet.
#'
#' @param unique_key unique value for this sample in the `key_colname` column
#' @param key_colname the column name containing the `unique_key`
#' @param samplesheet_path path to a tsv file containing a samplesheet
#' @param ... for col_classes (may be extended later - fix)
#'
#' @return metadata the metadata for the specified sample#'
#'
#' @family annotation functions
#' @importFrom cli cli_alert_danger col_green cli_alert_danger cli_text cli_h1
#' @importFrom tools toTitleCase
#' @importFrom english words
#' @export
read_metadata <- function(unique_key,
                              key_colname,
                              samplesheet_path,
                              ...) {

  args <- list(col_classes = NA)
  inargs <- list(...)
  args[names(inargs)] <- inargs

  cli::cli_h1("Reading sample metadata")

  if (!file.exists(samplesheet_path)) {
    stop(cli::cli_alert_danger("Samplesheet not found.\
                                Check that the path is correct."))
  } else {
    cat("Reading", cli::col_green(c(samplesheet_path, " \r\n")))
    samplesheet <- read.delim(
      samplesheet_path,
      colClasses = args$col_classes,
      stringsAsFactors = TRUE
      )
  }

  if (!(key_colname %in% colnames(samplesheet))) {
    stop(cli::cli_alert_danger(sprintf(
    "Column '%s' is not present in the samplesheet.", key_colname))
    )
  }

  if (!(unique_key %in% samplesheet[[key_colname]])) {
    stop(cli::cli_alert_danger(sprintf(
      "unique_key '%s' is not present in column %s of the samplesheet.",
      unique_key,
      key_colname))
    )
  }

  if (sum(unique_key %in% samplesheet[[key_colname]]) > 1) {
    stop(cli::cli_alert_danger("unique_key is not unique."))
  }

  metadata <- samplesheet[samplesheet[[key_colname]] == unique_key, ]

  if (dim(metadata)[[1]] != 1) {
    stop(cli::cli_alert_danger("Something went wrong, more than one match."))
  } else {
    cli::cli_alert_success(sprintf("%s metadata variables loaded for %s='%s'",
                           tools::toTitleCase(
                             english::words(dim(metadata)[[2]])),
                           key_colname,
                           unique_key)
    )
    for (var in names(metadata)) {
      if(!is.na(metadata[[var]][[1]]) & metadata[[var]][[1]] != "") {
        cli::cli_text(
          c("{.strong {var}}: {.envvar {as.character(metadata[[var]][[1]])}} ",
          "{.emph ({class(metadata[[var]][[1]])}})")
          )
      } else {
        metadata[[var]][[1]] <- NA
        #addNA(metadata[[var]])
        cli::cli_text(c("{.strong {var}}: {.envvar NA }",
                        "{.emph ({class(metadata[[var]][[1]])}})"))
      }
    }
  }

  return(metadata)

}
