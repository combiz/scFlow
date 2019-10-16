################################################################################
#' Return mappings for an ensembl gene id
#'
#' Returns a data frame with the provided ensembl gene ids and requested
#' mappings from biomaRt or a local file.  If the mappings are absent from the
#' local file, or the local file is not provided or missing, the mappings are
#' retrieved by querying biomaRt.
#'
#' @param ensembl_ids vector of ensembl_ids, if versioned the version will be
#'   stripped.
#' @param mappings the biomaRt attributes to be mapped or mappings file column
#'   names.  This should always include ensembl_gene_id.
#' @param mappings_filepath path to the mappings tsv file
#'
#' @return mapped_df a data.frame of the provided ensembl_id's with mappings.
#'
#' @examples
#' map_ensembl_gene_id("ENSG00000130707")
#'
#' @family annotation functions
#'
#' @import cli purrr biomaRt stringr
#' @importFrom utils read.delim read.csv
#'
#' @export

map_ensembl_gene_id <- function(ensembl_ids,
                                mappings = c("ensembl_gene_id",
                                             "gene_biotype",
                                             "external_gene_name"),
                                mappings_filepath = NULL) {

  # strip the . number suffix with version number
  ensembl_ids <- stringr::str_replace(ensembl_ids,
                                    pattern = "\\..*",
                                    replacement = "")

  if (!(is.null(mappings_filepath))) {
    if (!file.exists(mappings_filepath)) {
      stop(cli::cli_alert_danger(
        "File not found.  Specify a valid path.")
      )
    } else { # file provided and is found
      ensembl_mappings <- utils::read.delim(mappings_filepath)
      check_mappings_are_present <- purrr::map_lgl(
        mappings,
        function(x) x %in% colnames(ensembl_mappings)
      )
      if (!(Reduce("&", check_mappings_are_present))) {
        stop(cli::cli_alert_danger(
          "Mappings file is missing the requested mappings.")
        )
      } else { # mappings requested are present in the file
        if (!(Reduce("&", ensembl_ids %in% ensembl_mappings$ensembl_gene_id))) {
          stop(cli::cli_alert_danger(
            "The requested ensembl_gene_ids are absent from the mappings file.")
          )
        } else { # all ok -- mappings and genes are present
          #colnames(ensembl_mappings) <- c("ensembl_gene_id",
          #"gene_biotype", "external_gene_name")
          mapped_df <- ensembl_mappings[
            ensembl_mappings$ensembl_gene_id %in% ensembl_ids]
        }
      }
    }
  } else { # mappings file not provided, use biomaRt

    cli::cli_alert_info("Mappings file not provided, using biomaRt (slower).")

    ensembl <- biomaRt::useEnsembl(
      biomart = "ensembl",
      dataset = "hsapiens_gene_ensembl"
    )

    mapped_df <- biomaRt::getBM(
      attributes = mappings,
      filters = "ensembl_gene_id",
      #values = "",
      values = ensembl_ids,
      mart = ensembl
    )

  }

  if (dim(mapped_df)[[1]] != length(ensembl_ids)) {
    warning("Warning")
    cli::cli_alert_warning("Not all ensembl_ids were found.")
  }

  return(mapped_df)

}
