################################################################################
#' Return mappings for Ensembl Gene IDs
#'
#' Returns a data frame with the provided ensembl gene ids and requested
#' mappings from biomaRt or a local file.  If the mappings are absent from the
#' local file, or the local file is not provided or missing, the mappings are
#' retrieved by querying biomaRt.
#'
#' @param ensembl_ids vector of ensembl_ids, if versioned the version will be
#'   stripped.
#' @param mappings the biomaRt attributes to be mapped or mappings file column
#'   names.  ensembl_gene_id will be added if not included here.
#'   `c("ensembl_gene_id", "external_gene_name")`
#' @param species ensembl_ids are mapped to `human` or `mouse`
#' @param ensembl_mapping_file path to the mappings tsv file
#'
#' @return mapped_df a data.frame of the provided ensembl_id's with mappings.
#'
#' @family annotation functions
#'
#' @importFrom cli cli_alert_danger col_green cli_alert_info
#' @importFrom stringr str_replace
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom purrr map_lgl map_if
#' @importFrom utils read.delim read.csv
#'
#' @export

map_ensembl_gene_id <- function(ensembl_ids,
                                mappings = c("ensembl_gene_id",
                                             "gene_biotype",
                                             "external_gene_name",
                                             "percentage_gene_gc_content"),
                                species = getOption(
                                  "scflow_species",
                                  default = "human"),
                                ensembl_mapping_file = NULL) {

  if (!species %in% c("human", "mouse")) {
    stop("only human and mouse currently supported.")
  }

  # include the id for biomart
  mappings <- unique(c("ensembl_gene_id", mappings))

  # strip the . number suffix with version number
  ensembl_ids <- stringr::str_replace(ensembl_ids,
                                    pattern = "\\..*",
                                    replacement = "")

  if (!(is.null(ensembl_mapping_file))) {
    if (!file.exists(ensembl_mapping_file)) {
      stop(cli::cli_alert_danger(
        "File not found.  Specify a valid path.")
      )
    } else { # file provided and is found
        cat("Reading", cli::col_green(c(ensembl_mapping_file, " \r\n")))
        ensembl_mappings <- utils::read.delim(
          ensembl_mapping_file,
          stringsAsFactors = FALSE)

        check_mappings_are_present <- purrr::map_lgl(
          mappings,
          function(x) x %in% colnames(ensembl_mappings)
        )

        if (!all(check_mappings_are_present)) {
          stop(cli::cli_alert_danger(
            "Mappings file is missing the requested mappings.")
          )
        } else { # mappings requested are present in the file

          idx <- match(ensembl_ids, ensembl_mappings$ensembl_gene_id)

          mapped_df <- ensembl_mappings[idx, ]
        }
      }
  } else { # mappings file not provided, use biomaRt

    cli::cli_alert_info("Mappings file not provided, using biomaRt (slower).")

    species_dataset <- c(
      "human" = "hsapiens_gene_ensembl",
      "mouse" = "musculus_gene_ensembl"
    )

    ensembl <- biomaRt::useEnsembl(
      biomart = "ensembl",
      dataset = species_dataset[[species]]
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
    cli::cli_alert_info(c(
      "Note: not all ensembl_ids were found (",
      (length(ensembl_ids) - dim(mapped_df)[[1]]), "/",
      length(ensembl_ids), " not found).")
    )
  }

  mapped_df <- mapped_df %>%
    purrr::map_if(is.factor, as.character) %>%
    dplyr::as_tibble()

  return(mapped_df)

}
