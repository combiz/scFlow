#' Append a citation to a SingleCellExperiment
#'
#' @param sce a SingleCellExperiment object
#' @param key_var the key variable (e.g. "MENDELEY.TAGS")
#' @param key the value(s) for the key_var (e.g. "doubletfinder")
#'
#' @return sce a SingleCellExperiment object
#'
#' @family citation functions
#' @importFrom bib2df bib2df
#' @importFrom cli cli_text
#' @importFrom purrr map_lgl
#' @importFrom assertthat assert_that
#' @keywords internal
.append_citation_sce <- function(sce,
                                 key_var = "MENDELEY.TAGS",
                                 key,
                                 verbose = TRUE) {

  bib <- bib2df::bib2df(
    system.file("bibtex/references.bib", package = "scFlow")
    )

  assertthat::assert_that(
    all(purrr::map_lgl(key, ~ . %in% bib[[key_var]])),
    msg = "Invalid citation.")

  citations <- bib[bib[[key_var]] %in% key, ]
  sce@metadata$citations <- unique(rbind(sce@metadata$citations, citations))

  if (verbose) {
    cli::cli_text("Please consider citing: -")
    apply(citations, 1, .print_citation)
  }

  return(sce)

}

#' Print a single citation to stdout
#'
#' @param citation a single row from a bib2df tibble
#'
#' @family citation functions
#' @importFrom cli cli_text symbol
#' @keywords internal
.print_citation <- function(citation) {

  author_text <- paste0(unlist(citation$AUTHOR), collapse = ", ")

  cli::cli_text(c(
    "{cli::symbol$square_small_filled} ",
    "{author_text} ",
    "({citation$YEAR}). ",
    "{.strong {citation$TITLE}}. ",
    "{citation$JOURNAL}. ",
    "{citation$VOLUME}({citation$NUMBER}), ",
    "{citation$PAGES}. ",
    "{.emph {citation$DOI}}")
  )

}

#' Print a single citation to stdout as html
#'
#' @param citation a single row from a bib2df tibble
#'
#' @family citation functions
#' @importFrom cli symbol
#' @keywords internal
.print_citation_html <- function(citation) {

  author_text <- paste0(unlist(citation$AUTHOR), collapse = ", ")

  cat(paste0(
    "<p>",
    sprintf("%s ", author_text),
    sprintf("(%s). ", citation$YEAR),
    sprintf("<b>%s</b>. ", citation$TITLE),
    sprintf("%s. ", citation$JOURNAL),
    sprintf("%s(%s), ", citation$VOLUME, citation$NUMBER),
    sprintf("%s. ", citation$PAGES),
    sprintf("<a href=\"http://dx.doi.org/%s\"><strong>[DOI]</strong></a>",
            citation$DOI),
    "</p>"))

}
