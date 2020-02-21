################################################################################
#' Functional enrichment analysis using enrichR
#'
#' Performs impacted pathway analysis with a list of genes.
#'
#' @param gene_file A data frame or the path of a .tsv file containing
#' a list of genes, their fold-change, p-value and adjusted p-value.
#' Column names should be gene, logFC, pval and padj respectively.
#' @param enrichment_database Name of the database for enrichment. User can
#' specify one or more database names from [enrichR::listEnrichDbs()].
#' @param is_output If TRUE a folder will be created and results of enrichment
#' analysis will be saved otherwise a R list will be returned. Default FALSE
#' @param plot_output If TRUE plots will be saved as .png file otherwise
#' returned as a R list. Default FALSE.
#' @param output_dir Path for the output directory. Default is current dir.
#'
#' @return enrichment_result a list of data.frames containing enrichment output
#' and a list of plots of top 10 significant genesets.
#'
#' @family Impacted pathway analysis
#' @importFrom cli cli_alert_danger rule cli_alert_info cli_alert
#' @importFrom ggplot2 ggplot ggsave
#' @importFrom cowplot theme_cowplot background_grid
#' @importFrom stringr str_wrap
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @examples
#' set.seed(1234)
#' enrichment_result <- pathway_analysis_enrichr(
#'   gene_file = paste(system.file("extdata", package = "scFlowData"), "/",
#'     "de_result_table.tsv",
#'     sep = ""
#'   ),
#'   enrichment_database = "GO_Molecular_Function_2018",
#'   is_output = FALSE
#' )
pathway_analysis_enrichr <- function(gene_file = NULL,
                                     reference_file = NULL,
                                     enrichment_database = c(
                                       "GO_Molecular_Function_2018",
                                       "GO_Cellular_Component_2018",
                                       "GO_Biological_Process_2018",
                                       "MGI_Mammalian_Phenotype_2017",
                                       "MGI_Mammalian_Phenotype_Level_3",
                                       "MGI_Mammalian_Phenotype_Level_4",
                                       "ChEA_2016",
                                       "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
                                       "KEGG_2019_Human",
                                       "WikiPathways_2019_Human",
                                       "Reactome_2016",
                                       "Allen_Brain_Atlas_down",
                                       "Allen_Brain_Atlas_up"),
                                        is_output = FALSE,
                                        output_dir = ".") {

  dbs <- enrichR::listEnrichrDbs()
  assertthat::assert_that(
    all(enrichment_database %in% as.character(dbs$libraryName)),
    msg = "Invalid databases specified.  See enrichR::listEnrichDbs()."
    )

  assertthat::assert_that(
    !is.null(gene_file),
    msg = "No input gene list found!")

  if (is.data.frame(gene_file)) {
    interest_gene <- gene_file
  } else {
    assertthat::assert_that(file.exists(gene_file), msg = "File not found.")
    cli::cli_text("Reading: {.file {gene_file}}")
    interest_gene <- read.delim(gene_file, sep = "\t")
  }

  expected_cols <- c("gene")
  assertthat::assert_that(
    all(expected_cols %in% colnames(gene_file)),
    msg = sprintf(
      "One or more expected columns missing: %s",
      paste0(expected_cols, collapse = ","))
    )

  enrichr_res <- enrichR::enrichr(
    as.character(gene_file$gene),
    enrichment_database
    )

  enrichr_res <- purrr::map(enrichr_res, function(x) {
    x %>%
      dplyr::filter(Adjusted.P.value <= 0.05) %>%
      dplyr::mutate('-log(p-value)' = -log(Adjusted.P.value))})

  enrichr_res <- enrichr_res[lapply(enrichr_res, nrow) > 0]

  return(enrichr_res)


}
