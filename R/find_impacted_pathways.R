################################################################################
#' Functional enrichment analysis
#'
#' Performs gene ontology and impacted pathway enrichment analysis with a list
#' of gene names and their fold-change.
#'
#' @param gene_file A data frame or the path of a .tsv file containing
#' a list of genes, their fold-change, p-value and adjusted p-value. Column
#' names should be gene, logFC, pval and padj respectively.
#' @param reference_file A data frame containing all the genes that were used
#' as input for differential expression. Column name should be gene.
#' If not provided the human protein-coding genome will be used as background
#' genes.
#' @param enrichment_tool Enrichment tool to use. WebGestaltR, ROntoTools and
#' EnrichR is implemented.
#' @param enrichment_database Name of the database for enrichment.
#' If NULL then multiple databases will be used or user can specify one or more
#' database names.Default NULL.
#' @param is_output If TRUE a folder will be created and results of enrichment
#' analysis will be saved otherwise a R list will be returned. Default FALSE
#' @param output_dir Path for the output directory. Default is current
#' directory.
#'
#' @return enrichment_result a list of list containing enrichment output from
#' different enrichment tools and a list of plots.
#'
#' @family Functional enrichment and impacted pathway analysis
#'
#' @export
#'
#' @examples
#' enrichment_result <- find_impacted_pathways(
#'   gene_file = paste(system.file("extdata", package = "scFlowData"), "/",
#'     "de_result_table.tsv",
#'     sep = ""
#'   ),
#'   enrichment_method = "ORA",
#'   additional_enrichment_databse = FALSE
#' )
find_impacted_pathways <- function(gene_file = NULL,
                                   reference_file = NULL,
                                   enrichment_tool = c(
                                     "WebgestaltR", "ROntoTools"),
                                   enrichment_database = c(
                                     "pathway_KEGG",
                                     "pathway_Panther",
                                     "pathway_Reactome"
                                   ),
                                   is_output = FALSE,
                                   output_dir = ".",
                                   ...) {
  fargs <- list()
  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  res <- vector("list", length = length(enrichment_tool))
  names(res) <- enrichment_tool

  if ("WebgestaltR" %in% enrichment_tool) {
    res[["WebgestaltR"]] <- pathway_analysis_webgestaltr(
      gene_file = gene_file,
      reference_file = reference_file,
      enrichment_method = fargs$enrichment_method,
      enrichment_database = enrichment_database,
      additional_enrichment_databse = fargs$additional_enrichment_databse,
      is_output = is_output,
      output_dir = output_dir
    )
  }

  if ("ROntoTools" %in% enrichment_tool) {
    res[["ROntoTools"]] <- pathway_analysis_rontotools(
      gene_file = gene_file,
      reference_file = reference_file,
      enrichment_database = enrichment_database,
      is_output = is_output,
      output_dir = output_dir
    )
  }

  return(res)
}
