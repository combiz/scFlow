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
#' enrichR is implemented.
#' @param enrichment_database Name of the database for enrichment.
#' If NULL then multiple databases will be used or user can specify one or more
#' database names.Default NULL.
#' @param is_output If TRUE a folder will be created and results of enrichment
#' analysis will be saved otherwise a R list will be returned. Default FALSE
#' @param output_dir Path for the output directory. Default is current
#' directory.
#'
#' @return enrichment_result a list of list containing enrichment outputs from
#' different enrichment tools.
#'
#' @family Impacted pathway analysis
#'
#' @importFrom WebGestaltR listGeneSet
#' @importFrom enrichR listEnrichrDbs
#'
#' @export
#'
#' @examples
#' enrichment_result <- find_impacted_pathways(
#'   gene_file = paste(system.file("extdata", package = "scFlowData"), "/",
#'     "de_result_table.tsv",
#'     sep = ""
#'   ),
#'   enrichment_method = "ORA"
#' )
find_impacted_pathways <- function(gene_file = NULL,
                                   reference_file = NULL,
                                   enrichment_tool = c(
                                     "WebGestaltR", "ROntoTools", "enrichR"
                                   ),
                                   enrichment_database = c("KEGG",
                                                           "Panther",
                                                           "Reactome"),
                                   is_output = FALSE,
                                   output_dir = ".",
                                   ...) {
  fargs <- list()
  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  temp_dbs <- list_databases()

  res <- vector("list", length = length(enrichment_tool))
  names(res) <- enrichment_tool

  if ("WebGestaltR" %in% enrichment_tool) {
    if (is.null(fargs$enrichment_method)) {
      stop(cli::cli_alert_danger(
        "enrichment_method is not specified.
      Specify either {.strong ORA} or {.strong GSEA}.
      For details check the help page for pathway_analysis_webgestaltr."))
    } else {
      cli::cli_h2("Starting enrichment analysis by WebGestaltR")
      res[["WebGestaltR"]] <- pathway_analysis_webgestaltr(
        gene_file = gene_file,
        reference_file = reference_file,
        enrichment_method = fargs$enrichment_method,
        enrichment_database = temp_dbs$WebGestaltR$name[
          temp_dbs$WebGestaltR$db_alias %in% enrichment_database],
        is_output = is_output,
        output_dir = output_dir
      )
      res$WebGestaltR$metadata$enrichment_database_link <- temp_dbs$WebGestaltR$link[
        temp_dbs$WebGestaltR$db_alias %in% enrichment_database]
      names(res$WebGestaltR$metadata$enrichment_database_link) <- tolower(
        temp_dbs$WebGestaltR$db_alias[
          temp_dbs$WebGestaltR$db_alias %in% enrichment_database])
      cli_alert_success("WebGestaltR analysis completed")
    }
  }

  if ("ROntoTools" %in% enrichment_tool) {
    cli::cli_h2("Starting enrichment analysis by ROntoTools")
    res[["ROntoTools"]] <- pathway_analysis_rontotools(
      gene_file = gene_file,
      reference_file = reference_file,
      enrichment_database = temp_dbs$ROntoTools$name[
        temp_dbs$ROntoTools$db_alias %in% enrichment_database],
      is_output = is_output,
      output_dir = output_dir
    )
    res$ROntoTools$metadata$enrichment_database_link <- temp_dbs$ROntoTools$link[
      temp_dbs$ROntoTools$db_alias %in% enrichment_database]
    names(res$ROntoTools$metadata$enrichment_database_link) <- tolower(
      temp_dbs$ROntoTools$db_alias[
        temp_dbs$ROntoTools$db_alias %in% enrichment_database])
    cli_alert_success("ROntoTools analysis completed")
  }

  if ("enrichR" %in% enrichment_tool) {
    cli::cli_h2("Starting enrichment analysis by enrichR")
    res[["enrichR"]] <- pathway_analysis_enrichr(
      gene_file = gene_file,
      enrichment_database = temp_dbs$enrichR$libraryName[
        temp_dbs$enrichR$db_alias %in% enrichment_database],
      is_output = is_output,
      output_dir = output_dir
    )
    res$enrichR$metadata$enrichment_database_link <- temp_dbs$enrichR$link[
      temp_dbs$enrichR$db_alias %in% enrichment_database]
    names(res$enrichR$metadata$enrichment_database_link) <- tolower(
      temp_dbs$enrichR$db_alias[
        temp_dbs$enrichR$db_alias %in% enrichment_database])
    cli_alert_success("enrichR analysis completed")
  }

  return(res)
}

#' Check available databases and their alias for scFlow
#'
#' @return Returns a list containing available database names and their alias for scFlow
#'
#' @family Impacted pathway analysis
#'
#' @importFrom WebGestaltR listGeneSet
#' @importFrom enrichR listEnrichrDbs
#'
#' @export

list_databases <- function() {
  temp_dbs <- list()

  temp_dbs[["WebGestaltR"]] <- WebGestaltR::listGeneSet()
  temp_dbs$WebGestaltR <- temp_dbs$WebGestaltR[c(2, 4, 6:10), ]
  temp_dbs$WebGestaltR$db_alias <- c(
    "GO_Biological_Process",
    "GO_Cellular_Component",
    "GO_Molecular_Function",
    "KEGG",
    "Panther",
    "Reactome",
    "Wikipathway"
  )
  temp_dbs$WebGestaltR$link <- c(
    rep("http://amigo.geneontology.org/amigo/term/", 3),
    "https://www.genome.jp/dbget-bin/www_bget?pathway:",
    "http://www.pantherdb.org/pathway/pathwayDiagram.jsp?catAccession=",
    "https://reactome.org/content/detail/",
    "https://www.wikipathways.org/index.php/Pathway:"
  )

  library(enrichR)
  temp_dbs[["enrichR"]] <- enrichR::listEnrichrDbs()
  temp_dbs$enrichR <- temp_dbs$enrichR[
    c(130, 131, 132, 148, 102, 93, 145, 101, 104, 49, 53), ]
  temp_dbs$enrichR$db_alias <- c(
    "GO_Biological_Process",
    "GO_Cellular_Component",
    "GO_Molecular_Function",
    "KEGG",
    "Panther",
    "Reactome",
    "Wikipathway",
    "NCI",
    "ChEA_2016",
    "Allen_Brain_Atlas_up",
    "Allen_Brain_Atlas_down"
  )
  temp_dbs$enrichR$link <- c(
    rep("http://amigo.geneontology.org/amigo/term/", 3),
    "https://www.genome.jp/dbget-bin/www_bget?pathway:",
    "http://www.pantherdb.org/pathway/pathwayDiagram.jsp?catAccession=",
    "https://reactome.org/content/detail/",
    "https://www.wikipathways.org/index.php/Pathway:",
    rep(NA, 4)
  )

  temp_dbs[["ROntoTools"]] <- .listdb()
  temp_dbs$ROntoTools$db_alias <- c(
    "KEGG",
    "NCI",
    "Panther",
    "Reactome"
  )
  temp_dbs$ROntoTools$link <- c(
    "https://www.genome.jp/dbget-bin/www_bget?pathway:",
    NA,
    "http://www.pantherdb.org/pathway/pathwayDiagram.jsp?catAccession=",
    "https://reactome.org/content/detail/"
  )

  return(temp_dbs)
}
