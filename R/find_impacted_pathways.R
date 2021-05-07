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
#' @param organism Organism name. hsapiens for Homo sapiens,
#' mmusculus for Mus musculus
#' @param enrichment_tool Enrichment tool to use. WebGestaltR, ROntoTools and
#' enrichR is implemented. Use one or more of the tools.
#' @param enrichment_database Name of the database for enrichment. User can
#' specify one or more database names. Check [scFlow::list_databases()] for
#' available database alias.
#' @param is_output If TRUE a folder will be created and results of enrichment
#' analysis will be saved otherwise a R list will be returned. Default FALSE.
#' @param output_dir Path for the output directory. Default is current
#' directory.
#' @param ... Additional arguments
#'
#' @return enrichment_result a list of list containing enrichment outputs from
#' different enrichment tools.
#'
#' @family Impacted pathway analysis
#'
#' @importFrom WebGestaltR listGeneSet
#' @importFrom enrichR listEnrichrDbs
#' @importFrom dplyr %>% filter pull
#'
#' @export
find_impacted_pathways <- function(gene_file = NULL,
                                   reference_file = NULL,
                                   organism = c("hsapiens"),
                                   enrichment_tool = c(
                                     "WebGestaltR", "enrichR"
                                   ),
                                   enrichment_database = c(
                                     "GO_Biological_Process",
                                     "GO_Cellular_Component",
                                     "GO_Molecular_Function",
                                     "KEGG",
                                     "Reactome",
                                     "Wikipathway"
                                   ),
                                   is_output = FALSE,
                                   output_dir = ".",
                                   ...) {
  fargs <- list()
  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  temp_dbs <- list_databases()

  res <- vector("list", length = length(enrichment_tool))
  names(res) <- enrichment_tool

  if ("enrichR" %in% enrichment_tool) {
    cli::cli_h2("Starting enrichment analysis by enrichR")
    res[["enrichR"]] <- pathway_analysis_enrichr(
      gene_file = gene_file,
      enrichment_database = temp_dbs$enrichR %>%
        dplyr::filter(db_alias %in% enrichment_database) %>%
        dplyr::pull(libraryName) %>%
        as.character(),
      is_output = is_output,
      output_dir = output_dir
    )
    res$enrichR$metadata$enrichment_database_link <- temp_dbs$enrichR %>%
      dplyr::filter(db_alias %in% enrichment_database) %>%
      dplyr::pull(link) %>%
      as.character()
    names(res$enrichR$metadata$enrichment_database_link) <- tolower(
      temp_dbs$enrichR %>%
        dplyr::filter(db_alias %in% enrichment_database) %>%
        dplyr::pull(db_alias) %>%
        as.character()
    )
    if (length(setdiff(names(res$enrichR), "metadata")) == 0) {
      res$enrichR$metadata$result <- FALSE
    }
    cli::cli_alert_success("enrichR analysis completed")
  }

  if ("ROntoTools" %in% enrichment_tool) {
    cli::cli_h2("Starting enrichment analysis by ROntoTools")
    res[["ROntoTools"]] <- pathway_analysis_rontotools(
      gene_file = gene_file,
      reference_file = reference_file,
      enrichment_database = temp_dbs$ROntoTools %>%
        dplyr::filter(db_alias %in% enrichment_database) %>%
        dplyr::pull(name) %>%
        as.character(),
      is_output = is_output,
      output_dir = output_dir
    )
    res$ROntoTools$metadata$enrichment_database_link <- temp_dbs$ROntoTools %>%
      dplyr::filter(db_alias %in% enrichment_database) %>%
      dplyr::pull(link) %>%
      as.character()
    names(res$ROntoTools$metadata$enrichment_database_link) <- tolower(
      temp_dbs$ROntoTools %>%
        dplyr::filter(db_alias %in% enrichment_database) %>%
        dplyr::pull(db_alias) %>%
        as.character()
    )
    if (length(setdiff(names(res$ROntoTools), "metadata")) == 0) {
      res$ROntoTools$metadata$result <- FALSE
    }
    cli::cli_alert_success("ROntoTools analysis completed")
  }

  if ("WebGestaltR" %in% enrichment_tool) {
    if (is.null(fargs$enrichment_method)) {
      stop(cli::cli_alert_danger(
        "enrichment_method is not specified.
      Specify either {.strong ORA} or {.strong GSEA}.
      For details check the help page for pathway_analysis_webgestaltr."
      ))
    } else {
      cli::cli_h2("Starting enrichment analysis by WebGestaltR")
      res[["WebGestaltR"]] <- pathway_analysis_webgestaltr(
        gene_file = gene_file,
        reference_file = reference_file,
        organism = organism,
        enrichment_method = fargs$enrichment_method,
        enrichment_database = temp_dbs$WebGestaltR %>%
          dplyr::filter(db_alias %in% enrichment_database) %>%
          dplyr::pull(name) %>%
          as.character(),
        is_output = is_output,
        output_dir = output_dir
      )
      res$WebGestaltR$metadata$enrichment_database_link <- temp_dbs$WebGestaltR %>%
        dplyr::filter(db_alias %in% enrichment_database) %>%
        dplyr::pull(link) %>%
        as.character()
      names(res$WebGestaltR$metadata$enrichment_database_link) <- tolower(
        temp_dbs$WebGestaltR %>%
          dplyr::filter(db_alias %in% enrichment_database) %>%
          dplyr::pull(db_alias) %>%
          as.character()
      )
      if (length(setdiff(names(res$WebGestaltR), "metadata")) == 0) {
        res$WebGestaltR$metadata$result <- FALSE
      }
      cli::cli_alert_success("WebGestaltR analysis completed")
    }
  }
  return(res)
}

#' Check available databases and their alias for scFlow
#'
#' @return Returns a list containing available database names and
#' their alias for scFlow
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
  temp_dbs$WebGestaltR <- temp_dbs$WebGestaltR[c(1, 3, 5, 7:10), ]
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

  eval(parse(text = "enrichR:::.onAttach()")) # R CMD check workaround

  temp_dbs[["enrichR"]] <- enrichR::listEnrichrDbs()
  temp_dbs$enrichR <- temp_dbs$enrichR[
    c(130, 131, 132, 148, 102, 93, 145, 101, 104, 49, 53),
  ]
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
