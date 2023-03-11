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
#' @param enrichment_tool Enrichment tool to use. WebGestaltR and
#' enrichR is implemented. Use one or more of the tools.
#' @param enrichment_database Name of the database for enrichment. User can
#' specify one or more database names. Check [scFlow::list_databases()] for
#' available database alias.
#' @param ... Additional arguments
#'
#' @return enrichment_result a list of list containing enrichment outputs from
#' different enrichment tools.
#'
#' @family Impacted pathway analysis
#'
#' @importFrom WebGestaltR listGeneSet
#' @importFrom enrichR listEnrichrDbs
#' @importFrom dplyr filter pull
#' @importFrom magrittr %>%
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
                                   ...) {
  fargs <- list()
  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  temp_dbs <- list_databases()

  res <- vector("list", length = length(enrichment_tool))
  names(res) <- enrichment_tool

  if ("enrichR" %in% enrichment_tool) {
    cli::cli_h2("Starting enrichment analysis by enrichR")

    if( "GSEA" %in% fargs$enrichment_method) {

        assertthat::assert_that(
          all(c("padj_threshold", "logFC_threshold") %in% names(fargs)),
          msg = "Please provide values for padj_threshold and
          logFC_threshold in function call"
        )

      cli::cli_h3("Filtering gene_file for significant genes for enrichR")

      gene_file <- gene_file %>%
        dplyr::filter(padj <= fargs$padj_threshold,
                      abs(logFC) >= fargs$logFC_threshold)
    }


    res[["enrichR"]] <- pathway_analysis_enrichr(
      gene_file = gene_file,
      enrichment_database = temp_dbs$enrichR %>%
        dplyr::filter(db_alias %in% enrichment_database) %>%
        dplyr::pull(libraryName) %>%
        as.character()
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
          as.character()
      )
      res$WebGestaltR$metadata$enrichment_database_link <-
        temp_dbs$WebGestaltR %>%
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
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>%
#'
#' @export

list_databases <- function() {
  temp_dbs <- list()

  temp_dbs[["WebGestaltR"]] <- WebGestaltR::listGeneSet()
  temp_dbs$WebGestaltR <- temp_dbs$WebGestaltR %>%
    dplyr::filter(name %in% c("geneontology_Biological_Process_noRedundant",
                              "geneontology_Cellular_Component_noRedundant",
                              "geneontology_Molecular_Function_noRedundant",
                              "pathway_KEGG",
                              "pathway_Reactome",
                              "pathway_Wikipathway")) %>%
    dplyr::mutate(db_alias = c("GO_Biological_Process",
                               "GO_Cellular_Component",
                               "GO_Molecular_Function",
                               "KEGG",
                               "Reactome",
                               "Wikipathway"),
                  link = c(
                    rep("http://amigo.geneontology.org/amigo/term/", 3),
                    "https://www.genome.jp/dbget-bin/www_bget?pathway:",
                    "https://reactome.org/content/detail/",
                    "https://www.wikipathways.org/index.php/Pathway:")
                  )

  eval(parse(text = "enrichR:::.onAttach()")) # R CMD check workaround

  temp_dbs[["enrichR"]] <- enrichR::listEnrichrDbs()
  temp_dbs$enrichR <- temp_dbs$enrichR %>%
    dplyr::filter(libraryName %in% c("GO_Biological_Process_2021",
                                    "GO_Cellular_Component_2021",
                                    "GO_Molecular_Function_2021",
                                    "KEGG_2021_Human",
                                    "Reactome_2022",
                                    "WikiPathway_2021_Human")) %>%
    dplyr::mutate(db_alias = c("KEGG",
                               "Wikipathway",
                               "GO_Biological_Process",
                               "GO_Cellular_Component",
                               "GO_Molecular_Function",
                               "Reactome"),
                  link = c(
                    "https://www.genome.jp/dbget-bin/www_bget?pathway:",
                    "https://www.wikipathways.org/index.php/Pathway:",
                    rep("http://amigo.geneontology.org/amigo/term/", 3),
                    "https://reactome.org/content/detail/")
                  )

  return(temp_dbs)
}
