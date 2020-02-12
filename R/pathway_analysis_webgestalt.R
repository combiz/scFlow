################################################################################
#' Functional enrichment analysis using ROntoTools
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
#' @param enrichment_method Method of enrichment analysis.
#' Either over-representation analysis (ORA) or (Gene set enrichment analysis)
#' GSEA.
#' @param project_name Logical. If TRUE, the gene_file name will be used in
#' ouput. If FALSE the name will be generate with the timestamp. Default TRUE.
#' @param enrichment_database Name of the database for enrichment.
#' If NULL then multiple databases will be used or user can specify one or more
#' database names.Default NULL.
#' @param additional_enrichment_databse Logical. If TRUE additional databases
#' will be loaded. If FALSE only default databases will be used. Default TRUE.
#' @param is_output If TRUE a folder will be created and results of enrichment
#' analysis will be saved otherwise a R list will be returned. Default FALSE
#' @param output_dir Path for the output directory. Default is current
#' directory.
#'
#' @return enrichment_result a list of data.frames containing enrichment output
#' and a list of plots of top 10 significant genesets.
#'
#' @family Functional enrichment and impacted pathway analysis
#' @importfrom cli cli_alert_danger rule cli_alert_info
#' @importfrom dplyr filter
#' @importFrom WebGestaltR WebGestaltR
#' @importFrom ggplot2 ggplot ggsave
#' @importFrom cowplot theme_cowplot background_grid
#'
#' @export
#'
#' @examples
#' enrichment_result <- pathway_analysis_webgestalt(
#'   gene_file = paste(system.file("extdata", package = "scFlowData"), "/",
#'     "de_result_table.tsv",
#'     sep = ""
#'   ),
#'   enrichment_method = "ORA",
#'   project_name = TRUE,
#'   additional_enrichment_databse = FALSE,
#'   is_output = TRUE
#' )
pathway_analysis_webgestalt <- function(gene_file = NULL,
                                        reference_file = NULL,
                                        enrichment_method = "ORA",
                                        project_name = TRUE,
                                        enrichment_database = c(
                                "geneontology_Biological_Process_noRedundant",
                                "geneontology_Cellular_Component_noRedundant",
                                "geneontology_Molecular_Function_noRedundant",
                                "pathway_KEGG",
                                "pathway_Panther",
                                "pathway_Reactome",
                                "pathway_Wikipathway"
                                        ),
                                        additional_enrichment_databse = FALSE,
                                        is_output = FALSE,
                                        output_dir = ".") {
  if (is.null(gene_file)) {
    cli::cli_alert_danger("No input gene list found! \n")
  } else if (is.data.frame(gene_file)) {
    interest_gene <- gene_file
  } else {
    interest_gene <- read.delim(gene_file, sep = "\t")
  }


  if (enrichment_method == "ORA") {
    interest_gene <- as.vector(interest_gene$gene)
  } else if (enrichment_method == "GSEA") {
    interest_gene <- interest_gene[, c("gene", "logFC")]
  }

  if (is.null(reference_file)) {
    cat(
      cli::rule(
        "Using genome_protein-coding as background gene list",
        line = 2
      ),
      "\r\n"
    )
    reference_gene <- NULL
  } else if (is.data.frame(reference_file)) {
    reference_gene <- as.vector(reference_file$gene)
  }

  enrichment_database <- enrichment_database

  if (isTRUE(additional_enrichment_databse)) {
    enrichment_database_file <- list.files(
      path = paste(system.file(
        "extdata/pathway_database", package = "scFlowData"), "/", sep = ""),
      pattern = ".gmt", full.names = TRUE
    )[1]
  } else {
    enrichment_database_file <- NULL
  }

  enrichment_result <- WebGestaltR::WebGestaltR(
    enrichMethod = enrichment_method,
    organism = "hsapiens",
    enrichDatabase = enrichment_database,
    enrichDatabaseFile = enrichment_database_file,
    enrichDatabaseType = rep("genesymbol", length(enrichment_database_file)),
    enrichDatabaseDescriptionFile = NULL,
    interestGene = interest_gene,
    interestGeneType = "genesymbol",
    referenceGene = reference_gene,
    referenceGeneType = "genesymbol",
    referenceSet = "genome_protein-coding",
    networkConstructionMethod = "Network_Retrieval_Prioritization",
    isOutput = FALSE,
    outputDirectory = getwd(),
    projectName = project_name
  )


  enrichment_result$database <- factor(enrichment_result$database)
  enrichment_result <- split.data.frame(
    enrichment_result, enrichment_result$database
  )
  enrichment_result <- lapply(
    enrichment_result[names(enrichment_result) != "plot"],
    function(dt) {
      transform(
        dt,
        description = ifelse(description == "", geneSet, description)
      )
    }
  )

  if (enrichment_method == "ORA") {
    enrichment_result$plot <- lapply(
      enrichment_result,
      function(dt) .dotplot_ora(dt)
    )
  } else if (enrichment_method == "GSEA") {
    enrichment_result$plot <- lapply(
      enrichment_result,
      function(dt) .barplot_gsea(dt)
    )
  }

  project_name <- .generate_project_name(
    project_name = project_name,
    gene_file = gene_file,
    enrichment_method = enrichment_method
  )

  output_dir <- output_dir
  sub_dir <- "WebGestalt.Output"
  output_dir_path <- file.path(output_dir, sub_dir)
  project_dir <- file.path(output_dir_path, paste(project_name, sep = ""))

  if (isTRUE(is_output)) {
    dir.create(output_dir_path, showWarnings = FALSE)
    dir.create(project_dir, showWarnings = FALSE)
    lapply(
      names(enrichment_result)[names(enrichment_result) != "plot"],
      function(dt) {
        write.table(enrichment_result[dt],
          file = paste(project_dir, "/", dt, ".tsv", sep = ""),
          row.names = FALSE,
          col.names = gsub(
            dt, "", colnames(enrichment_result[[dt]])
          ), sep = "\t"
        )
      }
    )

    lapply(
      names(enrichment_result$plot),
      function(p) {
        ggplot2::ggsave(paste(project_dir, "/", p, ".png", sep = ""),
          enrichment_result$plot[[p]],
          device = "png", height = 8,
          width = 10, units = "in", dpi = 300
        )
      }
    )
  } else {
    cli::cli_alert_info("Outout is returned as a list!")
  }

  return(enrichment_result)
}



#' Generating project name
#' @keywords internal

.generate_project_name <- function(project_name = NULL,
                                   gene_file = NULL,
                                   enrichment_method = NULL) {
  if (isTRUE(project_name) && is.data.frame(gene_file)) {
    project_name <- paste(
      deparse(substitute(gene_file)), enrichment_method,
      sep = "_"
    )
  } else if (isTRUE(project_name) && !is.data.frame(gene_file)) {
    project_name <- paste(
      gsub("\\.tsv$", "", basename(gene_file)), enrichment_method,
      sep = "_"
    )
    project_name <- gsub("-", "_", project_name)
  } else {
    project_name <- as.character(as.integer(Sys.time()))
  }
}


#' dotplot for ORA. x axis enrichmentRatio, y axis description
#' @keywords internal


.dotplot_ora <- function(dt) {
  plot_title <- dt$database[1]
  dt <- dt[1:10, ]
  dt <- na.omit(dt)

  ggplot2::ggplot(dt, aes(
    x = enrichmentRatio,
    y = reorder(description, enrichmentRatio)
  )) +
    geom_point(aes(fill = FDR, size = overlap),
               shape = 21, alpha = 0.7, color = "black") +
    scale_size(name = "Size", range = c(3, 8)) +
    xlab("Enrichment Ratio") +
    ylab("") +
    ggtitle(plot_title) +
    scale_fill_gradient(
      low = "violetred", high = "navy", name = "FDR",
      guide = guide_colorbar(reverse = TRUE),
      limits = c(0, 0.05),
      aesthetics = c("fill")
    ) +
    guides(size = guide_legend(
      override.aes = list(fill = "violetred", color = "violetred")
    )) + theme_cowplot() +
    background_grid()
}


#' barplot for GSEA. x axis normalizedEnrichmentScore, y axis description
#' @keywords internal


.barplot_gsea <- function(dt) {
  plot_title <- dt$database[1]
  dt <- rbind(
    dplyr::filter(dt, normalizedEnrichmentScore > 0)[1:10, ],
    dplyr::filter(dt, normalizedEnrichmentScore < 0)[1:10, ]
  )
  dt <- na.omit(dt)

  ggplot2::ggplot(dt, aes(
    x = reorder(description, normalizedEnrichmentScore),
    y = normalizedEnrichmentScore, fill = FDR
  ), alpha = 0.7) +
    geom_bar(stat = "identity") +
    xlab("") +
    ylab("Normalised Enrichment Score") +
    ggtitle(plot_title) +
    scale_fill_continuous(
      low = "violetred", high = "navy", name = "FDR",
      guide = guide_colorbar(reverse = TRUE), limits = c(0, 0.05)
    ) +
    coord_flip() +
    theme_cowplot() +
    background_grid()
}
