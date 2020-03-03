################################################################################
#' Functional enrichment analysis using WebgestaltR
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
#' @param enrichment_database Name of the database for enrichment. If NULL
#' then multiple databases will be used or user can specify one or more
#' database names.Default NULL. User can specify one or more database names
#' from [WebGestaltR::listGeneSet()]
#' @param is_output If TRUE a folder will be created and results of enrichment
#' analysis will be saved otherwise a R list will be returned. Default FALSE
#' @param output_dir Path for the output directory. Default is current
#' directory.
#'
#' @return enrichment_result a list of data.frames containing enrichment output
#' and a list of plots of top 10 significant genesets.
#'
#' @family Impacted pathway analysis
#'
#' @importFrom cli cli_text cli_alert_info
#' @importFrom dplyr filter mutate
#' @importFrom WebGestaltR WebGestaltR listGeneSet
#' @importFrom ggplot2 ggplot ggsave
#' @importFrom cowplot theme_cowplot background_grid
#' @importFrom stringr str_wrap
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @examples
#' enrichment_result <- pathway_analysis_webgestaltr(
#'   gene_file = paste(system.file("extdata", package = "scFlowData"), "/",
#'     "de_result_table.tsv",
#'     sep = ""
#'   ),
#'   enrichment_method = "ORA",
#'   is_output = TRUE
#' )
pathway_analysis_webgestaltr <- function(gene_file = NULL,
                                         reference_file = NULL,
                                         enrichment_method = "ORA",
                                         enrichment_database = c(
                                "geneontology_Biological_Process_noRedundant",
                                "geneontology_Cellular_Component_noRedundant",
                                "geneontology_Molecular_Function_noRedundant",
                                           "pathway_KEGG",
                                           "pathway_Panther",
                                           "pathway_Reactome",
                                           "pathway_Wikipathway"
                                         ),
                                         is_output = FALSE,
                                         output_dir = ".") {
  assertthat::assert_that(
    !is.null(gene_file),
    msg = "No input gene list found!"
  )

  if (is.data.frame(gene_file)) {
    interest_gene <- gene_file
  } else {
    assertthat::assert_that(file.exists(gene_file), msg = "File not found.")
    cli::cli_text("Reading: {.file {gene_file}}")
    interest_gene <- read.delim(gene_file, sep = "\t")
  }


  if (enrichment_method == "ORA") {
    expected_cols <- c("gene")

    assertthat::assert_that(
      all(expected_cols %in% colnames(interest_gene)),
      msg = sprintf(
        "One or more expected columns missing: %s",
        paste0(expected_cols, collapse = ",")
      )
    )
    interest_gene <- as.vector(interest_gene$gene)
  } else if (enrichment_method == "GSEA") {
    expected_cols <- c("gene", "logFC")

    assertthat::assert_that(
      all(expected_cols %in% colnames(interest_gene)),
      msg = sprintf(
        "One or more expected columns missing: %s",
        paste0(expected_cols, collapse = ",")
      )
    )

    interest_gene <- interest_gene[, c("gene", "logFC")]
  }

  if (is.null(reference_file)) {
    cli::cli_text("Using genome_protein-coding as background gene list")
    reference_gene <- NULL
  } else if (is.data.frame(reference_file)) {
    cli::cli_text("Using custom background gene list provided by the user")
    reference_gene <- as.vector(reference_file$gene)
  }

  enrichment_database <- enrichment_database

  dbs <- WebGestaltR::listGeneSet()

  assertthat::assert_that(
    all(enrichment_database %in% as.character(dbs$name)),
    msg = "Invalid databases specified.  See WebGestaltR::listGeneSet()."
  )

  res <- WebGestaltR::WebGestaltR(
    enrichMethod = enrichment_method,
    organism = "hsapiens",
    enrichDatabase = enrichment_database,
    interestGene = interest_gene,
    interestGeneType = "genesymbol",
    referenceGene = reference_gene,
    referenceGeneType = "genesymbol",
    referenceSet = "genome_protein-coding",
    networkConstructionMethod = "Network_Retrieval_Prioritization",
    isOutput = FALSE,
    outputDirectory = getwd(),
    projectName = NULL
  )

  if (enrichment_method == "ORA") {
    enrichment_result <- .format_res_table_ORA(res, enrichment_database)
    enrichment_result$plot <- lapply(
      enrichment_result,
      function(dt) .dotplot_ora(dt)
    )
  } else if (enrichment_method == "GSEA") {
    enrichment_result <- .format_res_table_GSEA(res, enrichment_database)
    enrichment_result$plot <- lapply(
      enrichment_result,
      function(dt) .barplot_gsea(dt)
    )
  }

  if (is.data.frame(gene_file)) {
    project_name <- paste(
      deparse(substitute(gene_file)), enrichment_method,
      sep = "_"
    )
  } else if (!is.data.frame(gene_file)) {
    project_name <- paste(
      gsub("\\.tsv$", "", basename(gene_file)), enrichment_method,
      sep = "_"
    )
    project_name <- gsub("-", "_", project_name)
  }

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
    cli::cli_alert_info("Output is returned as a list!")
  }


  if (is.data.frame(gene_file)) {
    enrichment_result$metadata$gene_file <- deparse(substitute(gene_file))
  } else if (!is.data.frame(gene_file)) {
    enrichment_result$metadata$gene_file <- gsub(
      "\\.tsv$", "", basename(gene_file)
    )
  }
  enrichment_result$metadata$enrichment_method <- enrichment_method
  enrichment_result$metadata$enrichment_database <- enrichment_database


  return(enrichment_result)
}



#' Format result table
#' @keywords internal


.format_res_table_ORA <- function(res, enrichment_database) {
  res_table <- data.frame(
    geneset = res$geneSet,
    description = ifelse(res$description == "", res$geneSet, res$description),
    size = res$size,
    overlap = res$overlap,
    enrichment_ratio = round(res$enrichmentRatio, 2),
    pval = as.numeric(format(res$pValue, format = "e", digits = 2)),
    FDR = as.numeric(format(res$FDR, format = "e", digits = 2))
  )

  res_table <- res_table %>% dplyr::mutate("-Log10(FDR)" = as.numeric(
    format(-log10(FDR), format = "e", digits = 2)))

  res_table$genes <- res$userId

  res_table <- res_table[res_table$FDR <= 0.05, ]

  if (length(enrichment_database) == 1) {
    res_table$database <- rep(enrichment_database)
  } else {
    res_table$database <- res$database
  }

  res_table$database <- tolower(gsub("pathway_", "", res_table$database))

  res_table$database <- factor(res_table$database)

  res_table <- split.data.frame(
    res_table, res_table$database
  )

  return(res_table)
}



.format_res_table_GSEA <- function(res, enrichment_database) {
  res_table <- data.frame(
    geneset = res$geneSet,
    description = ifelse(res$description == "", res$geneSet, res$description),
    size = res$size,
    normalised_enrichment_ratio = round(res$normalizedEnrichmentScore, 2),
    pval = as.numeric(format(res$pValue, format = "e", digits = 2)),
    FDR = as.numeric(format(res$FDR, format = "e", digits = 2))
  )

  res_table <- res_table %>% dplyr::mutate("-Log10(FDR)" = -log10(FDR))

  res_table$genes <- res$userId

  res_table <- res_table[res_table$FDR <= 0.05, ]

  if (length(enrichment_database) == 1) {
    res_table$database <- rep(enrichment_database)
  } else {
    res_table$database <- res$database
  }

  res_table$database <- tolower(gsub("pathway_", "", res_table$database))

  res_table$database <- factor(res_table$database)

  res_table <- split.data.frame(
    res_table, res_table$database
  )

  return(res_table)
}

#' dotplot for ORA. x axis enrichment_ratio, y axis description
#' @keywords internal


.dotplot_ora <- function(dt) {
  dt <- dt[1:10, ]
  dt <- na.omit(dt)
  dt$description <- stringr::str_wrap(dt$description, 40)

  ggplot2::ggplot(dt, aes(
    x = enrichment_ratio,
    y = reorder(description, enrichment_ratio)
  )) +
    geom_point(aes(fill = FDR, size = overlap),
               shape = 21, alpha = 0.7, color = "black"
    ) +
    scale_size(name = "Overlap", range = c(3, 8)) +
    xlab("Enrichment Ratio") +
    ylab("") +
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
  dt <- rbind(
    dplyr::filter(dt, normalised_enrichment_ratio > 0)[1:10, ],
    dplyr::filter(dt, normalised_enrichment_ratio < 0)[1:10, ]
  )
  dt <- na.omit(dt)
  dt$description <- stringr::str_wrap(dt$description, 40)

  ggplot2::ggplot(dt, aes(
    x = reorder(description, normalised_enrichment_ratio),
    y = normalised_enrichment_ratio, fill = FDR
  ), alpha = 0.7) +
    geom_bar(stat = "identity") +
    xlab("") +
    ylab("Normalised Enrichment Score") +
    scale_fill_continuous(
      low = "violetred", high = "navy", name = "FDR",
      guide = guide_colorbar(reverse = TRUE), limits = c(0, 0.05)
    ) +
    coord_flip() +
    theme_cowplot() +
    background_grid()
}

