################################################################################
#' Functional enrichment analysis using WebgestaltR
#'
#' Performs gene ontology and impacted pathway enrichment analysis with a list
#' of gene names and their fold-change.
#'
#' @param gene_file For ORA, A data frame containing a list of significant genes
#' with column name `gene` or a vector of significant genes. For GSEA a
#' data frame containing a list of all genes in the analysis, their fold-change,
#' p-value and adjusted p-value. Column names should be gene, logFC, pval and
#' padj respectively.
#' @param reference_file A data frame containing all the genes that were used
#' as input for differential expression. Column name should be gene.
#' If not provided the human protein-coding genome will be used as background
#' genes.
#' @param organism default is human. From WebGestaltR supports 12 organisms,
#' common choices are "hsapiens" or "mmusculus". Users can use the function
#' [WebGestaltR::listOrganism()] to check available organisms. Users can also
#' input others to perform the enrichment analysis for other organisms not
#' supported by WebGestaltR. For other organisms, users need to provide the
#' functional categories, interesting list and reference list (for ORA method).
#' Because WebGestaltR does not perform the ID mapping for the other organisms,
#' the above data should have the same ID type.
#' @param enrichment_method Method of enrichment analysis.
#' Either over-representation analysis (ORA) or (Gene set enrichment analysis)
#' GSEA.
#' @param enrichment_database Name of the database for enrichment. If not
#' provided then multiple databases will be used or user can specify one or more
#' database names from [WebGestaltR::listGeneSet()]
#'
#' @return enrichment_result a list of data.frames containing enrichment output
#' and a list of plots of top 10 significant genesets.
#'
#' @family Impacted pathway analysis
#'
#' @importFrom cli cli_text cli_alert_info
#' @importFrom dplyr %>% filter mutate
#' @importFrom WebGestaltR WebGestaltR listGeneSet
#' @importFrom ggplot2 ggplot ggsave
#' @importFrom cowplot theme_cowplot background_grid
#' @importFrom stringr str_wrap
#' @importFrom assertthat assert_that
#' @importFrom stats reorder
#'
#' @export
pathway_analysis_webgestaltr <- function(gene_file = NULL,
                                         reference_file = NULL,
                                         organism =getOption(
                                           "scflow_species",
                                           default = "human"),
                                         enrichment_method = "ORA",
                                         enrichment_database = c(
                                           "geneontology_Biological_Process_noRedundant",
                                           "geneontology_Cellular_Component_noRedundant",
                                           "geneontology_Molecular_Function_noRedundant",
                                           "pathway_KEGG",
                                           "pathway_Reactome",
                                           "pathway_Wikipathway"
                                         )) {
  assertthat::assert_that(
    !is.null(gene_file),
    msg = "No input gene list found!"
  )

  #"human" to "hsapiens" and "mouse" to "mmusculus"
  if(organism %in% c("human","mouse")){
    if(organism=="human")
      organism <- "hsapiens"
    else
      organism <- "mmusculus"
  }

  if (enrichment_method == "ORA" & is.data.frame(gene_file)) {

    assertthat::assert_that(
      all("gene" %in% colnames(gene_file)),
      msg = sprintf(
        "Expected column missing: %s",
        paste0("gene", collapse = ",")
      )
    )
    interest_gene <- gene_file$gene

  } else if (enrichment_method == "ORA" & is.vector(gene_file)) {

    interest_gene <- gene_file

  } else if (enrichment_method == "GSEA") {

    assertthat::assert_that(
      all(c("gene", "logFC", "pval") %in% colnames(gene_file)),
      msg = sprintf(
        "One or more expected columns missing: %s",
        paste0(expected_cols, collapse = ",")
      )
    )

    interest_gene <- gene_file %>%
      dplyr::mutate(rank = sign(logFC)* -log10(pval)) %>%
      dplyr::select(c(gene, rank))
  }

  if (is.null(reference_file)) {
    cli::cli_text("Using genome_protein-coding as background gene list")
    reference_gene <- NULL
  } else if (is.data.frame(reference_file)) {
    cli::cli_text("Using custom background gene list provided by the user")
    reference_gene <- as.vector(reference_file$gene)
  }

  dbs <- WebGestaltR::listGeneSet()

  assertthat::assert_that(
    all(enrichment_database %in% as.character(dbs$name)),
    msg = "Invalid databases specified. See WebGestaltR::listGeneSet()."
  )


  res <- WebGestaltR::WebGestaltR(
    enrichMethod = enrichment_method,
    organism = organism,
    enrichDatabase = enrichment_database,
    interestGene = interest_gene,
    interestGeneType = "genesymbol",
    referenceGene = reference_gene,
    referenceGeneType = "genesymbol",
    referenceSet = "genome_protein-coding",
    minNum = 5,
    maxNum = 300,
    isOutput = FALSE,
    projectName = NULL
  )

  if (is.null(res)) {
    cli::cli_text(
      "{.strong No significant impacted pathways found at FDR <= 0.05! }"
    )
    enrichment_result <- NULL
  } else {
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

    cli::cli_alert_info("Output is returned as a list!")

    enrichment_result$metadata$gene_file <- deparse(substitute(gene_file))
    enrichment_result$metadata$enrichment_method <- enrichment_method
    enrichment_result$metadata$enrichment_database <- enrichment_database
  }

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
    format(-log10(FDR), format = "e", digits = 2)
  ))

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
  dt <- dt %>%
    dplyr::filter(!is.na(FDR)) %>%
    dplyr::top_n(., min(nrow(.), 10), -pval)
  dt$description <- stringr::str_wrap(dt$description, 40)
  p <- ggplot2::ggplot(dt, aes(
    x = enrichment_ratio,
    y = reorder(description, enrichment_ratio)
  )) +
    geom_point(aes(fill = FDR, size = size),
               shape = 21, alpha = 0.7, color = "black"
    ) +
    xlab("Enrichment Ratio") +
    ylab("") +
    scale_fill_gradient(
      low = "navy", high = "gold", name = "FDR",
      guide = guide_colorbar(reverse = TRUE),
      limits = c(0, 0.05),
      aesthetics = c("fill")
    ) +
    guides(size = guide_legend(
      override.aes = list(fill = "gold", color = "gold")
    )) + cowplot::theme_cowplot() +
    cowplot::background_grid()
  if ( nrow(dt) < 4 ){
    p <- p + scale_size(name = "size", range = c(3, 8))
  } else {
    p <- p +
      scale_size_binned(name = "Geneset size", range = c(3, 8),
                        n.breaks = 4, nice.breaks = TRUE)
  }
  return(p)
}

#' barplot for GSEA. x axis normalizedEnrichmentScore, y axis description
#' @importFrom stats reorder
#' @keywords internal


.barplot_gsea <- function(dt) {
  dt <- rbind(
    dt %>% dplyr::filter(normalised_enrichment_ratio > 0, !is.na(FDR)) %>%
      dplyr::top_n(., min(nrow(.), 10), -pval),
  dt %>% dplyr::filter(normalised_enrichment_ratio < 0, !is.na(FDR)) %>%
      dplyr::top_n(., min(nrow(.), 10), -pval)
  )
  dt$description <- stringr::str_wrap(dt$description, 40)

  ggplot2::ggplot(dt, aes(
    x = stats::reorder(description, normalised_enrichment_ratio),
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
    cowplot::theme_cowplot() +
    cowplot::background_grid()
}

