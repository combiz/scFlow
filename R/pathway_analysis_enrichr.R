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
#' @param output_dir Path for the output directory. Default is current dir.
#'
#' @return enrichment_result a list of data.frames containing enrichment output
#' and a list of plots of top 10 significant genesets.
#'
#' @family Impacted pathway analysis
#'
#' @importFrom enrichR enrichr listEnrichrDbs
#' @importFrom cli cli_alert_info cli_text
#' @importFrom ggplot2 ggplot ggsave
#' @importFrom cowplot theme_cowplot background_grid
#' @importFrom stringr str_wrap
#' @importFrom assertthat assert_that
#' @importFrom dplyr %>% mutate
#' @importFrom purrr map map_chr discard
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
                                     enrichment_database = c(
                                       "GO_Molecular_Function_2018",
                                       "GO_Cellular_Component_2018",
                                       "GO_Biological_Process_2018",
                                       "MGI_Mammalian_Phenotype_2017",
                                       "ChEA_2016",
                                       "KEGG_2019_Human",
                                       "WikiPathways_2019_Human",
                                       "Reactome_2016",
                                       "Allen_Brain_Atlas_down",
                                       "Allen_Brain_Atlas_up"
                                     ),
                                     is_output = FALSE,
                                     output_dir = ".") {
  library(enrichR)

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

  expected_cols <- c("gene")

  assertthat::assert_that(
    all(expected_cols %in% colnames(interest_gene)),
    msg = sprintf(
      "One or more expected columns missing: %s",
      paste0(expected_cols, collapse = ",")
    )
  )

  dbs <- enrichR::listEnrichrDbs()

  assertthat::assert_that(
    all(enrichment_database %in% as.character(dbs$libraryName)),
    msg = "Invalid databases specified. See enrichR::listEnrichrDbs()."
  )


  res <- enrichR::enrichr(
    genes = as.character(interest_gene$gene),
    databases = enrichment_database
  )

  enrichr_res <- purrr::map(
    res,
    ~ .format_res_table_enrichr(.)
  )

  enrichr_res <- purrr::discard(enrichr_res, function(x) {
    dim(x)[1] == 0
  })

  if (length(enrichr_res) == 0) {
    cli::cli_text(
      "{.strong No significant impacted pathways found at FDR <= 0.05! }"
    )
    enrichr_res <- NULL
  } else {
    enrichr_res$plot <- lapply(
      enrichr_res,
      function(dt) .dotplot_enrichr(dt)
    )

    if (is.data.frame(gene_file)) {
      project_name <- paste(deparse(substitute(gene_file)), sep = "")
    } else if (!is.data.frame(gene_file)) {
      project_name <- gsub("\\.tsv$", "", basename(gene_file))
      project_name <- gsub("-", "_", project_name)
    }

    output_dir <- output_dir
    sub_dir <- "enrichr_output"
    output_dir_path <- file.path(output_dir, sub_dir)
    project_dir <- file.path(output_dir_path, paste(project_name, sep = ""))

    if (isTRUE(is_output)) {
      dir.create(output_dir_path, showWarnings = FALSE)
      dir.create(project_dir, showWarnings = FALSE)
      lapply(
        names(enrichr_res)[names(enrichr_res) != "plot"],
        function(dt) {
          write.table(enrichr_res[dt],
            file = paste(project_dir, "/", dt, ".tsv", sep = ""),
            row.names = FALSE,
            col.names = gsub(
              dt, "", colnames(enrichr_res[[dt]])
            ), sep = "\t"
          )
        }
      )

      lapply(
        names(enrichr_res$plot),
        function(p) {
          ggplot2::ggsave(paste(project_dir, "/", p, ".png", sep = ""),
            enrichr_res$plot[[p]],
            device = "png", height = 8,
            width = 10, units = "in", dpi = 300
          )
        }
      )
    } else {
      cli::cli_alert_info("Output is returned as a list!")
    }


    if (is.data.frame(gene_file)) {
      enrichr_res$metadata$gene_file <- deparse(substitute(gene_file))
    } else if (!is.data.frame(gene_file)) {
      enrichr_res$metadata$gene_file <- gsub(
        "\\.tsv$", "", basename(gene_file)
      )
    }

    enrichr_res$metadata$enrichment_database <- enrichment_database
  }

  return(enrichr_res)
}


#' Format result table
#' @keywords internal

.format_res_table_enrichr <- function(res) {
  res_table <- data.frame(
    geneset = .get_geneset(res$Term),
    description = gsub(" *\\(.*?\\) *", "", res$Term),
    size = as.numeric(gsub(".*\\/", "", res$Overlap)),
    overlap = as.numeric(gsub("\\/.*", "", res$Overlap)),
    odds_ratio = round(res$Odds.Ratio, 2),
    pval = as.numeric(format(res$P.value, format = "e", digits = 2)),
    FDR = as.numeric(format(res$Adjusted.P.value, format = "e", digits = 2))
  )

  res_table$geneset <- ifelse(is.na(res_table$geneset),
    res_table$description,
    res_table$geneset
  )

  res_table <- res_table %>% dplyr::mutate("-Log10(FDR)" = as.numeric(
    format(-log10(FDR), format = "e", digits = 2)
  ))

  res_table$genes <- res$Genes

  res_table <- res_table[res_table$FDR <= 0.05, ]

  return(res_table)
}

.get_geneset <- function(term) {
  geneset <- purrr::map_chr(
    as.character(term),
    ~ strsplit(., "(", fixed = TRUE)[[1]][2]
  )
  geneset <- gsub(")", "", geneset, fixed = TRUE)
  return(geneset)
}

#' dotplot for ORA. x axis perturbation, y axis description
#' @keywords internal


.dotplot_enrichr <- function(dt) {
  dt <- dt[1:10, ]
  dt <- na.omit(dt)
  dt$description <- stringr::str_wrap(dt$description, 40)

  ggplot2::ggplot(dt, aes(
    x = odds_ratio,
    y = reorder(description, odds_ratio)
  )) +
    geom_point(aes(fill = FDR, size = overlap),
      shape = 21, alpha = 0.7, color = "black"
    ) +
    scale_size(name = "Overlap", range = c(3, 8)) +
    xlab("Total Odds Ratio") +
    ylab("") +
    scale_fill_gradient(
      low = "violetred", high = "navy", name = "FDR",
      guide = guide_colorbar(reverse = TRUE),
      limits = c(0, 0.05),
      aesthetics = c("fill")
    ) +
    guides(size = guide_legend(
      override.aes = list(fill = "violetred", color = "violetred")
    )) +
    cowplot::theme_cowplot() +
    cowplot::background_grid()
}
