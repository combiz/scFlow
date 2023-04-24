  ################################################################################
#' Functional enrichment analysis using enrichR
#'
#' Performs impacted pathway analysis with a list of genes.
#'
#' @param gene_file A data frame containing a list of significant genes with
#' column name gene or a vector of significant genes.
#' @param enrichment_database Name of the database for enrichment. User can
#' specify one or more database names from [enrichR::listEnrichrDbs()].
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
#' @importFrom dplyr %>% mutate top_n
#' @importFrom purrr map map_chr discard
#'
#' @export
pathway_analysis_enrichr <- function(gene_file = NULL,
                                     enrichment_database = c(
                                       "GO_Molecular_Function_2021",
                                       "GO_Cellular_Component_2021",
                                       "GO_Biological_Process_2021",
                                       "KEGG_2021_Human",
                                       "WikiPathways_2021_Human",
                                       "Reactome_2022")
                                     ) {
  assertthat::assert_that(
    !is.null(gene_file),
    msg = "No input gene list found!"
  )

  if (is.data.frame(gene_file)) {
    assertthat::assert_that(
    all("gene" %in% colnames(gene_file)),
    msg = sprintf(
      "Expected column missing: %s",
      paste0("gene", collapse = ",")
    )
  )
    interest_gene <- gene_file$gene

  } else {

    interest_gene <- gene_file
  }

  eval(parse(text = "enrichR:::.onAttach()")) # R CMD check workaround

  dbs <- enrichR::listEnrichrDbs()

  assertthat::assert_that(
    all(enrichment_database %in% as.character(dbs$libraryName)),
    msg = "Invalid databases specified. See enrichR::listEnrichrDbs()."
  )

  res <- enrichR::enrichr(
    genes = as.character(interest_gene),
    databases = enrichment_database
  )

  res <- purrr::discard(res, function(x) {
    dim(x)[1] == 0
  })

  enrichr_res <- lapply(names(res), function(x){
    dt <- res[[x]] %>%
    mutate(database = x)})

  names(enrichr_res) <- names(res)

  enrichr_res <- purrr::map(
    enrichr_res,
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
  }

    cli::cli_alert_info("Output is returned as a list!")

    enrichr_res$metadata$gene_file <- deparse(substitute(gene_file))

    enrichr_res$metadata$enrichment_database <- enrichment_database


  return(enrichr_res)
}


#' Format result table
#' @keywords internal


.format_res_table_enrichr <- function(res) {
  res_table <- res %>% as.data.frame() %>%
    dplyr::transmute(
      geneset = .get_geneset(Term),
      description = gsub("\\(GO:.*|Homo sapiens.R-HSA.*|WP.*" , "", Term),
      size = as.numeric(gsub(".*\\/", "", Overlap)),
      overlap = as.numeric(gsub("\\/.*", "", Overlap)),
      odds_ratio = round(Odds.Ratio, 2),
      pval = as.numeric(format(P.value, format = "e", digits = 2)),
      FDR = as.numeric(format(Adjusted.P.value, format = "e", digits = 2))
    )

  res_table$geneset <- ifelse(is.na(res_table$geneset),
    res_table$description,
    res_table$geneset
  )

  res_table$genes <- res$Genes

  res_table <- res_table[res_table$FDR <= 0.05, ]

  return(res_table)
}


#' @keywords internal

.get_geneset <- function(term) {

  geneset <- purrr::map_chr(as.character(term),
                            ~ stringr::str_extract(. , "GO:.*|R-HSA.*|WP.*"))
  geneset <- gsub("\\)|Homo sapiens", "", geneset)
  geneset <- as.character(geneset)
  return(geneset)
}


#' dotplot for ORA. x axis perturbation, y axis description
#' @importFrom stats reorder
#' @keywords internal
.dotplot_enrichr <- function(dt) {
  dt <- dt %>%
    dplyr::filter(!is.na(FDR)) %>%
    dplyr::top_n(., min(nrow(.), 10), -pval)
  dt$description <- stringr::str_wrap(dt$description, 40)
  p <- ggplot2::ggplot(dt, aes(
    x = odds_ratio,
    y = stats::reorder(description, odds_ratio)
  )) +
    geom_point(aes(fill = FDR, size = size),
               shape = 21, alpha = 0.7, color = "black"
    ) +
    xlab("Odds ratio") +
    ylab("") +
    scale_fill_gradient(
      low = "navy", high = "gold", name = "FDR",
      guide = guide_colorbar(reverse = TRUE),
      limits = c(0, 0.05),
      aesthetics = c("fill")
    ) +
    guides(size = guide_legend(
      override.aes = list(fill = "gold", color = "gold")
    )) +
    cowplot::theme_cowplot() +
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
