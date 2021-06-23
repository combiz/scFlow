################################################################################
#' Functional enrichment analysis using ROntoTools
#'
#' Performs impacted pathway analysis with a list of gene, their fold-change
#' and p-value. The main tool used here is pathway-express (pe).
#'
#' @param gene_file A data frame or the path of a .tsv file containing
#' a list of genes, their fold-change, p-value and adjusted p-value.
#' Column names should be gene, logFC, pval and padj respectively.
#' @param reference_file A .txt file or a data frame containing all the genes
#' that were used as input for differential expression.
#' column name should be gene. If not provided the human protein-coding genome
#' will be used as background genes.
#' @param enrichment_database Name of the database for enrichment. User can
#' specify one or more database names. Default kegg.
#' @param is_output If TRUE a folder will be created and results of enrichment
#' analysis will be saved otherwise a R list will be returned. Default FALSE
#' @param output_dir Path for the output directory. Default is current dir.
#'
#' @return enrichment_result a list of data.frames containing enrichment output
#' and a list of plots of top 10 significant genesets.
#'
#' @family Impacted pathway analysis
#'
#' @importFrom cli cli_text cli_alert_info
#' @importFrom ROntoTools setNodeWeights alphaMLG pe
#' @importFrom ggplot2 ggplot ggsave
#' @importFrom cowplot theme_cowplot background_grid
#' @importFrom stringr str_wrap
#' @importFrom assertthat assert_that
#' @importFrom dplyr mutate
#' @importFrom purrr discard
#' @importFrom methods Summary
#'
#' @export
pathway_analysis_rontotools <- function(gene_file = NULL,
                                        reference_file = NULL,
                                        enrichment_database = c(
                                          "kegg",
                                          "nci",
                                          "panther",
                                          "reactome"
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

  expected_cols <- c("gene", "logFC", "pval")

  assertthat::assert_that(
    all(expected_cols %in% colnames(interest_gene)),
    msg = sprintf(
      "One or more expected columns missing: %s",
      paste0(expected_cols, collapse = ",")
    )
  )

  fc <- interest_gene$logFC
  names(fc) <- interest_gene$gene
  names(fc) <- paste("SYMBOL:", names(fc), sep = "")

  pv <- interest_gene$pval
  names(pv) <- interest_gene$gene
  names(pv) <- paste("SYMBOL:", names(pv), sep = "")


  if (is.null(reference_file)) {
    cli::cli_text("Using genome_protein-coding as background gene list")
    reference_gene <- read.delim(
      file = paste(
        system.file("extdata", package = "scFlowData"),
        "/", "human_protein_coding_genome.tsv",
        sep = ""
      ),
      header = T
    )
    reference_gene <- paste("SYMBOL:", reference_gene$gene, sep = "")
  } else if (is.data.frame(reference_file)) {
    cli::cli_text("Using custom background gene list provided by the user")
    reference_gene <- paste("SYMBOL:", reference_gene$gene, sep = "")
  }

  enrichment_database <- tolower(gsub("pathway_", "", enrichment_database))

  dbs <- .listdb()

  assertthat::assert_that(
    all(enrichment_database %in% as.character(dbs$name)),
    msg = "Invalid databases specified.  See .listdb()."
  )

  enrichment_result <- vector("list", length = length(enrichment_database))
  names(enrichment_result) <- enrichment_database

  for (database_name in enrichment_database) {
    pathway_graph <- readRDS(
      file = paste(system.file(
        "extdata/pathway_database",
        package = "scFlowData"
      ),
      "/", database_name, "_graphNEL.rds",
      sep = ""
      )
    )

    pathway_graph <- ROntoTools::setNodeWeights(
      pathway_graph,
      weights = ROntoTools::alphaMLG(pv),
      defaultWeight = 1
    )

    pathway_id <- readRDS(
      file = paste(system.file(
        "extdata/pathway_database",
        package = "scFlowData"
      ),
      "/", database_name, "_pathway_id.rds",
      sep = ""
      )
    )
    pathway_id <- unlist(pathway_id)

    res <- ROntoTools::pe(
      x = fc,
      graphs = pathway_graph,
      ref = reference_gene,
      nboot = 200,
      verbose = FALSE
    )

    res_summary <- Summary(
      res,
      pathNames = pathway_id,
      totalAcc = FALSE,
      pAcc = FALSE,
      comb.pv = c("pPert", "pORA"),
      order.by = "pComb"
    )

    res_summary <- na.omit(res_summary)

    if (dim(res_summary)[1] == 0) {
      enrichment_result[[database_name]] <- NULL
    } else {
      res_table <- data.frame(
        geneset = gsub(":", "", res_summary$pathNames),
        description = rownames(res_summary),
        size = sapply(
          res@pathways, function(x) length(x@ref)
        )[rownames(res_summary)],
        overlap = sapply(
          res@pathways, function(x) length(x@input)
        )[rownames(res_summary)],
        perturbation = round(res_summary$totalPertNorm, 2),
        pval = as.numeric(format(
          res_summary$pComb,
          format = "e", digits = 2
        )),
        FDR = as.numeric(format(
          res_summary$pComb.fdr,
          format = "e", digits = 2
        ))
      )

      res_table <- res_table %>% dplyr::mutate("-Log10(FDR)" = as.numeric(
        format(-log10(FDR), format = "e", digits = 2)
      ))

      res_table$genes <- .get_overlap_id(res = res, res_summary = res_summary)

      res_table$database <- rep(database_name)

      res_table <- res_table[res_table$FDR <= 0.05, ]

      rownames(res_table) <- NULL

      enrichment_result[[database_name]] <- res_table
    }
  }

  enrichment_result <- purrr::discard(enrichment_result, function(x) {
    is.null(x)
  })

  if (length(enrichment_result) == 0) {
    cli::cli_text(
      "{.strong No significant impacted pathways found at FDR <= 0.05! }"
    )
    enrichment_result <- NULL
  } else {
    enrichment_result$plot <- lapply(
      enrichment_result, function(dt) .dotplot_pe(dt)
    )


    if (is.data.frame(gene_file)) {
      project_name <- paste(deparse(substitute(gene_file)), sep = "")
    } else if (!is.data.frame(gene_file)) {
      project_name <- gsub("\\.tsv$", "", basename(gene_file))
      project_name <- gsub("-", "_", project_name)
    }

    output_dir <- output_dir
    sub_dir <- "rontotools_output"
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

    enrichment_result$metadata$enrichment_database <- enrichment_database
  }


  return(enrichment_result)
}

#' Getting the list of available database for ROntoTools
#'
#' @return Returns a data frame of available databases for ROntoTools
#'
#' @family Impacted pathway analysis
#'
#' @export

.listdb <- function() {
  dbs <- data.frame(database = list.files(
    path = paste(system.file("extdata/pathway_database/",
      package = "scFlowData"
    ), sep = ""),
    pattern = "_graphNEL.rds"
  ))
  dbs$name <- gsub("_graphNEL.rds", "", dbs$database)
  return(dbs)
}


#' Getting the overlapping IDs after enrichment analysis is done
#'
#' @return Returns the overlapping gene symbols after enrichment analysis
#' is done
#'
#' @family Impacted pathway analysis
#'
#' @importFrom plyr ldply
#'
#' @keywords internal

.get_overlap_id <- function(res, res_summary) {
  overlapid_list <- lapply(
    res@pathways, function(x) names(x@input)
  )[rownames(res_summary)]
  overlapid_list <- plyr::ldply(overlapid_list, rbind)
  rownames(overlapid_list) <- overlapid_list$.id
  overlapid_list$.id <- NULL
  overlapid_list$collapsed <- apply(
    overlapid_list, 1, function(x) paste(x[!is.na(x)], collapse = ";")
  )
  overlapid_list$collapsed <- gsub("SYMBOL:", "", overlapid_list$collapsed)
  return(overlapid_list$collapsed)
}


#' dotplot for ORA. x axis perturbation, y axis description
#'
#' @family Impacted pathway analysis
#'
#' @importFrom stringr str_wrap
#' @importFrom ggplot2 ggplot
#' @importFrom cowplot theme_cowplot background_grid
#'
#' @keywords internal

.dotplot_pe <- function(dt) {
  dt <- na.omit(dt)
  dt <- dt %>%
    top_n(., 10, enrichment_ratio)
  dt$description <- stringr::str_wrap(dt$description, 40)

  ggplot2::ggplot(dt, aes(
    x = perturbation,
    y = reorder(description, perturbation)
  )) +
    geom_point(aes(fill = FDR, size = size),
      shape = 21, alpha = 0.7, color = "black"
    ) +
    scale_size(name = "Geneset size", range = c(3, 8)) +
    xlab("Total perturbation") +
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
}
