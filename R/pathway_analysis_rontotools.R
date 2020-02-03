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
#' @param project_name Logical. If TRUE, the gene_file name will be used in
#' ouput. If FALSE the name will be generate with the timestamp. Default TRUE.
#' @param enrichment_database Name of the database for enrichment. User can
#' specify one or more database names. Default kegg.
#' @param is_output If TRUE a folder will be created and results of enrichment
#' analysis will be saved otherwise a R list will be returned. Default FALSE
#' @param plot_output If TRUE plots will be saved as .png file otherwise
#' returned as a R list. Default FALSE.
#' @param output_dir Path for the output directory. Default is current dir.
#'
#' @return enrichment_result a list of data.frames containing enrichment output
#' and a list of plots of top 10 significant genesets.
#'
#' @family Functional enrichment and impacted pathway analysis
#' @importfrom cli cli_alert_danger rule cli_alert_info
#' @importFrom ROntoTools setNodeWeights alphaMLG pe Summary
#' @importFrom ggplot2 ggplot
#'
#' @export
#'
#' @examples
#' set.seed(1234)
#' enrichment_result <- pathway_analysis_rontotools(
#'   gene_file = "EN-L5-6-I_MASTZLM__DE_groupLow.tsv",
#'   project_name = TRUE,
#'   enrichment_database = "kegg",
#'   is_output = FALSE
#' )
pathway_analysis_rontotools <- function(gene_file = NULL,
                                        reference_file = NULL,
                                        project_name = TRUE,
                                        enrichment_database = "kegg",
                                        is_output = FALSE,
                                        output_dir = ".") {
  library(ggplot2)


  if (is.null(gene_file)) {
    cli::cli_alert_danger("No input gene list found! \n")
  } else if (is.data.frame(gene_file)) {
    interest_gene <- gene_file
  } else {
    interest_gene <- read.delim(gene_file, sep = "\t")
  }


  fc <- interest_gene$logFC
  names(fc) <- interest_gene$gene
  names(fc) <- paste("SYMBOL:", names(fc), sep = "")

  pv <- interest_gene$pval
  names(pv) <- interest_gene$gene
  names(pv) <- paste("SYMBOL:", names(pv), sep = "")


  if (is.null(reference_file)) {
    cat(cli::rule(
      "Using genome_protein-coding as background gene list", line = 2),
      "\r\n")
    reference_gene <- read.delim(
      file = paste(
        system.file("extdata", package = "scflow"),
        "/", "human_protein_coding_genome.tsv",
        sep = ""
      ),
      header = T
    )
    reference_gene <- paste("SYMBOL:", reference_gene$gene, sep = "")
  } else if (is.data.frame(reference_file)) {
    reference_gene <- paste("SYMBOL:", reference_gene$gene, sep = "")
  }


  enrichment_result <- vector("list", length = length(enrichment_database))
  names(enrichment_result) <- enrichment_database

  for (database_name in enrichment_database) {
    pathway_graph <- readRDS(
      file = paste(system.file("extdata", package = "scflow"),
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
      file = paste(system.file("extdata", package = "scflow"),
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
    res_summary <- ROntoTools::Summary(
      res,
      pathNames = pathway_id,
      totalAcc = FALSE,
      totalAccNorm = FALSE,
      pAcc = FALSE,
      pAcc.fdr = FALSE,
      comb.pv = NULL,
      order.by = "pPert"
    )

    res_summary <- na.omit(res_summary)


    res_table <- data.frame(
      geneSet = res_summary$pathNames,
      description = rownames(res_summary),
      size = sapply(
        res@pathways, function(x) length(x@ref)
      )[rownames(res_summary)],
      overlap = sapply(
        res@pathways, function(x) length(x@input)
      )[rownames(res_summary)],
      perturbation = res_summary$totalPert,
      pValue = res_summary$pPert,
      FDR = res_summary$pPert.fdr,
      pValue_ORA = res_summary$pORA,
      FDR_ORA = res_summary$pORA.fdr,
      overlapID = .get_overlap_id(res = res, res_summary = res_summary),
      database = rep(database_name)
    )

    rownames(res_table) <- NULL

    enrichment_result[[database_name]] <- res_table
  }

  enrichment_result$plot <- lapply(
    enrichment_result, function(dt) .dotplot_pe(dt)
  )


  output_dir <- output_dir
  sub_dir <- "ROntoTools.Output"
  output_dir_path <- file.path(output_dir, sub_dir)

  if (isTRUE(is_output)) {
    dir.create(output_dir_path, showWarnings = FALSE)
  } else {
    cli::cli_alert_info("Outout is returned as a list!")
  }

  if (isTRUE(project_name) && is.data.frame(gene_file)) {
    project_name <- paste(deparse(substitute(gene_file)), sep = "")
  } else if (isTRUE(project_name) && !is.data.frame(gene_file)) {
    project_name <- paste(gsub("\\.tsv$", "", basename(gene_file)), sep = "")
    project_name <- gsub("-", "_", project_name)
  } else {
    project_name <- as.character(as.integer(Sys.time()))
  }

  project_dir <- file.path(output_dir_path, paste(project_name, sep = ""))
  dir.create(project_dir, showWarnings = FALSE)

  if (isTRUE(is_output)) {
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
  }

  return(enrichment_result)
}


#' Getting the overlapping IDs after enrichment analysis is done
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
#' @keywords internal


.dotplot_pe <- function(dt) {
  plot_title <- dt$database[1]
  dt <- dt[1:10, ]
  dt <- na.omit(dt)

  ggplot2::ggplot(dt, aes(
    x = perturbation,
    y = reorder(description, perturbation)
  )) +
    geom_point(aes(color = FDR, fill = FDR, size = overlap), shape = 21) +
    scale_size(name = "Size", range = c(3, 8)) +
    xlab("Total perturbation") +
    ylab("") +
    ggtitle(plot_title) +
    scale_fill_gradient(
      low = "tomato", high = "blue", name = "FDR",
      guide = guide_colorbar(reverse = TRUE),
      limits = c(0, 0.05),
      aesthetics = c("fill", "colour")
    ) +
    guides(size = guide_legend(
      override.aes = list(fill = "tomato", color = "tomato")
    )) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(colour = "black", size = 12),
      legend.key = element_blank()
    )
}
