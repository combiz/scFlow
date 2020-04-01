################################################################################
#' Plot ipa result as network plot
#'
#' @param enrichment_result A dataframe containing pathway enrichment analysis.
#' @param show_category No or names of pathways to see in the netwrok plot.
#' @param de_table A dataframe containg the gene name and foldchange of the
#' genes. Ideally the input file for the impacted pathway analysis.
#'
#' @return A network plot object.
#'
#' @family Impacted pathway analysis
#'
#' @importFrom igraph V plot.igraph graph.data.frame
#'
#' @export
#'
#'

networkplot <- function(enrichment_result,
                        de_table = NULL,
                        show_category = 5,
                        ...) {
  fargs <- list()
  inargs <- list(...)
  fargs[names(inargs)] <- inargs

  genesets <- .extract_genesets(enrichment_result, show_category)

  g <- .list2graph(genesets)

  size <- sapply(genesets, length)

  igraph::V(g)$size <- min(size) * 2 / 3

  n <- length(genesets)

  igraph::V(g)$size[1:n] <- size

  if (!is.null(de_table)) {
    fc <- setNames(de_table$logFC, de_table$gene)
    fc <- fc[igraph::V(g)$name[(n + 1):length(igraph::V(g))]]
    igraph::V(g)$color <- NA
    igraph::V(g)$color[(n + 1):length(igraph::V(g))] <- fc
    igraph::V(g)$color[(n + 1):length(igraph::V(g))] <- ifelse(
      fc > 0, "tomato", "steelblue")
    igraph::V(g)$color[1:n] <- "papayawhip"

    p <- igraph::plot.igraph(g,
      layout = layout_with_kk,
      vertex.color = igraph::V(g)$color,
      vertex.size = igraph::V(g)$size,
      vertex.label.color = "black",
      vertex.label.dist = 2,
      vertex.label.cex = 0.8,
      vertex.label.degree = 0,
      edge.color = "black"
    )
  } else {
    print("Gene logFC was not provided. Using default colours for the nodes.")
    igraph::V(g)$color <- "aquamarine"
    igraph::V(g)$color[1:n] <- "papayawhip"
    p <- igraph::plot.igraph(g,
      layout = layout_with_kk,
      vertex.color = igraph::V(g)$color,
      vertex.size = igraph::V(g)$size,
      vertex.label.color = "black",
      vertex.label.dist = 2,
      vertex.label.cex = 0.8,
      vertex.label.degree = 0,
      edge.color = "black"
    )
  }

  return(p)
}

#' helper functions
#' @keywords internal

.extract_genesets <- function(enrichment_result, show_category) {
  show_category <- .update_category(enrichment_result, show_category)

  genesets <- setNames(
    strsplit(
      as.character(enrichment_result$genes), ";",
      fixed = TRUE
    ),
    enrichment_result$geneset
  )

  genesets <- genesets[enrichment_result$geneset]
  names(genesets) <- enrichment_result$description

  if (is.numeric(show_category)) {
    return(genesets[1:show_category])
  } else if (is.character(show_category)) {
    return(genesets[show_category])
  }
}


.update_category <- function(enrichment_result, show_category) {
  if (is.numeric(show_category)) {
    if (nrow(enrichment_result) < show_category) {
      show_category <- nrow(enrichment_result)
      cat("There is not enough category, returning", show_category, "category")
    }
  } else if (all(is.character(show_category))) {
    show_category <- show_category
  }
  return(show_category)
}




.list2graph <- function(genesets) {
  dt <- .list2df(genesets)
  g <- igraph::graph.data.frame(dt, directed = FALSE)
  return(g)
}


.list2df <- function(genesets) {
  ldf <- lapply(1:length(genesets), function(i) {
    data.frame(
      categoryID = rep(
        names(genesets[i]),
        length(genesets[[i]])
      ),
      Gene = genesets[[i]]
    )
  })

  do.call("rbind", ldf)
}
