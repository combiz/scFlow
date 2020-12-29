################################################################################
#' Convert pathway database into graphNEL format
#'
#' This function converts a pathway database into graphNEL format using R
#' package graphite. Currently available databases are kegg, reactome,
#' nci and panther.
#'
#' @param database_name Name of the database. i.e. "kegg"
#'
#' @return Saves the converted databases as rds in the extdata folder
#'
#' @family Impacted pathway analysis
#' @importFrom graph nodes edgeDataDefaults edgeData
#' @importFrom graphite pathways convertIdentifiers pathwayGraph pathwayId
#' @importFrom purrr map_chr
#'
#' @examples
#' pathway_to_graph("kegg")
#' @export

pathway_to_graph <- function(database_name) {

  # To download the pathway in graphite format.
  # The output is a list of all the pathways
  dat_pathway <- graphite::pathways("hsapiens", database_name)
  dat_pathway <- graphite::convertIdentifiers(dat_pathway, "SYMBOL")

  # To convert the list of pathways from graphite format to graphNEL format
  dat_graph <- lapply(dat_pathway, graphite::pathwayGraph, which = "protein")
  dat_pathway_id <- lapply(dat_pathway, graphite::pathwayId)

  # To check if all the list elements contain pathway data.
  node_num <- lapply(dat_graph, graph::nodes)
  zero_node <- which(lengths(node_num) == 0)

  # Keeping the pathways with nodes
  if (length(zero_node) == 0) {
    dat_graph <- dat_graph
    dat_pathway_id <- dat_pathway_id
  } else {
    dat_graph <- dat_graph[-zero_node]
    dat_pathway_id <- dat_pathway_id[-zero_node]
  }


  # Adding new attribute slot to the list of graphNEL objects
  dat_graph <- lapply(dat_graph, .add_attr)

  # Setting edge weights to the pathway

  edge_type <- levels(dat_pathway[[1]]@protEdges$type)
  sub_type <- purrr::map_chr(
    as.character(edge_type),
    ~ .clean_string(.)
  )

  edge_weight_type <- vector(
    mode = "list", length = length(grep("ACTIVATION|INHIBITION", sub_type))
  )
  names(edge_weight_type) <- sub_type[grep("ACTIVATION|INHIBITION", sub_type)]
  edge_weight_type[grep("ACTIVATION", names(edge_weight_type))] <- 1
  edge_weight_type[grep("INHIBITION", names(edge_weight_type))] <- -1

  dat_graph <- ROntoTools::setEdgeWeights(dat_graph,
    edgeTypeAttr = "subtype",
    edgeWeightByType = edge_weight_type,
    defaultWeight = 0
  )

  saveRDS(
    dat_graph,
    file = paste(
      system.file("extdata/pathway_database", package = "scFlowData"), "/",
      database_name, "_graphNEL.rds",
      sep = ""
    )
  )
  saveRDS(
    dat_pathway_id,
    file = paste(
      system.file("extdata/pathway_database", package = "scFlowData"), "/",
      database_name, "_pathway_id.rds",
      sep = ""
    )
  )
}


#' function to add new attributes to the edge
#' @keywords internal

.add_attr <- function(g) {
  graph::edgeDataDefaults(g, attr = "subtype") <- "undefined"

  edge_name <- names(graph::edgeData(g))

  from <- purrr::map_chr(
    as.character(edge_name),
    ~ strsplit(., "|", fixed = TRUE)[[1]][1]
  )

  to <- purrr::map_chr(
    as.character(edge_name),
    ~ strsplit(., "|", fixed = TRUE)[[1]][2]
  )

  edge_type <- unlist(graph::edgeData(g, attr = "edgeType"))

  sub_type <- purrr::map_chr(
    as.character(edge_type),
    ~ .clean_string(.)
  )

  graph::edgeData(g, from = from, to = to, attr = "subtype") <- sub_type
  return(g)
}


.clean_string <- function(x) {
  x <- gsub("Control[(]Out: ", "", x)
  x <- gsub("Control[(]In: ", "", x)
  x <- gsub("[(/]", "_", x)
  x <- gsub("[)]", "", x)
  x <- gsub(" ", "_", x)
  x <- purrr::map_chr(
    as.character(x),
    ~ strsplit(., ";", fixed = TRUE)[[1]][1]
  )
  return(x)
}
