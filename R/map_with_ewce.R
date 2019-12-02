
################################################################################
#' Map cluster celltypes with EWCE
#'
#' Currently works with Allan and Zeisel only
#'
#' @param sce a SingleCellExperiment
#' @param ctd_folder path to a folder with ctd RDS files
#' @param cells_to_sample total cells to consider for gene enrichment
#' @param clusters_colname the name of the colData column with cluster number
#'
#' @return sce a SingleCellExperiment object annotated with celltypes/metadata
#' @author Nathan Skene / Combiz Khozoie
#' @family clustering and dimensionality reduction
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom Matrix Matrix
#' @importFrom EWCE generate.celltype.data
#' @importFrom assertthat assert_that
#' @importFrom dplyr rename inner_join
#' @importFrom plyr adply
#' @importFrom magrittr %>%
#' @export
map_celltypes_sce <- function(sce,
                              ctd_folder,
                              cells_to_sample = 10000,
                              clusters_colname = "clusters") {

  assertthat::assert_that(dir.exists(ctd_folder))
  l_ctd <- .read_rds_files_to_list(ctd_folder)

  if(dim(sce)[[2]] > cells_to_sample) {
    message(sprintf("Subsetting %s cells", cells_to_sample))
    set.seed(42)
    sce_subset <- sce[, sample(dim(sce)[[2]], cells_to_sample)]
  } else {
    sce_subset <- sce
  }

  mat <- Matrix::Matrix(counts(sce_subset), sparse = TRUE)
  rownames(mat) <- as.character(
    SummarizedExperiment::rowData(sce_subset)$gene
  )

  annotLevels = list(level1class = sce_subset[[clusters_colname]])
  message("generating ctd with ewce")
  ctd <- EWCE::generate.celltype.data(exp = mat,
                                      annotLevels = annotLevels,
                                      groupName = "ctd")
  load(ctd[1])

  sce@metadata$ctd <- ctd

  message("mapping celltypes with ewce")
  mappings <- .map_celltypes_with_ewce(ctd, l_ctd)

  # generate lookup for merging into colData
  mappings_lookup <- mappings %>%
    dplyr::select(Cluster, Mapped, mapping) %>%
    tidyr::spread(mapping, Mapped) %>%
    dplyr::arrange(Cluster)

  # split long allan name into components
  allan_split <- purrr::map_df(
    as.character(mappings_lookup$Allan2019_ML2),
    function(x) data.frame(t(unlist(strsplit(x, split = " "))))
  )
  colnames(allan_split) <- c("allan_celltype", "allan_layer",
                             "allan_cluster_gene_1", "allan_cluster_gene_2")

  mappings_lookup <- cbind(mappings_lookup, allan_split)

  generate_ct_label <- plyr::adply(mappings_lookup, 1, function(x) {
    if (x$Allan2019_ML1 == "Non-neuronal") {
      x$allan_celltype
    } else if (x$allan_celltype == "Exc") {
      paste("EN", x$allan_layer, sep = "-")
    } else if (x$allan_celltype == "Inh") {
      paste("IN", x$allan_cluster_gene_1, sep = "-")
    } else x$allan_celltype
  }, .expand = FALSE, .id = "Cluster") %>%
    dplyr::rename(cluster_celltype = V1)

  mappings_lookup <- merge(mappings_lookup, generate_ct_label, on = "Cluster")
  mappings_lookup$Cluster <- as.factor(mappings_lookup$Cluster)

  # rename Cluster column
  mappings_lookup <- mappings_lookup %>%
    dplyr::rename(!!clusters_colname := Cluster)

  sce@metadata$mappings <- mappings_lookup

  #update sce
  new_coldata <- dplyr::inner_join(as.data.frame(colData(sce)),
                                   mappings_lookup, by = clusters_colname)

  sce <- SingleCellExperiment(
    assays = list(counts = counts(sce)),
    colData = new_coldata,
    rowData = as.data.frame(rowData(sce)),
    reducedDims = reducedDims(sce),
    metadata = sce@metadata
  )

  sce <- .append_celltype_plots_sce(sce)

  return(sce)


}


################################################################################
#' Map cluster celltypes with EWCE
#'
#' Works with Allan and Zeisel only
#'
#' @param ctd ctd for current expt
#' @param l_ctd a a list of reference ctd's
#'
#' @return sce a SingleCellExperiment object annotated with sample metadata
#' @author Nathan Skene / Combiz Khozoie
#' @family clustering and dimensionality reduction
#' @importFrom purrr pmap
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @keywords internal
.map_celltypes_with_ewce <- function(ctd, l_ctd) {

  # cortical types only
  cortical_types <- c(
    "ABC", "ACTE", "CHOR", "COP", "EPEN", "EPSC", "MFOL",
    "MGL", "MOL", "NFOL", "OPC", "PER", "PVM", "TEGLU",
    "TEINH", "VECA", "VECC", "VECV", "VEND", "VLMC",
    "VSMC", "VSMCA"
  )

  cortical_types <- paste(cortical_types, collapse = "|")

  idx <- grep(cortical_types, colnames(l_ctd$Zeisel2018[[5]]$specificity))
  l_ctd$Zeisel2018[[5]]$specificity <- l_ctd$Zeisel2018[[5]]$specificity[, idx]

  # set up the map_celltypes parameter combinations for pmap
  map_params <- list(
    ctd = c(
      "Zeisel2018", "Zeisel2018",
      "Allan2019", "Allan2019"
    ),
    species = c("mouse", "mouse", "human", "human"),
    mlevel = c(4, 5, 1, 2)
  ) # mapping levels

  l_mappings <- purrr::pmap(map_params, function(ctd_name, species, mlevel) {
    mapped_dt <- .map_celltypes(
      ctdToMap = ctd,
      ctdToMapAgainst = l_ctd[[ctd_name]],
      inputSpecies = "human",
      mapAgainstSpecies = species,
      annotLevel = 1,
      numTopMarkers = 50,
      mappingLevel = mlevel
    )
    mapped_dt$mapping <- sprintf("%s_ML%s", ctd_name, mlevel)
    mapped_dt <- mapped_dt %>%
      dplyr::rename(Cluster = Original)
    mapped_dt$Cluster <- as.numeric(as.character(mapped_dt$Cluster))
    return(mapped_dt)
  })

  mappings <- Reduce(rbind, l_mappings)
  return(mappings)

}

################################################################################
#' Main mapping function
#'
#' @param ctdToMap ctd for current expt
#' @param ctdToMapAgainst the ctd to map against
#' @param inputSpecies "human" or "mouse"
#' @param mapAgainstSpecies "human" or "mouse"
#' @param annotLevel level of annotation
#' @param numTopMarkers number of top markers to consider
#' @param mappingLevel level of mapping
#'
#' @return sce a SingleCellExperiment object annotated with sample metadata
#' @author Nathan Skene / Combiz Khozoie
#' @family clustering and dimensionality reduction
#' @importFrom EWCE bootstrap.enrichment.test
#' @importFrom purrr pmap
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @keywords internal
.map_celltypes <- function(ctdToMap,
                          ctdToMapAgainst = MAGMA.Celltyping::ctd_Tasic,
                          inputSpecies = "human",
                          mapAgainstSpecies = "mouse",
                          annotLevel = 2,
                          numTopMarkers = 500,
                          mappingLevel = 2) {

  ctd <- ctdToMapAgainst
  count <- 0
  if (!is.null(names(ctdToMap))) {
    numLevels <- sum(names(ctdToMap) == "")
  } else {
    numLevels <- length(ctdToMap)
  }

  for (ct in colnames(ctdToMap[[annotLevel]]$specificity)) {

    mostSpecificGenes <- .get_x_most_specific_genes(
      ct = ct,
      annotLevel = annotLevel,
      howMany = numTopMarkers,
      ctd = ctdToMap,
      exprPercentile = 0.9)$x_most_specific

    mostSpecificGenes <- mostSpecificGenes[!is.na(mostSpecificGenes)]
    full_results <- EWCE::bootstrap.enrichment.test(
      sct_data = ctd,
      hits = mostSpecificGenes,
      bg = rownames(ctdToMap[[1]]$specificity),
      reps = 1000,
      annotLevel = mappingLevel,
      genelistSpecies = inputSpecies,
      sctSpecies = mapAgainstSpecies
    )
    p <- full_results$results[order(
      full_results$results$sd_from_mean, decreasing = TRUE), ][1, ]$p
    z <- full_results$results[order(
      full_results$results$sd_from_mean, decreasing = TRUE), ][1, ]$sd_from_mean
    ctMapped <- data.frame(
      Original = ct,
      Mapped = as.character(
        full_results$results[order(
          full_results$results$sd_from_mean, decreasing = TRUE), ][1, ]$CellType),
      annotLevel = annotLevel,
      p = p,
      z = z
    )
    count <- count + 1
    if (count == 1) {
      allMapped <- ctMapped
    } else {
      allMapped <- rbind(ctMapped, allMapped)
    }
  }

  return(allMapped)

}

################################################################################
#' Find x most specific genes
#'
#' @param ct from colnames of ctdToMap specificity
#' @param annot_level annotation level
#' @param howMany number of genes to return
#' @param ctd the ctd to analyse
#' @param exprPercentile must be expressed in at least this proportion
#'
#' @return output the x most specific genes
#' @author Nathan Skene / Combiz Khozoie
#' @family clustering and dimensionality reduction
#' @importFrom stats quantile
#' @keywords internal
.get_x_most_specific_genes <- function(ct,
                                   annotLevel = 5,
                                   howMany = 20,
                                   ctd,
                                   exprPercentile = 0.9) {

  specificity_vector <- ctd[[annotLevel]]$specificity[, ct]
  expr <- ctd[[annotLevel]]$mean_exp[, ct]
  keep_genes <- expr > stats::quantile(expr, exprPercentile)
  specificity_vector <- specificity_vector[keep_genes]
  names(specificity_vector) <- rownames(
    ctd[[annotLevel]]$specificity[keep_genes, ]
  )
  mostSpec <- names(
    specificity_vector[order(specificity_vector, decreasing = TRUE)][1:howMany]
  )
  mostSpec_string <- paste(mostSpec, collapse = ", ")
  output <- list()
  output$x_most_specific <- mostSpec
  output$x_most_specific_string <- mostSpec_string
  output$x_with_expr <- expr[mostSpec]
  return(output)
}

################################################################################
#' Helper function to read rds files into a list
#'
#' @param folder_path path to folder with rds files
#'
#' @return rds_l list of objects imported with readRDS
#' @family helper functions
#' @importFrom purrr map
#' @importFrom tools file_path_sans_ext
#' @keywords internal
.read_rds_files_to_list <- function(folder_path){

  rds_l <- purrr::map(
    list.files(folder_path, full.names = TRUE),
    readRDS
  )
  names(rds_l) <- tools::file_path_sans_ext(list.files(folder_path))

  return(rds_l)
}

################################################################################
#' Helper function to plot celltypes in reduced dimensions
#'
#' @param sce SingleCellExperiment
#' @param celltype_dim the colData variable with the celltype annotation
#'
#' @return sce a SCE with plots appended to the metadata
#' @family helper functions
#' @importFrom SingleCellExperiment reducedDims reducedDim
#' @importFrom leaflet colorFactor
#' @importFrom threejs scatterplot3js
#' @importFrom plotly plot_ly
#' @keywords internal
.append_celltype_plots_sce <- function(sce,
                                       celltype_dim = "cluster_celltype") {

  # 2d
  sce@metadata$celltype_plots <- list()
  for (reddim in names(SingleCellExperiment::reducedDims(sce))) {
    sce@metadata$celltype_plots[[reddim]] <- plot_umap_with_feature(
      sce,
      feature_dim = celltype_dim,
      reduced_dim = reddim,
      label_clusters = TRUE,
      size = 1)
  }

  #3d
  if ("UMAP3D" %in% names(SingleCellExperiment::reducedDims(sce))) {

    umap_res <- SingleCellExperiment::reducedDim(sce, "UMAP3D")

    pal <- leaflet::colorFactor(palette = "Accent", domain = sce[[celltype_dim]])
    celltype_pal <- pal(sce[[celltype_dim]])

    sce@metadata$celltype_plots$umap3d_3js <- threejs::scatterplot3js(
      x = umap_res[, 1],
      y = umap_res[, 2],
      z = umap_res[, 3],
      color = celltype_pal,
      axis = FALSE,
      num.ticks = NULL,
      stroke = "#4682b4",
      size = .01
    )

    sce@metadata$celltype_plots$umap3d_plotly <- plotly::plot_ly(
      x = umap_res[, 1],
      y = umap_res[, 2],
      z = umap_res[, 3],
      color = celltype_pal,
      name = sce[[celltype_dim]],
      type = "scatter3d",
      mode = "markers",
      size = 0.1
    )
  }

  return(sce)
}


