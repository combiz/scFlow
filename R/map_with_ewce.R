
################################################################################
#' Map cluster celltypes with EWCE
#'
#' @param sce a SingleCellExperiment
#' @param ctd_folder path to a folder with ctd RDS file
#' @param cells_to_sample total cells to consider for gene enrichment
#' @param clusters_colname the name of the colData column with cluster number
#' @param species specify species of analysis dataset
#' @param save_path The directory where the generated ctd file produced by EWCE
#' is saved. Default is a temp directory
#' @param annotation_level The level(s) of annotation required.
#' Integer for single CTD/level. For multi-CTD/level a named list,
#' names are the CTD file names (without file extension) and
#' elements are vectors of levels. Levels must present in the provided CTD(s).
#' The first name and level will also be used as the primary annotation for
#' automated reporting and will be returned in 'cluster_celltype' column data.
#' @param reps Number of bootstrap repetitions for EWCE.
#' For publishable results set >=10000.
#' @param ctd_species Specify species used to build CTD if different from
#' species. If multiple CTDs are used with differing species specify as a
#' named list where names are the CTD file names (without file extension) and
#' elements are the species. Must be 'human', 'mouse' or listed in
#' EWCE::list_species()$id.
#' @param num_markers Number of cluster markers used to identify cell type
#' @param label_index The annotation_level index to use to name the
#' cluster_celltype. Default is 1
#'
#' @return sce a SingleCellExperiment object annotated with celltypes/metadata
#' @author Nathan Skene / Combiz Khozoie / Michael Thomas
#' @family Celltype annotation
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom Matrix Matrix
#' @importFrom EWCE generate_celltype_data
#' @importFrom assertthat assert_that
#' @importFrom dplyr rename inner_join
#' @importFrom plyr adply
#' @importFrom magrittr %>%
#' @importFrom purrr map_chr
#' @importFrom SingleCellExperiment counts
#' @export

map_celltypes_sce <- function(sce,
                              ctd_folder,
                              cells_to_sample = 10000,
                              clusters_colname = "clusters",
                              species = getOption(
                                "scflow_species",
                                default = "human"
                              ),
                              save_path = tempdir(),
                              annotation_level = 1,
                              reps = 1000,
                              ctd_species = NULL,
                              num_markers = 50,
                              label_index = 1) {
  assertthat::assert_that(dir.exists(ctd_folder))
  assertthat::assert_that(
    clusters_colname %in% names(SummarizedExperiment::colData(sce)),
    msg = "clusters_colname missing from colData"
  )
  assertthat::assert_that(
    species %in%
      c("human", "mouse", EWCE::list_species(verbose = FALSE)$id),
    msg = "species not supported for EWCE mapping"
  )
  if (!is.null(ctd_species)) {
    assertthat::assert_that(
      all(ctd_species %in%
        c("human", "mouse", EWCE::list_species(verbose = FALSE)$id)),
      msg = "ctd_species not supported for EWCE mapping"
    )
  } else {
    ctd_species <- species
  }
  assertthat::assert_that(
    num_markers >= 4,
    msg = "EWCE requires a minumum of 4 markers
           per cluster to indentify cell types"
  )

  assertthat::assert_that(
    if (class(annotation_level) == "list") {
      !is.null(names(annotation_level))
    } else if (class(annotation_level) == "numeric") {
      length(annotation_level) == 1
    } else {
      FALSE
    },
    msg = "Annotation level must be an integer for
           single CTD annotation or a named list for
           multi-level/CTD annotation"
  )

  l_ctd <- .read_rds_files_to_list(ctd_folder)


  if(!is.null(cells_to_sample)){
  set.seed(123)
  idx <- c()
  for(i in as.character(unique(sce[[clusters_colname]]))){
    idx_temp <- which( sce[[clusters_colname]] == i) %>%
      sample(size = min(cells_to_sample, length(.)), replace = FALSE )
    idx_temp <- idx_temp[order(idx_temp)]
    idx <- c(idx, idx_temp)}

    message(sprintf("Subsetting %s cells", cells_to_sample))

    sce_subset <- sce[, idx]

  } else {
    sce_subset <- sce
  }

  mat <- Matrix::Matrix(SingleCellExperiment::counts(sce_subset), sparse = TRUE)
  rownames(mat) <- as.character(
    SummarizedExperiment::rowData(sce_subset)$gene
  )

  annotLevels <- list(level1class = sce_subset[[clusters_colname]])
  message("generating ctd with ewce")
  ctd <- suppressMessages(EWCE::generate_celltype_data(
    exp = mat,
    annotLevels = annotLevels,
    groupName = "ctd",
    savePath = save_path
  ), classes = "message")
  load(ctd[1])

  sce@metadata$ctd <- ctd

  message("mapping celltypes with ewce")
  mappings <- .map_celltypes_with_ewce(ctd,
    l_ctd,
    inputSpecies = species,
    ctd_species = ctd_species,
    annotation_level = annotation_level,
    reps = reps,
    num_markers = num_markers
  )

  # generate lookup for merging into colData
  mappings_lookup <- mappings %>%
    dplyr::select(Cluster, Mapped, mapping) %>%
    tidyr::spread(mapping, Mapped) %>%
    dplyr::arrange(Cluster)

  mappings_lookup$Cluster <- as.factor(mappings_lookup$Cluster)

  # rename Cluster column
  mappings_df <- mappings_lookup %>%
    dplyr::rename(!!clusters_colname := Cluster)
  mappings_df[] <- lapply(mappings_df, as.character)

  if (class(annotation_level) == "list") {
    mapping_column <- paste0(
      names(annotation_level)[1],
      "_ML",
      annotation_level[[1]][label_index]
    )
  } else {
    mapping_column <- names(mappings_df)[2]
  }
  mappings_df <- mappings_df %>%
    dplyr::mutate(cluster_celltype = get(mapping_column),
                  cluster_celltype = gsub("_", "-", cluster_celltype))

  sce <- map_custom_celltypes(
    sce = sce,
    mappings = mappings_df
  )

  sce@metadata$mappings <- mappings_df

  sce <- .append_celltype_plots_sce(sce)

  return(sce)
}


################################################################################
#' Map cluster celltypes with EWCE
#'
#' Maps cluster celltypes using EWCE
#'
#' @param ctd ctd for current expt
#' @param l_ctd a list of reference ctd's
#' @param inputSpecies species of ctd
#' @param annotation_level level of annotation
#' @param ctd_species species of l_ctd
#' @param reps number pf repitions performed by EWCE boostrap
#'
#' @return sce a SingleCellExperiment object annotated with sample metadata
#' @author Nathan Skene / Combiz Khozoie
#' @family clustering and dimensionality reduction
#' @importFrom purrr pmap
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @keywords internal
.map_celltypes_with_ewce <- function(ctd,
                                     l_ctd,
                                     inputSpecies = "human",
                                     annotation_level = 1,
                                     ctd_species = "human",
                                     reps = 1000,
                                     num_markers = 50) {
  # set up the map_celltypes parameter combinations for pmap

  map_params <- list(
    ctd = names(l_ctd),
    ctd_species = ctd_species, # get species from ctd
    mlevel = annotation_level
    # need to include user specified annotaion levels
  ) # mapping levels

  l_mappings <- purrr::pmap(map_params, function(ctd_name,
                                                 ctd_species,
                                                 mlevel) {
    message(ctd_name)

    mapped_dt <- .map_celltypes(
      ctdToMap = ctd,
      ctdToMapAgainst = l_ctd[[ctd_name]],
      inputSpecies = inputSpecies,
      mapAgainstSpecies = ctd_species,
      annotLevel = 1,
      numTopMarkers = num_markers,
      mappingLevel = mlevel,
      reps = reps
    )
    mapped_dt$mapping <- sprintf("%s_ML%s", ctd_name, mapped_dt$mappingLevel)
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
#' @importFrom EWCE bootstrap_enrichment_test
#' @importFrom purrr pmap
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @keywords internal

### Allow variation in annotation level

.map_celltypes <- function(ctdToMap,
                           ctdToMapAgainst,
                           inputSpecies = "human",
                           mapAgainstSpecies = "human",
                           annotLevel = 1,
                           numTopMarkers = 50,
                           mappingLevel = 2,
                           reps = 1000) {
  ctd <- ctdToMapAgainst
  count <- 0
set.seed(123)
for (ml in mappingLevel) {
    message(" Level ", ml)
    for (ct in colnames(ctdToMap[[annotLevel]]$specificity)) {
      mostSpecificGenes <- .get_x_most_specific_genes(
        ct = ct,
        annotLevel = annotLevel, # Annotation level
        howMany = numTopMarkers,
        ctd = ctdToMap,
        exprPercentile = 0.9
      )$x_most_specific

      mostSpecificGenes <- mostSpecificGenes[!is.na(mostSpecificGenes)]

      # Handle check_ewce_inputs errors
      ewce_check <- tryCatch(
        EWCE::check_ewce_genelist_inputs(
          sct_data = ctd,
          hits = mostSpecificGenes,
          bg = rownames(ctdToMap[[1]]$specificity),
          genelistSpecies = inputSpecies,
          sctSpecies = mapAgainstSpecies,
          verbose = FALSE
        ),
        error = function(e) {
          message(paste(
            "Cannot map cluster",
            ct, "
                        to any cell type in",
            ctd_name
          ))
          message(paste("ERROR MESSAGE:", e))
          return(e)
        }
      )
      ewce_check <- !any(class(ewce_check) == "error")

      if (!ewce_check) {
        ctMapped <- data.frame(
          Original = ct,
          Mapped = "Not identified",
          annotLevel = ml,
          p = NA,
          z = NA
        )
      } else {
        full_results <- suppressMessages(EWCE::bootstrap_enrichment_test(
          sct_data = ctd,
          hits = mostSpecificGenes,
          annotLevel = ml,
          genelistSpecies = inputSpecies,
          sctSpecies = mapAgainstSpecies,
          verbose = FALSE
        ), classes = "message")
        p <- full_results$results[order(
          full_results$results$sd_from_mean,
          decreasing = TRUE
        ), ][1, ]$p
        z <- full_results$results[order(
          full_results$results$sd_from_mean,
          decreasing = TRUE
        ), ][1, ]$sd_from_mean
        ctMapped <- data.frame(
          Original = ct,
          Mapped = as.character(
            full_results$results[order(
              full_results$results$sd_from_mean,
              decreasing = TRUE
            ), ][1, ]$CellType
          ),
          mappingLevel = ml,
          p = p,
          z = z
        )
        cat(
          paste(c("  Cluster", "="), ctMapped[, c("Original", "Mapped")]),
          "\n"
        )
      }

      count <- count + 1
      if (count == 1) {
        allMapped <- ctMapped
      } else {
        allMapped <- rbind(ctMapped, allMapped)
      }
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


### Annotaiton level
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
.read_rds_files_to_list <- function(folder_path) {
  rds_l <- purrr::map(
    list.files(folder_path, pattern = "ctd_", full.names = TRUE),
    readRDS
  )
  names(rds_l) <- tools::file_path_sans_ext(list.files(folder_path, pattern = "ctd_"))

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
    sce@metadata$celltype_plots[[reddim]] <- plot_reduced_dim(
      sce,
      feature_dim = celltype_dim,
      reduced_dim = reddim,
      label_clusters = TRUE
    )
  }

  # 3d
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
