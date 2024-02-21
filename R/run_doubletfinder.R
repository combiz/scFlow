################################################################################
#' Find singlets in a SingleCellExperiment with DoubletFinder by McGinnis CS
#' et al
#'
#' Runs the doubletfinder algorithm on the SingleCellExperiment. Returns a
#' SingleCellExperiment annotated with is_singlet and tSNE coordinates.
#'
#' Remember to cite: -
#' DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using
#' Artificial Nearest Neighbors.
#' McGinnis CS, Murrow LM, Gartner ZJ.
#' Cell Syst. 2019 Apr 24;8(4):329-337.e4.
#' doi: \url{https://doi.org/10.1016/j.cels.2019.03.003}.
#'
#' @param sce a SingleCellExperiment object
#' @param pK a pK value to use in place of a parameter sweep
#' @param pca_dims the number of principal components to use
#' @param var_features the top n variable features to use
#' @param vars_to_regress_out the variables to regress out
#' @param doublet_rate either a fixed doublet rate (e.g. 0.075) or 0 to
#' @param dpk doublets per thousand cells increment if doublet_rate is 0.
#' @param num.cores the number of CPU cores to use
#'
#' @return sce a SingleCellExperiment object annotated for singlets
#'
#' @family annotation functions
#'
#' @importFrom cli cli_alert_danger cli_alert_success cli_h2 cli_h3
#' @importFrom cli cli_alert cli_text
#' @import ggplot2
#' @importFrom dplyr rename_all
#' @importFrom SingleCellExperiment counts reducedDim
#' @importFrom purrr map_int map_df
#' @importFrom DoubletFinder paramSweep summarizeSweep find.pK
#' @importFrom DoubletFinder modelHomotypic doubletFinder
#' @importFrom magrittr %>%
#' @importFrom Seurat FindVariableFeatures NormalizeData CreateSeuratObject
#' @importFrom Seurat ScaleData RunPCA RunTSNE RunUMAP
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom stats setNames na.omit
#' @importFrom sctransform vst get_residual_var get_residuals correct_counts
#' @importFrom future.apply future_lapply
#' @importFrom future availableCores
#'
#' @keywords internal

run_doubletfinder <- function(sce,
                              pK = NULL,
                              pca_dims = 10,
                              var_features = 2000,
                              vars_to_regress_out = "nCount_RNA",
                              doublet_rate = 0,
                              dpk = 8,
                              num_cores = max(1, future::availableCores() - 2)) {
  cli::cli_h2("Finding Singlets with DoubletFinder")
  cli::cli_h3("Pre-processing")

  # add citations and print
  sce <- .append_citation_sce(sce, key = c("doubletfinder", "seurat"))

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  cli::cli_alert("Selecting QC passed cells/genes")
  col_idx <- which(sce$qc_metric_passed)
  sce_ss <- filter_sce(
    sce,
    filter_genes = TRUE,
    filter_cells = TRUE
  )

  # calculate estimated doublet rate
  if (doublet_rate == 0) {
    abs_dpk <- (dim(sce_ss)[[2]] / 1000) * dpk
    doublet_rate <- abs_dpk / 1000
    cli::cli_text(
      "Estimating a doublet rate of {.val {doublet_rate}} ",
      "based on a doublets-per-thousand (dpk) increment of {.val {dpk}}"
    )
  } else {
    dpk <- NULL
  }

  sce@metadata$doubletfinder_params <- list(
    pK = pK,
    pca_dims = pca_dims,
    var_features = var_features,
    vars_to_regress_out = vars_to_regress_out,
    doublet_rate = doublet_rate,
    dpk = dpk
  )

  cli::cli_alert(c(
    "Selected {.val {dim(sce_ss)[[2]]}} cells ",
    "and {.val {dim(sce_ss)[[1]]}} genes"
  ))

  seu_metadata <- data.frame(sce_ss@colData) %>%
    droplevels()

  seu <- Seurat::CreateSeuratObject(
    counts = SingleCellExperiment::counts(sce_ss),
    meta.data = seu_metadata
  )

  cli::cli_alert("Finding variable features")

  seu <- Seurat::FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = var_features
  )

  cli::cli_alert("Normalizing data")
  seu <- Seurat::NormalizeData(seu)


  cli::cli_alert("Scaling data")

  seu <- Seurat::ScaleData(
    seu,
    features = NULL,
    vars.to.regress = vars_to_regress_out,
    model.use = "linear",
    use.umi = FALSE,
    do.scale = TRUE,
    do.center = TRUE,
    scale.max = Inf,
    block.size = 750,
    min.cells.to.block = 3000,
    verbose = TRUE
  )

  cli::cli_alert("Calculating PCA reduced dimensions")
  seu <- Seurat::RunPCA(seu)
  cli::cli_alert("Calculating tSNE reduced dimensions")
  seu <- Seurat::RunTSNE(seu, dims = 1:pca_dims)
  cli::cli_alert("Calculating UMAP reduced dimensions")
  seu <- Seurat::RunUMAP(seu, dims = 1:pca_dims)

  # pK Identification -------------------------------------------------------
  if (is.null(pK)) { # if not specified, use sweep
    cli::cli_h3("Performing optimal pK parameter sweep")
    sweep_res_list <- DoubletFinder::paramSweep(
      seu,
      PCs = 1:pca_dims, sct = FALSE,
      num.cores = num_cores
    )
    sweep_stats <- DoubletFinder::summarizeSweep(sweep_res_list, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep_stats)
    bcmvn$pK <- as.numeric(as.character(bcmvn$pK)) # oddly, pK are factors
    pK <- bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric), ]$pK
    sce <- .doublet_finder_plot_param_sweep(sce, bcmvn)
    # add the bcmvn dataframe to the sce metadata
    sce@metadata$qc_plot_data$doubletfinder_param_sweep <- bcmvn
    sce@metadata$doubletfinder_params$doubletfinder_sweep <- TRUE
  } else {
    cli::cli_text("Skipping parameter sweep and using pK={.value {pK}}.")
    sce@metadata$doubletfinder_params$doubletfinder_sweep <- FALSE
  }
  # add the pK value used to metadata
  sce@metadata$doubletfinder_params$pK <- pK

  # Homotypic Doublet Proportion Estimate -----------------------------------
  cli::cli_alert("Estimating homotypic doublet proportions")
  annotations <- seu@meta.data$seurat_clusters
  homotypic_prop <- DoubletFinder::modelHomotypic(annotations)
  n_exp_poi <- round(doublet_rate * length(seu@meta.data$orig.ident))
  n_exp_poi_adj <- round(n_exp_poi * (1 - homotypic_prop))
  doublet_vals <- data.frame(
    n_exp_poi = n_exp_poi,
    n_exp_poi_adj = n_exp_poi_adj
  )

  pann_hi <- sprintf("pANN_0.25_%s_%s", pK, doublet_vals$n_exp_poi)

  # Run DoubletFinder with varying classification stringencies -------------
  cli::cli_h3("Running DoubletFinder")
  seu <- DoubletFinder::doubletFinder(seu,
    PCs = 1:pca_dims, pN = 0.25, pK = as.numeric(pK),
    nExp = n_exp_poi,
    reuse.pANN = FALSE
  )

  seu <- DoubletFinder::doubletFinder(seu,
    pN = 0.25, pK = as.numeric(pK),
    nExp = doublet_vals$n_exp_poi_adj,
    reuse.pANN = pann_hi
  )

  # Annotate results -------------------------------------------------------
  cli::cli_alert("Annotating results")
  seu@meta.data[, "DF_hi.lo"] <- seu@meta.data[, sprintf(
    "DF.classifications_0.25_%s_%s",
    pK,
    doublet_vals$n_exp_poi
  )]

  seu@meta.data$DF_hi.lo[which(seu@meta.data$DF_hi.lo == "Doublet" &
    seu@meta.data[, sprintf(
      "DF.classifications_0.25_%s_%s",
      pK,
      doublet_vals$n_exp_poi_adj
    )] == "Singlet")] <- "Doublet_lo"

  seu@meta.data$DF_hi.lo[which(seu@meta.data$DF_hi.lo == "Doublet")] <-
    "Doublet_hi"

  # prepare data to return
  # sce$is_singlet <- seu@meta.data$DF_hi.lo == "Singlet"
  is_singlet <- rep(NA, dim(sce)[[2]])
  is_singlet[col_idx] <- seu@meta.data$DF_hi.lo == "Singlet"
  sce$is_singlet <- is_singlet
  sce@metadata$doubletfinder_params$singlets_found <-
    sum(sce$is_singlet, na.rm = TRUE)
  sce@metadata$doubletfinder_params$multiplets_found <-
    sum(!sce$is_singlet, na.rm = TRUE)

  # save dim_reductions and generate metadata plots
  for (dr_method in names(seu@reductions)) {
    dim_rd_name <- paste0(dr_method, "_by_individual")
    rd_mat <- as.matrix(seu@reductions[[dr_method]]@cell.embeddings)
    full_rd_mat <- matrix(
      data = NA,
      nrow = dim(sce)[[2]],
      ncol = dim(rd_mat)[[2]],
      dimnames = list(colnames(sce), colnames(rd_mat))
    )
    full_rd_mat[col_idx, ] <- rd_mat
    SingleCellExperiment::reducedDim(sce, dim_rd_name) <-
      full_rd_mat
    sce <- .doublet_finder_plot_dim_red(sce, dim_rd_name)
  }

  # append df to qc_summary metadata table

  doubletfinder_params <- sce@metadata$doubletfinder_params
  lidx <- which(purrr::map_int(doubletfinder_params, length) > 1)
  for (idx in lidx) {
    doubletfinder_params[[idx]] <- paste(doubletfinder_params[[idx]], collapse = ";")
  }

  doubletfinder_params <- purrr::map_df(
    unlist(doubletfinder_params), ~.
  ) %>%
    dplyr::rename_all(~ paste0("doubletfinder_", .))

  sce@metadata$qc_summary <- cbind(
    sce@metadata$qc_summary,
    doubletfinder_params
  )

  sce$qc_metric_passed <- sce$qc_metric_passed & sce$is_singlet

  cli::cli_alert_success("DoubletFinder completed succesfully")

  return(sce)
}

################################################################################
#' temp fn
#'
#' @param sce a SingleCellExperiment object
#'
#' @return seu a seurat object
#'
#' @family annotation functions
#'
#' @importFrom SingleCellExperiment counts
#' @import Seurat
#' @importFrom magrittr %>%
#'
#' @keywords internal
.preprocess_seurat_object <- function(sce,
                                      vars_to_regress_out,
                                      pca_dims,
                                      var_features) {
  cat("\r\n")
  cat(cli::rule("Creating SeuratObject", line = 1), "\r\n")

  seu_metadata <- data.frame(sce@colData) %>%
    droplevels()
  seu <- Seurat::CreateSeuratObject(
    counts = counts(sce),
    meta.data = seu_metadata
  )
  # doubletfinder not compatible yet, no slots error
  # seu <- Seurat::SCTransform(
  #  seu,
  #  vars.to.regress = vars_to_regress_out)
  cat(cli::rule("Normalizing data", line = 1), "\r\n")
  seu <- Seurat::NormalizeData(seu)
  cat(cli::rule("Scaling data", line = 1), "\r\n")
  cat(sprintf(
    "Regressing out %s, this may take a while..",
    vars_to_regress_out
  ))
  seu <- Seurat::ScaleData(
    seu,
    vars.to.regress = vars_to_regress_out
  )

  cat(cli::rule("Finding variable features", line = 1), "\r\n")
  seu <- Seurat::FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = var_features
  )
  cat(cli::rule("Calculating PCA reduced dimensions", line = 1), "\r\n")
  seu <- Seurat::RunPCA(seu)
  cat(cli::rule("Calculating tSNE reduced dimensions", line = 1), "\r\n")
  seu <- Seurat::RunTSNE(seu, dims = 1:pca_dims)
  cat(cli::rule("Calculating UMAP reduced dimensions", line = 1), "\r\n")
  seu <- Seurat::RunUMAP(seu, dims = 1:pca_dims)

  return(seu)
}


#' plot 2d dimensionality reduction with is_doublet colour
#' @export
#' @keywords internal
.doublet_finder_plot_dim_red <- function(sce,
                                         reduced_dim = NULL) {
  if (is.null(reduced_dim)) {
    stop("need a reduced dim slot, see reducedDims(sce).")
  }

  df <- data.frame(
    SingleCellExperiment::reducedDim(sce, reduced_dim)
  )
  colnames(df) <- paste0("dim_", seq_along(df))
  df$is_singlet <- sce$is_singlet
  df <- na.omit(df)

  p <- ggplot2::ggplot(data = df) +
    geom_point(aes(
      x = dim_1,
      y = dim_2,
      colour = is_singlet
    ), shape = 16, size = 1, alpha = .4) +
    scale_colour_manual(values = c("#E64B35", "#4DBBD5")) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      line = element_blank(),
      text = element_blank(),
      title = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 18, hjust = 0.5)
    )

  # p$plot_env <- rlang::new_environment()
  p <- .grobify_ggplot(p)
  sce@metadata$qc_plots$doublet_finder[[reduced_dim]] <- p

  return(sce)
}


#' pparam sweep plot, pK vs BCmetricx with BCmetric maxima highlighted
#' @export
#' @keywords internal
.doublet_finder_plot_param_sweep <- function(sce, bcmvn) {
  p <- ggplot2::ggplot(data = bcmvn, aes(x = pK, y = BCmetric)) +
    geom_line() +
    geom_point(size = 4) +
    geom_point(
      data = bcmvn[bcmvn$pK == sce@metadata$doubletfinder_params$pK, ],
      size = 4,
      colour = "red"
    ) +
    theme_bw() +
    theme(
      # panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      # line = element_blank(),
      text = element_text(size = 18, hjust = 0.5, colour = "black"),
      # title = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 18, hjust = 0.5)
    )

  p <- .grobify_ggplot(p)
  sce@metadata$qc_plots$doublet_finder$param_sweep <- p

  return(sce)
}
