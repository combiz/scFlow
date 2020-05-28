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
#' @param ... 'pK' passed here will override parameter sweep,
#' 'doublet_rate' passed here will override the default of 7.5\%
#' a doublet_rate of 0 is a special case that automates the estimated dr
#'
#' @return sce a SingleCellExperiment object annotated for singlets
#'
#' @family annotation functions
#' @import cli Matrix dplyr SingleCellExperiment
#' @import DoubletFinder ggplot2 purrr
#' @importFrom magrittr %>%
#' @import Seurat
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom stats setNames na.omit
#' @importFrom sctransform vst get_residual_var get_residuals correct_counts
#' @importFrom future.apply future_lapply
#' @importFrom future availableCores
#'
#' @keywords internal

run_doubletfinder <- function(sce, ...) {

  print(sys.nframe())
  print(sys.frames())

  cat(cli::rule("Finding Singlets with DoubletFinder", line = 1), "\r\n")

  # defaults
  fargs <- list(
    pK = NULL,
    pca_dims = 10,
    vars_to_regress_out = "nCount_RNA",
    var_features = 2000,
    doublet_rate = 0,
    dpk = 8, # estimated doublets per thousand cells
    num.cores = future::availableCores()
  )
  inargs <- list(...)
  fargs[names(inargs)] <- inargs #override defaults if provided



  # add citations and print
  sce <- .append_citation_sce(sce, key = c("doubletfinder", "seurat"))

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  cat(cli::rule("Selecting QC passed cells/genes", line = 1), "\r\n")
  col_idx <- which(sce$qc_metric_passed)
  sce_ss <- filter_sce(
    sce,
    filter_genes = TRUE,
    filter_cells = TRUE
  )

  # calculate estimated doublet rate
  if(fargs$doublet_rate == 0) {
    abs_dpk <- (dim(sce_ss)[[2]] / 1000) * fargs$dpk
    fargs$doublet_rate <-  abs_dpk / 1000
    cli::cli_text(
      "Estimating a doublet rate of {.var {fargs$doublet_rate}} ",
      "based on a doublets-per-thousand (dpk) increment of {.var {fargs$dpk}}")
  } else {
    fargs$dpk <- NULL
  }

  sce@metadata$doubletfinder_params <- fargs

  cat(cli::cli_text(
    "Selected {.var {dim(sce_ss)[[2]]}} cells ",
    "and {.var {dim(sce_ss)[[1]]}} genes"),
    "\r\n")

  # Pre-process Seurat object -----------------------------------------------
  #seu <- do.call(
  #  .preprocess_seurat_object,
  #  list(sce = sce_ss,
  #       vars_to_regress_out = fargs$vars_to_regress_out,
  #       pca_dims = fargs$pca_dims,
  #       var_features = fargs$var_features)
  #)
  #seu <- .preprocess_seurat_object(
  #  sce = sce_ss,
  #  vars_to_regress_out = fargs$vars_to_regress_out,
  #  pca_dims = fargs$pca_dims,
  #  var_features = fargs$var_features
  #)

  seu_metadata <- data.frame(sce_ss@colData) %>%
    droplevels()

  #seu <- do.call(
  #  Seurat::CreateSeuratObject,
  #  list(counts = SingleCellExperiment::counts(sce_ss),
  #       meta.data = seu_metadata)
  #)
  seu <- Seurat::CreateSeuratObject(
    counts = SingleCellExperiment::counts(sce_ss),
    meta.data = seu_metadata
  )

  cat(cli::rule("Finding variable features", line = 1), "\r\n")
  #seu <- do.call(
    #Seurat::FindVariableFeatures,
    #list(object = seu,
    #     selection.method = "vst",
    #     nfeatures = fargs$var_features
    #     )
    #)
  seu <- Seurat::FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = fargs$var_features
  )

  cat(cli::rule("Normalizing data", line = 1), "\r\n")
  seu <- Seurat::NormalizeData(seu)
  #seu <- do.call(Seurat::NormalizeData,
  #               list(object = seu)
  #)

  cat(cli::rule("Scaling data", line = 1), "\r\n")

  seu <- Seurat::ScaleData(
    seu,
    vars.to.regress = fargs$vars_to_regress_out
  )
  print(ls())
  #seu <- do.call(
    #Seurat::ScaleData,
    #list(object = seu,
    #     vars.to.regress = fargs$vars_to_regress_out)
    #)

  cat(cli::rule("Calculating PCA reduced dimensions", line = 1), "\r\n")
  seu <- Seurat::RunPCA(seu)
  cat(cli::rule("Calculating tSNE reduced dimensions", line = 1), "\r\n")
  seu <- Seurat::RunTSNE(seu, dims = 1:fargs$pca_dims)
  cat(cli::rule("Calculating UMAP reduced dimensions", line = 1), "\r\n")
  seu <- Seurat::RunUMAP(seu, dims = 1:fargs$pca_dims)

  #seu <- sce_to_seu(sce_ss)

  # pK Identification -------------------------------------------------------
  if (is.null(fargs$pK)) { # if not specified, use sweep
    cat(cli::rule(
      "Identifying optimal pK with parameter sweep", line = 1), "\r\n")
    sweep_res_list <- DoubletFinder::paramSweep_v3(
      seu,
      PCs = 1:fargs$pca_dims, sct = FALSE,
      num.cores = fargs$num.cores
    )
    sweep_stats <- DoubletFinder::summarizeSweep(sweep_res_list, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep_stats)
    bcmvn$pK <- as.numeric(as.character(bcmvn$pK)) # oddly, pK are factors
    fargs$pK <- bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric), ]$pK
    sce <- .doublet_finder_plot_param_sweep(sce, bcmvn)
    # add the bcmvn dataframe to the sce metadata
    sce@metadata$qc_plot_data$doubletfinder_param_sweep <- bcmvn

    sce@metadata$doubletfinder_params$doubletfinder_sweep <- TRUE
  } else {
    cli::cli_text("Skipping parameter sweep and using pK={.value {fargs$pK}}.")
    sce@metadata$doubletfinder_params$doubletfinder_sweep <- FALSE
  }
  # add the pK value used to metadata
  sce@metadata$doubletfinder_params$pK <- fargs$pK

  # Homotypic Doublet Proportion Estimate -----------------------------------
  cat(cli::rule("Estimating homotypic doublet proportions", line = 1), "\r\n")
  annotations <- seu@meta.data$seurat_clusters
  homotypic_prop <- DoubletFinder::modelHomotypic(annotations)
  n_exp_poi <- round(fargs$doublet_rate * length(seu@meta.data$orig.ident))
  n_exp_poi_adj <- round(n_exp_poi * (1 - homotypic_prop))
  doublet_vals <- data.frame(
    n_exp_poi = n_exp_poi,
    n_exp_poi_adj = n_exp_poi_adj
  )

  pann_hi <- sprintf("pANN_0.25_%s_%s", fargs$pK, doublet_vals$n_exp_poi)

  # Run DoubletFinder with varying classification stringencies -------------
  cat(cli::rule("Running DoubletFinder", line = 1), "\r\n")
  seu <- DoubletFinder::doubletFinder_v3(seu,
    PCs = 1:fargs$pca_dims, pN = 0.25, pK = as.numeric(fargs$pK),
    nExp = n_exp_poi,
    reuse.pANN = FALSE
  )

  seu <- DoubletFinder::doubletFinder_v3(seu,
    pN = 0.25, pK = as.numeric(fargs$pK),
    nExp = doublet_vals$n_exp_poi_adj,
    reuse.pANN = pann_hi
  )

  # Annotate results -------------------------------------------------------
  cat(cli::rule("Annotating results", line = 1), "\r\n")
  seu@meta.data[, "DF_hi.lo"] <- seu@meta.data[, sprintf(
    "DF.classifications_0.25_%s_%s",
    fargs$pK,
    doublet_vals$n_exp_poi
  )]

  seu@meta.data$DF_hi.lo[which(seu@meta.data$DF_hi.lo == "Doublet" &
    seu@meta.data[, sprintf(
      "DF.classifications_0.25_%s_%s",
      fargs$pK,
      doublet_vals$n_exp_poi_adj
    )] == "Singlet")] <- "Doublet_lo"

  seu@meta.data$DF_hi.lo[which(seu@meta.data$DF_hi.lo == "Doublet")] <-
    "Doublet_hi"

  # prepare data to return
  #sce$is_singlet <- seu@meta.data$DF_hi.lo == "Singlet"
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
  for (idx in lidx){
    doubletfinder_params[[idx]] <- paste(doubletfinder_params[[idx]], collapse = ";")
  }

  doubletfinder_params <- purrr::map_df(
    unlist(doubletfinder_params), ~ .) %>%
    dplyr::rename_all(~ paste0("doubletfinder_",.))

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
  #seu <- Seurat::SCTransform(
  #  seu,
  #  vars.to.regress = vars_to_regress_out)
  cat(cli::rule("Normalizing data", line = 1), "\r\n")
  seu <- Seurat::NormalizeData(seu)
  cat(cli::rule("Scaling data", line = 1), "\r\n")
  cat(sprintf(
    "Regressing out %s, this may take a while..",
    vars_to_regress_out))
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

  if (is.null(reduced_dim)){
    stop("need a reduced dim slot, see reducedDims(sce).")
  }

  df <- data.frame(
    SingleCellExperiment::reducedDim(sce, reduced_dim)
  )
  colnames(df) <- paste0("dim_", seq_along(df))
  df$is_singlet <- sce$is_singlet
  df <- na.omit(df)

  p <- ggplot2::ggplot(data = df)+
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

  p$plot_env <- rlang::new_environment()
  sce@metadata$qc_plots$doublet_finder[[reduced_dim]] <- p

  return(sce)

}


#' pparam sweep plot, pK vs BCmetricx with BCmetric maxima highlighted
#' @export
#' @keywords internal
.doublet_finder_plot_param_sweep <- function(sce, bcmvn) {

  p <- ggplot2::ggplot(data = bcmvn, aes(x = pK, y = BCmetric))+
    geom_line() +
    geom_point(size = 4) +
    geom_point(
      data = bcmvn[bcmvn$pK == sce@metadata$doubletfinder_params$pK, ],
      size = 4,
      colour = "red") +
    theme_bw() +
    theme(
      #panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      #line = element_blank(),
      text = element_text(size = 18, hjust = 0.5, colour = "black"),
      #title = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 18, hjust = 0.5))

  p$plot_env <- rlang::new_environment()
  sce@metadata$qc_plots$doublet_finder$param_sweep <- p

  return(sce)

}


