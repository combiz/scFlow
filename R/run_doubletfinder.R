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
#'
#' @return sce a SingleCellExperiment object annotated for singlets
#'
#' @family annotation functions
#' @import cli Matrix SummarizedExperiment dplyr SingleCellExperiment
#' @import Seurat DoubletFinder
#' @export

run_doubletfinder <- function(sce, ...) {

  cat(cli::rule("Finding Singlets with DoubletFinder", line = 1), "\r\n")

  # defaults
  args <- list(
    pK = NULL,
    doublet_rate = 0.075 # assume 7.5%
  )
  inargs <- list(...)
  args[names(inargs)] <- inargs #override defaults if provided

  cat(cli::boxx(c(
    "Remember to cite:",
    "DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing",
    "Data Using Artificial Nearest Neighbors.",
    "McGinnis CS, Murrow LM, Gartner ZJ",
    "Cell Syst. 2019 Apr 24;8(4):329-337.e4."),
    padding = 1, align = "center", float = "center"))

  if (class(sce) != "SingleCellExperiment") {
    stop(cli::cli_alert_danger("A SingleCellExperiment is required."))
  }

  mat <- counts(sce)
  colnames(mat) <- sce$barcode

  # Pre-process Seurat object -----------------------------------------------
  cat("\r\n")
  cat(cli::rule("Creating SeuratObject", line = 1), "\r\n")
  seu <- Seurat::CreateSeuratObject(mat)
  cat(cli::rule("Normalizing data", line = 1), "\r\n")
  seu <- Seurat::NormalizeData(seu)
  cat(cli::rule("Scaling data", line = 1), "\r\n")
  seu <- Seurat::ScaleData(seu)
  cat(cli::rule("Finding variable features", line = 1), "\r\n")
  seu <- Seurat::FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = 2000
  )
  cat(cli::rule("Calculating PCA reduced dimensions", line = 1), "\r\n")
  seu <- Seurat::RunPCA(seu, )
  cat(cli::rule("Calculating UMAP reduced dimensions", line = 1), "\r\n")
  seu <- Seurat::RunUMAP(seu, dims = 1:10)

  # pK Identification -------------------------------------------------------
  if (is.null(args$pK)) { # if not specified, use sweep
    cat(cli::rule(
      "Identifying optimal pK with parameter sweep", line = 1), "\r\n")
    sweep_res_list <- DoubletFinder::paramSweep_v3(seu, PCs = 1:10)
    sweep_stats <- DoubletFinder::summarizeSweep(sweep_res_list, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep_stats)
    bcmvn$pK <- as.numeric(as.character(bcmvn$pK)) # oddly, pK are factors
    args$pK <- bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric), ]$pK
  } else {
    cli::cli_text("Skipping parameter sweep and using pK={.value {args$pK}}.")
  }

  # Homotypic Doublet Proportion Estimate -----------------------------------
  cat(cli::rule("Estimating homotypic doublet proportions", line = 1), "\r\n")
  annotations <- seu@meta.data$seurat_clusters
  homotypic_prop <- DoubletFinder::modelHomotypic(annotations)
  n_exp_poi <- round(args$doublet_rate * length(seu@meta.data$orig.ident))
  n_exp_poi_adj <- round(n_exp_poi * (1 - homotypic_prop))
  doublet_vals <- data.frame(
    n_exp_poi = n_exp_poi,
    n_exp_poi_adj = n_exp_poi_adj
  )

  pann_hi <- sprintf("pANN_0.25_%s_%s", args$pK, doublet_vals$n_exp_poi)

  # Run DoubletFinder with varying classification stringencies -------------
  cat(cli::rule("Running DoubletFinder", line = 1), "\r\n")
  seu <- DoubletFinder::doubletFinder_v3(seu,
    PCs = 1:10, pN = 0.25, pK = as.numeric(args$pK),
    nExp = n_exp_poi,
    reuse.pANN = FALSE
  )

  seu <- DoubletFinder::doubletFinder_v3(seu,
    pN = 0.25, pK = as.numeric(args$pK),
    nExp = doublet_vals$n_exp_poi_adj,
    reuse.pANN = pann_hi
  )

  # Annotate results -------------------------------------------------------
  cat(cli::rule("Annotating results", line = 1), "\r\n")
  seu@meta.data[, "DF_hi.lo"] <- seu@meta.data[, sprintf(
    "DF.classifications_0.25_%s_%s",
    args$pK,
    doublet_vals$n_exp_poi
  )]

  seu@meta.data$DF_hi.lo[which(seu@meta.data$DF_hi.lo == "Doublet" &
    seu@meta.data[, sprintf(
      "DF.classifications_0.25_%s_%s",
      args$pK,
      doublet_vals$n_exp_poi_adj
    )] == "Singlet")] <- "Doublet_lo"

  seu@meta.data$DF_hi.lo[which(seu@meta.data$DF_hi.lo == "Doublet")] <-
    "Doublet_hi"

  # prepare data to return
  sce$is_singlet <- seu@meta.data$DF_hi.lo == "Singlet"
  SingleCellExperiment::reducedDim(sce, "seurat_umap_by_individual") <-
    as.matrix(seu@reductions$umap@cell.embeddings)

  cli::cli_alert_success("DoubletFinder completed succesfully")

  return(sce)

}
