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
#' @importFrom Seurat RunPCA RunTSNE RunUMAP CreateSeuratObject SCTransform
#' @importFrom magrittr %>%
#'
#' @keywords internal
sce_to_seu <- function(sce){

  seu_metadata <- data.frame(sce@colData) %>%
    droplevels()
  seu <- Seurat::CreateSeuratObject(
    counts = counts(sce),
    meta.data = seu_metadata
  )
  seu <- Seurat::SCTransform(
    seu,
    vars.to.regress = NULL)

  cat(cli::rule("Calculating PCA reduced dimensions", line = 1), "\r\n")
  seu <- Seurat::RunPCA(seu)
  cat(cli::rule("Calculating tSNE reduced dimensions", line = 1), "\r\n")
  seu <- Seurat::RunTSNE(seu, dims = 1:10)
  cat(cli::rule("Calculating UMAP reduced dimensions", line = 1), "\r\n")
  seu <- Seurat::RunUMAP(seu, dims = 1:10)

  return(seu)
}
