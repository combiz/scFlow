library(Seurat)
library(SeuratData)
data("pbmcsca")

my_fn <- function(...) {
  fnargs <- list(...)
  seu <- NormalizeData(fnargs$object)
  #seu <- ScaleData(seu)
  #seu <- Seurat::SCTransform(
  #  fargs$object,
  #  vars.to.regress = fargs$vars_to_regress_out,
  #  do.scale = FALSE)
  return(seu)
}

x <- my_fn(
  object = seu,
  vars_to_regress_out = NULL
)

pbmcsca@meta.data


seu <- Seurat::CreateSeuratObject(
  counts = fargs$mat,
  meta.data = fargs$seu_metadata
)

x <- NormalizeData(pbmcsca)

seufargs$mat
fargs$seu_metadata
