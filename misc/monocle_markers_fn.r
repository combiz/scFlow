library(scFlow)
sce <- read_sce("~/Documents/junk/MS_Custom_Mapped_SCE")
sce_all <- sce

idx <- as.numeric(caret::createDataPartition(sce$manifest, p = .05, list = FALSE))
sce <- sce[,idx]

