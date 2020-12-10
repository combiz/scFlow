# final
#sce <- read_sce("~/Documents/junk/enriched_with_ptau/")


library(parallel)
options(mc.cores = 10)
library(scFlow)
library(magrittr)
sce <- read_sce("~/Documents/junk/enriched/final_sce/")
new_metadata <- read.delim("~/Documents/junk/Enriched_Histology.csv", sep = ",")

# impute NAs
#preimpute <- new_metadata %>% dplyr::select(braak, amyloid_beta, p_tau, gfap, iba1, hla, gpnmb)
#preproc <- caret::preProcess(preimpute, method = "bagImpute")
#imputedata <- predict(preproc,preimpute)
#imputedata$individual <- new_metadata$individual
#imputedata$brain_region <- new_metadata$brain_region
#imputedata <- dplyr::filter(imputedata, individual %in% unique(sce$individual))
imputedata <- new_metadata
imputedata$individual <- gsub("/", "_", imputedata$individual, fixed = TRUE) # for seu only

old_df <- as.data.frame(SummarizedExperiment::colData(sce))
old_df$individual <- as.character(old_df$individual)
imputedata$individual <- as.character(imputedata$individual)
new_df <- dplyr::left_join(old_df, imputedata, by = c("individual", "brain_region"))
#table(new_df$braak, new_df$Braak_tangle_stage)

# check everything there
x <- paste0(sce$individual, sce$brain_region)
y <- paste0(imputedata$individual, imputedata$brain_region)
unique(x[!(x %in% y)])

# check
dim(new_df)

for (var in c("braak", "amyloid_beta", "p_tau", "gfap", "iba1", "hla", "gpnmb")) {
  sce[[var]] <- new_df[[var]]
}

write_sce(sce, "~/Documents/junk/enriched_with_ptau/")
