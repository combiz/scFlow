
#  with batch
var <- "total_features_by_counts"
x <- scater::isOutlier(
  sce[[var]],
  nmads = 5,
  batch = sce[["manifest"]],
  type = "both",
  log = FALSE)
x <- data.frame("higher" = as.numeric(attributes(x)$thresholds["higher",]), "id" = as.character(names(attributes(x)$thresholds["higher",])))
#sce@metadata$merged_plots[[var]][[var]] + geom_hline(yintercept=attributes(x)$thresholds[["higher"]], colour = "red")

#sce@metadata$merged_plots[[var]][[var]] + geom_(data = x, aes(x = id, y = higher), size = 12, colour = "red")
sce@metadata$merged_plots[[var]][[var]] + geom_segment(data = x, aes(x = as.numeric(id) - 0.5, y = higher, xend = as.numeric(id) + 0.5, yend = higher), colour = "red", size = 1)

var <- "total_counts"
y <- scater::isOutlier(
  sce[[var]],
  nmads = 5,
  batch = sce[["manifest"]],
  type = "both",
  log = FALSE)
y <- data.frame("higher" = as.numeric(attributes(y)$thresholds["higher",]), "id" = as.character(names(attributes(y)$thresholds["higher",])))
#sce@metadata$merged_plots[[var]][[var]] + geom_hline(yintercept=attributes(y)$thresholds[["higher"]], colour = "red")
sce@metadata$merged_plots[[var]][[var]] + geom_segment(data = y, aes(x = as.numeric(id) - 0.5, y = higher, xend = as.numeric(id) + 0.5, yend = higher), colour = "red", size = 1)

z <- scater::isOutlier(
  sce[[var]],
  nmads = 3,
  batch = sce[["cluster_celltype"]],
  type = "both",
  log = FALSE)

outfeats <- scater::isOutlier(
  sce[[var]],
  nmads = 3,
  batch = sce[["cluster_celltype"]],
  type = "both",
  log = FALSE)

#table(x,y)
badsce <- sce[, x | y]
badsce_feat <- sce[, x]
badsce_counts <- sce[, y]
badsce_ct <- sce[, z]
table(badsce$cluster_celltype)

table(badsce_ct$cluster_celltype)




data.frame(cluster_celltype = sce$cluster_celltype, pc_mito = sce$pc_mito) %>%
  group_by(cluster_celltype) %>%
  summarize(mean_percentage_mito_counts = mean(pc_mito * 100))


# without batch
var <- "total_features_by_counts"
x <- scater::isOutlier(sce[[var]], nmads = 5, batch = sce[["manifest"]], type = "higher")
x <- data.frame("higher" = as.numeric(attributes(x)$thresholds["higher",]), "id" = as.character(names(attributes(x)$thresholds["higher",])))
#sce@metadata$merged_plots[[var]][[var]] + geom_hline(yintercept=attributes(x)$thresholds[["higher"]], colour = "red")

library(DropletUtils)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DropletUtils")

set.seed(100)

ed_out <- emptyDrops(counts(sce))
summary(ed_out$FDR <= 0.001)

table(Sig=ed_out$FDR <= 0.001, Limited=ed_out$Limited)



limit <- 100
all_out <- emptyDrops(counts(sce), lower=limit, test.ambient=TRUE)

dt <- data.frame("pval" = all_out$PValue[all_out$Total <= limit & all_out$Total > 0])
brx <- pretty(range(dt$pval),
              n = nclass.Sturges(dt$pval),min.n = 1)


hist(all_out$PValue[all_out$Total <= limit & all_out$Total > 0],
     xlab="P-value", main="", col="grey80")
