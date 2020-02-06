library(dplyr)

dt <- as.data.frame(table(sce$manifest, sce$cluster_celltype) / rowSums((table(sce$manifest, sce$cluster_celltype))) * 100)

crit_z <- qnorm(1-.05/2)

dt_with_z <- dt %>%
  group_by(Var2) %>%
  mutate(outlier = scale(Freq) < -crit_z | scale(Freq) > crit_z)

dt_ct_summary <- dt %>%
  group_by(Var2) %>%
  summarize(mean = mean(Freq, na.rm = TRUE), sd = sd(Freq, na.rm = TRUE))

dt <- dplyr::left_join(dt_with_z, dt_ct_summary)

ggplot(dt) +
         geom_col(aes(fill = Var2, x = Var1, y = Freq), position = "fill", colour = "black") +
  theme_minimal()

ggplot(dt) +
  geom_col(aes(fill = Var2, x = Var1, y = Freq), position = "fill", colour = "black") +
  theme_minimal()


ggplot(dt) +
  geom_boxplot(aes(x = Var2, y = Freq)) +
  geom_text(data = dt[dt$outlier == TRUE,], aes(label = Var1, x = Var2, y = Freq), size = 6) +
  xlab("Cell Type")+
  ylab("Proportion (%)") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 18)
  )


ggplot(as.data.frame(table(dt[dt$outlier == TRUE,]$Var1))) +
  geom_col(aes(x = Var1, y = Freq), colour = "black") +
  xlab("") +
  ylab("Number of celltypes with outlying relative cell counts") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 18)
  )

dt[dt$outlier == TRUE,]

library(devtools)
install_github('theislab/kBET')
library(kBET)

x <- kBET(SingleCellExperiment::reducedDim(sce, "UMAP"), sce$manifest)


mat <- Matrix::Matrix(t(as.matrix(SingleCellExperiment::reducedDim(sce, "UMAP"))), sparse = FALSE)
mat <- as.matrix(mat)

x <- kBET(mat, sce$manifest)
