library(ggplot2)
library(ggpubr)
library(patchwork)

pca_mat <- SingleCellExperiment::reducedDim(sce, "PCA")
plot_dt <- cbind(pca_mat, as.data.frame(SummarizedExperiment::colData(sce)))
#sce$RIN <- as.numeric(as.character(sce$RIN))
plot_l <- list()
y_var <- "cngeneson"
for (x_var in c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")) {
  r_val <- cor(plot_dt[[x_var]], plot_dt[[y_var]], method = "pearson")
  plot_l[[x_var]] <- ggplot(plot_dt, aes(x=.data[[x_var]], y=.data[[y_var]])) +
    geom_point(size = 1) +
    #ggtitle("NR3C1 vs NEK7") +
    geom_smooth(method=lm, se=FALSE) +
    stat_cor(label.y = 4, size = 6)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
    stat_regline_equation(label.y = 3.5, size = 6)+ #this means at 30th unit regresion line equation will be shown
    #scale_x_continuous(name = "NR3C1", limits = c(5, 15), breaks = seq(5, 15, 2)) +
    #scale_y_continuous(name = "NEK7", limits = c(5,15), breaks = seq(5, 15, 2)) +
    #annotation_custom(grob3) +
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      axis.line = element_line(color="black"),
      axis.line.x = element_line(color="black"),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      )
}

p <- Reduce(`+`, plot_l)
p

###
r_dt <- data.frame()
y_var <- "cngeneson"
pcs <- paste0("PC", seq(1, 30))
for (x_var in pcs) {
  r_val <- cor(plot_dt[[x_var]], plot_dt[[y_var]])
  new_row <- data.frame(x_var = x_var, y_var = y_var, r = r_val)
  r_dt <- rbind(r_dt, new_row)
}
r_dt$x_var <- factor(r_dt$x_var, levels = rev(pcs))

ggplot(r_dt, aes(x = x_var, y = r)) +
  geom_col() +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "grey80")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60")+
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey40")+
  geom_hline(yintercept = -0.3, linetype = "dashed", color = "grey80")+
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "grey60")+
  geom_hline(yintercept = -0.7, linetype = "dashed", color = "grey40")+
  scale_y_continuous(limits = c(-1, 1))+
  ylab("Correlation coefficient (r)")+
  theme_bw()+
  ggtitle(y_var) +
  coord_flip() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    axis.line = element_line(color="black"),
    axis.line.x = element_line(color="black"),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_blank()
  )
