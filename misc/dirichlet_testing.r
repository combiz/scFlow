library(scFlow)
library(dplyr)
library(magrittr)
library(ggplot2)
#sce <- read_sce("~/Documents/nf-sc/results/celltype_mapped_sce/celltype_mapped_sce")
sce <- read_sce("~/Documents/junk/mergedsce")
results <- model_celltype_freqs(sce, unique_id_var = "individual")
report_celltype_model(results)

qs::qsave(results, "test.qs")
qs::qsave(results$dirichlet_plot, "test.qs")

library(microbenchmark)
microbenchmark(saveRDS(results, "test.rds"), qs::qsave(results, "test.qs"), times=1)

paletteer::paletteer_dynamic("cartography::red.pal", 4)
results <- model_celltype_freqs(
  sce,
  unique_id_var = "individual",
  var_order = c("Control", "Low", "High"),
  palette = c("#2574A0FF", "#EA3B36FF", "#7C000CFF")
  )

#x <- results$unique_id_plot_table %>% dplyr::filter(cluster_celltype == "Oligo", group == "Control")

x <- results$unique_id_plot_table %>% dplyr::filter(cluster_celltype == "Oligo")

x$cells_pc

beta_l <- list()
for (group in unique(x$group)) {
  x <- results$unique_id_plot_table %>% dplyr::filter(cluster_celltype == "Oligo")
  m <- MASS::fitdistr(x$cells_pc, dbeta,
                    start = list(shape1 = 1, shape2 = 10))
  alpha0 <- m$estimate[1] # parameter 1
  beta0 <- m$estimate[2] # parameter 2
  p_layer <- geom_histogram(data = x, aes(cells_pc, y = ..density..), binwidth = .005)
  beta_l[[group]] <- list(
    x = x, m = m, alpha0 = alpha0, beta0 = beta0, p_layer = p_layer
  )
}

p <- ggplot(x) +
  geom_histogram(aes(cells_pc, y = ..density..), binwidth = .005) +
  stat_function(fun = function(x) dbeta(x, alpha0, beta0), color = "red",
                size = 1) +
  xlab("p(celltype)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
    axis.title = element_text(size = 18),
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black", hjust = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.justification = "left",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 16),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black")
  )

#####

results <- model_celltype_freqs(sce, var_order = c("Control", "Low", "High"), palette = c("#3B4992FF", "#EE0000FF", "#7C000CFF"))

results <- model_celltype_freqs(sce)

dependent_var <- "group"
celltype_var <- "cluster_celltype"
var_order <- c("Control", "Low", "High")
unique_id_var <- "manifest"


x <- results$dirichlet_plot_table
results$prop_counts_mat
x <- results$DR_data_df
x$counts <- NULL


x$group <- as.factor(x$group)
x[[dependent_var]] <- factor(
  x[[dependent_var]],
  levels = var_order
  )
x <- x[order(x[[dependent_var]]), ]
x$manifest <- factor(x$manifest, levels = x$manifest)
x <- x %>% tidyr::pivot_longer(
  cols = setdiff(colnames(x),
                 c(dependent_var, unique_id_var)),
  names_to = celltype_var,
  values_to = "cells_pc"
  )


ggplot(x, aes(x = .data[[unique_id_var]], y = cells_pc)) +
  geom_col(aes(fill = .data[[dependent_var]]), colour = "black") +
  ylab("Relative Proportion") +
  xlab(NULL) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
    axis.title = element_text(size = 18),
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black", angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.justification = "left",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 16),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black")
  ) +
  facet_wrap(~ .data[[celltype_var]], ncol = 3)


m <- MASS::fitdistr(career_filtered$average, dbeta,
                    start = list(shape1 = 1, shape2 = 10))

alpha0 <- m$estimate[1] # parameter 1
beta0 <- m$estimate[2] # parameter 2
