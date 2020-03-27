library(scFlow)
library(dplyr)
library(magrittr)
library(ggplot2)
sce <- read_sce("~/Documents/nf-sc/results/celltype_mapped_sce/celltype_mapped_sce")

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


ggplot(results$unique_id_plots, aes(x = .data[[unique_id_var]], y = cells_pc)) +
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
