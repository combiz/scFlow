p <- plot_x_y(
     sce,
     group_var = "log_p_tau",
     subset_var = "cluster_celltype",
     subset_group = "Micro",
     gene = "ABCA1",
     size = 0.5
   )
p

p <- plot_x_y(
  sce,
  group_var = "log_p_tau",
  subset_var = "cluster_celltype",
  subset_group = "Micro",
  gene = "MEF2A", # up 3x
  size = 0.5
)
p

p <- plot_x_y(
  sce,
  group_var = "log_p_tau",
  subset_var = "cluster_celltype",
  subset_group = "Micro",
  gene = "PRKAG2", # down 20%
  size = 0.5
)
p

p <- plot_x_y(
  sce,
  group_var = "log_p_tau",
  subset_var = "cluster_celltype",
  subset_group = "Micro",
  gene = "CD163",
  size = 0.5
)
p

p <- plot_x_y(
  sce,
  group_var = "log_p_tau",
  subset_var = "cluster_celltype",
  subset_group = "Micro",
  gene = "CD74",
  size = 0.5
)
p
HPSE2
sce$log_hist <- sce$log_p_tau
sce$hist <- sce$p_tau

#for (gene in c("OR3A2","MAP3K8","TLN2","LIMK2","KLHL2","FAM20A","TBC1D8","JAK3","IL4R","SETBP1","MRC1","BCL6","ADGRG1","HIVEP2","CCSER2")) { #down
for (gene in c("DPYD","STARD13","HP1BP3","ATG16L2","ASAH1","MYO1E","SAMD12","ATG7","ITM2B","PDGFC","MITF","DLEU7","TANC2","PLCL1","ZFP36L2")) { #up
  direction <- "UP"
  p <- plot_x_y(
    sce,
    group_var = "log_p_tau",
    subset_var = "cluster_celltype",
    subset_group = "Micro",
    gene = gene,
    size = 0.5
  )
  p
  ggsave(sprintf("Micro_%s_%s.png", direction, gene), width = 4.5, height = 6.5)
}


#for (gene in c("OR3A2","MAP3K8","TLN2","LIMK2","KLHL2","FAM20A","TBC1D8","JAK3","IL4R","SETBP1","MRC1","BCL6","ADGRG1","HIVEP2","CCSER2")) { #down
for (gene in c("MAP3K8","TLN2","LIMK2","KLHL2","FAM20A","TBC1D8","JAK3","IL4R","SETBP1","MRC1","BCL6","ADGRG1","HIVEP2","CCSER2")) { #down
  direction <- "DOWN"
  p <- plot_expr_by_numeric_var(
    sce,
    numeric_var = "p_tau",
    subset_var = "cluster_celltype",
    subset_group = "Micro",
    gene = gene,
    size = 0.5
  )
  p
  ggsave(sprintf("Micro_%s_%s.png", direction, gene), width = 4.5, height = 6.5)
}

sce$log_amyloid_beta <- log10(sce$amyloid_beta + 1)
p <- plot_expr_by_numeric_var(
  sce,
  #group_var = "log_amyloid_beta",
  numeric_var = "amyloid_beta",
  subset_var = "cluster_celltype",
  subset_group = "Micro",
  gene = "TLN2",
  size = 0.5
)
p

p <- plot_expr_by_numeric_var(
  sce,
  numeric_var = "log_amyloid_beta",
  subset_var = "cluster_celltype",
  subset_group = "Micro",
  gene = "GPNMB",
  size = 0.5
)
p




