# ms 2020

sce <- read_sce("~/Documents/junk/integrated_sce", read_metadata = TRUE)

sce <- reduce_dims_sce(
  sce,
  input_reduced_dim = "Liger",
  reduction_methods = "tSNE",
  dims= 2,
  initial_dims= 50,
  perplexity= 100,
  theta= 0.5,
  stop_lying_iter= 250,
  mom_switch_iter= 250,
  max_iter= 1000,
  pca_center= TRUE,
  pca_scale= FALSE,
  normalize= TRUE,
  momentum= 0.5,
  final_momentum= 0.8,
  eta= 1000,
  exaggeration_factor= 12
  )


plot_reduced_dim(sce, feature_dim = "individual", reduced_dim = "tSNE_Liger")
sce <- cluster_sce(sce, reduction_method = "tSNE_Liger", k = 50, res = 0.0001)
plot_reduced_dim(sce, feature_dim = "clusters", reduced_dim = "tSNE_Liger", label_clusters = TRUE)
sce <- map_celltypes_sce(sce, ctd_folder = "~/Documents/junk/ctd/")



sce <- reduce_dims_sce(
  sce,
  input_reduced_dim = "Liger",
  reduction_methods = "UMAP",
  pca_dims= 30,
  n_neighbors= 35,
  n_components= 2,
  init= "spectral",
  metric= "euclidean",
  n_epochs= 200,
  learning_rate= 1,
  min_dist= 0.4,
  spread= 0.85,
  set_op_mix_ratio= 1,
  local_connectivity= 1,
  repulsion_strength= 1,
  negative_sample_rate= 5,
  fast_sgd= FALSE,
  n_threads = 10
)

plot_reduced_dim(sce, feature_dim = "individual", reduced_dim = "UMAP_Liger")
sce <- cluster_sce(sce, reduction_method = "UMAP_Liger", k = 50, res = 0.001)
plot_reduced_dim(sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = TRUE)
sce <- map_celltypes_sce(sce, ctd_folder = "~/Documents/junk/ctd/")
plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger", label_clusters = TRUE, size = .5)

plot_reduced_dim(sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = TRUE, size = 0.5)
plot_reduced_dim(sce, feature_dim = "individual", reduced_dim = "UMAP_Liger", label_clusters = TRUE, size = 0.5)

plot_reduced_dim_gene(sce, reduced_dim = "tSNE_Liger", gene = "PLP1", size = 0.5)

plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "PLP1", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "SYT1", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "RBFOX3", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "CUX2", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "RORB", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "DRD1", size = 0.5)

plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "ITGA1", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "THY1", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "CD74", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "CX3CR1", size = 0.5)
#CD137, CD82, and CD355
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "SKAP1", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "CCL5", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "PRF1", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "IL32", size = 0.5)
plot_reduced_dim_gene(sce, reduced_dim = "UMAP_Liger", gene = "IL7R", size = 0.5)

