bioc_pkgs<-c(
  'batchelor',
  'BiocGenerics',
  'BiocStyle',
  'biomaRt',
  'DelayedArray',
  'DelayedMatrixStats',
  'DropletUtils',
  'edgeR',
  'GenomicRanges',
  'graph',
  'IRanges',
  'limma',
  'MAST',
  'multtest',
  'preprocessCore',
  'rhdf5',
  'S4Vectors',
  'scater',
  'SingleCellExperiment',
  'SummarizedExperiment'
)

requireNamespace("BiocManager")
BiocManager::install(bioc_pkgs, ask=F, type = "source")
