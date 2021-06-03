bioc_pkgs<-c(
'biomaRt',
'SingleCellExperiment',
'SummarizedExperiment',
'MAST',
'limma',
'scater',
'ROntoTools',
'BiocGenerics',
'S4Vectors',
'IRanges',
'GenomicRanges',
'DelayedArray',
'graphite',
'batchelor',
'DelayedMatrixStats',
#'monocle3',#won't install correctly, known issue for docker 
'DropletUtils',
'multtest',
'graph',
'rhdf5',
'preprocessCore',
'BiocStyle'
)

requireNamespace("BiocManager")
BiocManager::install(bioc_pkgs,ask=F)
