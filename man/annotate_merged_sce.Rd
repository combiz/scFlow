% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_merged_sce.R
\name{annotate_merged_sce}
\alias{annotate_merged_sce}
\title{Annotate a post-merge SingleCellExperiment with plots}
\usage{
annotate_merged_sce(
  sce,
  plot_vars = c("total_features_by_counts", "total_counts", "pc_mito", "pc_ribo"),
  unique_id_var = "manifest",
  facet_vars = NULL,
  outlier_vars = c("total_features_by_counts", "total_counts")
)
}
\arguments{
\item{sce}{a SingleCellExperiment object}

\item{plot_vars}{the colData variable(s) to generate plots for}

\item{unique_id_var}{the colData variable identifying unique samples}

\item{facet_vars}{the colData variable(s) to facet/subset by}

\item{outlier_vars}{the colData variable(s) to apply adaptive thresholding}
}
\value{
sce a annotated SingleCellExperiment object
}
\description{
Annotate a post-merge SingleCellExperiment with plots
}
\seealso{
Other annotation functions: 
\code{\link{.preprocess_seurat_object}()},
\code{\link{annotate_celltype_metrics}()},
\code{\link{annotate_integrated_sce}()},
\code{\link{annotate_sce}()},
\code{\link{annotate_sce_cells}()},
\code{\link{annotate_sce_genes}()},
\code{\link{filter_sce}()},
\code{\link{find_cells}()},
\code{\link{find_singlets}()},
\code{\link{generate_sce}()},
\code{\link{map_ensembl_gene_id}()},
\code{\link{merge_sce}()},
\code{\link{read_metadata}()},
\code{\link{report_celltype_metrics}()},
\code{\link{report_celltype_model}()},
\code{\link{report_merged_sce}()},
\code{\link{report_qc_sce}()},
\code{\link{run_doubletfinder}()},
\code{\link{sce_to_seu}()}
}
\concept{annotation functions}
