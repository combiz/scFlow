% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_celltype_freqs.R
\name{.tally_cells}
\alias{.tally_cells}
\title{Tally the cells in each celltype for each unique sample}
\usage{
.tally_cells(
  sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  ...
)
}
\arguments{
\item{unique_id_var}{the unique identifier variable for each sample}

\item{celltype_var}{the colData variable specifying celltype or subtype}
}
\value{
mat a matrix sample x celltypes with absolute cell numbers
}
\description{
Tally the cells in each celltype for each unique sample
}
\seealso{
Other helper: 
\code{\link{.append_cell_metric_plots_sce}()},
\code{\link{.append_cell_prop_plots_sce}()},
\code{\link{.clean_ggplot_plot_env}()},
\code{\link{.clean_ggplot_quosures}()},
\code{\link{.clean_quosure_env}()},
\code{\link{.fancy_scientific}()},
\code{\link{.get_d_palette}()},
\code{\link{.grobify_ggplot}()},
\code{\link{.model_fisher_celltype}()},
\code{\link{.plot_dirichlet_results}()},
\code{\link{.plot_unique_id_var}()},
\code{\link{.prepare_dirichlet_plot_table}()},
\code{\link{.prepare_fisher_counts_table}()},
\code{\link{.prepare_unique_id_var_plot_table}()},
\code{\link{.process_dirichlet_fit}()},
\code{\link{.retrieve_covariates}()},
\code{\link{.uniftest}()}
}
\concept{helper}
\keyword{internal}
