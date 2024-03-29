% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_celltype_freqs.R
\name{.retrieve_covariates}
\alias{.retrieve_covariates}
\title{Retrieve the covariates from a SCE colData}
\usage{
.retrieve_covariates(
  sce,
  unique_id_var,
  dependent_var,
  confounding_vars,
  ref_class,
  ...
)
}
\arguments{
\item{sce}{a SingleCellExperiment object}

\item{unique_id_var}{the unique identifier variable for each sample}

\item{dependent_var}{the name of the colData variable for contrasts}

\item{ref_class}{the class of dependent_var used as reference}
}
\value{
mat a dataframe with covariates ordered by unique_id_var
}
\description{
Retrieve the covariates from a SCE colData
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
\code{\link{.tally_cells}()},
\code{\link{.uniftest}()}
}
\concept{helper}
\keyword{internal}
