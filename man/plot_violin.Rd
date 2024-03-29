% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_violin.R
\name{plot_violin}
\alias{plot_violin}
\title{Plot Gene Expression Violin Plots Stratified by Group}
\usage{
plot_violin(
  sce,
  group_var = "group",
  subset_var = "cluster_celltype",
  subset_group = "Oligo",
  gene = "PLP1",
  var_order = NULL,
  palette = NULL,
  alpha = 0.05,
  size = 0.01,
  label_angle = 0
)
}
\arguments{
\item{sce}{a SingleCellExperiment object}

\item{group_var}{The colData variable for x-axis groups}

\item{subset_var}{The colData variable to subset on}

\item{subset_group}{The specific subset_var group to subset}

\item{gene}{The gene of interest}

\item{var_order}{Optional re-ordering of subset_group factor levels}

\item{palette}{Optional custom palette}

\item{alpha}{point alpha}

\item{size}{point size}

\item{label_angle}{The angle of x-axis labels (e.g. 0, 45)}
}
\value{
p a ggplot object
}
\description{
Plots the cells in reduced dimensionality space with log10 expression values
for a specified gene.
}
\details{
Generates a 2d plot of cells with gene expression
}
\seealso{
Other plotting functions: 
\code{\link{plot_expr_by_numeric_var}()},
\code{\link{plot_reduced_dim}()},
\code{\link{plot_reduced_dim_gene}()}
}
\concept{plotting functions}
