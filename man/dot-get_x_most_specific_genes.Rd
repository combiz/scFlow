% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_with_ewce.R
\name{.get_x_most_specific_genes}
\alias{.get_x_most_specific_genes}
\title{Find x most specific genes}
\usage{
.get_x_most_specific_genes(
  ct,
  annotLevel = 5,
  howMany = 20,
  ctd,
  exprPercentile = 0.9
)
}
\arguments{
\item{ct}{from colnames of ctdToMap specificity}

\item{howMany}{number of genes to return}

\item{ctd}{the ctd to analyse}

\item{exprPercentile}{must be expressed in at least this proportion}

\item{annot_level}{annotation level}
}
\value{
output the x most specific genes
}
\description{
Find x most specific genes
}
\seealso{
Other clustering and dimensionality reduction: 
\code{\link{.map_celltypes}()},
\code{\link{.map_celltypes_with_ewce}()},
\code{\link{cluster_sce}()},
\code{\link{map_custom_celltypes}()},
\code{\link{reduce_dims_sce}()}
}
\author{
Nathan Skene / Combiz Khozoie
}
\concept{clustering and dimensionality reduction}
\keyword{internal}
