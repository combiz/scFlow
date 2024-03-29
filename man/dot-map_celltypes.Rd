% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_with_ewce.R
\name{.map_celltypes}
\alias{.map_celltypes}
\title{Main mapping function}
\usage{
.map_celltypes(
  ctdToMap,
  ctdToMapAgainst,
  inputSpecies = "human",
  mapAgainstSpecies = "human",
  annotLevel = 1,
  numTopMarkers = 50,
  mappingLevel = 2,
  reps = 1000
)
}
\arguments{
\item{ctdToMap}{ctd for current expt}

\item{ctdToMapAgainst}{the ctd to map against}

\item{inputSpecies}{"human" or "mouse"}

\item{mapAgainstSpecies}{"human" or "mouse"}

\item{annotLevel}{level of annotation}

\item{numTopMarkers}{number of top markers to consider}

\item{mappingLevel}{level of mapping}
}
\value{
sce a SingleCellExperiment object annotated with sample metadata
}
\description{
Main mapping function
}
\seealso{
Other clustering and dimensionality reduction: 
\code{\link{.get_x_most_specific_genes}()},
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
