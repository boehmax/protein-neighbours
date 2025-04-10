% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_functions.R
\name{getNeighborProteins}
\alias{getNeighborProteins}
\title{Get neighboring proteins}
\usage{
getNeighborProteins(gff.df, protein.id, bp = 300, n = 15)
}
\arguments{
\item{gff.df}{The GFF data frame.}

\item{protein.id}{The protein ID to search for neighbors.}

\item{bp}{Number of base pairs to consider for neighbors (default is 300).}

\item{n}{Number of neighbors to find (default is 15).}
}
\value{
A data frame with neighboring proteins.
}
\description{
Get neighboring proteins
}
