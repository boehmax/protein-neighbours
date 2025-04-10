% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_open.R
\name{read_representatives}
\alias{read_representatives}
\title{Read Representatives}
\usage{
read_representatives(
  PATH = "data",
  path = "representatives",
  ipg_file = "ipg_representative.txt",
  pdb_file = "pdb_representative.txt",
  cluster_file = "cluster_representative.txt"
)
}
\arguments{
\item{path}{The base path to the data directory (default is 'data/representatives').}

\item{ipg_file}{The name of the IPG alias file (default is 'ipg_representative.txt').}

\item{pdb_file}{The name of the PDB alias file (default is 'pdb_representative.txt').}

\item{cluster_file}{The name of the cluster alias file (default is 'cluster_representative.txt').}
}
\value{
A data frame with combined alias data.
}
\description{
This function reads the IPG, PDB, and cluster alias data from text files and combines them.
}
