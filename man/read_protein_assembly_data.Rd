% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_io.R
\name{read_protein_assembly_data}
\alias{read_protein_assembly_data}
\title{Input/Output Functions for Protein Neighbor Analysis}
\usage{
read_protein_assembly_data(
  protein_file = "proteins.csv",
  assembly_file = "assm_accs.csv",
  protein_assembly_file = "assm_accs_protein.csv",
  PATH = "data",
  interactive = TRUE,
  protein_of_interest = NULL
)
}
\arguments{
\item{protein_file}{The path to the protein file (CSV with protein IDs).}

\item{assembly_file}{The path to the assembly file (CSV with assembly IDs).}

\item{protein_assembly_file}{The path to the protein assembly file (maps proteins to assemblies).}

\item{PATH}{The base path to the data directory.}

\item{interactive}{Whether to prompt the user for input (default TRUE).}

\item{protein_of_interest}{Optional protein ID to use (overrides prompt).}
}
\value{
A list containing the protein, assembly, and protein assembly data frames.
}
\description{
This file contains functions for reading and writing data files,
including protein and assembly data, clades, and representative proteins.
}
\author{
Maximilian Böhm
}
