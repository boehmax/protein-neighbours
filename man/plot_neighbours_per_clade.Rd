% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05_plotting.R
\name{plot_neighbours_per_clade}
\alias{plot_neighbours_per_clade}
\title{Visualization Functions}
\usage{
plot_neighbours_per_clade(
  combined_data,
  exclude_unknown_clade = FALSE,
  exclude_unknown_cog = FALSE,
  output_path = NULL,
  plot_count_codh = FALSE,
  width = 30,
  height = 15
)
}
\arguments{
\item{combined_data}{A combined data frame with neighbors and clade information.}

\item{exclude_unknown_clade}{Logical indicating whether to exclude "unknown" clades. Default is FALSE.}

\item{exclude_unknown_cog}{Logical indicating whether to exclude "unknown" COG_category. Default is FALSE.}

\item{output_path}{A string representing the output path. Default is NULL.}

\item{plot_count_codh}{Logical indicating whether to plot the count of CODH per clade or count of
Neighbours per clade. Default is FALSE.}

\item{width}{Plot width in cm (default is 30).}

\item{height}{Plot height in cm (default is 15).}
}
\value{
Invisible NULL, called for side effects.
}
\description{
This file contains functions for creating visualizations of protein neighbor data,
including bar plots, correlation matrices, and histograms.
}
\author{
Maximilian Böhm
}
