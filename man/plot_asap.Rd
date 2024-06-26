% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_ASAP.R
\name{plot_asap}
\alias{plot_asap}
\title{plot_asap}
\usage{
plot_asap(asap_result, output_name, type_ggplot = NULL)
}
\arguments{
\item{asap_result}{The matrix returned by asap()}

\item{output_name}{output name for asap pdf plot}

\item{type_ggplot}{NULL by default or 'YES' to plot with ggplot2}
}
\value{
Returns an ASAP_plot.pdf file in the working directory
}
\description{
plot_asap allows to plot ASAP assignations as barplot, using either R base (default) or ggplot2 (type_ggplot = 'YES')
}
