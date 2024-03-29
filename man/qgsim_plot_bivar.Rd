% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qgsim_plot_bivar.R
\name{qgsim_plot_bivar}
\alias{qgsim_plot_bivar}
\title{Generate a bivariate plot}
\usage{
qgsim_plot_bivar(qgsim_res, populations, vars, type = "p", pch = 16)
}
\arguments{
\item{qgsim_res}{the results of a qgsim function call.}

\item{populations}{either a single number indicating a population to plot, or a vector with two values indicating populations to plot against each other.}

\item{vars}{either a single value or a vector of with two values indicating the z/theta values to plot against each other (can be z1, z2, theta1, or theta2).}

\item{type}{same as for base R plot() function}

\item{pch}{same as for base R plot() function}
}
\value{
the output of an r plot() function displaying the required bivariate plot
}
\description{
This function takes the output of qgsim and some specifiers and generates a bivariate plot
}
