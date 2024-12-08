% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gmm_k3n1000o10}
\alias{gmm_k3n1000o10}
\title{Data set from \code{simulate_gmm}. 3 components, 1000 true observations, 10
outliers.}
\format{
\subsection{\code{gmm_k3n1000o10}}{

A data frame with 1010 rows and 3 columns:
\describe{
\item{X1, X2}{Continuous variables.}
\item{G}{Component label: 0 for outliers; 1, 2, or 3 for true points.}
}
}
}
\usage{
gmm_k3n1000o10
}
\description{
Simulated data set with 1000 observations from 3 Gaussian components and 10
outliers from a uniform distribution with a rejection step.
}
\keyword{datasets}
