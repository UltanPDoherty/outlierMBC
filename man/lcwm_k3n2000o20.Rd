% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{lcwm_k3n2000o20}
\alias{lcwm_k3n2000o20}
\title{Simulated data set consisting of 2000 observations from 3 Gaussian
components and 20 outliers.}
\format{
\subsection{\code{lcwm_k3n2000o20}}{

A data frame with 2020 rows and 3 columns:
\describe{
\item{X1}{Continuous explanatory variable.}
\item{Y}{Continuous response variable.}
\item{G}{Component label: 0 for outliers; 1, 2, or 3 for true points.}
}
}
}
\source{
For simulation code, see \code{lcwm_k3n2000o20.R} in \code{data-raw} folder at
\url{https://github.com/UltanPDoherty/outlierMBC}.
}
\usage{
lcwm_k3n2000o20
}
\description{
This data set was simulated using \code{simulate_lcwm}. There are 600 observations
in Component 1, 600 observations in Component 2, and 800 observations in
Component 3.
}
\keyword{datasets}
