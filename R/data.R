#' @title Simulated data set consisting of 1000 observations from 3 Gaussian
#'        components and 10 outliers.
#'
#' @description
#' This data set was simulated using `simulate_gmm`. There are 500 observations
#' in Component 1, 250 observations in Component 2, and 250 observations in
#' Component 3
#'
#' @format ## `gmm_k3n1000o10`
#' A data frame with 1010 rows and 3 columns:
#' \describe{
#'   \item{X1, X2}{Continuous variables.}
#'   \item{G}{Component label: 0 for outliers; 1, 2, or 3 for true points.}
#' }
#'
#' @source For simulation code, see `gmm_k3n1000o10.R` in `data-raw` folder at
#'         <https://github.com/UltanPDoherty/outlierMBC>.
"gmm_k3n1000o10"

#' @title Simulated data set consisting of 2000 observations from 3 Gaussian
#'        components and 20 outliers.
#'
#' @description
#' This data set was simulated using `simulate_gmm`. There are 1000 observations
#' in Component 1, 500 observations in Component 2, and 500 observations in
#' Component 3.
#'
#' @format ## `gmm_k3n2000o20`
#' A data frame with 2020 rows and 3 columns:
#' \describe{
#'   \item{X1, X2}{Continuous variables.}
#'   \item{G}{Component label: 0 for outliers; 1, 2, or 3 for true points.}
#' }
#'
#' @source For simulation code, see `gmm_k3n2000o20.R` in `data-raw` folder at
#'         <https://github.com/UltanPDoherty/outlierMBC>.
"gmm_k3n2000o20"

#' @title Simulated data set consisting of 4000 observations from 3 Gaussian
#'        components and 40 outliers.
#'
#' @description
#' This data set was simulated using `simulate_gmm`. There are 2000 observations
#' in Component 1, 1000 observations in Component 2, and 1000 observations in
#' Component 3.
#'
#' @format ## `gmm_k3n4000o40`
#' A data frame with 4040 rows and 3 columns:
#' \describe{
#'   \item{X1, X2}{Continuous variables.}
#'   \item{G}{Component label: 0 for outliers; 1, 2, or 3 for true points.}
#' }
#'
#' @source For simulation code, see `gmm_k3n4000o40.R` in `data-raw` folder at
#'         <https://github.com/UltanPDoherty/outlierMBC>.
#'
"gmm_k3n4000o40"

#' @title Simulated data set consisting of 1000 observations from 3 Gaussian
#'        components and 10 outliers.
#'
#' @description
#' This data set was simulated using `simulate_lcwm`. There are 300 observations
#' in Component 1, 300 observations in Component 2, and 400 observations in
#' Component 3
#'
#' @format ## `lcwm_k3n1000o10`
#' A data frame with 1010 rows and 3 columns:
#' \describe{
#'   \item{X1}{Continuous explanatory variable.}
#'   \item{Y}{Continuous response variable.}
#'   \item{G}{Component label: 0 for outliers; 1, 2, or 3 for true points.}
#' }
#'
#' @source For simulation code, see `lcwm_k3n1000o10.R` in `data-raw` folder at
#'         <https://github.com/UltanPDoherty/outlierMBC>.
"lcwm_k3n1000o10"

#' @title Simulated data set consisting of 2000 observations from 3 Gaussian
#'        components and 20 outliers.
#'
#' @description
#' This data set was simulated using `simulate_lcwm`. There are 600 observations
#' in Component 1, 600 observations in Component 2, and 800 observations in
#' Component 3.
#'
#' @format ## `lcwm_k3n2000o20`
#' A data frame with 2020 rows and 3 columns:
#' \describe{
#'   \item{X1}{Continuous explanatory variable.}
#'   \item{Y}{Continuous response variable.}
#'   \item{G}{Component label: 0 for outliers; 1, 2, or 3 for true points.}
#' }
#'
#' @source For simulation code, see `lcwm_k3n2000o20.R` in `data-raw` folder at
#'         <https://github.com/UltanPDoherty/outlierMBC>.
"lcwm_k3n2000o20"

#' @title Simulated data set consisting of 4000 observations from 3 Gaussian
#'        components and 40 outliers.
#'
#' @description
#' This data set was simulated using `simulate_lcwm`. There are 1200
#' observations in Component 1, 1200 observations in Component 2, and 1600
#' observations in Component 3.
#'
#' @format ## `lcwm_k3n4000o40`
#' A data frame with 4040 rows and 3 columns:
#' \describe{
#'   \item{X1}{Continuous explanatory variable.}
#'   \item{Y}{Continuous response variable.}
#'   \item{G}{Component label: 0 for outliers; 1, 2, or 3 for true points.}
#' }
#'
#' @source For simulation code, see `lcwm_k3n4000o40.R` in `data-raw` folder at
#'         <https://github.com/UltanPDoherty/outlierMBC>.
"lcwm_k3n4000o40"
