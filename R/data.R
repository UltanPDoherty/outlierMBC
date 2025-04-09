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
#' @source <https://github.com/UltanPDoherty/outlierMBC/blob/main/data-raw/gmm_k3n1000o10.R>
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
#' @source <https://github.com/UltanPDoherty/outlierMBC/blob/main/data-raw/gmm_k3n2000o20.R>
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
#' @source <https://github.com/UltanPDoherty/outlierMBC/blob/main/data-raw/gmm_k3n4000o40.R>
"gmm_k3n4000o40"

#' @title Simulated data set consisting of 1000 observations from 3 Gaussian
#'        components and 10 outliers.
#'
#' @description
#' This data set was simulated using `simulate_lcwm`. There are 500 observations
#' in Component 1, 250 observations in Component 2, and 250 observations in
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
#' @source <https://github.com/UltanPDoherty/outlierMBC/blob/main/data-raw/lcwm_k3n1000o10.R>
"lcwm_k3n1000o10"

#' @title Simulated data set consisting of 2000 observations from 3 Gaussian
#'        components and 20 outliers.
#'
#' @description
#' This data set was simulated using `simulate_lcwm`. There are 1000 observations
#' in Component 1, 500 observations in Component 2, and 500 observations in
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
#' @source <https://github.com/UltanPDoherty/outlierMBC/blob/main/data-raw/lcwm_k3n2000o20.R>
"lcwm_k3n2000o20"

#' @title Simulated data set consisting of 4000 observations from 3 Gaussian
#'        components and 40 outliers.
#'
#' @description
#' This data set was simulated using `simulate_lcwm`. There are 2000 observations
#' in Component 1, 1000 observations in Component 2, and 1000 observations in
#' Component 3.
#'
#' @format ## `lcwm_k3n4000o40`
#' A data frame with 4040 rows and 3 columns:
#' \describe{
#'   \item{X1}{Continuous explanatory variable.}
#'   \item{Y}{Continuous response variable.}
#'   \item{G}{Component label: 0 for outliers; 1, 2, or 3 for true points.}
#' }
#'
#' @source <https://github.com/UltanPDoherty/outlierMBC/blob/main/data-raw/lcwm_k3n4000o40.R>
"lcwm_k3n4000o40"
