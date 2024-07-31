#' @keywords internal
#' @description
#'
#' Exported Functions:
#' * `ombc_gmm` - Identify multivariate outliers while clustering the data with
#'                a Gaussian mixture model.
#' * `ombc_mlr` - Identify response variable outliers while fitting a multiple
#'                linear regression model to the data.
#' * `simulate_noisy_gmm` - Simulate data from a Gaussian mixture model with
#'                          multivariate outliers.
#' * `simulate_noisy_mlr` - Simulate data from a multiple linear regression
#'                          model with response variable outliers.
"_PACKAGE"

## usethis namespace: start
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom flexCWM cwm
#' @importFrom mixture gpcm
#' @importFrom mvtnorm dmvnorm
#' @importFrom spatstat.univar ewcdf
#' @importFrom stats dnorm
#' @importFrom stats quantile
## usethis namespace: end
NULL
