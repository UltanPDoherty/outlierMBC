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
#' @importFrom changepoint cpt.meanvar
#' @importFrom clue solve_LSAP
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom cpop cpop
#' @importFrom dbscan kNNdist
#' @importFrom dplyr mutate
#' @importFrom flexCWM cwm
#' @importFrom ggpubr ggarrange
#' @importFrom mixture gpcm
#' @importFrom mvtnorm dmvnorm
#' @importFrom spatstat.univar ewcdf
#' @importFrom stats as.dist
#' @importFrom stats dist
#' @importFrom stats dnorm
#' @importFrom stats lm
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom tidyr pivot_longer
## usethis namespace: end
NULL
