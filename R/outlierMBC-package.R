#' @keywords internal
#' @description
#'
#' This package has the following exported functions:
#' \describe{
#'   \item{`ombc_gmm`}{Sequentially identify outliers while fitting a Gaussian
#'                     mixture model.}
#'   \item{`ombc_lcwm`}{Sequentially identify outliers while fitting a linear
#'                      cluster-weighted model.}
#'   \item{`simulate_gmm`}{Simulate data from a Gaussian mixture model with
#'                         outliers.}
#'   \item{`simulate_lcwm`}{Simulate data from a linear cluster-weighted model
#'                          with outliers.}
#'   \item{`find_gross`}{Find gross outliers.}
#'   \item{`backtrack`}{Move backwards from the minimum to a more conservative
#'                      solution.}
#'   \item{`backtrack_gmm`}{Fit a Gaussian mixture model to the backtrack
#'                          solution.}
#'   \item{`backtrack_lcwm`}{Fit a linear cluster-weighted model to the
#'                           backtrack solution.}
#'   \item{`plot_curve`}{Plot the dissimilarity curve.}
#'   \item{`plot_backtrack`}{Plot the dissimilarity curve showing the backtrack
#'                           solution.}
#'   \item{`plot_comparison`}{Plot multiple dissimilarity curves.}
#'   \item{`plot_selection`}{Plot dissimilarity values for multiple solutions.}
#' }
#'
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
