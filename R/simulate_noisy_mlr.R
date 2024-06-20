#' simulate_noisy_mlr
#'
#' @description
#' Simulate a multiple linear regression model with response variable outliers.
#'
#' @param n Number of true observations.
#' @param mu Covariate mean vector.
#' @param sigma Covariate covariance matrix.
#' @param beta Regression coefficient vector.
#' @param error_sd Regression error standard deviations.
#' @param outlier_num Desired number of outliers.
#' @param seed Seed.
#' @param crit_val Critical value for uniform sample rejection.
#' @param unif_range_multiplier How much greater should the range of the Uniform
#'                              samples be than the range of the Normal samples?
#' @param print_interval How frequently the iteration count is printed.
#'
#' @return List:
#' * covariates: matrix of samples from a multivariate Normal distribution.
#' * responses: vector of regression dependent variable values.
#' * labels: vector of 1s (true observations) and 0s (outliers).
#' @export
#'
#' @examples
#' n <- 1000
#' mu <- c(+1)
#' sigma <- as.matrix(0.1)
#' beta <- c(1, 1)
#' error_sd <- 0.5
#' noisy_mlr_p1 <- simulate_noisy_mlr(
#'   n, mu, sigma, beta,
#'   error_sd,
#'   outlier_num = 20, seed = 123,
#'   crit_val = 0.9999
#' )
#' plot(
#'   x = noisy_mlr_p1$covariates[, 1], y = noisy_mlr_p1$responses,
#'   col = 1 + noisy_mlr_p1$labels, pch = 1 + noisy_mlr_p1$labels
#' )
simulate_noisy_mlr <- function(
    n, mu, sigma, beta, error_sd,
    outlier_num, seed = 123, crit_val = 0.9999, unif_range_multiplier = 1.5,
    print_interval = Inf) {
  var_num <- length(mu)

  set.seed(seed)
  samp <- mvtnorm::rmvnorm(n, mu, sigma)
  colnames(samp) <- paste0("V", seq_len(var_num))

  errors <- stats::rnorm(n, 0, sd = error_sd)
  resp <- errors + beta[1] + samp %*% beta[-1]

  err_width <- diff(range(errors))

  set.seed(seed)
  count <- 0
  attempts <- 0
  out_norm <- matrix(nrow = outlier_num, ncol = var_num)
  err_unif <- rep(NA, outlier_num)
  out_unif <- rep(NA, outlier_num)
  while (count < outlier_num) {
    out_norm[count + 1, ] <- mvtnorm::rmvnorm(1, mu, sigma)

    err_unif[count + 1] <- stats::runif(
      1,
      0 - (unif_range_multiplier / 2) * err_width,
      0 + (unif_range_multiplier / 2) * err_width
    )

    out_unif[count + 1] <-
      beta[1] + sum(out_norm[count + 1, ] * beta[-1]) + err_unif[count + 1]

    check <- stats::pnorm(abs(err_unif[count + 1]), 0, error_sd) > crit_val

    count <- count + check

    attempts <- attempts + 1
    if (attempts %% print_interval == 0) {
      cat(paste0("attempts = ", attempts, ", count = ", count, "\n"))
    }
  }

  labels <- c(rep(1, n), rep(0, outlier_num))

  return(list(
    covariates = rbind(samp, out_norm),
    responses = c(resp, out_unif),
    labels = labels
  ))
}
