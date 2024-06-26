#' simulate_lcwm
#'
#' @description
#' Simulate a multiple linear regression model with response variable outliers.
#'
#' @param n Vector of component sizes.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
#' @param beta List of component regression coefficient vectors.
#' @param error_sd Vector of component regression error standard deivations.
#' @param outlier_num Desired number of outliers.
#' @param outlier_type .
#' @param seed Seed.
#' @param crit_val Critical value for uniform sample rejection.
#' @param range_multipliers .
#' @param print_interval How frequently the iteration count is printed.
#'
#' @return `data.frame`:
#' * $X1: first covariate
#' * ...
#' * $Y: response
#' * $G: label
#'
#' @export
#'
#' @examples
#' lcwm <- simulate_lcwm(
#'   n = c(1000, 1000),
#'   mu = list(c(-1), c(+1)),
#'   sigma = list(as.matrix(0.2), as.matrix(0.2)),
#'   beta = list(c(1, 0), c(1, 3)),
#'   error_sd = c(1, 1),
#'   outlier_num = c(25, 25),
#'   outlier_type = "x_and_y",
#'   seed = 123,
#'   crit_val = 0.9999,
#'   range_multipliers = c(1.5, 2)
#' )
#'
#' plot(lcwm[, c("X1", "Y")], col = lcwm$G + 1, pch = lcwm$G + 1)
simulate_lcwm <- function(
    n,
    mu,
    sigma,
    beta,
    error_sd,
    outlier_num,
    outlier_type = c("x_and_y", "x_only", "y_only"),
    seed = 123,
    crit_val = 0.9999,
    range_multipliers = c(1.5, 1.5),
    print_interval = Inf) {
  outlier_type <- match.arg(outlier_type)

  var_num <- length(mu[[1]])
  comp_num <- length(n)

  set.seed(seed)
  observations <- covariates <- errors <- responses <- uniform_spans <- list()
  for (g in seq_len(comp_num)) {
    covariates[[g]] <- mvtnorm::rmvnorm(n[g], mu[[g]], sigma[[g]])

    errors[[g]] <- stats::rnorm(n[g], 0, sd = error_sd[g])

    responses[[g]] <-
      errors[[g]] + beta[[g]][1] + covariates[[g]] %*% beta[[g]][-1]

    observations[[g]] <- cbind(covariates[[g]], responses[[g]], rep(g, n[g]))

    uniform_spans[[g]] <- uniform_spans_lcwm(
      range_multipliers, covariates[[g]], errors[[g]]
    )
  }

  set.seed(seed)
  outliers <- lapply(outlier_num, function(x) matrix(NA, x, var_num + 2))
  for (g in seq_len(comp_num)) {
    for (j in seq_len(outlier_num[g])) {
      outliers[[g]][j, ] <- uniform_outlier_ombc(
        outlier_type, mu, sigma, beta, error_sd, g, uniform_spans, crit_val
      )
    }
  }

  lcwm <- as.data.frame(rbind(
    Reduce(rbind, observations),
    Reduce(rbind, outliers)
  ))
  colnames(lcwm) <- c(paste0("X", seq_len(var_num)), "Y", "G")

  return(lcwm)
}

# ==============================================================================

#' Obtain the span of the observations for each component.
#'
#' @inheritParams simulate_lcwm
#' @param covariates_g Covariate values of the sampled observations.
#' @param errors_g Response errors of the sampled observations.
#'
#' @return `matrix`: minimum and maximum columns, rows for each covariate and
#'                   one for the response errors.
#'
#' @export
uniform_spans_lcwm <- function(range_multipliers, covariates_g, errors_g) {
  ranges_x <- apply(covariates_g, 2, range)
  centres_x <- colMeans(ranges_x)
  widths_x <- apply(ranges_x, 2, diff)
  mins_x <- centres_x - ((range_multipliers[1]) / 2) * widths_x
  maxs_x <- centres_x + ((range_multipliers[1]) / 2) * widths_x

  range_err <- range(errors_g)
  centre_err <- mean(range_err)
  width_err <- diff(range_err)
  min_err <- centre_err - (range_multipliers[2] / 2) * width_err
  max_err <- centre_err + (range_multipliers[2] / 2) * width_err

  spans <- cbind(c(mins_x, min_err), c(maxs_x, max_err))

  return(spans)
}

# ==============================================================================

#' Produce a single sample that passes the outlier checks.
#'
#' @inheritParams simulate_lcwm
#' @param g Component index.
#' @param uniform_spans Covariate and response error spans.
#'
#' @return Vector consisting of covariate values, response value, and label 0.
#'
#' @export
uniform_outlier_ombc <- function(
    outlier_type,
    mu, sigma, beta, error_sd, g,
    uniform_spans, crit_val) {
  test <- FALSE

  while (!test) {
    uniform_sample <- uniform_sample_lcwm(
      outlier_type,
      mu[[g]], sigma[[g]],
      beta[[g]], error_sd[g],
      uniform_spans[[g]]
    )

    test <- test_outlier_ombc(
      outlier_type,
      mu, sigma, beta, error_sd,
      uniform_sample$x, uniform_sample$y, crit_val
    )
  }

  return(c(uniform_sample$x, uniform_sample$y, 0))
}

# ==============================================================================

#' Sample a potential outlier.
#'
#' @inheritParams simulate_lcwm
#' @param mu_g Covariate mean vector for one component.
#' @param sigma_g Covariate covariance matrix for one component.
#' @param beta_g Regression coefficient vector for one component.
#' @param error_sd_g Regression error standard deviation for one component.
#' @param uniform_spans_g Covariate and response error ranges for one component.
#'
#' @return List:
#' * x = Covariate sample
#' * y = Response sample
#'
#' @export
uniform_sample_lcwm <- function(
    outlier_type, mu_g, sigma_g, beta_g, error_sd_g, uniform_spans_g) {
  var_num <- length(beta_g) - 1

  if (outlier_type == "y_only") {
    outlier_x_g <- mvtnorm::rmvnorm(1, mu_g, sigma_g)
  } else {
    outlier_x_g <- rep(NA, var_num)
    for (p in seq_len(var_num)) {
      outlier_x_g[p] <- stats::runif(
        1, uniform_spans_g[p, 1], uniform_spans_g[p, 2]
      )
    }
  }

  if (outlier_type == "x_only") {
    outlier_err <- stats::rnorm(1, 0, error_sd_g)
  } else {
    outlier_err <- stats::runif(
      1, uniform_spans_g[var_num + 1, 1], uniform_spans_g[var_num + 1, 2]
    )
  }

  outlier_fitted <- beta_g[1] + outlier_x_g %*% beta_g[-1]

  outlier_y_g <- outlier_fitted + outlier_err

  return(list(x = outlier_x_g, y = outlier_y_g))
}

# ==============================================================================

#' Check whether a new sample is an outlier for each component.
#'
#' @inheritParams simulate_lcwm
#' @param x_sample New covariate sample.
#' @param y_sample New response sample.
#'
#' @return Logical: TRUE if the new sample is an outlier for each component.
#'
#' @export
test_outlier_ombc <- function(
    outlier_type,
    mu, sigma, beta, error_sd,
    x_sample, y_sample, crit_val) {
  comp_num <- length(mu)
  var_num <- length(mu[[1]])

  prob_x <- prob_y <- temp_outliers_fitted <- checks <- rep(NA, comp_num)

  for (h in seq_len(comp_num)) {
    prob_x[h] <- stats::pchisq(
      stats::mahalanobis(x_sample, mu[[h]], sigma[[h]]),
      df = var_num, lower.tail = FALSE
    )

    temp_outliers_fitted[h] <-
      (beta[[h]][1] + x_sample %*% beta[[h]][-1])

    prob_y[h] <- stats::pnorm(
      abs(y_sample - temp_outliers_fitted[h]),
      mean = 0, sd = error_sd[h], lower.tail = FALSE
    )

    checks[h] <- switch(outlier_type,
      x_and_y = (prob_x[h] * prob_y[h]) < 1 - crit_val,
      x_only = prob_x[h] < 1 - crit_val,
      y_only = prob_y[h] < 1 - crit_val
    )
  }

  return(all(checks))
}
