#' @title Simulate data from a linear cluster-weighted model with outliers.
#'
#' @description
#' Simulates data from a linear cluster-weighted model, then simulates outliers
#' from a region around each mixture component, with a rejection step to
#' control how unlikely the outliers are under the model.
#'
#' @details
#' `simulate_lcwm` samples a user-defined number of outliers for each component.
#' However, even though an outlier may be associated with one component, it must
#' be outlying with respect to every component.
#'
#' The covariate values of the simulated outliers for a given component `g` are
#' sampled from a Uniform distribution over a hyper-rectangle which is specific
#' to that component. For each covariate dimension, the hyper-rectangle is
#' centred at the midpoint between the maximum and minimum values for that
#' variable from all of the Gaussian observations from component `g`. Its width
#' in that dimension is the distance between the minimum and maximum values for
#' that variable multiplied by the value of `range_multiplier[1]`.
#'
#' The response values of the simulated outliers for a given component `g` are
#' obtained by sampling random errors from a Uniform distribution over a
#' univariate interval, simulating covariate values as discussed above,
#' computing the mean response value for those covariate values, then adding
#' this simulated error to the response. The error sampling interval is centred
#' at the midpoint between the maximum and minimum errors for that variable from
#' all of the Gaussian observations from component `g`. Its width is the
#' distance between the minimum and maximum errors multiplied by the value of
#' `range_multiplier[2]`.
#'
#' A proposed outlier for component `g` is rejected if the probability of
#' sampling a more extreme point from any of the components is greater than
#' `prob_range[2]` or if the probability of sampling a less extreme point from
#' component `g` is less than `prob_range[1]`. This can be visualised as a pair
#' of inner and outer envelopes around each component. To be accepted, a
#' proposed outlier must lie inside the outer envelope for its component and
#' outside the inner envelopes of all components. Setting `prob_range[1] = 0`
#' will eliminate the outer envelope, while setting `prob_range[2] = 0` will
#' eliminate the inner envelope.
#'
#' By setting `outlier_type` = `"x_only"` and giving arbitrary values to
#' `error_sd` (e.g. a zero vector) and `beta` (e.g. a list of zero vectors),
#' then ignoring the simulated `Y` variable, `simulate_lcwm` can be used to
#' simulate a Gaussian mixture model. Since `simulate_lcwm` simulates
#' component-specific outliers from sampling regions around each component,
#' rather than a single sampling region around all of the components, this will
#' not be equivalent to [simulate_gmm]. `simulate_lcwm` also allows the user to
#' set an upper bound on how unlikely an outlier is, as well as a lower bound,
#' whereas [simulate_gmm] only sets a lower bound.
#'
#' @param n Vector of component sizes.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
#' @param beta List of component regression coefficient vectors.
#' @param error_sd Vector of component regression error standard deivations.
#' @param outlier_num Desired number of outliers.
#' @param outlier_type Character string governing whether the outliers are
#'                     outlying with respect to the explanatory variable only
#'                     (`"x_only"`), the response variable only (`"y_only"`), or
#'                     both (`"x_and_y"`). `"x_and_y"` is the default value.
#' @param seed Seed.
#' @param prob_range Values for uniform sample rejection.
#' @param range_multipliers For every explanatory variable, the sampling region
#' The sampling region for the Uniform distribution
#'                          used to simulate proposed outliers is
#'                          controlled by multiplying the component widths by
#'                          these values.
#' @param more_extreme Whether to return a column in the data frame consisting
#'                     of the probabilities of sampling more extreme true
#'                     observations than the simulated outliers.
#'
#' @returns
#' `simulate_lcwm` returns a `data.frame` with continuous variables
#' `X1`, `X2`, ..., followed by a continuous response variable, `Y`, and a
#' mixture component label vector `G` with outliers denoted by `0`. The
#' optional variable `more_extreme` may be included, if specified by the
#' corresponding argument.
#'
#' @export
#'
#' @examples
#' lcwm_k3n1000o10 <- simulate_lcwm(
#'   n = c(300, 300, 400),
#'   mu = list(c(3), c(6), c(3)),
#'   sigma = list(as.matrix(1), as.matrix(0.1), as.matrix(1)),
#'   beta = list(c(0, 0), c(-75, 15), c(0, 5)),
#'   error_sd = c(1, 1, 1),
#'   outlier_num = c(3, 3, 4),
#'   outlier_type = "x_and_y",
#'   seed = 123,
#'   prob_range = c(1e-8, 1e-6),
#'   range_multipliers = c(1, 2)
#' )
#'
#' plot(
#'   lcwm_k3n1000o10[, c("X1", "Y")],
#'   col = lcwm_k3n1000o10$G + 1,
#'   pch = lcwm_k3n1000o10$G + 1
#' )
simulate_lcwm <- function(
    n,
    mu,
    sigma,
    beta,
    error_sd,
    outlier_num,
    outlier_type = c("x_and_y", "x_only", "y_only"),
    seed = NULL,
    prob_range = c(1e-8, 1e-6),
    range_multipliers = c(3, 3),
    more_extreme = FALSE) {
  outlier_type <- match.arg(outlier_type)

  var_num <- length(mu[[1]])
  comp_num <- length(n)

  if (!is.null(seed)) {
    set.seed(seed)
  }
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

  if (!is.null(seed)) {
    set.seed(seed)
  }
  outliers <- lapply(outlier_num, function(x) matrix(NA, x, var_num + 2))
  more_extreme_prob <- list()
  for (g in seq_len(comp_num)) {
    more_extreme_prob[[g]] <- double(outlier_num[g])
    for (j in seq_len(outlier_num[g])) {
      out <- uniform_outlier_ombc(
        outlier_type, mu, sigma, beta, error_sd, g, uniform_spans, prob_range
      )
      outliers[[g]][j, ] <- out[seq_len(var_num + 2)]

      more_extreme_prob[[g]][j] <- out[var_num + 3]
    }
  }

  lcwm <- as.data.frame(rbind(
    Reduce(rbind, observations),
    Reduce(rbind, outliers)
  ))
  colnames(lcwm) <- c(paste0("X", seq_len(var_num)), "Y", "G")

  if (more_extreme) {
    lcwm$more_extreme <- c(rep(NA, sum(n)), Reduce(c, more_extreme_prob))
  }

  lcwm
}

# ==============================================================================

#' @title Obtain the span of the observations for each component.
#'
#' @description
#' Determine the minimum and maximum values for each covariate / explanatory
#' variable and for the response errors from all Gaussian observations.
#'
#' @inheritParams simulate_lcwm
#' @param covariates_g Covariate values of the sampled observations.
#' @param errors_g Response errors of the sampled observations.
#'
#' @returns
#' `uniform_spans_lcwm` returns a 2-column matrix. The final row contains the
#' minimum and maximum values of the response errors, while the previous rows
#' contain the minimum and maximum values for each covariate.
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

  cbind(c(mins_x, min_err), c(maxs_x, max_err))
}

# ==============================================================================

#' @title Produce a single sample that passes the outlier checks.
#'
#' @description
#' This function calls [uniform_sample_lcwm] to sample a proposed outlier and
#' then calls [test_outlier_ombc] to check if it satisfies the required
#' criteria.
#'
#' @inheritParams simulate_lcwm
#' @param g Component index.
#' @param uniform_spans Covariate and response error spans.
#'
#' @returns
#' `uniform_outlier_ombc` returns a simulated outlier as a vector containing its
#' covariate values, response value, and its component label `0`. This vector's
#' final element is the probability of sampling a more extreme Gaussian point
#' from this outlier's associated component.
uniform_outlier_ombc <- function(
    outlier_type,
    mu, sigma, beta, error_sd, g,
    uniform_spans, prob_range) {
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
      uniform_sample$x, uniform_sample$y, prob_range, g
    )
    prob_g <- test[2]
    test <- test[1]
  }

  # This zero is the outlier label.
  c(uniform_sample$x, uniform_sample$y, 0, prob_g)
}

# ==============================================================================

#' @title Sample a potential outlier.
#'
#' @description
#' If `outlier_type = "x_and_y"`, then both the covariate values and response
#' error of the outlier proposed by this function will be Uniformly distributed.
#' If `outlier_type = "x_only"`, then the covariate values will be Uniformly
#' distributed but the response error will be Normally distributed. If
#' `outlier_type = "y_only"`, then the response error will be Uniformly
#' distributed but the covariate values will be Normally distributed.
#'
#' @inheritParams simulate_lcwm
#' @param mu_g Covariate mean vector for component `g`.
#' @param sigma_g Covariate covariance matrix for component `g`.
#' @param beta_g Regression coefficient vector for component `g`.
#' @param error_sd_g Regression error standard deviation for component `g`.
#' @param uniform_spans_g Covariate and response error ranges for component `g`.
#'
#' @returns
#' `uniform_sample_lcwm` returns a list with the following elements:
#' \describe{
#'   \item{`x`}{Vector of covariate values.}
#'   \item{`y`}{Response value.}
#' }
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

  list(x = outlier_x_g, y = outlier_y_g)
}

# ==============================================================================

#' @title Check if a new sample satisfies the outlier criteria.
#'
#' @description
#' This function checks whether a given sample is an acceptable outlier with
#' respect to `prob_range` and also computes the probability of sampling a more
#' extreme point from component `g`.
#'
#' @inheritParams simulate_lcwm
#' @param x_sample New covariate sample.
#' @param y_sample New response sample.
#' @param g Component number.
#'
#' @returns
#' `test_outlier_ombc` returns a vector consisting of a logical value indicating
#' whether the new sample satisfies the outlier checks, and a numeric value
#' giving the probability of sampling a more extreme point from component `g`.
test_outlier_ombc <- function(
    outlier_type,
    mu, sigma, beta, error_sd,
    x_sample, y_sample, prob_range, g) {
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
      x_and_y = (prob_x[h] * prob_y[h]) < prob_range[2],
      x_only = prob_x[h] < prob_range[2],
      y_only = prob_y[h] < prob_range[2]
    )
  }

  check2 <- switch(outlier_type,
    x_and_y = (prob_x[g] * prob_y[g]) > prob_range[1],
    x_only = prob_x[g] > prob_range[1],
    y_only = prob_y[g] > prob_range[1]
  )

  prob_g <- switch(outlier_type,
    x_and_y = min(prob_x[g] * prob_y[g]),
    x_only = min(prob_x[g]),
    y_only = min(prob_y[g])
  )

  c(all(checks) && check2, prob_g)
}
