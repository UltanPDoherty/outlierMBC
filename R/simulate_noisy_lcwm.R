#' simulate_noisy_lcwm
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
#' @param seed Seed.
#' @param crit_val Critical value for uniform sample rejection.
#' @param range_multipliers .
#' @param print_interval How frequently the iteration count is printed.
#'
#' @return Matrix with a label column.
#' @export
#'
#' @examples
#' noisy_lcwm_p1 <- simulate_noisy_lcwm(
#'   n = c(1000, 1000),
#'   mu = list(c(-1), c(+1)),
#'   sigma = list(as.matrix(0.1), as.matrix(0.1)),
#'   beta = list(c(-1, 1), c(1, 1)),
#'   error_sd = c(0.5, 0.5),
#'   outlier_num = c(10, 10),
#'   seed = 123,
#'   crit_val = 0.9999
#' )
#' plot(
#'   x = noisy_lcwm_p1$covariates[, 1], y = noisy_lcwm_p1$responses,
#'   col = 1 + noisy_lcwm_p1$labels, pch = 1 + noisy_lcwm_p1$labels
#' )
simulate_noisy_lcwm <- function(
    n,
    mu,
    sigma,
    beta,
    error_sd,
    outlier_num,
    seed = 123,
    crit_val = 0.9999,
    range_multipliers = c(1.5, 1.5),
    print_interval = Inf) {

  var_num <- length(mu[[1]])
  comp_num <- length(n)

  set.seed(seed)
  covariates <- list()
  responses <- list()
  errors <- list()
  range_mats <- rep(list(matrix(nrow = var_num, ncol = 2)), comp_num)
  dim_means <- matrix(nrow = comp_num, ncol = var_num)
  dim_widths <- matrix(nrow = comp_num, ncol = var_num)
  err_width <- rep(NA, comp_num)
  for (g in seq_len(comp_num)) {
    covariates[[g]] <- mvtnorm::rmvnorm(n[g], mu[[g]], sigma[[g]])

    errors[[g]] <- stats::rnorm(n[g], 0, sd = error_sd[g])

    responses[[g]] <- errors[[g]] + beta[[g]][1] + covariates[[g]] %*% beta[[g]][-1]

    range_mats[[g]] <- apply(covariates[[g]], 2, range)
    dim_means[g, ] <- colMeans(range_mats[[g]])
    dim_widths[g, ] <- range_mats[[g]][2, ] - range_mats[[g]][1, ]

    err_width[g] <- diff(range(errors[[g]]))
  }

  set.seed(seed)
  prob_x <- prob_y <- checks <- temp_outliers_fitted <- rep(NA, comp_num)
  outliers_x <- lapply(outlier_num, function(x) matrix(nrow = x, ncol = var_num))
  outliers_err <- outliers_fitted <- outliers_y <- lapply(outlier_num, function(x) rep(NA, x))
  for (g in seq_len(comp_num)) {
    count <- attempts <- 0
    while (count < outlier_num[g]) {
      count <- count + 1

      for (p in seq_len(var_num)) {
        outliers_x[[g]][count, p] <- stats::runif(
          1,
          dim_means[g, p] - (range_multipliers[1] / 2) * dim_widths[g, p],
          dim_means[g, p] + (range_multipliers[1] / 2) * dim_widths[g, p]
        )
      }

      outliers_err[[g]][count] <- stats::runif(
        1,
        min = 0 - (range_multipliers[2] / 2) * err_width[g],
        max = 0 + (range_multipliers[2] / 2) * err_width[g]
      )

      outliers_fitted[[g]][count] <-
        (beta[[g]][1] + outliers_x[[g]][count, ] %*% beta[[g]][-1])

      outliers_y[[g]][count] <-
        (outliers_fitted[[g]][count] + outliers_err[[g]][count])

      for (h in seq_len(comp_num)) {
        prob_x[h] <- stats::pchisq(
          stats::mahalanobis(outliers_x[[g]][count, ], mu[[h]], sigma[[h]]),
          df = var_num, lower.tail = FALSE
        )

        temp_outliers_fitted[h] <-
          (beta[[h]][1] + outliers_x[[g]][count, ] %*% beta[[h]][-1])

        prob_y[h] <- stats::pnorm(
          abs(outliers_y[[g]][count] - temp_outliers_fitted[h]),
          mean = 0, sd = error_sd[h], lower.tail = FALSE
        )

        checks[h] <- (prob_x[h] * prob_y[h]) < 1 - crit_val
      }

      count <- count - any(!checks)

      attempts <- attempts + 1
      if (attempts %% print_interval == 0) {
        cat(paste0(
          "component ", g, ": ",
          "attempts = ", attempts, ", ",
          "count = ", count,
          "\n"
        ))
      }
    }
  }

  labels <- rep(seq_len(comp_num), n)
  labels <- c(labels, rep(0, sum(outlier_num)))

  covariates <- rbind(Reduce(rbind, covariates), Reduce(rbind, outliers_x))
  colnames(covariates) <- paste0("V", seq_len(var_num))

  responses <- c(Reduce(c, responses), Reduce(c, outliers_y))

  return(list(
    covariates = covariates,
    responses = responses,
    labels = labels
  ))
}
