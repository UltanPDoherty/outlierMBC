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
#' @param unif_range_multipliers .
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
    unif_range_multipliers = c(1.5, 1.5),
    print_interval = Inf) {
  var_num <- length(mu[[1]])
  comp_num <- length(n)

  set.seed(seed)
  comps <- list()
  responses <- list()
  errors <- list()
  range_mats <- rep(list(matrix(nrow = var_num, ncol = 2)), comp_num)
  dim_means <- matrix(nrow = comp_num, ncol = var_num)
  dim_widths <- matrix(nrow = comp_num, ncol = var_num)
  for (g in seq_len(comp_num)) {
    comps[[g]] <- mvtnorm::rmvnorm(
      n[g],
      mu[[g]],
      sigma[[g]]
    )

    errors[[g]] <- stats::rnorm(n[g], 0, sd = error_sd[g])

    responses[[g]] <- errors[[g]] + beta[[g]][1] + comps[[g]] %*% beta[[g]][-1]

    range_mats[[g]] <- apply(comps[[g]], 2, range)
    dim_means[g, ] <- colMeans(range_mats[[g]])
    dim_widths[g, ] <- range_mats[[g]][2, ] - range_mats[[g]][1, ]
  }

  samp <- Reduce(rbind, comps)
  colnames(samp) <- paste0("V", seq_len(var_num))
  resp <- Reduce(c, responses)

  err_width <- vapply(errors, function(x) diff(range(x)), double(1L))

  set.seed(seed)
  count <- rep(0, comp_num)
  mahala_probs <- rep(NA, comp_num)
  resp_probs <- rep(NA, comp_num)
  checks <- rep(NA, comp_num)
  attempts <- rep(0, comp_num)
  unif_samp <- lapply(outlier_num, function(x) matrix(nrow = x, ncol = var_num))
  err_unif <- lapply(outlier_num, function(x) rep(NA, x))
  pred_unif <- lapply(outlier_num, function(x) rep(NA, x))
  temp_pred_unif <- rep(NA, comp_num)
  out_unif <- lapply(outlier_num, function(x) rep(NA, x))
  for (g in seq_len(comp_num)) {
    while (count[g] < outlier_num[g]) {
      for (p in seq_len(var_num)) {
        unif_samp[[g]][count[g] + 1, p] <- stats::runif(
          1,
          dim_means[g, p] - (unif_range_multipliers[1] / 2) * dim_widths[g, p],
          dim_means[g, p] + (unif_range_multipliers[1] / 2) * dim_widths[g, p]
        )
      }

      err_unif[[g]][count[g] + 1] <- stats::runif(
        1,
        0 - (unif_range_multipliers[2] / 2) * err_width[g],
        0 + (unif_range_multipliers[2] / 2) * err_width[g]
      )

      pred_unif[[g]][count[g] + 1] <-
        (beta[[g]][1] + sum(unif_samp[[g]][count[g] + 1, ] * beta[[g]][-1]))

      out_unif[[g]][count[g] + 1] <-
        (pred_unif[[g]][count[g] + 1] + err_unif[[g]][count[g] + 1])

      for (h in seq_len(comp_num)) {
        unif_mahala <- stats::mahalanobis(
          unif_samp[[g]][count[g] + 1, ],
          mu[[h]], sigma[[h]]
        )
        mahala_probs[h] <- stats::pchisq(
          unif_mahala,
          df = var_num, lower.tail = FALSE
        )

        temp_pred_unif[h] <-
          (beta[[h]][1] + sum(unif_samp[[g]][count[g] + 1, ] * beta[[h]][-1]))

        resp_probs[h] <- stats::pnorm(
          abs(out_unif[[g]][count[g] + 1] - temp_pred_unif[h]),
          mean = 0, sd = error_sd[h], lower.tail = FALSE
        )

        checks[h] <- (mahala_probs[h] * resp_probs[h]) < 1 - crit_val
      }

      count[g] <- count[g] + all(checks)

      attempts[g] <- attempts[g] + 1
      if (attempts[g] %% print_interval == 0) {
        cat(paste0(
          "attempts[", g, "] = ", attempts[g], ", count[", g, "] = ", count[g],
          "\n"
        ))
      }
    }
  }

  labels <- rep(seq_len(comp_num), n)
  labels <- c(labels, rep(0, sum(outlier_num)))

  return(list(
    covariates = rbind(samp, Reduce(rbind, unif_samp)),
    responses = c(resp, Reduce(c, out_unif)),
    labels = labels
  ))
}
