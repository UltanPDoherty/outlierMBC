#' simulate_noisy_mlr
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
#' @param unif_range_multiplier How much greater should the range of the Uniform
#'                              samples be than the range of the Normal samples?
#' @param print_interval How frequently the iteration count is printed.
#'
#' @return Matrix with a label column.
#' @export
#'
#' @examples
#' n_vec <- c(1000)
#' mu_list <- list(+1)
#' sigma_list <- list(as.matrix(0.1))
#' beta_list <- list(c(1, 1))
#' error_sd_vec <- c(0.1)
#' noisy_mlr_p1 <- simulate_noisy_mlr(n_vec, mu_list, sigma_list, beta_list,
#'                                    error_sd_vec,
#'                                    outlier_num = 20, seed = 123,
#'                                    crit_val = 0.9999)
#' plot(x = noisy_mlr_p1$covariates[, 1], y = noisy_mlr_p1$responses,
#'      col = 1 + noisy_mlr_p1$labels, pch = 1 + noisy_mlr_p1$labels)
simulate_noisy_mlr <- function(
  n, mu, sigma, beta, error_sd,
  outlier_num, seed = 123, crit_val = 0.9999, unif_range_multiplier = 1.5,
  print_interval = Inf
) {
  var_num <- length(mu[[1]])
  comp_num <- length(n)

  set.seed(seed)
  comps <- list()
  responses <- list()
  errors <- list()
  for (g in seq_len(comp_num)) {
    comps[[g]] <- mvtnorm::rmvnorm(
      n[g],
      mu[[g]],
      sigma[[g]]
    )

    errors[[g]] <- stats::rnorm(n[[g]], 0, sd = error_sd[g])

    responses[[g]] <- errors[[g]] + beta[[g]][1] + comps[[g]] %*% beta[[g]][-1]
  }
  samp <- Reduce(rbind, comps)
  colnames(samp) <- paste0("V", seq_len(var_num))
  resp <- Reduce(c, responses)

  err_width <- diff(range(Reduce(c, errors)))

  set.seed(123)
  count <- 0
  out_error <- rep(NA, comp_num)
  out_pred <- rep(NA, comp_num)
  checks <- rep(NA, comp_num)
  attempts <- 0
  out_norm <- matrix(nrow = outlier_num, ncol = var_num)
  err_unif <- rep(NA, outlier_num)
  out_unif <- rep(NA, outlier_num)
  while (count < outlier_num) {
    out_g <- sample(seq_len(comp_num), 1)

    out_norm[count + 1, ] <- mvtnorm::rmvnorm(1, mu[[out_g]], sigma[[out_g]])

    err_unif[count + 1] <- stats::runif(
      1,
      0 - (unif_range_multiplier / 2) * err_width,
      0 + (unif_range_multiplier / 2) * err_width
    )

    out_unif[count + 1] <- (beta[[g]][1]
                            + sum(out_norm[count + 1, ] * beta[[g]][-1])
                            + err_unif[count + 1])

    for (g in seq_len(comp_num)) {
      checks[g] <- (stats::pnorm(abs(err_unif[count + 1]), 0, error_sd[g])
                    > crit_val)
    }

    count <- count + all(checks)

    attempts <- attempts + 1
    if (attempts %% print_interval == 0) {
      cat(paste0("attempts = ", attempts, ", count = ", count, "\n"))
    }
  }

  labels <- rep(seq_len(comp_num), n)
  labels <- c(labels, rep(0, outlier_num))

  return(list(covariates = rbind(samp, out_norm),
              responses = c(resp, out_unif),
              labels = labels))
}
