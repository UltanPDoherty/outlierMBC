#' @title Simulate data from a Gaussian mixture model with outliers.
#'
#' @description
#' Simulates data from a Gaussian mixture model, then simulates outliers from a
#' hyper-rectangle, with a rejection step to ensure that the outliers are
#' sufficiently unlikely under the model.
#'
#' @details
#' The simulated outliers are sampled from a Uniform distribution over a
#' hyper-rectangle. For each dimension, the hyper-rectangle is centred at the
#' midpoint between the maximum and minimum values for that variable from all of
#' the Gaussian observations. Its width in that dimension is the distance
#' between the minimum and maximum values for that variable multiplied by the
#' value of `range_multiplier`. If `range_multiplier = 1`, then this
#' hyper-rectangle is the axis-aligned minimum bounding box for all of the
#' Gaussian data points in this data set.
#'
#' The `crit_val` ensures that it would have been sufficiently unlikely for a
#' simulated outlier to have been sampled from any of the Gaussian components.
#' The Mahalanobis distances of a proposed outlier from each component's mean
#' vector with respect to that component's covariance matrix are computed. If
#' any of these Mahalanobis distances are smaller than the critical value of the
#' appropriate Chi-squared distribution, then the proposed outlier is rejected.
#' In summary, for a Uniform sample to be accepted, it must be sufficiently far
#' from each component in terms of Mahalanobis distance.
#'
#' @param n Vector of component sizes.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
#' @param outlier_num Desired number of outliers.
#' @param seed Seed.
#' @param crit_val Critical value for uniform sample rejection.
#' @param range_multiplier How much greater should the range of the Uniform
#'                         samples be than the range of the Normal samples?
#' @param print_interval How frequently the iteration count is printed.
#'
#' @returns
#' `simulate_gmm` return a `data.frame` with continuous variables
#' `X1`, `X2`, ..., followed by a mixture component label vector `G` with
#' outliers denoted by `0`.
#'
#' @export
#'
#' @examples
#' gmm_k3n1000o10 <- simulate_gmm(
#'   n = c(500, 250, 250),
#'   mu = list(c(-1, 0), c(+1, -1), c(+1, +1)),
#'   sigma = list(diag(c(0.2, 4 * 0.2)), diag(c(0.2, 0.2)), diag(c(0.2, 0.2))),
#'   outlier_num = 10,
#'   seed = 123,
#'   crit_val = 0.9999,
#'   range_multiplier = 1.5
#' )
#'
#' plot(
#'   gmm_k3n1000o10[, c("X1", "X2")],
#'   col = gmm_k3n1000o10$G + 1, pch = gmm_k3n1000o10$G + 1
#' )
simulate_gmm <- function(
    n,
    mu,
    sigma,
    outlier_num,
    seed = NULL,
    crit_val = 0.9999,
    range_multiplier = 1.5,
    print_interval = Inf) {
  var_num <- length(mu[[1]])
  comp_num <- length(n)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  comps <- list()
  for (g in seq_len(comp_num)) {
    comps[[g]] <- mvtnorm::rmvnorm(
      n[g],
      mu[[g]],
      sigma[[g]]
    )
  }
  samp <- Reduce(rbind, comps)
  colnames(samp) <- paste0("X", seq_len(var_num))

  range_mat <- matrix(nrow = var_num, ncol = 2)
  for (p in seq_len(var_num)) {
    range_mat[p, ] <- range(samp[, p])
  }
  dim_means <- rowMeans(range_mat)
  dim_widths <- range_mat[, 2] - range_mat[, 1]

  if (!is.null(seed)) {
    set.seed(seed)
  }
  count <- 0
  checks <- rep(NA, comp_num)
  attempts <- 0
  unif_samp <- matrix(nrow = outlier_num, ncol = var_num)
  while (count < outlier_num) {
    for (p in seq_len(var_num)) {
      unif_samp[count + 1, p] <- stats::runif(
        1,
        dim_means[p] - (range_multiplier / 2) * dim_widths[p],
        dim_means[p] + (range_multiplier / 2) * dim_widths[p]
      )
    }

    for (g in seq_len(comp_num)) {
      unif_mahala <- stats::mahalanobis(
        unif_samp[count + 1, ],
        mu[[g]], sigma[[g]]
      )
      checks[g] <- stats::pchisq(unif_mahala, df = var_num) > crit_val
    }

    count <- count + all(checks)

    attempts <- attempts + 1
    if (attempts %% print_interval == 0) {
      cat(paste0("attempts = ", attempts, ", count = ", count, "\n"))
    }
  }

  labels <- rep(seq_len(comp_num), n)
  labels <- c(labels, rep(0, outlier_num))

  data.frame(
    rbind(samp, unif_samp),
    G = labels
  )
}
