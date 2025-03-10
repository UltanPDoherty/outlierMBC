#' simulate_gmm
#'
#' @description
#' Simulate a Gaussian mixture model with multivariate outliers.
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
#' @return `data.frame` with columns:
#' * $X1: first covariate
#' * ...
#' * $G: label
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
    seed = 123,
    crit_val = 0.9999,
    range_multiplier = 1.5,
    print_interval = Inf) {
  var_num <- length(mu[[1]])
  comp_num <- length(n)

  set.seed(seed)
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

  set.seed(123)
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
