#' simulate_noisy_gmm
#'
#' @param n Vector of component sizes.
#' @param mu List of component mean vectors.
#' @param Sigma List of component covariance matrices.
#' @param outlier_num Desired number of outliers.
#' @param seed Seed.
#' @param crit_val Critical value for uniform sample rejection.
#' @param unif_range_multiplier How much greater should the range of the Uniform
#'                              samples be than the range of the Normal samples?
#'
#' @return Matrix with a label column.
#' @export
#'
#' @examples
#' n_vec <- c(100, 100)
#' mu_list <- list(c(-1, -1), c(1, 1))
#' Sigma_list <- list(diag(c(0.1, 0.1)), diag(c(0.1, 0.1)))
#' noisy_gmm <- simulate_noisy_gmm(n_vec, mu_list, Sigma_list, 123, 20)
#' plot(noisy_gmm[, 1:2], col = 1 + noisy_gmm[, 3])
simulate_noisy_gmm <- function(
    n, mu, Sigma,
    outlier_num, seed = 123, crit_val = 0.999, unif_range_multiplier = 1.5
) {
  var_num <- length(mu[[1]])
  comp_num <- length(n)

  set.seed(seed)
  comps <- list()
  for (g in seq_len(comp_num)) {
    comps[[g]] <- mvtnorm::rmvnorm(
      n[g],
      mu[[g]],
      Sigma[[g]]
    )
  }
  samp <- Reduce(rbind, comps)
  colnames(samp) <- paste0("V", seq_len(var_num))

  range_mat <- matrix(nrow = var_num, ncol = 2)
  for (p in seq_len(var_num)) {
    range_mat[p, ] <- range(samp[, p])
  }
  dim_means <- rowMeans(range_mat)
  dim_widths <- range_mat[, 2] - range_mat[, 1]

  chisq_crit <- stats::qchisq(crit_val, df = var_num)

  set.seed(123)
  count <- 0
  checks <- rep(NA, comp_num)
  attempts <- 0
  unif_samp <- matrix(nrow = outlier_num, ncol = var_num)
  while (count < outlier_num) {
    for (p in seq_len(var_num)) {
      unif_samp[count + 1, p] <- stats::runif(
        1,
        dim_means[p] - (unif_range_multiplier / 2) * dim_widths[p],
        dim_means[p] + (unif_range_multiplier / 2) * dim_widths[p])
    }

    for (g in seq_len(comp_num)) {
      unif_mahala <- stats::mahalanobis(unif_samp[count + 1, ],
                                        mu[[g]], Sigma[[g]])
      checks[g] <- unif_mahala > chisq_crit
    }

    count <- count + all(checks)

    attempts <- attempts + 1
    cat(paste0("attempts = ", attempts, ", count = ", count, "\n"))
  }

  labels <- rep(seq_len(comp_num), n)
  labels <- c(labels, rep(0, outlier_num))

  return(cbind(rbind(samp, unif_samp), labels))
}
