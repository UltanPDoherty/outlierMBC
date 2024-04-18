#' idio_gmm
#'
#' @description
#' Iterative Detection & Identification of Outliers for a Gaussian Mixture Model
#'
#'
#' @param x Data.
#' @param comp_num Number of components.
#' @param max_out Maximum number of outliers.
#' @param mnames Model names for mixture::gpcm.
#' @param seed Seed.
#' @param print_interval How frequently the iteration count is printed.
#'
#' @return List of
#' * distrib_diffs
#' * outlier_bool
#' * outlier_num
#' * outlier_rank
#' * gmm_labels
#' @export
#'
#' @examples
#' n_vec <- c(2000, 1000, 1000)
#' mu_list <- list(c(-1, 0), c(+1, -1), c(+1, +1))
#' sigma_list <- list(diag(c(0.2, 4 * 0.2)),
#'                    diag(c(0.2, 0.2)),
#'                    diag(c(0.2, 0.2)))
#' noisy_gmm_p2g3 <- simulate_noisy_gmm(
#'   n_vec, mu_list, sigma_list,
#'   outlier_num = 40, seed = 123, crit_val = 0.9999,
#'   unif_range_multiplier = 1.5
#' )
#' idio_gmm_p2g3 <- idio_gmm(noisy_gmm_p2g3[, 1:2], comp_num = 3, max_out = 100,
#'                           print_interval = 10)
#' par(mfrow = c(1, 2))
#' plot(0:100, idio_gmm_p2g3$distrib_diffs, type = "l")
#' abline(v = idio_gmm_p2g3$outlier_num)
#' plot(noisy_gmm_p2g3[, 1:2], col = idio_gmm_p2g3$gmm_labels,
#'      pch = 1 + noisy_gmm_p2g3[, 3])
#' par(mfrow = c(1, 1))
idio_gmm <- function(x, comp_num, max_out, mnames = "VVV", seed = 123,
                     print_interval = Inf) {

  x <- as.matrix(x)
  x0 <- x

  z0 <- 2

  var_num <- ncol(x)

  distrib_diffs <- c()
  outlier_rank <- rep(0, nrow(x))
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    set.seed(seed)
    mix <- mixture::gpcm(x, G = comp_num, mnames = mnames, start = z0)

    if (any(colSums(mix$z) < var_num + 1)) {
      warning(paste0("One of the components became too small after removing ",
                     i - 1, " outliers."))
      break()
    }

    out <- distrib_diff_gmm(
      x,
      mix$z,
      mix$best_model$model_obj[[1]]$mu,
      mix$best_model$model_obj[[1]]$sigs
    )

    distrib_diffs[i] <- out$distrib_diff

    outlier_rank[!outlier_rank][out$choice_id] <- i
    x <- x[-out$choice_id, ]
    z0 <- mix$z[-out$choice_id, ]
  }

  outlier_num <- which.min(distrib_diffs) - 1

  outlier_bool <- outlier_rank <= outlier_num & outlier_rank != 0

  set.seed(seed)
  mix <- mixture::gpcm(x0[!outlier_bool, ], G = comp_num, mnames = mnames)

  gmm_labels <- rep(1, nrow(x0))
  gmm_labels[!outlier_bool] <- 1 + mix$map

  return(list(distrib_diffs = distrib_diffs,
              outlier_bool = outlier_bool,
              outlier_num = outlier_num,
              outlier_rank = outlier_rank,
              gmm_labels = gmm_labels))
}
