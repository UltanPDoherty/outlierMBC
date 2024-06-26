#' ombc_gmm
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
#' * labels
#' * final_gmm
#' @export
#'
#' @examples
#' gmm <- simulate_gmm(
#'   n = c(2000, 1000, 1000),
#'   mu = list(c(-1, 0), c(+1, -1), c(+1, +1)),
#'   sigma = list(diag(c(0.2, 4 * 0.2)), diag(c(0.2, 0.2)), diag(c(0.2, 0.2))),
#'   outlier_num = 40,
#'   seed = 123,
#'   crit_val = 0.9999,
#'   range_multiplier = 1.5
#' )
#'
#' ombc_gmm <- ombc_gmm(gmm[, 1:2], comp_num = 3, max_out = 80)
#'
#' par(mfrow = c(1, 2))
#' plot(0:80, ombc_gmm$distrib_diffs, type = "l")
#' abline(v = ombc_gmm$outlier_num)
#' plot(gmm[, "X1", "X2"], col = ombc_gmm$labels + 1, pch = gmm$G + 1)
#' par(mfrow = c(1, 1))
ombc_gmm <- function(x, comp_num, max_out, mnames = "VVV", seed = 123,
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
      warning(paste0(
        "One of the components became too small after removing ",
        i - 1, " outliers."
      ))
      break()
    }

    out <- distrib_diff_gmm(
      x,
      mix$z,
      mix$best_model$model_obj[[1]]$pi_gs,
      mix$best_model$model_obj[[1]]$mu,
      mix$best_model$model_obj[[1]]$sigs
    )

    distrib_diffs[i] <- out$distrib_diff

    outlier_rank[!outlier_rank][out$choice_id] <- i
    x <- x[-out$choice_id, , drop = FALSE]
    z0 <- mix$z[-out$choice_id, , drop = FALSE]
  }

  outlier_num <- which.min(distrib_diffs) - 1

  outlier_bool <- outlier_rank <= outlier_num & outlier_rank != 0

  set.seed(seed)
  mix <- mixture::gpcm(x0[!outlier_bool, ], G = comp_num, mnames = mnames)

  labels <- rep(0, nrow(x0))
  labels[!outlier_bool] <- mix$map

  return(list(
    distrib_diffs = distrib_diffs,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    labels = labels,
    final_gmm = mix
  ))
}
