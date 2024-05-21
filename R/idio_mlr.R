#' idio_mlr
#'
#' @description
#' Iterative Detection & Identification of Outliers for Multiple Linear
#' Regression
#'
#'
#' @param x Covariate data.
#' @param y Response variable.
#' @param max_out Maximum number of outliers.
#' @param print_interval How frequently the iteration count is printed.
#'
#' @return List of
#' * distrib_diffs
#' * outlier_bool
#' * outlier_num
#' * outlier_rank
#' * mod
#' @export
#'
#' @examples
#' n_vec <- c(1000)
#' mu_list <- list(+1)
#' sigma_list <- list(as.matrix(0.1))
#' beta_list <- list(c(1, 1))
#' error_sd_vec <- c(0.5)
#' noisy_mlr_p1 <- simulate_noisy_mlr(n_vec, mu_list, sigma_list, beta_list,
#'   error_sd_vec,
#'   outlier_num = 20, seed = 123,
#'   crit_val = 0.9999
#' )
#' idio_mlr_p1 <- idio_mlr(noisy_mlr_p1$covariates, noisy_mlr_p1$responses,
#'   max_out = 40
#' )
#' # par(mfrow = c(1, 2))
#' # plot(0:40, idio_mlr_p1$distrib_diffs, type = "l")
#' # abline(v = idio_mlr_p1$outlier_num)
#' # plot(x = noisy_mlr_p1$covariates[, 1], y = noisy_mlr_p1$responses,
#' #      pch = 1 + noisy_mlr_p1$labels, col = 1 + idio_mlr_p1$outlier_bool)
#' # par(mfrow = c(1, 1))
idio_mlr <- function(x, y, max_out, print_interval = Inf) {
  x <- as.matrix(x)
  x0 <- x
  y0 <- y

  distrib_diffs <- c()
  outlier_rank <- rep(0, nrow(x))
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    mod <- stats::lm(y ~ x)

    out <- distrib_diff_mlr(mod)

    distrib_diffs[i] <- out$distrib_diff

    outlier_rank[!outlier_rank][out$choice_id] <- i

    x <- x[-out$choice_id, , drop = FALSE]
    y <- y[-out$choice_id]
  }

  outlier_num <- which.min(distrib_diffs) - 1

  outlier_bool <- outlier_rank <= outlier_num & outlier_rank != 0

  mod <- stats::lm(y0 ~ x0,
    data = list(x0 = x0, y0 = y0),
    subset = !outlier_bool
  )

  return(list(
    distrib_diffs = distrib_diffs,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    mod = mod
  ))
}
