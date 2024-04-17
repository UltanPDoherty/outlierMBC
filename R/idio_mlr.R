#' idio_mlr
#'
#' @description
#' Iterative Detection & Identification of Outliers while Clustering
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
#' faithful_plus <- rbind(faithful, c(4.5, 60), c(4.5, 55), c(4.5, 50), c(4.5, 45), c(4.5, 40))
#' faithful_idio_mlr <- idio_mlr(faithful_plus[, 1], faithful_plus[, 2],
#'                               max_out = 20)
#' par(mfrow = c(1, 2))
#' plot(0:20, faithful_idio_mlr$distrib_diffs, type = "l")
#' abline(v = faithful_idio_mlr$outlier_num)
#' plot(faithful_plus, col = 1 + faithful_idio_mlr$outlier_bool)
#' abline(faithful_idio_mlr$mod)
#' par(mfrow = c(1, 1))
idio_mlr <- function(x, y, max_out, print_interval = Inf) {

  x <- as.matrix(x)
  x0 <- x
  y0 <- y

  obs_num <- nrow(x)
  var_num <- ncol(x)

  distrib_diffs <- c()
  outlier_rank <- rep(0, nrow(x))
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    mod <- stats::lm(y ~ x)

    out <- distrib_diff_mlr(
      stats::rstudent(mod),
      obs_num - var_num
    )

    distrib_diffs[i] <- out$distrib_diff

    outlier_rank[!outlier_rank][out$choice_id] <- i

    x <- x[-out$choice_id, , drop = FALSE]
    y <- y[-out$choice_id]
  }

  outlier_num <- which.min(distrib_diffs) - 1

  outlier_bool <- outlier_rank <= outlier_num & outlier_rank != 0

  mod <- stats::lm(y0[!outlier_bool] ~ x0[!outlier_bool, , drop = FALSE])

  return(list(distrib_diffs = distrib_diffs,
              outlier_bool = outlier_bool,
              outlier_num = outlier_num,
              outlier_rank = outlier_rank,
              mod = mod))
}
