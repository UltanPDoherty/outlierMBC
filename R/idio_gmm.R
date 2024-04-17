#' idio_gmm
#'
#' @description
#' Iterative Detection & Identification of Outliers while Clustering
#'
#'
#' @param x Data.
#' @param G Number of components.
#' @param max_out Maximum number of outliers.
#' @param mnames Model names for mixture::gpcm.
#' @param seed Seed.
#' @param print_progress Should iteration count be printed?
#'
#' @return List of
#' * distrib_diffs
#' * outlier_bool
#' * outlier_num
#' * outlier_rank
#' * oGMM_labels
#' @export
#'
#' @examples
#' faithful_idio_gmm <- idio_gmm(faithful, G = 2, max_out = 20, seed = 123)
#' par(mfrow = c(1, 2))
#' plot(0:20, faithful_idio_gmm$distrib_diffs, type = "l")
#' abline(v = faithful_idio_gmm$outlier_num)
#' plot(faithful, col = faithful_idio_gmm$oGMM_labels)
#' par(mfrow = c(1, 1))
idio_gmm <- function(x, G, max_out, mnames = "VVV", seed = 123,
                      print_progress = FALSE) {

  x <- as.matrix(x)

  x0 <- x

  distrib_diffs <- c()
  outlier_rank <- rep(0, nrow(x))
  for (i in seq_len(max_out + 1)) {
    if (print_progress) cat("i = ", i, "\n")

    set.seed(seed)

    mix <- mixture::gpcm(x, G = G, mnames = mnames)

    out <- distrib_diff_gmm(
      x,
      mix$z,
      mix$best_model$model_obj[[1]]$mu,
      mix$best_model$model_obj[[1]]$sigs
    )

    distrib_diffs[i] <- out$distrib_diff

    outlier_rank[!outlier_rank][out$choice_id] <- i
    x <- x[-out$choice_id, ]
  }

  outlier_num <- which.min(distrib_diffs) - 1

  outlier_bool <- outlier_rank <= outlier_num & outlier_rank != 0

  oGMM_labels <- rep(1, nrow(x0))
  set.seed(seed)
  mix <- mixture::gpcm(x0[!outlier_bool, ], G = G, mnames = mnames)
  oGMM_labels[!outlier_bool] <- 1 + mix$map

  return(list(distrib_diffs = distrib_diffs,
              outlier_bool = outlier_bool,
              outlier_num = outlier_num,
              outlier_rank = outlier_rank,
              oGMM_labels = oGMM_labels))
}
