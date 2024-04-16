#' idioclust
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
#' * comparisons
#' * outlier_bool
#' * outlier_num
#' * outlier_rank
#' * oGMM_labels
#' @export
#'
#' @examples
#' faithful_idio <- idioclust(faithful, G = 2, max_out = 20, seed = 123)
#' par(mfrow = c(1, 2))
#' plot(faithful_idio$comparisons, type = "l")
#' abline(v = faithful_idio$outlier_num)
#' plot(faithful, col = faithful_idio$oGMM_labels)
#' par(mfrow = c(1, 1))
idioclust <- function(x, G, max_out, mnames = "VVV", seed = 123,
                      print_progress = FALSE) {

  x <- as.matrix(x)

  x0 <- x

  comparisons <- c()
  outlier_rank <- rep(0, nrow(x))
  for (i in seq_len(max_out)) {
    if (print_progress) cat("i = ", i, "\n")

    set.seed(seed)

    mix <- mixture::gpcm(x, G = G, mnames = mnames)

    out <- compare_mahalas(x, mix$z,
                           mix$best_model$model_obj[[1]]$mu,
                           mix$best_model$model_obj[[1]]$sigs)

    comparisons[i] <- out$comparison

    outlier_rank[!outlier_rank][out$choice_id] <- i
    x <- x[-out$choice_id, ]
  }

  outlier_num <- which.min(comparisons)

  outlier_bool <- outlier_rank <= outlier_num & outlier_rank != 0

  oGMM_labels <- rep(1, nrow(x0))
  set.seed(seed)
  mix <- mixture::gpcm(x0[!outlier_bool, ], G = G, mnames = mnames)
  oGMM_labels[!outlier_bool] <- 1 + mix$map

  return(list(comparisons = comparisons,
              outlier_bool = outlier_bool,
              outlier_num = outlier_num,
              outlier_rank = outlier_rank,
              oGMM_labels = oGMM_labels))
}
