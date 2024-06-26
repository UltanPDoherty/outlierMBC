#' distrib_diff_gmm
#'
#' @param x Data.
#' @param z Component assignment probability matrix.
#' @param prop Vector of component proportions.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
#'
#' @return List of
#' * distrib_diff
#' * distrib_diff_vec
#' * choice_id
distrib_diff_gmm <- function(x, z, prop, mu, sigma) {
  obs_num <- nrow(x)
  comp_num <- ncol(z)

  distrib_diff_vec <- c()
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    out_g <- distrib_diff_mahalanobis(x, z[, g], mu[[g]], sigma[[g]])
    distrib_diff_vec[g] <- out_g$diff
    dens_mat[, g] <- out_g$dens
  }

  mix_dens <- dens_mat %*% t(prop)

  choice_id <- which.min(mix_dens)

  distrib_diff <- sum(prop * distrib_diff_vec)

  return(list(
    distrib_diff = distrib_diff,
    distrib_diff_vec = distrib_diff_vec,
    choice_id = choice_id
  ))
}
