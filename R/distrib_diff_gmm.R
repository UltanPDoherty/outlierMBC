#' distrib_diff_gmm
#'
#' @param x Data.
#' @param z Component assignment probability matrix.
#' @param mu List of component mean vectors.
#' @param Sigma List of component covariance matrices.
#'
#' @return List of
#' * distrib_diff
#' * distrib_diff_vec
#' * choice_id
distrib_diff_gmm <- function(x, z, mu, Sigma) {

  obs_num <- nrow(x)
  comp_num <- ncol(z)
  prop <- colMeans(z)

  distrib_diff_vec <- c()
  scaled_mahalas <- matrix(nrow = obs_num, ncol = comp_num)
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    out_g <- distrib_diff_gmm_g(x, z[, g], mu[[g]], Sigma[[g]])
    distrib_diff_vec[g] <- out_g$distrib_diff_g
    scaled_mahalas[, g] <- out_g$scaled_mahalas_g
    dens_mat[, g] <- out_g$dens_g
  }

  mix_scaled_mahalas <- rowSums(z * scaled_mahalas)
  mix_dens <- rowSums(z * dens_mat)

  choice_id <- which.min(mix_dens)

  distrib_diff <- sum(prop * distrib_diff_vec)

  return(list(distrib_diff = distrib_diff,
              distrib_diff_vec = distrib_diff_vec,
              choice_id = choice_id))
}

#' distrib_diff_gmm_g
#'
#' @param x Data.
#' @param z_g Component assignment probability vector.
#' @param mu_g Component mean vector.
#' @param Sigma_g Component covariance matrix.
#'
#' @return List of
#' * distrib_diff_g
#' * scaled_mahalas_g
#' * dens_g
distrib_diff_gmm_g <- function(x, z_g, mu_g, Sigma_g) {
  var_num <- ncol(x)
  obs_num <- nrow(x)

  n_g <- sum(z_g)
  mahalas_g <- stats::mahalanobis(x, mu_g, (n_g / (n_g - 1)) * Sigma_g)

  w_g <- z_g / n_g

  scaling_g <- (n_g) / (n_g - 1)^2

  scaled_mahalas_g <- scaling_g * mahalas_g

  mahala_ewcdf_g_func <- spatstat.geom::ewcdf(scaled_mahalas_g, w_g)

  # checkpoints <- sort(scaled_mahalas_g)
  checkpoints <- stats::qbeta(seq(0.001, 0.999, 0.001),
                              var_num / 2, (n_g - var_num - 1) / 2)

  mahala_ewcdf_g <- mahala_ewcdf_g_func(checkpoints)
  beta_cdf_g <- stats::pbeta(checkpoints, var_num / 2, (n_g - var_num - 1) / 2)

  distrib_diff_g <- sum(abs(mahala_ewcdf_g - beta_cdf_g)) / obs_num
  # distrib_diff_g <- max(abs(mahala_ewcdf_g - beta_cdf_g))

  dens_g_factor <- (2 * pi)^(- 0.5 * var_num) * det(Sigma_g)^(-0.5)
  dens_g <- dens_g_factor * exp(- 0.5 * mahalas_g)

  return(list(distrib_diff_g = distrib_diff_g,
              scaled_mahalas_g = scaled_mahalas_g,
              dens_g = dens_g))
}
