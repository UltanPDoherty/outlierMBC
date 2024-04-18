#' distrib_diff_gmm
#'
#' @param x Data.
#' @param z Component assignment probability matrix.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
#'
#' @return List of
#' * distrib_diff
#' * distrib_diff_vec
#' * choice_id
distrib_diff_gmm <- function(x, z, mu, sigma) {

  obs_num <- nrow(x)
  comp_num <- ncol(z)
  prop <- colMeans(z)

  distrib_diff_vec <- c()
  scaled_mahalas <- matrix(nrow = obs_num, ncol = comp_num)
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    out_g <- distrib_diff_gmm_g(x, z[, g], mu[[g]], sigma[[g]])
    distrib_diff_vec[g] <- out_g$distrib_diff_g
    scaled_mahalas[, g] <- out_g$scaled_mahalas_g
    dens_mat[, g] <- out_g$dens_g
  }

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
#' @param sigma_g Component covariance matrix.
#'
#' @return List of
#' * distrib_diff_g
#' * scaled_mahalas_g
#' * dens_g
distrib_diff_gmm_g <- function(x, z_g, mu_g, sigma_g) {
  var_num <- ncol(x)
  n_g <- sum(z_g)

  check_seq <- seq(0.001, 0.999, 0.001)
  checkpoints <- stats::qbeta(check_seq, var_num / 2, (n_g - var_num - 1) / 2)

  mahalas_g <- stats::mahalanobis(x, mu_g, (n_g / (n_g - 1)) * sigma_g)
  scaled_mahalas_g <- ((n_g) / (n_g - 1)^2) * mahalas_g
  mahala_ewcdf_g_func <- spatstat.geom::ewcdf(scaled_mahalas_g, z_g / n_g)

  mahala_ewcdf_g <- mahala_ewcdf_g_func(checkpoints)
  beta_cdf_g <- stats::pbeta(checkpoints, var_num / 2, (n_g - var_num - 1) / 2)
  distrib_diff_g <- mean(abs(mahala_ewcdf_g - beta_cdf_g))

  dens_g <- (2 * pi)^(-var_num / 2) * det(sigma_g)^(-0.5) * exp(-mahalas_g / 2)

  return(list(distrib_diff_g = distrib_diff_g,
              scaled_mahalas_g = scaled_mahalas_g,
              dens_g = dens_g))
}
