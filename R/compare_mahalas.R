#' compare_mahalas
#'
#' @param x Data.
#' @param z Component assignment probability matrix.
#' @param mu List of component mean vectors.
#' @param Sigma List of component covariance matrices.
#'
#' @return List of
#' * mahala_comparison
#' * mahala_comparison_vec
#' * choice_id
compare_mahalas <- function(x, z, mu, Sigma) {

  obs_num <- nrow(x)
  comp_num <- ncol(z)
  prop <- colMeans(z)

  mahala_comparison_vec <- c()
  scaled_mahalas <- matrix(nrow = obs_num, ncol = comp_num)
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    out_g <- compare_mahalas_g(x, z[, g], mu[[g]], Sigma[[g]])
    mahala_comparison_vec[g] <- out_g[[1]]
    scaled_mahalas[, g] <- out_g[[2]]
    dens_mat[, g] <- out_g[[3]]
  }

  mix_scaled_mahalas <- rowSums(z * scaled_mahalas)
  mix_dens <- rowSums(z * dens_mat)

  choice_id <- which.min(mix_dens)

  mahala_comparison <- sum(prop * mahala_comparison_vec)

  return(list(comparison = mahala_comparison,
              mahala_comparison_vec = mahala_comparison_vec,
              choice_id = choice_id))
}

#' compare_mahalas_g
#'
#' @param x Data.
#' @param z_g Component assignment probability vector.
#' @param mu_g Component mean vector.
#' @param Sigma_g Component covariance matrix.
#'
#' @return List of
#' * mahala_comparison_g
#' * scaled_mahalas_g
#' * dens_g
compare_mahalas_g <- function(x, z_g, mu_g, Sigma_g) {
  var_num <- ncol(x)
  obs_num <- nrow(x)

  n_g <- sum(z_g)
  mahalas_g <- stats::mahalanobis(x, mu_g, (n_g / (n_g - 1)) * Sigma_g)

  w_g <- z_g / n_g

  scaling_g <- (n_g) / (n_g - 1)^2

  scaled_mahalas_g <- scaling_g * mahalas_g
  sorted_scaled_mahalas_g <- sort(scaled_mahalas_g)

  mahala_ewcdf_g_func <- spatstat.geom::ewcdf(scaled_mahalas_g, w_g)
  mahala_ewcdf_g <- mahala_ewcdf_g_func(sorted_scaled_mahalas_g)

  beta_cdf_g <- stats::pbeta(sorted_scaled_mahalas_g,
                             var_num / 2, (n_g - var_num - 1) / 2)

  mahala_comparison_g <- sum(abs(mahala_ewcdf_g - beta_cdf_g)) / obs_num
  # mahala_comparison_g <- max(abs(mahala_ewcdf_g - beta_cdf_g))

  dens_g_factor <- (2 * pi)^(- 0.5 * var_num) * det(Sigma_g)^(-0.5)
  dens_g <- dens_g_factor * exp(- 0.5 * mahalas_g)

  return(list(mahala_comparison_g = mahala_comparison_g,
              scaled_mahalas_g = scaled_mahalas_g,
              dens_g = dens_g))
}
