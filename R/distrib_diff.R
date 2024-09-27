#' distrib_diff_gmm
#'
#' @inheritParams ombc_gmm
#' @param z Component assignment probability matrix.
#' @param prop Vector of component proportions.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
#' @param logdet Vector of log-determinants for covariance matrices.
#'
#' @return List of
#' * distrib_diff
#' * distrib_diff_vec
#' * choice_id
#' * removal_dens
distrib_diff_gmm <- function(
    x, z, prop, mu, sigma, logdet,
    tail_probs = c(0, 0.5, 0.75, 0.9, 0.95, 0.99)) {
  obs_num <- nrow(x)
  comp_num <- ncol(z)
  track_num <- length(tail_probs)

  distrib_diff_mat <- matrix(nrow = comp_num, ncol = track_num)
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  mahala_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    dd_g <- distrib_diff_mahalanobis(
      x, z[, g], mu[[g]], sigma[[g]], logdet[g], tail_probs
    )
    distrib_diff_mat[g, ] <- dd_g$diff
    dens_mat[, g] <- dd_g$dens
    mahala_mat[, g] <- dd_g$mahalas
  }

  mix_dens <- dens_mat %*% t(prop)

  choice_id <- which.min(mix_dens)
  removal_dens <- mix_dens[choice_id]

  distrib_diff_vec <- as.numeric(sqrt(prop %*% (distrib_diff_mat^2)))

  return(list(
    distrib_diff_mat = distrib_diff_mat,
    distrib_diff_vec = distrib_diff_vec,
    choice_id = choice_id,
    removal_dens = removal_dens
  ))
}


# ==============================================================================

#' distrib_diff_mahalanobis
#'
#' @inheritParams distrib_diff_gmm
#' @param z_g Assignment probability vector for component g.
#' @param mu_g Mean vector for component g.
#' @param sigma_g Covariance matrix for component g.
#' @param logdet_g Log-determinants of covariance matrix for component g.
#'
#' @return List of
#' * diff
#' * dens
distrib_diff_mahalanobis <- function(
    x,
    z_g,
    mu_g,
    sigma_g,
    logdet_g,
    tail_probs = c(0, 0.5, 0.75, 0.9, 0.95, 0.99)) {
  var_num <- ncol(x)
  n_g <- sum(z_g)
  stopifnot("A cluster has become too small (< 4 points).\n" = n_g > 3)

  param1 <- var_num / 2
  param2 <- (n_g - var_num - 1) / 2

  mahalas_g <- stats::mahalanobis(x, mu_g, (n_g / (n_g - 1)) * sigma_g)
  scaled_mahalas_g <- ((n_g) / (n_g - 1)^2) * mahalas_g
  mahala_ewcdf_g_func <- spatstat.univar::ewcdf(scaled_mahalas_g, z_g / n_g)

  ewcdf_tail <- spatstat.univar::quantile.ewcdf(mahala_ewcdf_g_func, tail_probs)

  check_seq <- seq(0, max(scaled_mahalas_g), length.out = 1e5)
  mahala_ewcdf_g <- mahala_ewcdf_g_func(check_seq)
  beta_cdf_g <- stats::pbeta(check_seq, param1, param2)

  cdf_diffs <- beta_cdf_g - mahala_ewcdf_g

  distrib_diff_g_x <- double(length(tail_probs))
  for (i in seq_along(tail_probs)) {
    tail_subset <- check_seq >= ewcdf_tail[i]

    distrib_diff_g_x[i] <- sqrt(mean(cdf_diffs[tail_subset]^2))
  }

  dens_g_x <- exp(
    -0.5 * (var_num * log(2 * pi) + logdet_g + mahalas_g)
  )

  return(list(
    diff = distrib_diff_g_x,
    dens = dens_g_x,
    mahalas = scaled_mahalas_g
  ))
}
