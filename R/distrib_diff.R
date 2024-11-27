#' distrib_diff_gmm
#'
#' @inheritParams ombc_gmm
#' @param z Component assignment probability matrix.
#' @param prop Vector of component proportions.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
#' @param logdet Vector of log-determinants for covariance matrices.
#' @param tail_prop .
#'
#' @return List of
#' * distrib_diff
#' * distrib_diff_vec
#' * choice_id
#' * removal_dens
distrib_diff_gmm <- function(
    x, z, prop, mu, sigma, logdet,
    tail_prop) {
  obs_num <- nrow(x)
  comp_num <- ncol(z)
  track_num <- 2

  bin_z_vec <- apply(z, 1, which.max)
  bin_z_mat <- matrix(0, nrow = obs_num, ncol = comp_num)

  distrib_diff_mat <- matrix(nrow = comp_num, ncol = track_num)
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  mahala_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    bin_z_mat[, g] <- bin_z_vec == g

    dd_g <- distrib_diff_mahalanobis(
      x, z[, g], mu[[g]], sigma[[g]], logdet[g], tail_prop, as.logical(bin_z_mat[, g])
    )
    distrib_diff_mat[g, ] <- dd_g$diff
    dens_mat[, g] <- dd_g$dens
    mahala_mat[, g] <- dd_g$mahalas
  }

  mix_dens <- dens_mat %*% t(prop)

  choice_id <- which.min(mix_dens)
  removal_dens <- mix_dens[choice_id]

  distrib_diff_vec <- c()
  distrib_diff_vec[1] <- sum(prop * distrib_diff_mat[, 1])
  distrib_diff_vec[2] <- sum(distrib_diff_mat[, 2])

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
    tail_prop,
    bin_z_g) {
  var_num <- ncol(x)
  n_g <- sum(z_g)
  stopifnot("A cluster has become too small (< 4 points).\n" = n_g > 3)

  param1 <- var_num / 2
  param2 <- (n_g - var_num - 1) / 2

  mahalas_g <- stats::mahalanobis(x, mu_g, (n_g / (n_g - 1)) * sigma_g)
  scaled_mahalas_g <- ((n_g) / (n_g - 1)^2) * mahalas_g
  mahala_ewcdf_g_func <- spatstat.univar::ewcdf(scaled_mahalas_g, z_g / n_g)

  tail_quant <- stats::qbeta(1 - tail_prop, param1, param2)
  # mahala_tail_prop <- 1 - mahala_ewcdf_g_func(tail_quant)

  eps <- 1e-5
  check_seq <- seq(eps, 1, by = eps)
  mahala_ewcdf_g_vals <- mahala_ewcdf_g_func(check_seq)
  beta_cdf_g_vals <- stats::pbeta(check_seq, param1, param2)
  cdf_diffs <- mahala_ewcdf_g_vals - beta_cdf_g_vals

  distrib_diff_g_x <- c(
    mean(abs(cdf_diffs)),
    # n_g * mahala_tail_prop
    sum(scaled_mahalas_g[bin_z_g] > tail_quant)
  )

  dens_g_x <- exp(
    -0.5 * (var_num * log(2 * pi) + logdet_g + mahalas_g)
  )

  return(list(
    diff = distrib_diff_g_x,
    dens = dens_g_x,
    mahalas = scaled_mahalas_g
  ))
}
