#' distrib_diff_gmm
#'
#' @inheritParams ombc1_gmm
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
distrib_diff_gmm <- function(x, z, prop, mu, sigma, logdet, p_range = c(1, 2)) {
  obs_num <- nrow(x)
  comp_num <- ncol(z)
  track_num <- 10

  distrib_diff_mat <- matrix(nrow = comp_num, ncol = track_num)
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  mahala_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    dd_g <- distrib_diff_mahalanobis(
      x, z[, g], mu[[g]], sigma[[g]], logdet[g], p_range
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
    p_range = c(1, 2)) {
  var_num <- ncol(x)
  n_g <- sum(z_g)
  stopifnot("A cluster has become too small (< 4 points).\n" = n_g > 3)

  param1 <- var_num / 2
  param2 <- (n_g - var_num - 1) / 2

  mahalas_g <- stats::mahalanobis(x, mu_g, (n_g / (n_g - 1)) * sigma_g)
  scaled_mahalas_g <- ((n_g) / (n_g - 1)^2) * mahalas_g
  mahala_ewcdf_g_func <- spatstat.univar::ewcdf(scaled_mahalas_g, z_g / n_g)

  eps <- 1e-4
  check_seq <- seq(eps, 1, eps)

  mahala_ewcdf_g <- mahala_ewcdf_g_func(check_seq)
  beta_cdf_g <- stats::pbeta(check_seq, param1, param2)

  cdf_diffs <- beta_cdf_g - mahala_ewcdf_g
  abs_cdf_diffs <- abs(cdf_diffs)
  pos_cdf_diffs <- pmax(rep(0, length(check_seq)), cdf_diffs)
  neg_cdf_diffs <- pmax(rep(0, length(check_seq)), -cdf_diffs)

  p_vals <- round(seq(p_range[1], p_range[2], length.out = 10), 2)
  pos_cdf_pmeans <- vapply(
    p_vals,
    function(y) (mean(pos_cdf_diffs^(y)))^(1 / y),
    double(1L)
  )

  distrib_diff_g_x <- c(
    pos_cdf_pmeans
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
