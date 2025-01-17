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

  bin_z <- apply(z, 1, which.max)

  distrib_diff_mat <- matrix(nrow = comp_num, ncol = track_num)
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  mahala_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    bin_z_g <- bin_z == g
    dd_g <- distrib_diff_mahalanobis(
      x, z[, g], mu[[g]], sigma[[g]], logdet[g], tail_prop, bin_z_g
    )
    distrib_diff_mat[g, ] <- dd_g$diff
    dens_mat[, g] <- dd_g$dens
    mahala_mat[, g] <- dd_g$mahalas
  }

  mix_dens <- dens_mat %*% t(prop)

  choice_id <- which.min(mix_dens)
  removal_dens <- mix_dens[choice_id]

  distrib_diff_vec <- c()
  distrib_diff_vec[1] <- sqrt(sum(prop * distrib_diff_mat[, 1]^2))
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
#' @param bin_z_g .
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
  mahala_tail_prop <- 1 - mahala_ewcdf_g_func(tail_quant)

  eps <- 1e-5
  check_seq <- seq(eps, 1, by = eps)
  mahala_ewcdf_g_vals <- mahala_ewcdf_g_func(check_seq)
  beta_cdf_g_vals <- stats::pbeta(check_seq, param1, param2)
  cdf_diffs <- mahala_ewcdf_g_vals - beta_cdf_g_vals

  distrib_diff_g_x <- c(
    mean(abs(cdf_diffs)),
    n_g * mahala_tail_prop
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

# ==============================================================================


#' distrib_diff_gmm
#'
#' @inheritParams ombc_gmm
#' @param bin_z Binary version of component assignment probability matrix.
#' @param prop Vector of component proportions.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
#' @param logdet Vector of log-determinants for covariance matrices.
#' @param dens .
#'
#' @return List of
#' * distrib_diff
#' * distrib_diff_vec
#' * choice_id
#' * removal_dens
count_extremes <- function(
    x, bin_z, prop, mu, sigma, logdet, dens) {
  obs_num <- nrow(x)
  comp_num <- ncol(bin_z)

  outlier_bool <- rep(FALSE, obs_num)
  outlier_num <- 0

  sort_dens <- sort(dens)

  extreme_count <- Inf
  while (extreme_count > 1) {
    extreme_counts <- c()

    for (g in seq_len(comp_num)) {
      extreme_counts[g] <- count_extremes_g(
        x[!outlier_bool, ], bin_z[!outlier_bool, g],
        mu[[g]], sigma[[g]], logdet[g],
        1 / (obs_num - outlier_num)
      )
    }
    extreme_count <- sum(extreme_counts)

    temp_outlier_bool <- dens < sort_dens[extreme_count]
    temp_outlier_num <- sum(temp_outlier_bool)

    outlier_bool[!outlier_bool][temp_outlier_bool] <- TRUE
    outlier_num <- sum(outlier_bool)

    dens <- dens[!temp_outlier_bool]
    sort_dens <- sort_dens[-seq_len(temp_outlier_num)]
  }

  return(outlier_bool)
}

# ==============================================================================

#' distrib_diff_mahalanobis
#'
#' @inheritParams distrib_diff_gmm
#' @param bin_z_g Binary version of assignment probability vector for component
#'                g.
#' @param mu_g Mean vector for component g.
#' @param sigma_g Covariance matrix for component g.
#' @param logdet_g Log-determinants of covariance matrix for component g.
#'
#' @return List of
#' * diff
#' * dens
count_extremes_g <- function(
    x,
    bin_z_g,
    mu_g,
    sigma_g,
    logdet_g,
    tail_prop) {
  bin_z_g <- as.logical(bin_z_g)

  var_num <- ncol(x)
  n_g <- sum(bin_z_g)
  stopifnot("A cluster has become too small (< 4 points).\n" = n_g > 3)

  param1 <- var_num / 2
  param2 <- (n_g - var_num - 1) / 2

  mahalas_g <- stats::mahalanobis(
    x[bin_z_g, ], mu_g, (n_g / (n_g - 1)) * sigma_g
  )
  scaled_mahalas_g <- ((n_g) / (n_g - 1)^2) * mahalas_g

  tail_quant <- stats::qbeta(1 - tail_prop, param1, param2)
  extreme_count_g <- sum(scaled_mahalas_g > tail_quant)

  return(extreme_count_g)
}
