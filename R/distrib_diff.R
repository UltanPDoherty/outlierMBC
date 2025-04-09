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
#' * distrib_diff_vec: Vector containing dissimilarity value for each component.
#' * distrib_diff: Aggregated dissimilarity across components.
#' * choice_id: Index of observation with lowest mixture density.
#' * removal_dens: Value of the lowest mixture density.
distrib_diff_gmm <- function(
    x, z, prop, mu, sigma, logdet) {
  obs_num <- nrow(x)
  comp_num <- ncol(z)

  distrib_diff_vec <- double(comp_num)
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  mahala_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    dd_g <- distrib_diff_mahalanobis(
      x, z[, g], mu[[g]], sigma[[g]], logdet[g]
    )
    distrib_diff_vec[g] <- dd_g$diff
    dens_mat[, g] <- dd_g$dens
    mahala_mat[, g] <- dd_g$mahalas
  }

  mix_dens <- dens_mat %*% t(prop)

  choice_id <- which.min(mix_dens)
  removal_dens <- mix_dens[choice_id]

  distrib_diff <- sqrt(sum(prop * distrib_diff_vec^2))

  list(
    distrib_diff_vec = distrib_diff_vec,
    distrib_diff = distrib_diff,
    choice_id = choice_id,
    removal_dens = removal_dens
  )
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
#' * diff: Dissimilarity value for this component.
#' * dens: Gaussian density of all observations for this component.
#' * mahalas: Scaled squared sample Mahalanobis distances for all observations
#'            with respect to this component.
distrib_diff_mahalanobis <- function(
    x,
    z_g,
    mu_g,
    sigma_g,
    logdet_g) {
  var_num <- ncol(x)
  n_g <- sum(z_g)
  stopifnot("A cluster has become too small (< 4 points).\n" = n_g > 3)

  param1 <- var_num / 2
  param2 <- (n_g - var_num - 1) / 2

  mahalas_g <- stats::mahalanobis(x, mu_g, (n_g / (n_g - 1)) * sigma_g)
  scaled_mahalas_g <- ((n_g) / (n_g - 1)^2) * mahalas_g
  mahala_ewcdf_g_func <- spatstat.univar::ewcdf(scaled_mahalas_g, z_g / n_g)

  eps <- 1e-5
  check_seq <- seq(eps, 1, by = eps)
  mahala_ewcdf_g_vals <- mahala_ewcdf_g_func(check_seq)
  beta_cdf_g_vals <- stats::pbeta(check_seq, param1, param2)
  cdf_diffs <- mahala_ewcdf_g_vals - beta_cdf_g_vals

  distrib_diff_g_x <- mean(abs(cdf_diffs))

  dens_g_x <- exp(
    -0.5 * (var_num * log(2 * pi) + logdet_g + mahalas_g)
  )

  list(
    diff = distrib_diff_g_x,
    dens = dens_g_x,
    mahalas = scaled_mahalas_g
  )
}

# ==============================================================================

#' distrib_diff_lcwm
#'
#' @inheritParams ombc_lcwm
#' @inheritParams distrib_diff_gmm
#' @param mu Matrix of component mean vectors.
#' @param sigma Array of component covariance matrices.
#' @param mod_list List of component regression models.
#' @param y_sigma Vector of component regression standard deviations.
#'
#' @return List of
#' * distrib_diff
#' * distrib_diff_vec
#' * choice_id
distrib_diff_lcwm <- function(
    x,
    z,
    prop,
    mu,
    sigma,
    mod_list,
    y_sigma,
    dd_weight = 0.5,
    dens_power = 0.5) {
  obs_num <- nrow(x)
  comp_num <- ncol(z)

  dens_arr <- array(dim = c(obs_num, comp_num, 2))
  distrib_diff_mat <- matrix(nrow = comp_num, ncol = 2)
  distrib_diff_vec <- c()
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    dd_g <- distrib_diff_lcwm_g(
      x, z[, g],
      mu[, g, drop = FALSE],
      as.matrix(sigma[, , g]),
      mod_list[[g]],
      y_sigma[g],
      dd_weight,
      dens_power
    )
    distrib_diff_vec[g] <- dd_g$diff
    dens_mat[, g] <- dd_g$dens

    distrib_diff_mat[g, ] <- c(dd_g$diff_x, dd_g$diff_y)
    dens_arr[, g, 1] <- dd_g$dens_x
    dens_arr[, g, 2] <- dd_g$dens_y
  }

  mix_dens <- dens_mat %*% prop

  choice_id <- which.min(mix_dens)

  removal_dens <- mix_dens[choice_id]

  distrib_diff <- sqrt(sum(prop * distrib_diff_vec^2))

  list(
    distrib_diff = distrib_diff,
    distrib_diff_vec = distrib_diff_vec,
    choice_id = choice_id,
    removal_dens = removal_dens,
    distrib_diff_mat = distrib_diff_mat
  )
}

# ------------------------------------------------------------------------------

#' distrib_diff_lcwm_g
#'
#' @inheritParams distrib_diff_lcwm
#' @param z_g Component assignment probability vector.
#' @param mu_g Component mean vector.
#' @param sigma_g Component covariance matrix.
#' @param mod_g Component regression model.
#' @param y_sigma_g Component regression standard deviation.
#'
#' @return List of
#' * diff
#' * dens
distrib_diff_lcwm_g <- function(
    x,
    z_g,
    mu_g,
    sigma_g,
    mod_g,
    y_sigma_g,
    dd_weight = 0.5,
    dens_power = 0.5) {
  if (dd_weight == 0 && dens_power == 0) {
    diff_g_x <- 0
    dens_g_x <- 1
  } else {
    dd_g_x <- distrib_diff_mahalanobis(x, z_g, mu_g, sigma_g, log(det(sigma_g)))
    diff_g_x <- dd_g_x$diff
    dens_g_x <- dd_g_x$dens
  }

  if (dd_weight == 1 && dens_power == 1) {
    diff_g_y <- 0
    dens_g_y <- 1
  } else {
    dd_g_y <- distrib_diff_residual(x, z_g, mod_g, y_sigma_g)
    diff_g_y <- dd_g_y$diff
    dens_g_y <- dd_g_y$dens
  }

  diff_g <- sqrt(dd_weight * diff_g_x^2 + (1 - dd_weight) * diff_g_y^2)
  dens_g <- dens_g_x^dens_power * dens_g_y^(1 - dens_power)

  list(
    diff = diff_g,
    dens = dens_g,
    diff_x = diff_g_x,
    diff_y = diff_g_y,
    dens_x = dens_g_x,
    dens_y = dens_g_y
  )
}

# ------------------------------------------------------------------------------

#' distrib_diff_residual
#'
#' @inheritParams distrib_diff_lcwm_g
#'
#' @return List of
#' * diff
#' * dens
distrib_diff_residual <- function(
    x,
    z_g,
    mod_g,
    y_sigma_g) {
  # Let rss_zg = sqrt(sum(z_g * (mod_g$residuals^2)))
  # flexCWM uses y_sigma_g = rss_zg / sqrt(n_g)
  # stats::sigma uses rss_zg / sqrt(nrow(x) - 2))

  n_g <- sum(z_g)

  var_num <- ncol(x)
  reg_param_num <- var_num + 1
  df_g <- n_g - reg_param_num

  param1 <- 1 / 2
  param2 <- (df_g - 1) / 2

  hat_g <- stats::hatvalues(mod_g)

  student_resids_g <- mod_g$residuals / (y_sigma_g * sqrt(1 - hat_g))
  scsqst_res_g <- student_resids_g^2 / df_g
  scsqst_res_ewcdf_g_func <- spatstat.univar::ewcdf(scsqst_res_g, z_g / n_g)

  eps <- 1e-5
  check_seq <- seq(eps, 1, by = eps)
  scsqst_res_ewcdf_g_vals <- scsqst_res_ewcdf_g_func(check_seq)
  beta_cdf_g_vals <- stats::pbeta(check_seq, param1, param2)
  cdf_diffs <- scsqst_res_ewcdf_g_vals - beta_cdf_g_vals

  distrib_diff_g_y <- mean(abs(cdf_diffs))

  dens_g_y <- stats::dnorm(mod_g$residuals, mean = 0, sd = y_sigma_g)

  list(
    diff = distrib_diff_g_y,
    dens = dens_g_y
  )
}
