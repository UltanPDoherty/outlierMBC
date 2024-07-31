#' distrib_diff_gmm
#'
#' @inheritParams ombc_lcwm
#' @param z Component assignment probability matrix.
#' @param prop Vector of component proportions.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
#'
#' @return List of
#' * distrib_diff
#' * distrib_diff_vec
#' * choice_id
#' * min_dens
distrib_diff_gmm <- function(x, z, prop, mu, sigma) {
  obs_num <- nrow(x)
  comp_num <- ncol(z)
  track_num <- 4

  distrib_diff_mat <- matrix(nrow = comp_num, ncol = track_num)
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    dd_g <- distrib_diff_mahalanobis(x, z[, g], mu[[g]], sigma[[g]])
    distrib_diff_mat[g, ] <- dd_g$diff
    dens_mat[, g] <- dd_g$dens
  }

  mix_dens <- dens_mat %*% t(prop)

  choice_id <- which.min(mix_dens)
  min_dens <- mix_dens[choice_id]

  distrib_diff_vec <- as.numeric(prop %*% distrib_diff_mat)

  return(list(
    distrib_diff_mat = distrib_diff_mat,
    distrib_diff_vec = distrib_diff_vec,
    choice_id = choice_id,
    min_dens = min_dens
  ))
}


# ==============================================================================

#' distrib_diff_mahalanobis
#'
#' @inheritParams distrib_diff_lcwm_g
#'
#' @return List of
#' * diff
#' * dens
distrib_diff_mahalanobis <- function(
    x,
    z_g,
    mu_g,
    sigma_g) {
  var_num <- ncol(x)
  n_g <- sum(z_g)

  eps <- 1 / 1000
  check_seq <- seq(eps, 1 - eps, eps)

  checkpoints_x <- stats::qbeta(check_seq, var_num / 2, (n_g - var_num - 1) / 2)

  mahalas_g <- stats::mahalanobis(x, mu_g, (n_g / (n_g - 1)) * sigma_g)
  scaled_mahalas_g <- ((n_g) / (n_g - 1)^2) * mahalas_g
  mahala_ewcdf_g_func <- spatstat.univar::ewcdf(scaled_mahalas_g, z_g / n_g)

  mahala_ewcdf_g <- mahala_ewcdf_g_func(checkpoints_x)

  abs_cdf_diffs <- abs(mahala_ewcdf_g - check_seq)

  distrib_diff_g_x <- c(
    mean(abs_cdf_diffs),
    stats::quantile(abs_cdf_diffs, c(0.5, 0.75, 1))
  )

  dens_g_x <-
    (2 * pi)^(-var_num / 2) * det(sigma_g)^(-0.5) * exp(-mahalas_g / 2)

  return(list(
    diff = distrib_diff_g_x,
    dens = dens_g_x
  ))
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
    alpha = 0.5,
    outlier_type = c("x_and_y", "x_only", "y_only", "hybrid")) {
  outlier_type <- match.arg(outlier_type)

  obs_num <- nrow(x)
  comp_num <- ncol(z)

  distrib_diff_vec <- c()
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    dd_g <- distrib_diff_lcwm_g(
      x, z[, g],
      mu[, g, drop = FALSE],
      as.matrix(sigma[, , g]),
      mod_list[[g]],
      y_sigma[g],
      alpha,
      outlier_type
    )
    distrib_diff_vec[g] <- dd_g$diff
    dens_mat[, g] <- dd_g$dens
  }

  mix_dens <- dens_mat %*% prop

  choice_id <- which.min(mix_dens)

  distrib_diff <- sum(prop * distrib_diff_vec)

  return(list(
    distrib_diff = distrib_diff,
    distrib_diff_vec = distrib_diff_vec,
    choice_id = choice_id
  ))
}

# ==============================================================================

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
    alpha = 0.5,
    outlier_type = c("x_and_y", "x_only", "y_only", "hybrid")) {
  outlier_type <- match.arg(outlier_type)

  if (outlier_type == "x_and_y") {
    dd_g_x <- distrib_diff_mahalanobis(x, z_g, mu_g, sigma_g)
    dd_g_y <- distrib_diff_residual(x, z_g, mod_g, y_sigma_g)

    diff_g <- sqrt(alpha * dd_g_x$diff^2 + (1 - alpha) * dd_g_y$diff^2)
    dens_g <- dd_g_x$dens * dd_g_y$dens
  } else if (outlier_type == "x_only") {
    dd_g_x <- distrib_diff_mahalanobis(x, z_g, mu_g, sigma_g)

    diff_g <- dd_g_x$diff
    dens_g <- dd_g_x$dens
  } else if (outlier_type == "y_only") {
    dd_g_y <- distrib_diff_residual(x, z_g, mod_g, y_sigma_g)

    diff_g <- dd_g_y$diff
    dens_g <- dd_g_y$dens
  } else {
    dd_g_x <- distrib_diff_mahalanobis(x, z_g, mu_g, sigma_g)
    dd_g_y <- distrib_diff_residual(x, z_g, mod_g, y_sigma_g)

    diff_g <- dd_g_x$diff
    dens_g <- dd_g_x$dens * dd_g_y$dens
  }

  return(list(
    diff = diff_g,
    dens = dens_g
  ))
}

# ==============================================================================

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
  var_num <- ncol(x)
  n_g <- sum(z_g)

  eps <- 0.001
  check_seq <- seq(eps, 1 - eps, eps)

  df_g <- round(n_g) - (var_num + 1)

  hat_g <- stats::hatvalues(mod_g)

  student_resids_g <- mod_g$residuals / (y_sigma_g * sqrt(1 - hat_g))
  scsqst_res_g <- student_resids_g^2 / df_g

  checkpoints_y <- stats::qbeta(check_seq, 1 / 2, (df_g - 1) / 2)

  scsqst_res_ecdf_g_func <- spatstat.univar::ewcdf(scsqst_res_g, z_g / n_g)

  scsqst_res_ecdf_g <- scsqst_res_ecdf_g_func(checkpoints_y)
  beta_cdf_g_y <- stats::pbeta(checkpoints_y, 1 / 2, (df_g - 1) / 2)
  distrib_diff_g_y <- mean(abs(scsqst_res_ecdf_g - beta_cdf_g_y))

  dens_g_y <- stats::dnorm(mod_g$residuals, mean = 0, sd = y_sigma_g)

  return(list(
    diff = distrib_diff_g_y,
    dens = dens_g_y
  ))
}
