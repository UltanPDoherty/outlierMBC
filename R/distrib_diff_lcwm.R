#' distrib_diff_lcwm
#'
#' @inheritParams ombc_lcwm
#' @param z Component assignment probability matrix.
#' @param prop Vector of component proportions.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
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
    alpha = 0.5) {
  obs_num <- nrow(x)
  comp_num <- ncol(z)

  distrib_diff_vec <- c()
  scaled_mahalas <- matrix(nrow = obs_num, ncol = comp_num)
  dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    out_g <- distrib_diff_lcwm_g(
      x, z[, g],
      mu[, g, drop = FALSE],
      as.matrix(sigma[, , g]),
      mod_list[[g]],
      y_sigma[g],
      alpha
    )
    distrib_diff_vec[g] <- out_g$distrib_diff_g
    scaled_mahalas[, g] <- out_g$scaled_mahalas_g
    dens_mat[, g] <- out_g$dens_g
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

#' distrib_diff_lcwm_g
#'
#' @inheritParams distrib_diff_lcwm
#' @param z_g Component assignment probability vector.
#' @param mu_g Component mean vector.
#' @param sigma_g Component covariance matrix.
#' @param mod_g Component regression model.
#' @param y_sigma_g Component regression standard deviation.
#'
#'
#' @return List of
#' * distrib_diff_g
#' * scaled_mahalas_g
#' * dens_x_g
distrib_diff_lcwm_g <- function(
    x,
    z_g,
    mu_g,
    sigma_g,
    mod_g,
    y_sigma_g,
    alpha = 0.5) {
  var_num <- ncol(x)
  n_g <- sum(z_g)

  eps <- 0.001
  check_seq <- seq(eps, 1 - eps, eps)

  # --------

  df_g <- round(n_g) - 2

  hat_g <- stats::hatvalues(mod_g)

  student_resids_g <- mod_g$residuals / (y_sigma_g * sqrt(1 - hat_g))
  scsqst_res_g <- student_resids_g^2 / df_g

  checkpoints_y <- stats::qbeta(check_seq, 1 / 2, (df_g - 1) / 2)

  scsqst_res_ecdf_g_func <- spatstat.geom::ewcdf(scsqst_res_g, z_g / n_g)

  scsqst_res_ecdf_g <- scsqst_res_ecdf_g_func(checkpoints_y)
  beta_cdf_y_g <- stats::pbeta(checkpoints_y, 1 / 2, (df_g - 1) / 2)
  distrib_diff_y_g <- mean(abs(scsqst_res_ecdf_g - beta_cdf_y_g))

  dens_y_g <- stats::dnorm(mod_g$residuals, mean = 0, sd = y_sigma_g)

  # ----------

  checkpoints_x <- stats::qbeta(check_seq, var_num / 2, (n_g - var_num - 1) / 2)

  mahalas_g <- stats::mahalanobis(x, mu_g, (n_g / (n_g - 1)) * sigma_g)
  scaled_mahalas_g <- ((n_g) / (n_g - 1)^2) * mahalas_g
  mahala_ewcdf_g_func <- spatstat.geom::ewcdf(scaled_mahalas_g, z_g / n_g)

  mahala_ewcdf_g <- mahala_ewcdf_g_func(checkpoints_x)
  beta_cdf_x_g <- stats::pbeta(
    checkpoints_x, var_num / 2, (n_g - var_num - 1) / 2
  )
  distrib_diff_x_g <- mean(abs(mahala_ewcdf_g - beta_cdf_x_g))

  dens_x_g <-
    (2 * pi)^(-var_num / 2) * det(sigma_g)^(-0.5) * exp(-mahalas_g / 2)

  # ------------

  distrib_diff_g <- sqrt(
    alpha * distrib_diff_x_g^2 + (1 - alpha) * distrib_diff_y_g^2
  )

  dens_g <- dens_x_g * dens_y_g

  return(list(
    distrib_diff_g = distrib_diff_g,
    scaled_mahalas_g = scaled_mahalas_g,
    dens_g = dens_g
  ))
}
