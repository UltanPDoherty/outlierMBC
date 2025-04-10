#' @title Compute the dissimilarity for a Gaussian mixture model and identify
#' the lowest density observation.
#'
#' @description
#' At each iteration of [ombc_gmm], `distrib_diff_gmm` computes the
#' dissimilarity value of the current Gaussian mixture model. It also
#' identifies the observation with the lowest mixture density.
#'
#' @inheritParams ombc_gmm
#' @param z Component assignment probability matrix.
#' @param prop Vector of component proportions.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
#' @param logdet Vector of log-determinants for covariance matrices.
#'
#' @returns
#' `distrib_diff_gmm` returns a list with the following elements:
#' \describe{
#'   \item{`distrib_diff`}{Aggregated dissimilarity across components.}
#'   \item{`distrib_diff_vec`}{Vector containing dissimilarity value for each
#'                             component.}
#'   \item{`choice_id`}{Index of observation with lowest mixture density.}
#'   \item{`removal_dens`}{Value of the lowest mixture density.}
#' }
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
    distrib_diff = distrib_diff,
    distrib_diff_vec = distrib_diff_vec,
    choice_id = choice_id,
    removal_dens = removal_dens
  )
}

# ==============================================================================

#' @title Compute the dissimilarity for a single multivariate Gaussian
#' distribution.
#'
#' @description
#' Compute the dissimilarity value and observation densities for a single
#' multivariate Gaussian distribution. This could be a whole component in a
#' Gaussian mixture model or the covariate part of a component in a Linear CWM.
#'
#' @inheritParams distrib_diff_gmm
#' @param z_g Assignment probability vector for component g.
#' @param mu_g Mean vector for component g.
#' @param sigma_g Covariance matrix for component g.
#' @param logdet_g Log-determinants of covariance matrix for component g.
#'
#' @returns
#' `distrib_diff_mahalanobis` returns a list with the following elements:
#' \describe{
#'   \item{`diff`}{Dissimilarity value for this component.}
#'   \item{`dens`}{Gaussian density of all observations for this component.}
#'   \item{`mahalas`}{Scaled squared sample Mahalanobis distances for all
#'                    observations with respect to this component.}
#' }
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

#' @title Compute the dissimilarity for a linear cluster-weighted model and
#' identify the lowest density observation.
#'
#' @description
#' At each iteration of [ombc_lcwm], `distrib_diff_lcwm` computes the
#' dissimilarity value of the current linear cluster-weighted model. It also
#' identifies the observation with the lowest mixture density.
#'
#' @inheritParams ombc_lcwm
#' @inheritParams distrib_diff_gmm
#' @param mu Matrix of component mean vectors.
#' @param sigma Array of component covariance matrices.
#' @param mod_list List of component regression models.
#' @param y_sigma Vector of component regression standard deviations.
#'
#' @returns
#' `distrib_diff_lcwm_lcwm` returns a list with the following elements:
#' \describe{
#'   \item{`distrib_diff`}{Aggregated dissimilarity across components.}
#'   \item{`distrib_diff_vec`}{Vector containing dissimilarity value for each
#'                           component.}
#'   \item{`choice_id`}{Index of observation with lowest mixture density.}
#'   \item{`removal_dens`}{Value of the lowest mixture density.}
#'   \item{`distrib_diff_mat`}{Two-column matrix containing response and
#'                             covariate dissimilarities across components.}
#' }
distrib_diff_lcwm <- function(
    x,
    z,
    prop,
    mu,
    sigma,
    mod_list,
    y_sigma,
    dd_weight = 0.5) {
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
      dd_weight
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

# ==============================================================================

#' @title Compute the dissimilarity for a single component of a Linear CWM.
#'
#' @description
#' Computes the covariate dissimilarity value, the response dissimilarity value,
#' and their aggregated dissimilarity value. It also obtains the covariate,
#' response, and joint densities for every observation.
#'
#' @inheritParams distrib_diff_lcwm
#' @param z_g Component assignment probability vector.
#' @param mu_g Component mean vector for the covariates.
#' @param sigma_g Component covariance matrix for the covariates.
#' @param mod_g Component regression model.
#' @param y_sigma_g Component regression standard deviation for the response.
#'
#' @returns
#' `distrib_diff_lcwm_lcwm_g` returns a list with the following elements:
#' \describe{
#'   \item{`diff`}{Aggregated dissimilarity value for this component.}
#'   \item{`dens`}{Joint (covariate & response) density of all observations for
#'                 this component.}
#'   \item{`diff_x`}{Covariate dissimilarity value for this component.}
#'   \item{`diff_y`}{Response dissimilarity value for this component.}
#'   \item{`dens_x`}{Covariate density of all observations for this component.}
#'   \item{`dens_y`}{Response density of all observations for this component.}
#' }
distrib_diff_lcwm_g <- function(
    x,
    z_g,
    mu_g,
    sigma_g,
    mod_g,
    y_sigma_g,
    dd_weight = 0.5) {
  dd_g_x <- distrib_diff_mahalanobis(x, z_g, mu_g, sigma_g, log(det(sigma_g)))
  diff_g_x <- dd_g_x$diff
  dens_g_x <- dd_g_x$dens

  dd_g_y <- distrib_diff_residual(x, z_g, mod_g, y_sigma_g)
  diff_g_y <- dd_g_y$diff
  dens_g_y <- dd_g_y$dens

  diff_g <- sqrt(dd_weight * diff_g_x^2 + (1 - dd_weight) * diff_g_y^2)
  dens_g <- dens_g_x * dens_g_y

  list(
    diff = diff_g,
    dens = dens_g,
    diff_x = diff_g_x,
    diff_y = diff_g_y,
    dens_x = dens_g_x,
    dens_y = dens_g_y
  )
}

# ==============================================================================

#' @title Compute the response dissimilarity for a single component of a Linear
#' CWM.
#'
#' @description
#' Computes the response dissimilarity value and the response density for every
#' observation.
#'
#' @inheritParams distrib_diff_lcwm_g
#'
#' @returns
#' `distrib_diff_lcwm_residual` returns a list with the following elements:
#' \describe{
#'   \item{`diff`}{Response dissimilarity value for this component.}
#'   \item{`dens`}{Response density of all observations for this component.}
#' }
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
