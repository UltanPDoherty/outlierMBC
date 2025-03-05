#' ombc_lcwm
#'
#' @description
#' Iterative Detection & Identification of Outliers for a Linear
#' Cluster-Weighted Model
#'
#' @inheritParams ombc_gmm
#' @param xy `data.frame` containing covariates and response.
#' @param x Covariate data only.
#' @param y_formula Regression formula.
#' @param dd_weight .
#' @param dens_power .
#'
#' @return List of
#' *labels = labels,
#' * outlier_bool = outlier_bool,
#' * outlier_num = outlier_num,
#' * outlier_rank = outlier_rank,
#' * gross_outs = gross_outs,
#' * lcwm = lcwm,
#' * loglike = loglike,
#' * removal_dens = removal_dens,
#' * distrib_diff_vec = distrib_diff_vec,
#' * distrib_diff_mat = distrib_diff_mat,
#' * call = this_call,
#' * version = ombc_version,
#' * conv_status = conv_status
#' @export
#'
#' @examples
#'
#' gross_lcwm_k3n1000o10 <- find_gross(lcwm_k3n1000o10, 20)
#'
#' ombc_lcwm_k3n1000o10 <- ombc_lcwm(
#'   xy = lcwm_k3n1000o10[, c("X1", "Y")],
#'   x = lcwm_k3n1000o10$X1,
#'   y_formula = Y ~ X1,
#'   comp_num = 3,
#'   max_out = 20,
#'   mnames = "V",
#'   gross_outs = gross_lcwm_k3n1000o10$gross_bool
#' )
#'
ombc_lcwm <- function(
    xy,
    x,
    y_formula,
    comp_num,
    max_out,
    gross_outs = rep(FALSE, nrow(x)),
    init_scheme = c("update", "reinit", "reuse"),
    mnames = "VVV",
    nmax = 1000,
    atol = 1e-8,
    init_z = NULL,
    init_method = c("hc", "kmpp"),
    init_scaling = TRUE,
    kmpp_seed = 123,
    print_interval = Inf,
    dd_weight = 0.5,
    dens_power = 0.5) {
  init_method <- match.arg(init_method)
  init_scheme <- match.arg(init_scheme)
  init_model <- NULL

  this_call <- call(
    "ombc_gmm",
    "xy" = substitute(xy), "x" = substitute(x), "y_formula" = y_formula,
    "comp_num" = comp_num, "max_out" = max_out,
    "gross_outs" = substitute(gross_outs), "init_scheme" = init_scheme,
    "mnames" = mnames, "nmax" = nmax, "atol" = atol,
    "init_z" = substitute(init_z),
    "init_method" = init_method, "init_scaling" = init_scaling,
    "kmpp_seed" = kmpp_seed, "print_interval" = print_interval,
    "dd_weight" = dd_weight, "dens_power" = dens_power
  )

  ombc_version <- utils::packageVersion("outlierMBC")

  xy0 <- xy
  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x)

  xy1 <- scale(xy0, center = init_scaling, scale = init_scaling)
  dist_mat0 <- as.matrix(stats::dist(xy1))
  dist_mat <- dist_mat0

  gross_num <- sum(gross_outs)
  xy <- xy[!gross_outs, , drop = FALSE]
  x <- x[!gross_outs, , drop = FALSE]
  max_out <- max_out - gross_num
  dist_mat <- dist_mat[!gross_outs, !gross_outs]

  if (!is.null(init_model) && !is.null(init_z)) {
    stop("Only one of init_model and init_z may be provided.")
  } else if (!is.null(init_z)) {
    z <- init_z
  } else {
    z <- get_init_z(
      comp_num = comp_num, dist_mat = dist_mat, x = xy,
      init_method = init_method, kmpp_seed = kmpp_seed
    )
  }

  dd_min <- Inf

  conv_status <- c()
  loglike <- c()
  removal_dens <- c()
  distrib_diff_arr <- array(dim = c(max_out + 1, comp_num, 2))
  distrib_diff_mat <- matrix(nrow = max_out + 1, ncol = comp_num)
  distrib_diff_vec <- double(max_out + 1)
  outlier_rank_temp <- rep(0, obs_num - gross_num)
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    if (init_scheme %in% c("update", "reuse")) {
      invisible(utils::capture.output(lcwm <- flexCWM::cwm(
        formulaY = y_formula,
        familyY = stats::gaussian(link = "identity"),
        data = xy,
        Xnorm = x,
        modelXnorm = mnames,
        k = comp_num,
        initialization = "manual",
        start.z = z,
        iter.max = nmax,
        threshold = atol
      )))
    } else {
      reinit_z <- get_init_z(
        comp_num = comp_num, dist_mat = dist_mat, x = xy,
        init_method = init_method, kmpp_seed = kmpp_seed
      )
      invisible(utils::capture.output(lcwm <- flexCWM::cwm(
        formulaY = y_formula,
        familyY = stats::gaussian(link = "identity"),
        data = xy,
        Xnorm = x,
        modelXnorm = mnames,
        k = comp_num,
        initialization = "manual",
        start.z = reinit_z,
        iter.max = nmax,
        threshold = atol
      )))
    }

    loglike[i] <- lcwm$models[[1]]$logLik
    conv_status[i] <- lcwm$models[[1]]$converged

    mod_list <- lapply(lcwm$models[[1]]$GLModel, function(x) x$model)
    names(mod_list) <- paste0("comp.", seq_along(mod_list))
    y_sigma <- vapply(lcwm$models[[1]]$GLModel, function(x) x$sigma, double(1L))

    dd <- distrib_diff_lcwm(
      x,
      lcwm$models[[1]]$posterior,
      lcwm$models[[1]]$prior,
      lcwm$models[[1]]$concomitant$normal.mu,
      lcwm$models[[1]]$concomitant$normal.Sigma,
      mod_list,
      y_sigma,
      dd_weight,
      dens_power
    )

    distrib_diff_arr[i, , ] <- dd$distrib_diff_mat
    distrib_diff_mat[i, ] <- dd$distrib_diff_vec
    distrib_diff_vec[i] <- dd$distrib_diff
    removal_dens[i] <- dd$removal_dens

    outlier_rank_temp[!outlier_rank_temp][dd$choice_id] <- i
    x <- x[-dd$choice_id, , drop = FALSE]
    xy <- xy[-dd$choice_id, , drop = FALSE]
    dist_mat <- dist_mat[-dd$choice_id, -dd$choice_id]

    if (dd$distrib_diff < dd_min) {
      dd_min <- dd$distrib_diff
      best_z <- lcwm$models[[1]]$posterior
    }

    if (init_scheme %in% c("update")) {
      z <- lcwm$models[[1]]$posterior[-dd$choice_id, , drop = FALSE]
    } else if (init_scheme == "reuse") {
      z <- z[-dd$choice_id, , drop = FALSE]
    }
  }

  outlier_rank <- double(length(gross_outs))
  outlier_rank[gross_outs] <- 1
  outlier_rank[!gross_outs] <- outlier_rank_temp +
    gross_num * (outlier_rank_temp != 0)

  outlier_num <- which.min(distrib_diff_vec)
  outlier_num <- outlier_num - 1 + gross_num

  outlier_bool <- logical(obs_num)
  labels <- integer(obs_num)

  outlier_bool <- outlier_rank <= outlier_num & outlier_rank != 0

  invisible(utils::capture.output(lcwm <- flexCWM::cwm(
    formulaY = y_formula,
    familyY = stats::gaussian(link = "identity"),
    data = xy0[!outlier_bool, , drop = FALSE],
    Xnorm = x0[!outlier_bool, , drop = FALSE],
    modelXnorm = mnames,
    k = comp_num,
    initialization = "manual",
    start.z = best_z,
    iter.max = nmax,
    threshold = atol,
    pwarning = TRUE
  )))

  labels[!outlier_bool] <- lcwm$models[[1]]$cluster

  colnames(distrib_diff_mat) <- paste0("k", seq_len(comp_num))

  return(list(
    labels = labels,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    gross_outs = gross_outs,
    lcwm = lcwm,
    loglike = loglike,
    removal_dens = removal_dens,
    distrib_diff_vec = distrib_diff_vec,
    distrib_diff_mat = distrib_diff_mat,
    distrib_diff_arr = distrib_diff_arr,
    call = this_call,
    version = ombc_version,
    conv_status = conv_status
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

  return(list(
    distrib_diff = distrib_diff,
    distrib_diff_vec = distrib_diff_vec,
    choice_id = choice_id,
    removal_dens = removal_dens,
    distrib_diff_mat = distrib_diff_mat
  ))
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

  return(list(
    diff = diff_g,
    dens = dens_g,
    diff_x = diff_g_x,
    diff_y = diff_g_y,
    dens_x = dens_g_x,
    dens_y = dens_g_y
  ))
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

  return(list(
    diff = distrib_diff_g_y,
    dens = dens_g_y
  ))
}

# ==============================================================================

#' simulate_lcwm
#'
#' @description
#' Simulate a multiple linear regression model with response variable outliers.
#'
#' @param n Vector of component sizes.
#' @param mu List of component mean vectors.
#' @param sigma List of component covariance matrices.
#' @param beta List of component regression coefficient vectors.
#' @param error_sd Vector of component regression error standard deivations.
#' @param outlier_num Desired number of outliers.
#' @param outlier_type .
#' @param seed Seed.
#' @param crit_val Critical value for uniform sample rejection.
#' @param range_multipliers .
#' @param print_interval How frequently the iteration count is printed.
#'
#' @return `data.frame`:
#' * $X1: first covariate
#' * ...
#' * $Y: response
#' * $G: label
#'
#' @export
#'
#' @examples
#' lcwm_k3n1000o10 <- simulate_lcwm(
#'   n = c(300, 300, 400),
#'   mu = list(c(3), c(6), c(3)),
#'   sigma = list(as.matrix(1), as.matrix(0.1), as.matrix(1)),
#'   beta = list(c(0, 0), c(-75, 15), c(0, 5)),
#'   error_sd = c(1, 1, 1),
#'   outlier_num = c(3, 3, 4),
#'   outlier_type = "x_and_y",
#'   seed = 123,
#'   crit_val = 1 - 1e-6,
#'   range_multipliers = c(1, 2)
#' )
#'
#' plot(
#'   lcwm_k3n1000o10[, c("X1", "Y")],
#'   col = lcwm_k3n1000o10$G + 1,
#'   pch = lcwm_k3n1000o10$G + 1
#' )
simulate_lcwm <- function(
    n,
    mu,
    sigma,
    beta,
    error_sd,
    outlier_num,
    outlier_type = c("x_and_y", "x_only", "y_only"),
    seed = 123,
    crit_val = 0.9999,
    range_multipliers = c(1.5, 1.5),
    print_interval = Inf) {
  outlier_type <- match.arg(outlier_type)

  var_num <- length(mu[[1]])
  comp_num <- length(n)

  set.seed(seed)
  observations <- covariates <- errors <- responses <- uniform_spans <- list()
  for (g in seq_len(comp_num)) {
    covariates[[g]] <- mvtnorm::rmvnorm(n[g], mu[[g]], sigma[[g]])

    errors[[g]] <- stats::rnorm(n[g], 0, sd = error_sd[g])

    responses[[g]] <-
      errors[[g]] + beta[[g]][1] + covariates[[g]] %*% beta[[g]][-1]

    observations[[g]] <- cbind(covariates[[g]], responses[[g]], rep(g, n[g]))

    uniform_spans[[g]] <- uniform_spans_lcwm(
      range_multipliers, covariates[[g]], errors[[g]]
    )
  }

  set.seed(seed)
  outliers <- lapply(outlier_num, function(x) matrix(NA, x, var_num + 2))
  for (g in seq_len(comp_num)) {
    for (j in seq_len(outlier_num[g])) {
      outliers[[g]][j, ] <- uniform_outlier_ombc(
        outlier_type, mu, sigma, beta, error_sd, g, uniform_spans, crit_val
      )
    }
  }

  lcwm <- as.data.frame(rbind(
    Reduce(rbind, observations),
    Reduce(rbind, outliers)
  ))
  colnames(lcwm) <- c(paste0("X", seq_len(var_num)), "Y", "G")

  return(lcwm)
}

# ------------------------------------------------------------------------------

#' Obtain the span of the observations for each component.
#'
#' @inheritParams simulate_lcwm
#' @param covariates_g Covariate values of the sampled observations.
#' @param errors_g Response errors of the sampled observations.
#'
#' @return `matrix`: minimum and maximum columns, rows for each covariate and
#'                   one for the response errors.
uniform_spans_lcwm <- function(range_multipliers, covariates_g, errors_g) {
  ranges_x <- apply(covariates_g, 2, range)
  centres_x <- colMeans(ranges_x)
  widths_x <- apply(ranges_x, 2, diff)
  mins_x <- centres_x - ((range_multipliers[1]) / 2) * widths_x
  maxs_x <- centres_x + ((range_multipliers[1]) / 2) * widths_x

  range_err <- range(errors_g)
  centre_err <- mean(range_err)
  width_err <- diff(range_err)
  min_err <- centre_err - (range_multipliers[2] / 2) * width_err
  max_err <- centre_err + (range_multipliers[2] / 2) * width_err

  spans <- cbind(c(mins_x, min_err), c(maxs_x, max_err))

  return(spans)
}

# ------------------------------------------------------------------------------

#' Produce a single sample that passes the outlier checks.
#'
#' @inheritParams simulate_lcwm
#' @param g Component index.
#' @param uniform_spans Covariate and response error spans.
#'
#' @return Vector consisting of covariate values, response value, and label 0.
uniform_outlier_ombc <- function(
    outlier_type,
    mu, sigma, beta, error_sd, g,
    uniform_spans, crit_val) {
  test <- FALSE

  while (!test) {
    uniform_sample <- uniform_sample_lcwm(
      outlier_type,
      mu[[g]], sigma[[g]],
      beta[[g]], error_sd[g],
      uniform_spans[[g]]
    )

    test <- test_outlier_ombc(
      outlier_type,
      mu, sigma, beta, error_sd,
      uniform_sample$x, uniform_sample$y, crit_val
    )
  }

  return(c(uniform_sample$x, uniform_sample$y, 0))
}

# ------------------------------------------------------------------------------

#' Sample a potential outlier.
#'
#' @inheritParams simulate_lcwm
#' @param mu_g Covariate mean vector for one component.
#' @param sigma_g Covariate covariance matrix for one component.
#' @param beta_g Regression coefficient vector for one component.
#' @param error_sd_g Regression error standard deviation for one component.
#' @param uniform_spans_g Covariate and response error ranges for one component.
#'
#' @return List:
#' * x = Covariate sample
#' * y = Response sample
uniform_sample_lcwm <- function(
    outlier_type, mu_g, sigma_g, beta_g, error_sd_g, uniform_spans_g) {
  var_num <- length(beta_g) - 1

  if (outlier_type == "y_only") {
    outlier_x_g <- mvtnorm::rmvnorm(1, mu_g, sigma_g)
  } else {
    outlier_x_g <- rep(NA, var_num)
    for (p in seq_len(var_num)) {
      outlier_x_g[p] <- stats::runif(
        1, uniform_spans_g[p, 1], uniform_spans_g[p, 2]
      )
    }
  }

  if (outlier_type == "x_only") {
    outlier_err <- stats::rnorm(1, 0, error_sd_g)
  } else {
    outlier_err <- stats::runif(
      1, uniform_spans_g[var_num + 1, 1], uniform_spans_g[var_num + 1, 2]
    )
  }

  outlier_fitted <- beta_g[1] + outlier_x_g %*% beta_g[-1]

  outlier_y_g <- outlier_fitted + outlier_err

  return(list(x = outlier_x_g, y = outlier_y_g))
}

# ------------------------------------------------------------------------------

#' Check whether a new sample is an outlier for each component.
#'
#' @inheritParams simulate_lcwm
#' @param x_sample New covariate sample.
#' @param y_sample New response sample.
#'
#' @return Logical: TRUE if the new sample is an outlier for each component.
test_outlier_ombc <- function(
    outlier_type,
    mu, sigma, beta, error_sd,
    x_sample, y_sample, crit_val) {
  comp_num <- length(mu)
  var_num <- length(mu[[1]])

  prob_x <- prob_y <- temp_outliers_fitted <- checks <- rep(NA, comp_num)

  for (h in seq_len(comp_num)) {
    prob_x[h] <- stats::pchisq(
      stats::mahalanobis(x_sample, mu[[h]], sigma[[h]]),
      df = var_num, lower.tail = FALSE
    )

    temp_outliers_fitted[h] <-
      (beta[[h]][1] + x_sample %*% beta[[h]][-1])

    prob_y[h] <- stats::pnorm(
      abs(y_sample - temp_outliers_fitted[h]),
      mean = 0, sd = error_sd[h], lower.tail = FALSE
    )

    checks[h] <- switch(outlier_type,
      x_and_y = (prob_x[h] * prob_y[h]) < 1 - crit_val,
      x_only = prob_x[h] < 1 - crit_val,
      y_only = prob_y[h] < 1 - crit_val
    )
  }

  return(all(checks))
}

# ==============================================================================

#' Fit a LCWM to the backtrack solution.
#'
#' @inheritParams ombc_lcwm
#' @inheritParams backtrack
#' @inheritParams backtrack_gmm
#' @param ombc_lcwm_out Output from ombc_lcwm.
#'
#' @returns List:
#' * labels
#' * outlier_bool
#' * outlier_num
#' * lcwm
#' * call
#'
#' @export
#'
#' @examples
#'
#' gross_lcwm_k3n1000o10 <- find_gross(lcwm_k3n1000o10, 20)
#'
#' ombc_lcwm_k3n1000o10 <- ombc_lcwm(
#'   xy = lcwm_k3n1000o10[, c("X1", "Y")],
#'   x = lcwm_k3n1000o10$X1,
#'   y_formula = Y ~ X1,
#'   comp_num = 2,
#'   max_out = 20,
#'   mnames = "V",
#'   gross_outs = gross_lcwm_k3n1000o10$gross_bool
#' )
#'
#' backtrack_lcwm_k3n1000o10 <- backtrack_lcwm(
#'   xy = lcwm_k3n1000o10[, c("X1", "Y")],
#'   x = lcwm_k3n1000o10$X1,
#'   ombc_lcwm_out = ombc_lcwm_k3n1000o10
#' )
backtrack_lcwm <- function(
    xy, x, ombc_lcwm_out,
    max_total_rise = 0.1, max_step_rise = 0.05,
    init_z = NULL, manual_outlier_num = NULL) {
  this_call <- call(
    "backtrack_lcwm",
    "xy" = substitute(xy), "x" = substitute(x),
    "ombc_lcwm_out" = substitute(ombc_lcwm_out),
    "max_total_rise" = max_total_rise, "max_step_rise" = max_step_rise,
    "init_z" = substitute(init_z)
  )

  x <- as.matrix(x)
  x0 <- x
  xy0 <- xy

  gross_num <- sum(ombc_lcwm_out$gross_outs)

  if (is.null(manual_outlier_num)) {
    backtrack_out <- backtrack(
      ombc_lcwm_out$distrib_diff_vec, max_total_rise, max_step_rise
    )

    if (backtrack_out$backtrack$ind == backtrack_out$minimum$ind) {
      cat(paste0(
        "backtrack stayed at the minimum.",
        "backtrack_lcwm will return ombc_lcwm results directly.\n"
      ))

      return(list(
        "labels" = ombc_lcwm_out$labels,
        "outlier_bool" = ombc_lcwm_out$outlier_bool,
        "outlier_num" = ombc_lcwm_out$outlier_num,
        "lcwm" = ombc_lcwm_out$lcwm,
        "call" = this_call
      ))
    }

    outlier_num <- backtrack_out$backtrack$ind - 1 + gross_num
  } else {
    outlier_num <- manual_outlier_num
  }

  outlier_bool <-
    ombc_lcwm_out$outlier_rank <= outlier_num & ombc_lcwm_out$outlier_rank != 0

  init_scheme <- ombc_lcwm_out$call$init_scheme
  init_scaling <- ombc_lcwm_out$call$init_scaling

  stopifnot(
    "init_scheme must be 'update' or 'reinit' or 'reuse'." = (
      init_scheme %in% c("update", "reinit", "reuse")
    )
  )
  stopifnot(
    "init_z cannot be used with the 'reinit' init_scheme." = (
      (init_scheme != "reinit") || is.null(init_z)
    )
  )

  if (!is.null(init_z)) {
    z0 <- init_z
  } else if (init_scheme != "reinit") {
    xy1 <- scale(xy0, center = init_scaling, scale = init_scaling)
    z0 <- get_init_z(
      comp_num = ombc_lcwm_out$call$comp_num,
      dist_mat = as.matrix(
        stats::dist(xy1[!ombc_lcwm_out$gross_outs, ])
      ),
      x = x0[!ombc_lcwm_out$gross_outs, , drop = FALSE],
      init_method = ombc_lcwm_out$call$init_method,
      kmpp_seed = ombc_lcwm_out$call$kmpp_seed
    )
  }

  if (init_scheme == "reinit") {
    xy1 <- scale(xy0, center = init_scaling, scale = init_scaling)
    z <- get_init_z(
      comp_num = ombc_lcwm_out$call$comp_num,
      dist_mat = as.matrix(
        stats::dist(xy1[!outlier_bool, , drop = FALSE])
      ),
      x = xy0[!outlier_bool, , drop = FALSE],
      init_method = ombc_lcwm_out$call$init_method,
      kmpp_seed = ombc_lcwm_out$call$kmpp_seed
    )

    invisible(utils::capture.output(lcwm <- flexCWM::cwm(
      formulaY = ombc_lcwm_out$call$y_formula,
      familyY = stats::gaussian(link = "identity"),
      data = xy0[!outlier_bool, , drop = FALSE],
      Xnorm = x0[!outlier_bool, , drop = FALSE],
      modelXnorm = ombc_lcwm_out$call$mnames,
      k = ombc_lcwm_out$call$comp_num,
      initialization = "manual",
      start.z = z,
      iter.max = ombc_lcwm_out$call$nmax,
      threshold = ombc_lcwm_out$call$atol
    )))
  } else if (init_scheme == "reuse") {
    short_outlier_bool <- outlier_bool[!ombc_lcwm_out$gross_outs]
    z <- z0[!short_outlier_bool, !short_outlier_bool]

    invisible(utils::capture.output(lcwm <- flexCWM::cwm(
      formulaY = ombc_lcwm_out$call$y_formula,
      familyY = stats::gaussian(link = "identity"),
      data = xy0[!outlier_bool, , drop = FALSE],
      Xnorm = x0[!outlier_bool, , drop = FALSE],
      modelXnorm = ombc_lcwm_out$call$mnames,
      k = ombc_lcwm_out$call$comp_num,
      initialization = "manual",
      start.z = z,
      iter.max = ombc_lcwm_out$call$nmax,
      threshold = ombc_lcwm_out$call$atol
    )))
  } else {
    xy <- xy0[!ombc_lcwm_out$gross_outs, , drop = FALSE]
    x <- x0[!ombc_lcwm_out$gross_outs, , drop = FALSE]
    z <- z0

    cat("Fitting backtrack model:\n")
    if (outlier_num > gross_num) {
      temp_outlier_rank <- ombc_lcwm_out$outlier_rank[!ombc_lcwm_out$gross_outs]

      prog_bar <- utils::txtProgressBar(
        gross_num, outlier_num,
        style = 3
      )
      removals <- c()
      for (i in seq(gross_num, outlier_num)) {
        utils::setTxtProgressBar(prog_bar, i)

        invisible(utils::capture.output(lcwm <- flexCWM::cwm(
          formulaY = ombc_lcwm_out$call$y_formula,
          familyY = stats::gaussian(link = "identity"),
          data = xy,
          Xnorm = x,
          modelXnorm = ombc_lcwm_out$call$mnames,
          k = ombc_lcwm_out$call$comp_num,
          initialization = "manual",
          start.z = z,
          iter.max = ombc_lcwm_out$call$nmax,
          threshold = ombc_lcwm_out$call$atol
        )))

        next_removal <- which(temp_outlier_rank == i + 1)
        removals <- append(removals, next_removal)
        xy <- xy[-next_removal, , drop = FALSE]
        x <- x[-next_removal, , drop = FALSE]
        z <- lcwm$models[[1]]$posterior[-next_removal, -next_removal]
        temp_outlier_rank <- temp_outlier_rank[-next_removal]
      }
      close(prog_bar)
    } else {
      invisible(utils::capture.output(lcwm <- flexCWM::cwm(
        formulaY = ombc_lcwm_out$call$y_formula,
        familyY = stats::gaussian(link = "identity"),
        data = xy,
        Xnorm = x,
        modelXnorm = ombc_lcwm_out$call$mnames,
        k = ombc_lcwm_out$call$comp_num,
        initialization = "manual",
        start.z = z,
        iter.max = ombc_lcwm_out$call$nmax,
        threshold = ombc_lcwm_out$call$atol
      )))
    }
    cat("backtrack model fitted.\n")
  }

  labels <- integer(nrow(x))
  labels[outlier_bool] <- 0
  labels[!outlier_bool] <- lcwm$models[[1]]$cluster

  return(list(
    "labels" = labels,
    "outlier_bool" = outlier_bool,
    "outlier_num" = outlier_num,
    "lcwm" = lcwm,
    "call" = this_call
  ))
}
