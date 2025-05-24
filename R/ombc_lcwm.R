#' @title Sequentially identify outliers while fitting a linear cluster-weighted
#'        model.
#'
#' @description
#' This function performs model-based clustering, clusterwise regression, and
#' outlier identification. It does so by iteratively fitting a linear
#' cluster-weighted model and removing the observation that is least likely
#' under the model. Its procedure is summarised below:
#'
#' 1. Fit a linear cluster-weighted model to the data.
#' 2. Compute a dissimilarity between the theoretical and observed distributions
#'    of the scaled squared sample Mahalanobis distances for each mixture
#'    component.
#' 3. Compute a dissimilarity between the theoretical and observed distributions
#'    of the scaled squared studentised residuals for each mixture component.
#' 4. Aggregate these two dissimilarities to obtain one dissimilarity value
#'    for each component.
#' 5. Aggregate across the components to obtain a single dissimilarity value.
#' 6. Remove the observation  with the lowest mixture density.
#' 7. Repeat Steps 1-6 until `max_out` observations have been removed.
#' 8. Identify the number of outliers which minimised the aggregated
#'    dissimilarity, remove only those observations, and fit a linear
#'    cluster-weighted model to the remaining data.
#'
#' @inheritParams ombc_gmm
#' @param xy `data.frame` containing covariates and response.
#' @param x Covariate data only.
#' @param y_formula Regression formula.
#' @param nmax Maximum number of iterations for `flexCWM::cwm`.
#' @param atol EM convergence threshold for `flexCWM::cwm`.
#' @param dd_weight A value between `0` and `1` which controls the weighting of
#'                  the response and covariate dissimilarities when aggregating.
#'
#' @returns
#' `ombc_lcwm` returns an object of class "outliermbc_lcwm", which is
#' essentially a list with the following elements:
#' \describe{
#'   \item{`labels`}{Vector of mixture component labels with outliers denoted by
#'                   0.}
#'   \item{`outlier_bool`}{Logical vector indicating if an observation has been
#'                         classified as an outlier.}
#'   \item{`outlier_num`}{Number of observations classified as outliers.}
#'   \item{`outlier_rank`}{Order in which observations are removed from the data
#'                         set. Observations which were provisionally removed,
#'                         including those that were eventually not classified
#'                         as outliers, are ranked from `1` to `max_out`. All
#'                         gross outliers have rank `1`. If there are
#'                         `gross_num` gross outliers, then the observations
#'                         removed during the main algorithm itself will be
#'                         numbered from `gross_num + 1` to `max_out`.
#'                         Observations that were ever removed have rank `0`.}
#'   \item{`gross_outs`}{Logical vector identifying the gross outliers. This is
#'                       identical to the `gross_outs` vector passed to this
#'                       function as an argument / input.}
#'   \item{`lcwm`}{Output from `flexCWM::cwm` fitted to the non-outlier
#'                 observations.}
#'   \item{`loglike`}{Vector of log-likelihood values for each iteration.}
#'   \item{`removal_dens`}{Vector of mixture densities for the removed
#'                         observations. These are the lowest mixture densities
#'                         at each iteration.}
#'   \item{`distrib_diff_vec`}{Vector of aggregated cross-component
#'                             dissimilarity values for each iteration.}
#'   \item{`distrib_diff_mat`}{Matrix of component-specific dissimilarity values
#'                             for each iteration.}
#'   \item{`distrib_diff_arr`}{Array of component-specific response and
#'                             covariate dissimilarity values for each
#'                             iteration.}
#'   \item{`call`}{Arguments / parameter values used in this function call.}
#'   \item{`version`}{Version of `outlierMBC` used in this function call.}
#'   \item{`conv_status`}{Logical vector indicating which iterations' mixture
#'                        models reached convergence during model-fitting.}
#' }
#'
#' @export
#'
#' @examples
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
    verbose = TRUE,
    dd_weight = 0.5) {
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
    "kmpp_seed" = kmpp_seed, "verbose" = verbose,
    "dd_weight" = dd_weight
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
    if (init_scheme %in% c("update", "reuse")) {
      suppressWarnings(invisible(utils::capture.output(lcwm <- flexCWM::cwm(
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
      ))))
    } else {
      reinit_z <- get_init_z(
        comp_num = comp_num, dist_mat = dist_mat, x = xy,
        init_method = init_method, kmpp_seed = kmpp_seed
      )
      suppressWarnings(invisible(utils::capture.output(lcwm <- flexCWM::cwm(
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
      ))))
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
      dd_weight
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

    if (verbose) {
      if (i %% 10 == 0) {
        message("*: ", i, " provisional outliers.")
      } else {
        message("*", appendLF = FALSE)
      }
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

  new_outliermbc_lcwm(list(
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
