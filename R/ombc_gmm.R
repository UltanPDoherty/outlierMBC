#' Sequentially identify outliers while fitting a Gaussian mixture model.
#'
#' @description
#' This function performs model-based clustering and outlier identification. It
#' does so by iteratively fitting a Gaussian mixture model and removing the
#' observation that is least likely under the model. Its procedure is summarised
#' below:
#'
#' 1. Fit a Gaussian mixture model to the data.
#' 2. Compute a dissimilarity between the theoretical and observed distributions
#'   of the scaled squared sample Mahalanobis distances for each mixture
#'   component.
#' 3. Aggregate across the components to obtain a single dissimilarity value.
#' 4. Remove the observation  with the lowest mixture density.
#' 5. Repeat Steps 1-4 until `max_out` observations have been removed.
#' 6. Identify the number of outliers which minimised the aggregated
#'    dissimilarity, remove only those observations, and fit a Gaussian mixture
#'    model to the remaining data.
#'
#' @param x Data.
#' @param comp_num Number of mixture components.
#' @param max_out Maximum number of outliers.
#' @param gross_outs Logical vector identifying gross outliers.
#' @param init_scheme Which initialisation scheme to use.
#' @param mnames Model names for mixture::gpcm.
#' @param nmax Maximum number of iterations for mixture::gpcm.
#' @param atol EM convergence tolerance threshold for mixture::gpcm.
#' @param init_z Initial component assignment probability matrix.
#' @param init_model Initial mixture model (mixture::gpcm best_model).
#' @param init_method Method used to initialise each mixture model.
#' @param init_scaling Logical value controlling whether the data should be
#'                     scaled for initialisation.
#' @param kmpp_seed Optional seed for k-means++ initialisation.
#' @param print_interval How frequently the iteration count is printed.
#'
#' @returns
#' `ombc_gmm` returns a list with the following elements:
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
#'                         numbered from gross_num + 1` to `max_out`.
#'                         Observations that were ever removed have rank `0`.}
#'   \item{`gross_outs`}{Logical vector identifying the gross outliers. This is
#'                       identical to the `gross_outs` vector passed to this
#'                       function as an argument / input.}
#'   \item{`mix`}{Output from mixture::gpcm fitted to the non-outlier
#'                observations.}
#'   \item{`loglike`}{Vector of log-likelihood values for each iteration.}
#'   \item{`removal_dens`}{Vector of mixture densities for the removed
#'                         observations. These are the lowest mixture densities
#'                         at each iteration.}
#'   \item{`distrib_diff_vec`}{Vector of aggregated cross-component
#'                             dissimilarity values for each iteration.}
#'   \item{`distrib_diff_mat`}{Matrix of component-specific dissimilarity values
#'                             for each iteration.}
#'   \item{`call`}{Arguments / parameter values used in this function call.}
#'   \item{`version`}{Version of `outlierMBC` used in this function call.}
#'   \item{`conv_status`}{Logical vector indicating which iterations' mixture
#'                        models reached convergence during model-fitting.}
#' }
#'
#' @export
#'
#' @examples
#'
#' ombc_gmm_k3n1000o10 <- ombc_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   comp_num = 3, max_out = 20
#' )
#'
#' plot_curve(ombc_gmm_k3n1000o10)
#'
ombc_gmm <- function(
    x,
    comp_num,
    max_out,
    gross_outs = rep(FALSE, nrow(x)),
    init_scheme = c("update", "reinit", "reuse"),
    mnames = "VVV",
    nmax = 1000,
    atol = 1e-8,
    init_z = NULL,
    init_model = NULL,
    init_method = c("hc", "kmpp"),
    init_scaling = FALSE,
    kmpp_seed = 123,
    print_interval = Inf) {
  init_method <- match.arg(init_method)
  init_scheme <- match.arg(init_scheme)

  this_call <- call(
    "ombc_gmm",
    "x" = substitute(x), "comp_num" = comp_num, "max_out" = max_out,
    "gross_outs" = substitute(gross_outs), "init_scheme" = init_scheme,
    "mnames" = mnames, "nmax" = nmax, "atol" = atol,
    "init_z" = substitute(init_z), "init_model" = substitute(init_model),
    "init_method" = init_method, "init_scaling" = init_scaling,
    "kmpp_seed" = kmpp_seed, "print_interval" = print_interval
  )

  ombc_version <- utils::packageVersion("outlierMBC")

  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x)

  x1 <- scale(x0, center = init_scaling, scale = init_scaling)
  dist_mat0 <- as.matrix(stats::dist(x1))
  dist_mat <- dist_mat0

  gross_num <- sum(gross_outs)
  x <- x[!gross_outs, ]
  max_out <- max_out - gross_num
  dist_mat <- dist_mat[!gross_outs, !gross_outs]

  if (!is.null(init_model) && !is.null(init_z)) {
    stop("Only one of init_model and init_z may be provided.")
  } else if (!is.null(init_z)) {
    z <- init_z
  } else if (!is.null(init_model)) {
    z <- mixture::e_step(x, init_model)$z
  } else {
    z <- get_init_z(
      comp_num = comp_num, dist_mat = dist_mat, x = x,
      init_method = init_method, kmpp_seed = kmpp_seed
    )
  }

  dd_min <- Inf

  conv_status <- c()
  loglike <- c()
  removal_dens <- c()
  distrib_diff_mat <- matrix(nrow = max_out + 1, ncol = comp_num)
  distrib_diff_vec <- double(max_out + 1)
  outlier_rank_temp <- rep(0, obs_num - gross_num)
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    if (init_scheme %in% c("update", "reuse")) {
      mix <- try_mixture_gpcm(x, comp_num, mnames, z, nmax, atol)
    } else {
      reinit_z <- get_init_z(
        comp_num = comp_num, dist_mat = dist_mat, x = x,
        init_method = init_method, kmpp_seed = kmpp_seed
      )
      mix <- try_mixture_gpcm(x, comp_num, mnames, reinit_z, nmax, atol)
    }

    loglike[i] <- mix$best_model$loglik
    conv_status[i] <- mix$best_model$status

    dd <- distrib_diff_gmm(
      x,
      mix$z,
      mix$best_model$model_obj[[1]]$pi_gs,
      mix$best_model$model_obj[[1]]$mu,
      mix$best_model$model_obj[[1]]$sigs,
      mix$best_model$model_obj[[1]]$log_dets
    )

    distrib_diff_mat[i, ] <- dd$distrib_diff_vec
    distrib_diff_vec[i] <- dd$distrib_diff
    removal_dens[i] <- dd$removal_dens

    outlier_rank_temp[!outlier_rank_temp][dd$choice_id] <- i
    x <- x[-dd$choice_id, , drop = FALSE]
    dist_mat <- dist_mat[-dd$choice_id, -dd$choice_id]

    if (dd$distrib_diff < dd_min) {
      dd_min <- dd$distrib_diff
      best_z <- mix$z
    }

    if (init_scheme %in% c("update")) {
      z <- mix$z[-dd$choice_id, , drop = FALSE]
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

  mix <- try_mixture_gpcm(
    x0[!outlier_bool, ], comp_num, mnames, best_z, nmax, atol
  )

  labels[!outlier_bool] <- mix$map

  colnames(distrib_diff_mat) <- paste0("k", seq_len(comp_num))

  list(
    labels = labels,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    gross_outs = gross_outs,
    mix = mix,
    loglike = loglike,
    removal_dens = removal_dens,
    distrib_diff_vec = distrib_diff_vec,
    distrib_diff_mat = distrib_diff_mat,
    call = this_call,
    version = ombc_version,
    conv_status = conv_status
  )
}

# ------------------------------------------------------------------------------

get_init_z <- function(
    comp_num,
    dist_mat = NULL, x = NULL,
    init_method = c("hc", "kmpp"), kmpp_seed = NULL) {
  init_method <- match.arg(init_method)

  if (init_method == "hc") {
    hc <- stats::hclust(stats::as.dist(dist_mat), method = "ward.D2")
    init <- stats::cutree(hc, k = comp_num)
  } else if (init_method == "kmpp") {
    kmpp <- ClusterR::KMeans_rcpp(
      x,
      clusters = comp_num, num_init = 10, seed = kmpp_seed
    )
    init <- kmpp$clusters
  }

  z <- matrix(nrow = nrow(dist_mat), ncol = comp_num)
  for (k in seq_len(comp_num)) {
    z[, k] <- as.integer(init == k)
  }

  z
}

# ------------------------------------------------------------------------------

try_mixture_gpcm <- function(x, comp_num, mnames, z, nmax, atol) {
  mix <- NULL
  try(mix <- mixture::gpcm(
    x,
    G = comp_num, mnames = mnames, start = z, nmax = nmax, atol = atol
  ))
  if (is.null(mix)) {
    cat(paste0("Trying alternative covariance structures.\n"))
    try(mix <- mixture::gpcm(
      x,
      G = comp_num, start = z, nmax = nmax, atol = atol,
      mnames = setdiff(
        c(
          "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE",
          "EEV", "VEV", "VVV", "EVE", "VVE", "VEE", "EVV"
        ),
        mnames
      )
    ))
    if (!is.null(mix)) {
      cat(paste0(
        mix$best_model$cov_type, " covariance structure implemented.\n"
      ))
    } else {
      cat(paste0("Trying default gpcm k-means initialisation.\n"))
      try(
        mix <- mixture::gpcm(
          x,
          G = comp_num, start = 2, nmax = nmax, atol = atol
        )
      )
      if (!is.null(mix)) {
        cat(paste0(
          mix$best_model$cov_type, " covariance structure implemented.\n"
        ))
      }
    }
  }

  mix
}
