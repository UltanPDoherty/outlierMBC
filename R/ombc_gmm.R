#' ombc_gmm
#'
#' @description
#' Iterative Detection & Identification of Outliers for a Gaussian Mixture Model
#'
#' @param x Data.
#' @param comp_num Number of components.
#' @param max_out Maximum number of outliers.
#' @param gross_outs Logical vector identifying gross outliers.
#' @param init_scheme Which initialisation scheme to use.
#' @param mnames Model names for mixture::gpcm.
#' @param nmax Maximum number of iterations for mixture::gpcm.
#' @param atol EM convergence tolerance threshold for mixture::gpcm.
#' @param init_z Initial z matrix.
#' @param init_model mixture::gpcm best_model
#' @param init_method Method used to initialise each mixture model.
#' @param kmpp_seed Optional seed for k-means++ initialisation. Default is
#'                  hierarchical clustering.
#' @param print_interval How frequently the iteration count is printed.
#'
#' @return List of
#' * distrib_diffs
#' * distrib_diff_mat
#' * outlier_bool
#' * outlier_num
#' * outlier_rank
#' * labels
#' * final_gmm
#' * loglike
#' * removal_dens
#' @export
#'
#' @examples
#'
#' ombc_gmm_k3n1000o10 <- ombc_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   comp_num = 3, max_out = 20
#' )
#'
#' ombc_gmm_k3n1000o10$plot_tail_curve
#' ombc_gmm_k3n1000o10$plot_full_curve
#'
ombc_gmm <- function(
    x,
    comp_num,
    max_out,
    gross_outs = rep(FALSE, nrow(x)),
    init_scheme = c("choice", "reinit", "update", "reuse"),
    mnames = "VVV",
    nmax = 1000,
    atol = 1e-8,
    init_z = NULL,
    init_model = NULL,
    init_method = c("hc", "kmpp"),
    kmpp_seed = 123,
    print_interval = Inf) {
  init_method <- match.arg(init_method)
  init_scheme <- match.arg(init_scheme)

  this_call <- call(
    "ombc_gmm",
    "x" = substitute(x), "comp_num" = comp_num, "max_out" = max_out,
    "gross_outs" = substitute(gross_outs), "resets" = substitute(resets),
    "mnames" = mnames, "nmax" = nmax,
    "atol" = atol, "init_method" = init_method, "kmpp_seed" = kmpp_seed,
    "print_interval" = print_interval
  )

  ombc_version <- utils::packageVersion("outlierMBC")

  expect_num <- 1
  accept_num <- expect_num + 1
  reject_num <- accept_num

  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x)

  dist_mat0 <- as.matrix(stats::dist(x0))
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

  track_num <- 2
  tail_props <- expect_num / (seq(obs_num, obs_num - max_out) - gross_num)

  dd_min <- Inf
  best_z <- list()
  tail_accepted <- FALSE

  z_choice <- c()
  conv_status <- c()
  loglike <- c()
  removal_dens <- c()
  distrib_diff_arr <- array(dim = c(comp_num, max_out + 1, track_num))
  distrib_diff_mat <- matrix(nrow = max_out + 1, ncol = track_num)
  outlier_rank_temp <- rep(0, obs_num - gross_num)
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    if (init_scheme %in% c("update", "reuse")) {
      mix <- try_mixture_gpcm(x, comp_num, mnames, z, nmax, atol)
    } else if (init_scheme == "reinit") {
      reinit_z <- get_init_z(
        comp_num = comp_num, dist_mat = dist_mat, x = x,
        init_method = init_method, kmpp_seed = kmpp_seed
      )
      mix <- try_mixture_gpcm(x, comp_num, mnames, reinit_z, nmax, atol)
    } else {
      update_mix <- try_mixture_gpcm(x, comp_num, mnames, z, nmax, atol)

      reinit_z <- get_init_z(
        comp_num = comp_num, dist_mat = dist_mat, x = x,
        init_method = init_method, kmpp_seed = kmpp_seed
      )
      reinit_mix <- try_mixture_gpcm(x, comp_num, mnames, reinit_z, nmax, atol)

      if (reinit_mix$best_model$loglik > update_mix$best_model$loglik) {
        mix <- reinit_mix
        z_choice[i] <- "reinit"
      } else {
        mix <- update_mix
        z_choice[i] <- "update"
      }
    }

    loglike[i] <- mix$best_model$loglik
    conv_status[i] <- mix$best_model$status

    dd <- distrib_diff_gmm(
      x,
      mix$z,
      mix$best_model$model_obj[[1]]$pi_gs,
      mix$best_model$model_obj[[1]]$mu,
      mix$best_model$model_obj[[1]]$sigs,
      mix$best_model$model_obj[[1]]$log_dets,
      tail_props[i]
    )

    distrib_diff_arr[, i, ] <- dd$distrib_diff_mat
    distrib_diff_mat[i, ] <- dd$distrib_diff_vec
    removal_dens[i] <- dd$removal_dens

    outlier_rank_temp[!outlier_rank_temp][dd$choice_id] <- i
    x <- x[-dd$choice_id, , drop = FALSE]
    dist_mat <- dist_mat[-dd$choice_id, -dd$choice_id]

    if (dd$distrib_diff_vec[1] < dd_min) {
      dd_min <- dd$distrib_diff_vec[1]
      best_z[[1]] <- mix$z
    }
    if (!tail_accepted && dd$distrib_diff_vec[2] < accept_num) {
      tail_accepted <- TRUE
      best_z[[2]] <- mix$z
    } else if (tail_accepted && dd$distrib_diff_vec[2] > reject_num) {
      tail_accepted <- FALSE
      best_z[[2]] <- NULL
    }

    if (init_scheme %in% c("update", "choice")) {
      z <- mix$z[-dd$choice_id, , drop = FALSE]
    } else if (init_scheme == "reuse") {
      z <- z[-dd$choice_id, , drop = FALSE]
    }
  }

  outlier_rank <- double(length(gross_outs))
  outlier_rank[gross_outs] <- 1
  outlier_rank[!gross_outs] <- outlier_rank_temp +
    gross_num * (outlier_rank_temp != 0)

  outlier_num <- integer(track_num)

  outlier_num[1] <- which.min(distrib_diff_mat[, 1])

  reject_bool <- distrib_diff_mat[, 2] > reject_num
  final_reject <- max(which(c(TRUE, reject_bool))) - 1
  after_final_reject <- seq_len(max_out + 1) > final_reject
  accept_bools <- distrib_diff_mat[, 2] < accept_num
  outlier_num[2] <- which.max(accept_bools & after_final_reject)

  outlier_num <- outlier_num - 1 + gross_num

  outlier_bool <- matrix(nrow = obs_num, ncol = track_num)
  mix <- list()
  labels <- matrix(0, nrow = obs_num, ncol = track_num)
  for (j in seq_len(track_num)) {
    outlier_bool[, j] <- outlier_rank <= outlier_num[j] & outlier_rank != 0

    mix[[j]] <- try_mixture_gpcm(
      x0[!outlier_bool[, j], ], comp_num, mnames, best_z[[j]], nmax, atol
    )

    labels[!outlier_bool[, j], j] <- mix[[j]]$map
  }

  ombc_names <- c("full", "tail")
  colnames(distrib_diff_mat) <- ombc_names
  colnames(outlier_bool) <- ombc_names
  colnames(labels) <- ombc_names
  names(outlier_num) <- ombc_names
  dimnames(distrib_diff_arr) <- list(
    paste0("k", seq_len(comp_num)), NULL, ombc_names
  )
  names(mix) <- ombc_names

  labels <- as.data.frame(labels)
  outlier_bool <- as.data.frame(outlier_bool)

  quick_tail <- count_extremes(x0[!gross_outs, ], mix$full$best_model)

  return(list(
    labels = labels,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    gross_outs = gross_outs,
    mix = mix,
    loglike = loglike,
    removal_dens = removal_dens,
    distrib_diff_mat = distrib_diff_mat,
    distrib_diff_arr = distrib_diff_arr,
    call = this_call,
    version = ombc_version,
    quick_tail = quick_tail,
    z_choice = z_choice,
    conv_status = conv_status
  ))
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

  return(z)
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

  return(mix)
}
