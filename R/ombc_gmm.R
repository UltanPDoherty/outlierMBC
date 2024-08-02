#' ombc_gmm
#'
#' @description
#' Iterative Detection & Identification of Outliers for a Gaussian Mixture Model
#'
#' @param x Data.
#' @param comp_num Number of components.
#' @param max_out Maximum number of outliers.
#' @param mnames Model names for mixture::gpcm.
#' @param seed Seed.
#' @param reinit_interval How frequently to reinitialise.
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
#' * min_dens
#' @export
#'
#' @examples
#' gmm <- simulate_gmm(
#'   n = c(2000, 1000, 1000),
#'   mu = list(c(-1, 0), c(+1, -1), c(+1, +1)),
#'   sigma = list(diag(c(0.2, 4 * 0.2)), diag(c(0.2, 0.2)), diag(c(0.2, 0.2))),
#'   outlier_num = 40,
#'   seed = 123,
#'   crit_val = 0.9999,
#'   range_multiplier = 1.5
#' )
#'
#' ombc_gmm <- ombc_gmm(gmm[, 1:2], comp_num = 3, max_out = 80)
#'
#' plot(0:80, ombc_gmm$distrib_diffs, type = "l")
#' abline(v = ombc_gmm$outlier_num)
#'
#' plot(gmm[, c("X1", "X2")], col = ombc_gmm$labels + 1, pch = gmm$G + 1)
ombc_gmm <- function(
    x,
    comp_num,
    max_out,
    mnames = "VVV",
    seed = 123,
    reinit_interval = Inf,
    print_interval = Inf) {
  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x0)
  track_num <- 6

  z <- init_kmpp(x, comp_num, seed)

  var_num <- ncol(x)

  min_diff <- rep(Inf, track_num)
  min_diff_z <- list()
  loglike <- c()
  min_dens <- c()
  distrib_diff_arr <- array(dim = c(comp_num, max_out + 1, track_num))
  distrib_diff_mat <- matrix(nrow = max_out + 1, ncol = track_num)
  outlier_rank <- rep(0, obs_num)
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    set.seed(seed)
    mix <- mixture::gpcm(
      x,
      G = comp_num, mnames = mnames,
      start = z, seed = seed
    )

    if (i %% reinit_interval == 0) {
      alt_z <- init_kmpp(x, comp_num, seed)
      alt_mix <- mixture::gpcm(
        x,
        G = comp_num, mnames = mnames,
        start = alt_z, seed = seed
      )

      if (alt_mix$best_model$loglik > mix$best_model$loglik) {
        cat(paste0("Iteration ", i, ": k-means++ reinitialisation accepted.\n"))
        mix <- alt_mix
      }
    }

    if (any(colSums(mix$z) < var_num + 1)) {
      warning(paste0(
        "One of the components became too small after removing ",
        i - 1, " outliers."
      ))
      break()
    }

    dd <- distrib_diff_gmm(
      x,
      mix$z,
      mix$best_model$model_obj[[1]]$pi_gs,
      mix$best_model$model_obj[[1]]$mu,
      mix$best_model$model_obj[[1]]$sigs,
      mix$best_model$model_obj[[1]]$log_dets
    )

    distrib_diff_arr[, i, ] <- dd$distrib_diff_mat
    distrib_diff_mat[i, ] <- dd$distrib_diff_vec

    for (j in 1:track_num) {
      if (distrib_diff_mat[i, j] < min_diff[j]) {
        min_diff[j] <- distrib_diff_mat[i, j]
        min_diff_z[[j]] <- mix$z
      }
    }

    loglike[i] <- mix$best_model$loglik
    min_dens[i] <- dd$min_dens

    outlier_rank[!outlier_rank][dd$choice_id] <- i
    x <- x[-dd$choice_id, , drop = FALSE]
    z <- mix$z[-dd$choice_id, , drop = FALSE]
  }

  outlier_num <- apply(distrib_diff_mat, 2, which.min) - 1

  outlier_bool <- matrix(nrow = obs_num, ncol = track_num)
  for (j in 1:track_num) {
    outlier_bool[, j] <- outlier_rank <= outlier_num[j] & outlier_rank != 0
  }

  mix <- list()
  for (j in 1:track_num) {
    mix[[j]] <- mixture::gpcm(
      x0[!outlier_bool[, j], ],
      G = comp_num, mnames = mnames,
      start = min_diff_z[[j]], seed = seed
    )

    alt_z <- init_kmpp(x0[!outlier_bool[, j], ], comp_num, seed)
    alt_mix <- mixture::gpcm(
      x0[!outlier_bool[, j], ],
      G = comp_num, mnames = mnames,
      start = alt_z, seed = seed
    )

    if (alt_mix$best_model$loglik > mix[[j]]$best_model$loglik) {
      cat(paste0("Final k-means++ reinitialisation, ", j, " accepted.\n"))
      mix[[j]] <- alt_mix
    }
  }

  labels <- matrix(0, nrow = obs_num, ncol = track_num)
  for (j in 1:track_num) {
    labels[!outlier_bool[, j], j] <- mix[[j]]$map
  }

  return(list(
    distrib_diff_arr = distrib_diff_arr,
    distrib_diff_mat = distrib_diff_mat,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    labels = labels,
    final_gmm = mix,
    loglike = loglike,
    min_dens = min_dens
  ))
}

# ------------------------------------------------------------------------------

init_kmpp <- function(x, comp_num, seed) {
  init <- ClusterR::KMeans_rcpp(x, comp_num, 10, seed = seed)$clusters

  z <- matrix(nrow = nrow(x), ncol = comp_num)
  for (k in seq_len(comp_num)) {
    z[, k] <- as.integer(init == k)
  }

  return(z)
}

# ------------------------------------------------------------------------------

ombc_gmm_reverse <- function(
    ombc_gmm_out,
    track_choice,
    x,
    comp_num,
    mnames = "VVV",
    seed = 123,
    print_interval = Inf) {
  outlier_bool <- ombc_gmm_out$outlier_bool[, track_choice]
  outlier_rank <- ombc_gmm_out$outlier_rank
  outlier_num <- ombc_gmm_out$outlier_num[track_choice]
  mix <- ombc_gmm_out$final_gmm[[track_choice]]

  x0 <- as.matrix(x)
  x <- x0[!outlier_bool, ]

  obs_num <- nrow(x0)
  var_num <- ncol(x0)
  track_num <- 4

  min_diff <- rep(Inf, track_num)
  min_diff_z <- list()
  loglike <- c()
  min_dens <- c()
  distrib_diff_arr <- array(dim = c(comp_num, outlier_num + 1, track_num))
  distrib_diff_mat <- matrix(nrow = outlier_num + 1, ncol = track_num)
  for (i in seq_len(outlier_num + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    dd <- distrib_diff_gmm(
      x,
      mix$z,
      mix$best_model$model_obj[[1]]$pi_gs,
      mix$best_model$model_obj[[1]]$mu,
      mix$best_model$model_obj[[1]]$sigs,
      mix$best_model$model_obj[[1]]$log_dets
    )

    distrib_diff_arr[, i, ] <- dd$distrib_diff_mat
    distrib_diff_mat[i, ] <- dd$distrib_diff_vec

    for (j in 1:track_num) {
      if (distrib_diff_mat[i, j] < min_diff[j]) {
        min_diff[j] <- distrib_diff_mat[i, j]
        min_diff_z[[j]] <- mix$z
      }
    }

    loglike[i] <- mix$best_model$loglik
    min_dens[i] <- dd$min_dens

    x <- x0[!outlier_bool | outlier_rank > outlier_num - i, ]
    z <- mixture::e_step(x, mix$best_model)$z

    set.seed(seed)
    mix <- mixture::gpcm(
      x,
      G = comp_num, mnames = mnames,
      start = z, seed = seed
    )

    if (any(colSums(mix$z) < var_num + 1)) {
      warning(paste0(
        "One of the components became too small after removing ",
        i - 1, " outliers."
      ))
      break()
    }
  }

  outlier_num <- outlier_num - (apply(distrib_diff_mat, 2, which.min) - 1)

  outlier_bool <- matrix(nrow = obs_num, ncol = track_num)
  for (j in 1:track_num) {
    outlier_bool[, j] <- outlier_rank <= outlier_num[j] & outlier_rank != 0
  }

  mix <- list()
  for (j in 1:track_num) {
    mix[[j]] <- mixture::gpcm(
      x0[!outlier_bool[, j], ],
      G = comp_num, mnames = mnames,
      start = min_diff_z[[j]], seed = seed
    )

    alt_z <- init_kmpp(x0[!outlier_bool[, j], ], comp_num, seed)
    alt_mix <- mixture::gpcm(
      x0[!outlier_bool[, j], ],
      G = comp_num, mnames = mnames,
      start = alt_z, seed = seed
    )

    if (alt_mix$best_model$loglik > mix[[j]]$best_model$loglik) {
      cat(paste0("Final k-means++ reinitialisation, ", j, " accepted.\n"))
      mix[[j]] <- alt_mix
    }
  }

  labels <- matrix(0, nrow = obs_num, ncol = track_num)
  for (j in 1:track_num) {
    labels[!outlier_bool[, j], j] <- mix[[j]]$map
  }

  return(list(
    distrib_diff_arr = distrib_diff_arr,
    distrib_diff_mat = distrib_diff_mat,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    labels = labels,
    final_gmm = mix,
    loglike = loglike,
    min_dens = min_dens
  ))
}
