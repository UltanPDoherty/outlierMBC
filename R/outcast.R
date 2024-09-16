#' @title Alternate between trimming the lowest density point and refitting a
#' Gaussian mixture model.
#'
#' @description
#' Record the order and mixture density of the trimmed points.
#'
#' @param x Data set.
#' @param comp_num Number of mixture components.
#' @param max_out Maximum number of outliers.
#' @param mnames Gaussian mixture model covariance structure.
#' @param seed Seed.
#' @param reinit_interval Number of iterations between proposed
#'                        reinitialisations.
#' @param print_interval Number of iterations between print statements.
#'
#' @return List:
#' * $densities
#' * $outlier_rank
#' * $loglike
#'
#' @export
outcast_gmm <- function(
    x,
    comp_num,
    max_out,
    mnames = "VVV",
    seed = 123,
    reinit_interval = Inf,
    print_interval = Inf) {
  x0 <- as.matrix(x)
  x <- x0

  obs_num <- nrow(x)
  var_num <- ncol(x)

  z <- init_kmpp(x, comp_num, seed)

  removal_density <- double(max_out)
  loglike <- double(max_out)
  outlier_rank <- rep(0, obs_num)
  for (i in seq_len(max_out)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    mix <- mixture::gpcm(x, G = comp_num, mnames = mnames, start = z)

    if (i %% reinit_interval == 0) {
      alt_z <- init_kmpp(x, comp_num, seed)
      alt_mix <- mixture::gpcm(x, G = comp_num, mnames = mnames, start = alt_z)

      if (alt_mix$best_model$loglik > mix$best_model$loglik) {
        cat(paste0("Iteration ", i, ": k-means++ reinitialisation accepted.\n"))
        mix <- alt_mix
      }
    }

    if (any(colSums(mix$z) < var_num + 1)) {
      message(paste0(
        "One of the components became too small after removing ",
        i - 1, " outliers.\n"
      ))

      alt_z <- init_kmpp(x, comp_num, seed)
      mix <- mixture::gpcm(x, G = comp_num, mnames = mnames, start = alt_z)

      if (any(colSums(mix$z) < var_num + 1)) {
        stop("Emergency reinitialisation unsuccessful.\n")
      } else {
        message("Emergency reinitialisation successful.\n")
      }
    }

    dens_mat <- matrix(nrow = nrow(x), ncol = comp_num)
    for (k in 1:comp_num) {
      dens_mat[, k] <- mvtnorm::dmvnorm(
        x,
        mean = mix$best_model$model_obj[[1]]$mu[[k]],
        sigma = mix$best_model$model_obj[[1]]$sigs[[k]]
      )
    }
    dens_vec <- as.numeric(
      dens_mat %*% t(mix$best_model$model_obj[[1]]$pi_gs)
    )

    rem_id <- which.min(dens_vec)
    removal_density[i] <- dens_vec[rem_id]
    loglike[i] <- mix$best_model$loglik

    outlier_rank[!outlier_rank][rem_id] <- i
    x <- x[-rem_id, , drop = FALSE]
    z <- mix$z[-rem_id, , drop = FALSE]
  }

  return(list(
    densities = removal_density,
    outlier_rank = outlier_rank,
    loglike = loglike
  ))
}

# ------------------------------------------------------------------------------

#' @title Alternate between returning the highest density trimmed point and
#' refitting a Gaussian mixture model.
#'
#' @inheritParams outcast_gmm
#' @param outlier_rank Order in which points were trimmed / removed.
#' @param turning_point Point at which to start returning points.
#'
#' @return List:
#' * $densities
#' * $outlier_rank
#' * $loglike
#'
#' @export
outback_gmm <- function(
    x,
    outlier_rank,
    comp_num,
    turning_point = NULL,
    mnames = "VVV",
    seed = 123,
    print_interval = Inf) {
  x0 <- as.matrix(x)
  obs_num <- nrow(x0)

  if (is.null(turning_point)) {
    turning_point <- max(outlier_rank)
  }

  subset_bool <- outlier_rank == 0 | outlier_rank > turning_point
  return_bool <- rep(FALSE, obs_num)

  z <- init_kmpp(x0[subset_bool, ], comp_num, seed)

  loglike <- double(turning_point)
  replace_density <- double(turning_point)
  outlier_rank <- rep(0, obs_num)
  for (i in seq_len(turning_point)) {
    if (((turning_point + 1 - i) %% print_interval) == 0) {
      cat("turning_point + 1 - i = ", turning_point + 1 - i, "\n")
    }

    mix <- mixture::gpcm(
      x0[subset_bool, ],
      G = comp_num,
      mnames = mnames,
      start = z
    )
    loglike[i] <- mix$best_model$loglik

    dens_mat <- matrix(nrow = obs_num, ncol = comp_num)
    for (k in 1:comp_num) {
      dens_mat[, k] <- dmvnorm(
        x0,
        mean = mix$best_model$model_obj[[1]]$mu[[k]],
        sigma = mix$best_model$model_obj[[1]]$sigs[[k]]
      )
    }
    prop_dens_mat <- sweep(
      dens_mat, 2, mix$best_model$model_obj[[1]]$pi_gs, "*"
    )
    dens_vec <- rowSums(prop_dens_mat)

    return_bool <- rep(FALSE, obs_num)
    return_bool[!subset_bool] <-
      seq_len(sum(!subset_bool)) == which.max(dens_vec[!subset_bool])
    replace_density[i] <- dens_vec[return_bool]
    outlier_rank[return_bool] <- turning_point - i + 1

    subset_bool <- subset_bool | return_bool
    z <- (prop_dens_mat / dens_vec)[subset_bool, ]
  }

  return(list(
    densities = rev(replace_density),
    outlier_rank = outlier_rank,
    loglike = rev(loglike)
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
