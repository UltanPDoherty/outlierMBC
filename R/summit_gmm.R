#' @export
outcast_gmm_forward <- function(
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
    removal_density = removal_density,
    outlier_rank = outlier_rank,
    loglike = loglike,
    z = z,
    x0 = x0,
    mnames = mnames
  ))
}

# ------------------------------------------------------------------------------

#' @export
outcast_gmm_backward <- function(
    forward,
    print_interval = Inf
) {
  z <- forward$z

  comp_num <- ncol(z)
  max_out <- max(forward$outlier_rank)
  mnames <- forward$mnames
  x0 <- forward$x0
  obs_num <- nrow(x0)
  var_num <- ncol(x0)

  subset_bool <- forward$outlier_rank == 0
  return_bool <- rep(FALSE, obs_num)

  loglike <- double(max_out)
  replace_density <- double(max_out)
  outlier_rank <- rep(0, obs_num)
  for (i in seq_len(max_out)) {
    if (((max_out + 1 - i) %% print_interval) == 0) {
      cat("max_out + 1 - i = ", max_out + 1 - i, "\n")
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
    outlier_rank[return_bool] <- max_out - i + 1

    subset_bool <- subset_bool | return_bool
    z <- (prop_dens_mat / dens_vec) [subset_bool, ]
  }

  return(list(
    replace_density = rev(replace_density),
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
