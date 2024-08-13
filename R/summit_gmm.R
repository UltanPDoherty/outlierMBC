#' @export
summit_gmm <- function(
    x,
    comp_num,
    max_out,
    mnames = "VVV",
    seed = 123,
    reinit_interval = Inf,
    print_interval = Inf) {
  forward <- summit_gmm_forward(
    x, comp_num, max_out, mnames, seed, reinit_interval, print_interval
  )
  backward <- summit_gmm_backward(forward, print_interval)

  return(list(forward = forward, backward = backward))
}

# ------------------------------------------------------------------------------

#' @export
summit_gmm_forward <- function(
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

  rem_dens <- double(max_out)
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
      warning(paste0(
        "One of the components became too small after removing ",
        i - 1, " outliers."
      ))
      break()
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
    rem_dens[i] <- dens_vec[rem_id]

    outlier_rank[!outlier_rank][rem_id] <- i
    x <- x[-rem_id, , drop = FALSE]
    z <- mix$z[-rem_id, , drop = FALSE]
  }

  return(list(
    outlier_rank = outlier_rank,
    rem_dens = rem_dens,
    z = z,
    x0 = x0,
    mnames = mnames
  ))
}

# ------------------------------------------------------------------------------

#' @export
summit_gmm_backward <- function(
    forward,
    print_interval = Inf) {
  outlier_rank <- forward$outlier_rank
  z <- forward$z

  comp_num <- ncol(z)
  max_out <- max(outlier_rank)
  mnames <- forward$mnames
  x <- forward$x

  x0 <- as.matrix(x)
  x <- x0[(outlier_rank == 0), ]

  mix <- mixture::gpcm(x, G = comp_num, mnames = mnames, start = z)

  obs_num <- nrow(x0)
  var_num <- ncol(x0)

  rem_dens <- double(max_out)
  for (i in seq_len(max_out)) {
    if (i %% print_interval == 0) cat("max_out + 1 - i = ", max_out + 1 - i, "\n")

    x <- x0[(outlier_rank == 0) | outlier_rank > max_out - i, ]
    z <- mixture::e_step(x, mix$best_model)$z

    mix <- mixture::gpcm(x, G = comp_num, mnames = mnames, start = z)

    if (any(colSums(mix$z) < var_num + 1)) {
      warning(paste0(
        "One of the components became too small after removing ",
        i - 1, " outliers."
      ))
      break()
    }

    z <- mix$z

    dens_mat <- matrix(nrow = nrow(x), ncol = comp_num)
    for (k in 1:comp_num) {
      dens_mat[, k] <- dmvnorm(
        x,
        mean = mix$best_model$model_obj[[1]]$mu[[k]],
        sigma = mix$best_model$model_obj[[1]]$sigs[[k]]
      )
    }
    dens_vec <- as.numeric(
      dens_mat %*% t(mix$best_model$model_obj[[1]]$pi_gs)
    )

    rem_id <- which.min(dens_vec)
    rem_dens[i] <- dens_vec[rem_id]
  }

  return(list(
    rem_dens = rem_dens
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
