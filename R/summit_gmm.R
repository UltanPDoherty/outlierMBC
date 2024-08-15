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

#' @export
outback_gmm <- function(
    x,
    outlier_rank,
    comp_num,
    turning_point = NULL,
    mnames = "VVV",
    seed = 123,
    print_interval = Inf
) {
  x0 <- as.matrix(x)
  obs_num <- nrow(x0)
  var_num <- ncol(x0)

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
    z <- (prop_dens_mat / dens_vec) [subset_bool, ]
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
# ------------------------------------------------------------------------------

#' @export
use_cpop <- function(y, search_centre) {
  y_len <- length(y)

  search_radius <- min(c(
    floor((y_len - search_centre - 1) / 3),
    floor((y_len - 1) / 3)
  ))

  upper <- search_centre + search_radius
  lower <- search_centre - search_radius
  search_interval <- c(lower, upper)

  sd1 <- sqrt(mean(diff(diff(y[1:search_centre]))^2)/6)
  sd2 <- sqrt(mean(diff(diff(y[(search_centre + 1):y_len]))^2)/6)

  cpop_out <- suppressMessages(cpop::cpop(
    y, seq_len(y_len),
    grid = lower:upper, minseglen = 2 * search_radius + 1,
    sd = c(rep(sd1, search_centre), rep(sd2, y_len - search_centre))
  ))

  stopifnot(length(cpop_out@changepoints) == 3)

  cat("search centre = ", search_centre)
  cat(", search radius = ", search_radius)
  cat(", search interval = ", search_interval)
  cat(", search choice = ", cpop_out@changepoints[2], "\n")

  return(list(
    choice = cpop_out@changepoints[2],
    search_interval = search_interval,
    cpop_out = cpop_out
