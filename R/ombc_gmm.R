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
#' @param gross_outs logical vector identifying gross outliers to be removed.
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
    print_interval = Inf,
    gross_outs = NULL) {

  if (!is.null(gross_outs)) {
    gross_num <- sum(gross_outs)
    x <- x[!gross_outs, ]
    max_out <- max_out - gross_num
  }

  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x0)
  track_num <- 10

  dist_mat0 <- as.matrix(stats::dist(x0))
  dist_mat <- dist_mat0
  # z <- init_kmpp(x, comp_num, seed)
  z <- init_hc(dist_mat, comp_num)

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
    # mix <- mixture::gpcm(
    #   x,
    #   G = comp_num, mnames = mnames,
    #   start = z, seed = seed
    # )

    # if (i %% reinit_interval == 0) {
      # alt_z <- init_kmpp(x, comp_num, seed)
      alt_z <- init_hc(dist_mat, comp_num)
      alt_mix <- mixture::gpcm(
        x,
        G = comp_num, mnames = mnames,
        start = alt_z, seed = seed
      )

      # if (alt_mix$best_model$loglik > mix$best_model$loglik) {
        # cat(paste0("Iteration ", i, ": reinitialisation accepted.\n"))
        mix <- alt_mix
      # }
    # }

    # if (any(colSums(mix$z) < var_num + 1)) {
    #   message(paste0(
    #     "One of the components became too small after removing ",
    #     i - 1, " outliers.\n"
    #   ))
    #
    #   # alt_z <- init_kmpp(x, comp_num, seed)
    #   alt_z <- init_hc(dist_mat, comp_num)
    #   mix <- mixture::gpcm(x, G = comp_num, mnames = mnames, start = alt_z)
    #
    #   if (any(colSums(mix$z) < var_num + 1)) {
    #     stop("Emergency reinitialisation unsuccessful.\n")
    #   } else {
    #     message("Emergency reinitialisation successful.\n")
    #   }
    # }

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

    dist_mat <- dist_mat[-dd$choice_id, -dd$choice_id]
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

    # # alt_z <- init_kmpp(x0[!outlier_bool[, j], ], comp_num, seed)
    # alt_z <- init_hc(dist_mat0[!outlier_bool[, j], !outlier_bool[, j]], comp_num)
    # alt_mix <- mixture::gpcm(
    #   x0[!outlier_bool[, j], ],
    #   G = comp_num, mnames = mnames,
    #   start = alt_z, seed = seed
    # )
    #
    # if (alt_mix$best_model$loglik > mix[[j]]$best_model$loglik) {
    #   cat(paste0("Final reinitialisation, ", j, " accepted.\n"))
    #   mix[[j]] <- alt_mix
    # }
  }

  labels <- matrix(0, nrow = obs_num, ncol = track_num)
  for (j in 1:track_num) {
    labels[!outlier_bool[, j], j] <- mix[[j]]$map
  }

  if (!is.null(gross_outs)) {
    outlier_bool0 <- outlier_bool
    outlier_rank0 <- outlier_rank
    outlier_num0 <- outlier_num
    labels0 <- labels

    outlier_bool <- matrix(nrow = length(gross_outs), ncol = track_num)
    outlier_rank <- double(length(gross_outs))
    labels <- matrix(nrow = length(gross_outs), ncol = track_num)

    for (j in 1:track_num) {
      outlier_bool[gross_outs, j] <- TRUE
      outlier_bool[!gross_outs, j] <- outlier_bool0[, j]

      outlier_rank[gross_outs] <- 1
      outlier_rank[!gross_outs] <- outlier_rank0 + (outlier_rank0 != 0)

      outlier_num <- outlier_num0 + gross_num

      labels[gross_outs, j] <- 0
      labels[!gross_outs, j] <- labels0[, j]
    }
  }

  p_vals <- round(seq(1, 2, length.out = 10), 2)
  out_nums <- seq(0, max_out)
  if (!is.null(gross_outs)) { out_nums <- out_nums + gross_num }
  gg_curves_list <- list()
  for (j in 1:10) {
    gg_curves_list[[j]] <- as.data.frame(cbind(
      "outlier_number" = out_nums,
      "distrib_diff" = distrib_diff_mat[, j]
    )) |>
      ggplot2::ggplot(ggplot2::aes(x = outlier_number, y = distrib_diff)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = outlier_num[j]) +
      ggplot2::labs(
        title = paste0(
          j, ": ", outlier_num[j]
        ),
        x = "Outlier Number",
        y = "Distibutional Difference"
      ) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y.left = ggplot2::element_blank()
      )
  }
  gg_curves <- ggpubr::ggarrange(plotlist = gg_curves_list, nrow = 2, ncol = 5)

  gg_choice <- as.data.frame(cbind(
    "outlier_number" = out_nums,
    "removal_density" = min_dens
    )) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_number, y = removal_density)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = outlier_num) +
    ggplot2::labs(
      title = paste0(1:10, ": ", outlier_num, collapse = ", "),
      x = "Outlier Number",
      y = "Removal Density"
    )

  return(list(
    distrib_diff_arr = distrib_diff_arr,
    distrib_diff_mat = distrib_diff_mat,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    labels = labels,
    final_gmm = mix,
    loglike = loglike,
    min_dens = min_dens,
    plot_curves = gg_curves,
    plot_choice = gg_choice
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

init_hc <- function(dist_mat, comp_num) {
  hc <- stats::hclust(stats::as.dist(dist_mat), method = "ward.D2")
  init <- stats::cutree(hc, k = comp_num)

  z <- matrix(nrow = nrow(dist_mat), ncol = comp_num)
  for (k in seq_len(comp_num)) {
    z[, k] <- as.integer(init == k)
  }

  return(z)
}

# ------------------------------------------------------------------------------

ombc_gmm_reverse <- function(
    ombc_gmm_out,
    track_choice,
    fixed = FALSE,
    x,
    comp_num,
    mnames = "VVV",
    seed = 123,
    print_interval = Inf) {
  outlier_bool <- ombc_gmm_out$outlier_bool[, track_choice]
  outlier_rank <- ombc_gmm_out$outlier_rank
  outlier_num <- ombc_gmm_out$outlier_num[track_choice]
  mix <- ombc_gmm_out$final_gmm[[track_choice]]
  z <- mix$z

  x0 <- as.matrix(x)
  x <- x0[!outlier_bool, ]

  obs_num <- nrow(x0)
  var_num <- ncol(x0)
  track_num <- 6

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
      z,
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
        min_diff_z[[j]] <- z
      }
    }

    loglike[i] <- mix$best_model$loglik
    min_dens[i] <- dd$min_dens

    x <- x0[!outlier_bool | outlier_rank > outlier_num - i, ]
    z <- mixture::e_step(x, mix$best_model)$z

    if (!fixed) {
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

      z <- mix$z
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
