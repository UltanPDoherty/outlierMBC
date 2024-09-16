#' ombc_gmm
#'
#' @description
#' Iterative Detection & Identification of Outliers for a Gaussian Mixture Model
#'
#' @param x Data.
#' @param comp_num Number of components.
#' @param max_out Maximum number of outliers.
#' @param gross_outs logical vector identifying gross outliers to be removed.
#' @param p_range Range for power mean parameter, p, when summarising CDF
#'                differences.
#' @param mnames Model names for mixture::gpcm.
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
#'
#' ombc_gmm_k3n1000o10 <- ombc_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   comp_num = 3, max_out = 20
#' )
#'
#' ombc_gmm_k3n1000o10$plot_curves
#'
#' ombc_gmm_k3n1000o10$plot_choice
#'
ombc_gmm <- function(
    x,
    comp_num,
    max_out,
    gross_outs = NULL,
    p_range = c(1, 2),
    mnames = "VVV",
    print_interval = Inf) {
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
  z <- init_hc(dist_mat, comp_num)

  loglike <- c()
  min_dens <- c()
  distrib_diff_arr <- array(dim = c(comp_num, max_out + 1, track_num))
  distrib_diff_mat <- matrix(nrow = max_out + 1, ncol = track_num)
  outlier_rank <- rep(0, obs_num)
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    z <- init_hc(dist_mat, comp_num)
    mix <- mixture::gpcm(x, G = comp_num, mnames = mnames, start = z)

    if (i > 1 && mix$best_model$loglik < loglike[i - 1]) {
      alt_mix <- mixture::gpcm(x, G = comp_num, mnames = mnames, start = prev_z)

      if (alt_mix$best_model$loglik > mix$best_model$loglik) {
        cat("Previous z matrix carried forward at i = ", i, ".\n")
        mix <- alt_mix
      }
    }
    loglike[i] <- mix$best_model$loglik

    dd <- distrib_diff_gmm(
      x,
      mix$z,
      mix$best_model$model_obj[[1]]$pi_gs,
      mix$best_model$model_obj[[1]]$mu,
      mix$best_model$model_obj[[1]]$sigs,
      mix$best_model$model_obj[[1]]$log_dets,
      p_range
    )

    distrib_diff_arr[, i, ] <- dd$distrib_diff_mat
    distrib_diff_mat[i, ] <- dd$distrib_diff_vec
    min_dens[i] <- dd$min_dens

    outlier_rank[!outlier_rank][dd$choice_id] <- i
    x <- x[-dd$choice_id, , drop = FALSE]
    prev_z <- mix$z[-dd$choice_id, , drop = FALSE]

    dist_mat <- dist_mat[-dd$choice_id, -dd$choice_id]
  }

  outlier_num <- apply(distrib_diff_mat, 2, which.min) - 1

  outlier_bool <- matrix(nrow = obs_num, ncol = track_num)
  mix <- list()
  labels <- matrix(0, nrow = obs_num, ncol = track_num)
  for (j in 1:track_num) {
    outlier_bool[, j] <- outlier_rank <= outlier_num[j] & outlier_rank != 0

    z <- init_hc(dist_mat0[!outlier_bool[, j], !outlier_bool[, j]], comp_num)
    mix[[j]] <- mixture::gpcm(
      x0[!outlier_bool[, j], ],
      G = comp_num, mnames = mnames, start = z
    )

    labels[!outlier_bool[, j], j] <- mix[[j]]$map
  }

  outlier_seq <- seq(0, max_out)

  if (!is.null(gross_outs)) {
    outlier_bool0 <- outlier_bool
    outlier_rank0 <- outlier_rank
    outlier_num0 <- outlier_num
    labels0 <- labels

    outlier_bool <- matrix(nrow = length(gross_outs), ncol = track_num)
    outlier_rank <- double(length(gross_outs))
    labels <- matrix(nrow = length(gross_outs), ncol = track_num)

    outlier_rank[gross_outs] <- 1
    outlier_rank[!gross_outs] <- outlier_rank0 + (outlier_rank0 != 0)

    outlier_num <- outlier_num0 + gross_num

    outlier_bool[gross_outs, ] <- TRUE
    outlier_bool[!gross_outs, ] <- outlier_bool0

    labels[gross_outs, ] <- 0
    labels[!gross_outs, ] <- labels0

    outlier_seq <- outlier_seq + gross_num
  }

  p_vals <- round(seq(p_range[1], p_range[2], length.out = 10), 2)

  gg_curves_list <- list()
  for (j in 1:10) {
    distrib_diff_j <- distrib_diff_mat[, j]
    gg_curves_list[[j]] <- data.frame(outlier_seq, distrib_diff_j) |>
      ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = distrib_diff_j)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = outlier_num[j]) +
      ggplot2::labs(
        title = paste0(j, ": ", outlier_num[j], " (p = ", p_vals[j], ")"),
        x = "Outlier Number",
        y = "Distibutional Difference"
      ) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y.left = ggplot2::element_blank()
      )
  }
  gg_curves <- ggpubr::ggarrange(plotlist = gg_curves_list, nrow = 2, ncol = 5)

  gg_choice <- data.frame(outlier_seq, min_dens) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = min_dens)) +
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

init_hc <- function(dist_mat, comp_num) {
  hc <- stats::hclust(stats::as.dist(dist_mat), method = "ward.D2")
  init <- stats::cutree(hc, k = comp_num)

  z <- matrix(nrow = nrow(dist_mat), ncol = comp_num)
  for (k in seq_len(comp_num)) {
    z[, k] <- as.integer(init == k)
  }

  return(z)
}
