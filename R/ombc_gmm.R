#' ombc_gmm
#'
#' @description
#' Iterative Detection & Identification of Outliers for a Gaussian Mixture Model
#'
#' @param x Data.
#' @param comp_num Number of components.
#' @param max_out Maximum number of outliers.
#' @param p_range Range for power mean parameter, p, when summarising CDF
#'                differences.
#' @param mnames Model names for mixture::gpcm.
#' @param nmax Maximum number of iterations for mixture::gpcm.
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
#' ombc1_gmm_k3n1000o10 <- ombc1_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   comp_num = 3, max_out = 20
#' )
#'
#' ombc1_gmm_k3n1000o10$plot_curves
#' ombc1_gmm_k3n1000o10$plot_stacked
#' ombc1_gmm_k3n1000o10$plot_changes
#' ombc1_gmm_k3n1000o10$plot_removal
#'
ombc1_gmm <- function(
    x,
    comp_num,
    max_out,
    p_range = c(1, 2),
    mnames = "VVV",
    nmax = 10,
    print_interval = Inf) {
  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x0)
  track_num <- 10

  dist_mat0 <- as.matrix(stats::dist(x0))
  dist_mat <- dist_mat0

  loglike <- c()
  removal_dens <- c()
  distrib_diff_arr <- array(dim = c(comp_num, max_out + 1, track_num))
  distrib_diff_mat <- matrix(nrow = max_out + 1, ncol = track_num)
  outlier_rank <- rep(0, obs_num)
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    z <- init_hc(dist_mat, comp_num)
    mix <- mixture::gpcm(
      x,
      G = comp_num, mnames = mnames, start = z, nmax = nmax
    )
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
    removal_dens[i] <- dd$removal_dens

    outlier_rank[!outlier_rank][dd$choice_id] <- i
    x <- x[-dd$choice_id, , drop = FALSE]

    dist_mat <- dist_mat[-dd$choice_id, -dd$choice_id]
  }

  outlier_seq <- seq(0, max_out)
  p_vals <- round(seq(p_range[1], p_range[2], length.out = 10), 2)

  gg_curves_list <- list()
  for (j in 1:10) {
    distrib_diff_j <- distrib_diff_mat[, j]
    gg_curves_list[[j]] <- data.frame(outlier_seq, distrib_diff_j) |>
      ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = distrib_diff_j)) +
      ggplot2::geom_line() +
      ggplot2::labs(
        title = paste0(j, " (p = ", p_vals[j], ")"),
        x = "Outlier Number",
        y = "Distributional Difference"
      ) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y.left = ggplot2::element_blank()
      )
  }
  gg_curves <- ggpubr::ggarrange(plotlist = gg_curves_list, nrow = 2, ncol = 5)

  option <- diffs <- NULL
  gg_stacked <- as.data.frame(scale(distrib_diff_mat)) |>
    dplyr::mutate("outlier_seq" = outlier_seq) |>
    tidyr::pivot_longer(
      cols = !outlier_seq, names_to = "option", values_to = "diffs"
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = diffs, group = option)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Outlier_Number", y = "Scaled DD Values",
      title = "Scaled Distributional Differences"
    )

  gg_changes <- as.data.frame(diff(scale(distrib_diff_mat))) |>
    dplyr::mutate("outlier_seq" = outlier_seq[-1]) |>
    tidyr::pivot_longer(
      cols = !outlier_seq, names_to = "option", values_to = "diffs"
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = diffs, group = option)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Outlier_Number", y = "Changes in Scaled DD Values",
      title = "Changes in Scaled Distributional Differences"
    )

  gg_removal <- data.frame(outlier_seq, removal_dens) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = removal_dens)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Outlier Number",
      y = "Removal Density"
    )

  params <- list(
    "comp_num" = comp_num,
    "max_out" = max_out,
    "p_range" = p_range,
    "mnames" = mnames,
    "nmax" = nmax
  )

  return(list(
    distrib_diff_arr = distrib_diff_arr,
    distrib_diff_mat = distrib_diff_mat,
    outlier_rank = outlier_rank,
    loglike = loglike,
    removal_dens = removal_dens,
    plot_curves = gg_curves,
    plot_stacked = gg_stacked,
    plot_changes = gg_changes,
    plot_removal = gg_removal,
    params = params
  ))
}

# ------------------------------------------------------------------------------

#' ombc_gmm
#'
#' @description
#' Iterative Detection & Identification of Outliers for a Gaussian Mixture Model
#'
#' @param ombc1 Output from `ombc1_gmm1`.
#' @param pre_removed_num Number of points to be trimmed, chosen based on ombc1.
#' @inheritParams ombc1_gmm
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
#' ombc1_gmm_k3n1000o10 <- ombc1_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   comp_num = 3, max_out = 20
#' )
#'
#' ombc2_gmm_k3n1000o10 <- ombc2_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   ombc1_gmm_k3n1000o10,
#'   pre_removed_num = 0
#' )
#'
#' ombc2_gmm_k3n1000o10$plot_curves
#' ombc2_gmm_k3n1000o10$plot_removal
#'
ombc2_gmm <- function(
    x,
    ombc1,
    pre_removed_num = 0) {
  comp_num <- ombc1$params$comp_num
  max_out <- ombc1$params$max_out
  p_range <- ombc1$params$p_range
  mnames <- ombc1$params$mnames
  nmax <- ombc1$params$nmax
  outlier_rank <- ombc1$outlier_rank
  removal_dens <- ombc1$removal_dens
  distrib_diff_mat <- ombc1$distrib_diff_mat

  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x0)
  track_num <- 10

  dist_mat0 <- as.matrix(stats::dist(x0))

  if (pre_removed_num > 0) {
    distrib_diff_mat <- distrib_diff_mat[-seq_len(pre_removed_num), ]
  }

  outlier_num <- apply(distrib_diff_mat, 2, which.min) - 1 + pre_removed_num

  outlier_bool <- matrix(nrow = obs_num, ncol = track_num)
  mix <- list()
  labels <- matrix(0, nrow = obs_num, ncol = track_num)
  for (j in 1:track_num) {
    outlier_bool[, j] <- outlier_rank <= outlier_num[j] & outlier_rank != 0

    z <- init_hc(dist_mat0[!outlier_bool[, j], !outlier_bool[, j]], comp_num)
    mix[[j]] <- mixture::gpcm(
      x0[!outlier_bool[, j], ],
      G = comp_num, mnames = mnames, start = z, nmax = nmax
    )

    labels[!outlier_bool[, j], j] <- mix[[j]]$map
  }

  outlier_seq <- seq(0, max_out)
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
        y = "Distributional Difference"
      ) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y.left = ggplot2::element_blank()
      )
  }
  gg_curves <- ggpubr::ggarrange(plotlist = gg_curves_list, nrow = 2, ncol = 5)

  gg_removal <- data.frame(outlier_seq, removal_dens) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = removal_dens)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = outlier_num) +
    ggplot2::labs(
      title = paste0(1:10, ": ", outlier_num, collapse = ", "),
      x = "Outlier Number",
      y = "Removal Density"
    )

  return(list(
    distrib_diff_mat = distrib_diff_mat,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    labels = labels,
    plot_curves = gg_curves,
    plot_removal = gg_removal
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
