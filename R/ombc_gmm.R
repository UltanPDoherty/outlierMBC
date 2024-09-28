#' ombc_gmm
#'
#' @description
#' Iterative Detection & Identification of Outliers for a Gaussian Mixture Model
#'
#' @param x Data.
#' @param comp_num Number of components.
#' @param max_out Maximum number of outliers.
#' @param gross_outs Logical vector identifying gross outliers.
#' @param tail_probs Values for power mean parameter, p, when summarising CDF
#'               differences.
#' @param target_threshold The accepted number of outliers must have a tail
#'                         proportion difference less than this level.
#' @param reset_threshold If the tail proportion difference passes above this
#'                        level at some point, then no outlier number before
#'                        that point can be accepted.
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
#' ombc_gmm_k3n1000o10 <- ombc_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   comp_num = 3, max_out = 20
#' )
#'
#' ombc_gmm_k3n1000o10$plot_curves
#' ombc_gmm_k3n1000o10$plot_removal
#'
ombc_gmm <- function(
    x,
    comp_num,
    max_out,
    gross_outs = NULL,
    tail_probs = c(0.999, 0.9999),
    target_threshold = 1 - tail_probs,
    reset_threshold = 2 * target_threshold,
    mnames = "VVV",
    nmax = 10,
    print_interval = Inf) {
  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x)

  dist_mat0 <- as.matrix(stats::dist(x0))
  dist_mat <- dist_mat0

  if (!is.null(gross_outs)) {
    gross_num <- sum(gross_outs)
    x <- x[!gross_outs, ]
    max_out <- max_out - gross_num
    dist_mat <- dist_mat[!gross_outs, !gross_outs]
  } else {
    gross_num <- 0
  }

  track_num <- length(tail_probs)

  loglike <- c()
  removal_dens <- c()
  distrib_diff_arr <- array(dim = c(comp_num, max_out + 1, track_num))
  distrib_diff_mat <- matrix(nrow = max_out + 1, ncol = track_num)
  outlier_rank <- rep(0, obs_num - gross_num)
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    z <- init_hc(dist_mat, comp_num)
    mix <- try_mixture_gpcm(x, comp_num, mnames, z, nmax)

    loglike[i] <- mix$best_model$loglik

    dd <- distrib_diff_gmm(
      x,
      mix$z,
      mix$best_model$model_obj[[1]]$pi_gs,
      mix$best_model$model_obj[[1]]$mu,
      mix$best_model$model_obj[[1]]$sigs,
      mix$best_model$model_obj[[1]]$log_dets,
      tail_probs
    )

    distrib_diff_arr[, i, ] <- dd$distrib_diff_mat
    distrib_diff_mat[i, ] <- dd$distrib_diff_vec
    removal_dens[i] <- dd$removal_dens

    outlier_rank[!outlier_rank][dd$choice_id] <- i
    x <- x[-dd$choice_id, , drop = FALSE]
    dist_mat <- dist_mat[-dd$choice_id, -dd$choice_id]
  }

  if (!is.null(gross_outs)) {
    outlier_rank0 <- outlier_rank

    outlier_rank <- double(length(gross_outs))

    outlier_rank[gross_outs] <- 1
    outlier_rank[!gross_outs] <- outlier_rank0 +
      gross_num * (outlier_rank0 != 0)
  }

  outlier_num <- integer(track_num)
  last_reset <- integer(track_num)
  target_bools <- matrix(nrow = max_out + 1, ncol = track_num)
  last_reset_bools <- matrix(nrow = max_out + 1, ncol = track_num)
  for (j in seq_len(track_num)) {
    reset_bool <- distrib_diff_mat[, j] > reset_threshold[j]
    if (any(reset_bool)) {
      last_reset[j] <- max(which(reset_bool))
    }
    target_bools[, j] <- distrib_diff_mat[, j] < target_threshold[j]
    last_reset_bools[, j] <- seq_len(max_out + 1) > last_reset[j]
    outlier_num[j] <- which.max(target_bools[, j] & last_reset_bools[, j])
  }

  consensus_vec <- apply(
    target_bools & last_reset_bools, 1, all
  )
  if (any(consensus_vec)) {
    track_num_plus <- track_num + 1
    outlier_num[track_num_plus] <- which.max(consensus_vec)
  } else {
    cat(paste0("No consensus achieved.\n"))
  }

  outlier_num <- outlier_num - 1 + gross_num

  outlier_bool <- matrix(nrow = obs_num, ncol = track_num_plus)
  mix <- list()
  labels <- matrix(0, nrow = obs_num, ncol = track_num_plus)
  for (j in seq_len(track_num_plus)) {
    outlier_bool[, j] <- outlier_rank <= outlier_num[j] & outlier_rank != 0

    z <- init_hc(dist_mat0[!outlier_bool[, j], !outlier_bool[, j]], comp_num)

    mix[[j]] <- try_mixture_gpcm(
      x0[!outlier_bool[, j], ], comp_num, mnames, z, nmax
    )

    labels[!outlier_bool[, j], j] <- mix[[j]]$map
  }

  outlier_seq <- seq(gross_num, max_out + gross_num)

  gg_curves_list <- list()
  reset <- target <- choice <- consensus <- NULL
  point_size <- 1 / ceiling(max_out / 50)
  for (j in seq_len(track_num)) {
    lines_j <- data.frame(
      "reset" = reset_threshold[j], "target" = target_threshold[j],
      "choice" = outlier_num[j], "consensus" = outlier_num[track_num_plus]
    )
    distrib_diff_j <- distrib_diff_mat[, j]
    gg_curves_list[[j]] <- data.frame(outlier_seq, distrib_diff_j) |>
      ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = distrib_diff_j)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = point_size) +
      ggplot2::geom_vline(
        data = lines_j,
        ggplot2::aes(xintercept = choice, colour = "choice"),
        linetype = "solid", linewidth = 0.75
      ) +
      ggplot2::geom_vline(
        data = lines_j,
        ggplot2::aes(xintercept = consensus, colour = "consensus"),
        linetype = "solid", linewidth = 0.75
      ) +
      ggplot2::geom_hline(
        data = lines_j,
        ggplot2::aes(yintercept = reset, colour = "reset"),
        linetype = "dashed", linewidth = 0.75
      ) +
      ggplot2::geom_hline(
        data = lines_j,
        ggplot2::aes(yintercept = target, colour = "target"),
        linetype = "dashed", linewidth = 0.75
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          reset = "#D55E00", target = "#009E73",
          choice = "#CC79A7", consensus = "#0072B2"
        )
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
      ggplot2::labs(
        title = paste0(
          j, ": No. of Outliers = ", outlier_num[j],
          " (tail = ", tail_probs[j], ")"
        ),
        x = "Outlier Number",
        y = "Tail Proportion Difference",
        colour = "Thresholds:"
      ) +
      ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
      ggplot2::theme(
        legend.text = ggplot2::element_text(size = 11),
        legend.title = ggplot2::element_text(size = 11)
      )
  }
  gg_curves <- ggpubr::ggarrange(
    plotlist = gg_curves_list,
    nrow = min(2, track_num), ncol = ceiling(track_num / 2),
    common.legend = TRUE, legend = "bottom"
  )

  tp_names <- paste0("tp", tail_probs)
  colnames(distrib_diff_mat) <- tp_names
  colnames(outlier_bool) <- c(tp_names, "consensus")
  colnames(labels) <- c(tp_names, "consensus")
  names(outlier_num) <- c(tp_names, "consensus")
  dimnames(distrib_diff_arr) <- list(
    paste0("k", seq_len(comp_num)), NULL, tp_names
  )

  return(list(
    distrib_diff_mat = distrib_diff_mat,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    labels = labels,
    plot_curves = gg_curves,
    loglike = loglike,
    removal_dens = removal_dens,
    distrib_diff_arr = distrib_diff_arr
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

# ------------------------------------------------------------------------------

try_mixture_gpcm <- function(x, comp_num, mnames, z, nmax) {
  mix <- NULL
  try(mix <- mixture::gpcm(
    x,
    G = comp_num, mnames = mnames, start = z, nmax = nmax
  ))
  if (is.null(mix)) {
    try(mix <- mixture::gpcm(
      x,
      G = comp_num, start = z, nmax = nmax,
      mnames = setdiff(
        mnames,
        c(
          "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE",
          "EEV", "VEV", "VVV", "EVE", "VVE", "VEE", "EVV"
        )
      )
    ))
    if (!is.null(mix)) {
      cat(paste0(
        mix$best_model$cov_type, " covariance structure implemented.\n"
      ))
    }
  }

  return(mix)
}

# ------------------------------------------------------------------------------

two_set_dist <- function(x, y) {
  x_num <- nrow(x)
  y_num <- nrow(y)

  dist_mat <- matrix(nrow = x_num, ncol = y_num)
  for (i in seq_len(x_num)) {
    for (j in seq_len(y_num)) {
      dist_mat[i, j] <- as.numeric(dist(rbind(x[i, ], y[j, ])))
    }
  }

  lsap <- clue::solve_LSAP(dist_mat)

  mean_dist <- mean(dist_mat[cbind(seq_along(lsap), lsap)])

  return(mean_dist)
}
