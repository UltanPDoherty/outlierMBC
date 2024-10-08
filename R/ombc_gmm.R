#' ombc_gmm
#'
#' @description
#' Iterative Detection & Identification of Outliers for a Gaussian Mixture Model
#'
#' @param x Data.
#' @param comp_num Number of components.
#' @param max_out Maximum number of outliers.
#' @param gross_outs Logical vector identifying gross outliers.
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
    gross_outs = rep(FALSE, nrow(x)),
    mnames = "VVV",
    nmax = 10,
    print_interval = Inf) {
  expect_num <- 1
  accept_num <- expect_num + 1
  reject_num <- accept_num + 2

  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x)

  dist_mat0 <- as.matrix(stats::dist(x0))
  dist_mat <- dist_mat0

  gross_num <- sum(gross_outs)
  x <- x[!gross_outs, ]
  max_out <- max_out - gross_num
  dist_mat <- dist_mat[!gross_outs, !gross_outs]

  track_num <- length(expect_num)
  tail_props <- outer(
    1 / seq(obs_num - gross_num, obs_num - gross_num - max_out),
    expect_num
  )
  en_names <- paste0("en", expect_num)

  loglike <- c()
  removal_dens <- c()
  distrib_diff_arr <- array(dim = c(comp_num, max_out + 1, track_num))
  distrib_diff_mat <- matrix(nrow = max_out + 1, ncol = track_num)
  outlier_rank_temp <- rep(0, obs_num - gross_num)
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
      tail_props[i, ]
    )

    distrib_diff_arr[, i, ] <- dd$distrib_diff_mat
    distrib_diff_mat[i, ] <- dd$distrib_diff_vec
    removal_dens[i] <- dd$removal_dens

    outlier_rank_temp[!outlier_rank_temp][dd$choice_id] <- i
    x <- x[-dd$choice_id, , drop = FALSE]
    dist_mat <- dist_mat[-dd$choice_id, -dd$choice_id]
  }

  outlier_rank <- double(length(gross_outs))
  outlier_rank[gross_outs] <- 1
  outlier_rank[!gross_outs] <- outlier_rank_temp +
    gross_num * (outlier_rank_temp != 0)

  outlier_num <- integer(track_num)
  final_reject <- integer(track_num)
  after_final_reject <- matrix(nrow = max_out + 1, ncol = track_num)
  accept_bools <- matrix(nrow = max_out + 1, ncol = track_num)
  for (j in seq_len(track_num)) {
    reject_bool <- distrib_diff_mat[, j] > reject_num[j]
    final_reject[j] <- max(which(c(TRUE, reject_bool))) - 1
    after_final_reject[, j] <- seq_len(max_out + 1) > final_reject[j]

    accept_bools[, j] <- distrib_diff_mat[, j] < accept_num[j]
    outlier_num[j] <- which.max(accept_bools[, j] & after_final_reject[, j])
  }
  outlier_num <- outlier_num - 1 + gross_num

  consensus_vec <- apply(accept_bools & after_final_reject, 1, all)
  if (any(consensus_vec) && track_num > 1) {
    track_num_plus <- track_num + 1
    consensus_choice <- which.max(consensus_vec) - 1 + gross_num
    outlier_num[track_num_plus] <- consensus_choice
    en_names_plus <- c(en_names, "consensus")
  } else {
    consensus_choice <- NULL
    track_num_plus <- track_num
    en_names_plus <- en_names
  }

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
  observed <- rejection <- acceptance <- expected <- choice <- consensus <- NULL
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))
  for (j in seq_len(track_num)) {
    df_j <- data.frame(
      "outlier_seq" = outlier_seq,
      "observed" = distrib_diff_mat[, j],
      "expected" = expect_num[j],
      "acceptance" = accept_num[j],
      "rejection" = reject_num[j],
      "choice" = outlier_num[j]
    )

    gg_curves_list[[j]] <- df_j |>
      ggplot2::ggplot(ggplot2::aes(x = outlier_seq)) +
      ggplot2::geom_line(ggplot2::aes(y = observed, colour = "observed")) +
      ggplot2::geom_point(
        ggplot2::aes(y = observed, colour = "observed"),
        size = point_size
      ) +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = choice, colour = "choice"),
        linetype = "solid", linewidth = 0.75
      ) +
      ggplot2::geom_hline(
        ggplot2::aes(yintercept = rejection, colour = "rejection"),
        linetype = "dashed", linewidth = 0.75
      ) +
      ggplot2::geom_hline(
        ggplot2::aes(yintercept = acceptance, colour = "acceptance"),
        linetype = "dashed", linewidth = 0.75
      ) +
      ggplot2::geom_hline(
        ggplot2::aes(yintercept = expected, colour = "expected"),
        linetype = "dotted", linewidth = 0.75
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          observed = "#000000", expected = "#0072B2", acceptance = "#009E73",
          choice = "#CC79A7", rejection = "#D55E00", consensus = "#F0E442"
        )
      ) +
      ggplot2::labs(
        title = paste0(j, ": Number of Outliers = ", outlier_num[j]),
        x = "Outlier Number",
        y = "Number of Extreme Points",
        colour = ""
      ) +
      ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
      ggplot2::theme(
        legend.text = ggplot2::element_text(size = 11),
        legend.title = ggplot2::element_text(size = 11)
      ) +
      ggplot2::expand_limits(y = 0)

    if (track_num_plus > track_num) {
      df_j <- data.frame(
        "expected" = expect_num[j],
        "acceptance" = accept_num[j],
        "rejection" = reject_num[j],
        "choice" = outlier_num[j],
        "consensus" = consensus_choice
      )

      gg_curves_list[[j]] <- gg_curves_list[[j]] +
        ggplot2::geom_vline(
          data = df_j,
          ggplot2::aes(xintercept = consensus, colour = "consensus"),
          linetype = "solid", linewidth = 0.75
        )
    }
  }
  gg_curves <- ggpubr::ggarrange(
    plotlist = gg_curves_list,
    nrow = 1, ncol = track_num,
    common.legend = TRUE, legend = "bottom"
  )

  colnames(distrib_diff_mat) <- en_names
  colnames(outlier_bool) <- en_names_plus
  colnames(labels) <- en_names_plus
  names(outlier_num) <- en_names_plus
  dimnames(distrib_diff_arr) <- list(
    paste0("k", seq_len(comp_num)), NULL, en_names
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
    cat(paste0("Trying alternative covariance structures.\n"))
    try(mix <- mixture::gpcm(
      x,
      G = comp_num, start = z, nmax = nmax,
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
      cat(paste0("Trying defaul gpcm k-means initialisation.\n"))
      try(mix <- mixture::gpcm(x, G = comp_num, start = 2, nmax = nmax))
      if (!is.null(mix)) {
        cat(paste0(
          mix$best_model$cov_type, " covariance structure implemented.\n"
        ))
      }
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
