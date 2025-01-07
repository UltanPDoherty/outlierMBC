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
#' @param atol EM convergence tolerance threshold for mixture::gpcm.
#' @param init_method Method used to initialise each mixture model.
#' @param kmpp_seed Optional seed for k-means++ initialisation. Default is
#'                  hierarchical clustering.
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
#' ombc_gmm_k3n1000o10$plot_tail_curve
#' ombc_gmm_k3n1000o10$plot_full_curve
#'
ombc_gmm <- function(
    x,
    comp_num,
    max_out,
    gross_outs = rep(FALSE, nrow(x)),
    mnames = "VVV",
    nmax = 10,
    atol = 1e-8,
    init_method = c("hc", "kmpp"),
    kmpp_seed = 123,
    print_interval = Inf) {
  init_method <- match.arg(init_method)

  expect_num <- 1
  accept_num <- expect_num + 1
  reject_num <- accept_num

  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x)

  dist_mat0 <- as.matrix(stats::dist(x0))
  dist_mat <- dist_mat0

  gross_num <- sum(gross_outs)
  x <- x[!gross_outs, ]
  max_out <- max_out - gross_num
  dist_mat <- dist_mat[!gross_outs, !gross_outs]

  track_num <- 2 + 1
  tail_props <- expect_num / (seq(obs_num, obs_num - max_out) - gross_num)

  loglike <- c()
  removal_dens <- c()
  distrib_diff_arr <- array(dim = c(comp_num, max_out + 1, track_num))
  distrib_diff_mat <- matrix(nrow = max_out + 1, ncol = track_num)
  outlier_rank_temp <- rep(0, obs_num - gross_num)
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    z <- get_init_z(
      comp_num,
      dist_mat = dist_mat, x = x,
      init_method = init_method, kmpp_seed = kmpp_seed
    )

    mix <- try_mixture_gpcm(x, comp_num, mnames, z, nmax, atol)

    loglike[i] <- mix$best_model$loglik

    dd <- distrib_diff_gmm(
      x,
      mix$z,
      mix$best_model$model_obj[[1]]$pi_gs,
      mix$best_model$model_obj[[1]]$mu,
      mix$best_model$model_obj[[1]]$sigs,
      mix$best_model$model_obj[[1]]$log_dets,
      tail_props[i]
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

  outlier_num[1] <- which.min(distrib_diff_mat[, 1])

  reject_bool <- distrib_diff_mat[, 2] > reject_num
  final_reject <- max(which(c(TRUE, reject_bool))) - 1
  after_final_reject <- seq_len(max_out + 1) > final_reject
  accept_bools <- distrib_diff_mat[, 2] < accept_num
  outlier_num[2] <- which.max(accept_bools & after_final_reject)

  outlier_num <- outlier_num - 1 + gross_num

  outlier_bool <- matrix(nrow = obs_num, ncol = track_num)
  mix <- list()
  labels <- matrix(0, nrow = obs_num, ncol = track_num)
  for (j in seq_len(track_num)) {
    outlier_bool[, j] <- outlier_rank <= outlier_num[j] & outlier_rank != 0

    z <- get_init_z(
      comp_num,
      dist_mat = dist_mat0[!outlier_bool[, j], !outlier_bool[, j]],
      x = x0[!outlier_bool[, j], ],
      init_method = init_method, kmpp_seed = kmpp_seed
    )

    mix[[j]] <- try_mixture_gpcm(
      x0[!outlier_bool[, j], ], comp_num, mnames, z, nmax, atol
    )

    labels[!outlier_bool[, j], j] <- mix[[j]]$map
  }

  outlier_seq <- seq(gross_num, max_out + gross_num)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  full <- minimum <- NULL
  full_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "minimum" = outlier_num[1],
    "full" = distrib_diff_mat[, 1]
  )
  full_curve <- full_curve_df |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = full)) +
    ggplot2::geom_line(
      ggplot2::aes(colour = "full"),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = "full"),
      size = point_size, show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = minimum, colour = "minimum"),
      linetype = "solid", linewidth = 0.75, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(full = "#000000", minimum = "#CC79A7")
    ) +
    ggplot2::labs(
      title = paste0("outlierMBC: Number of Outliers = ", outlier_num[1]),
      x = "Outlier Number",
      y = "Mean Absolute CDF Difference",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq))

  observed <- acceptance <- expected <- choice <- NULL
  tail_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "observed" = distrib_diff_mat[, 2],
    "expected" = expect_num,
    "acceptance" = accept_num,
    "rejection" = reject_num,
    "choice" = outlier_num[2]
  )
  tail_curve <- tail_curve_df |>
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
        choice = "#CC79A7", rejection = "#D55E00"
      )
    ) +
    ggplot2::labs(
      title = paste0("outlierMBC-tail: Number of Outliers = ", outlier_num[2]),
      x = "Outlier Number",
      y = "Number of Extreme Points",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 11),
      legend.position = "bottom"
    ) +
    ggplot2::expand_limits(y = 0)

  ombc_names <- c("full", "tail", "bin_tail")
  colnames(distrib_diff_mat) <- ombc_names
  colnames(outlier_bool) <- ombc_names
  colnames(labels) <- ombc_names
  names(outlier_num) <- ombc_names
  dimnames(distrib_diff_arr) <- list(
    paste0("k", seq_len(comp_num)), NULL, ombc_names
  )
  names(mix) <- ombc_names

  outlier_class <- rep("normal", obs_num)
  outlier_class[outlier_bool[, 1] & outlier_bool[, 2]] <- "out_full_&_tail"
  outlier_class[outlier_bool[, 1] & !outlier_bool[, 2]] <- "out_full_only"
  outlier_class[!outlier_bool[, 1] & outlier_bool[, 2]] <- "out_tail_only"
  outlier_class[gross_outs] <- "out_gross"
  outlier_class <- as.factor(outlier_class)

  labels <- as.data.frame(labels)
  outlier_bool <- as.data.frame(outlier_bool)

  return(list(
    labels = labels,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    outlier_class = outlier_class,
    curve_plot = full_curve,
    tail_curve_plot = tail_curve,
    mix = mix,
    loglike = loglike,
    removal_dens = removal_dens,
    distrib_diff_mat = distrib_diff_mat,
    distrib_diff_arr = distrib_diff_arr
  ))
}

# ------------------------------------------------------------------------------

get_init_z <- function(
    comp_num,
    dist_mat = NULL, x = NULL,
    init_method = c("hc", "kmpp"), kmpp_seed = NULL) {
  init_method <- match.arg(init_method)

  if (init_method == "hc") {
    hc <- stats::hclust(stats::as.dist(dist_mat), method = "ward.D2")
    init <- stats::cutree(hc, k = comp_num)
  } else if (init_method == "kmpp") {
    kmpp <- ClusterR::KMeans_rcpp(
      x,
      clusters = comp_num, num_init = 10, seed = kmpp_seed
    )
    init <- kmpp$clusters
  }

  z <- matrix(nrow = nrow(dist_mat), ncol = comp_num)
  for (k in seq_len(comp_num)) {
    z[, k] <- as.integer(init == k)
  }

  return(z)
}

# ------------------------------------------------------------------------------

try_mixture_gpcm <- function(x, comp_num, mnames, z, nmax, atol) {
  mix <- NULL
  try(mix <- mixture::gpcm(
    x,
    G = comp_num, mnames = mnames, start = z, nmax = nmax, atol = atol
  ))
  if (is.null(mix)) {
    cat(paste0("Trying alternative covariance structures.\n"))
    try(mix <- mixture::gpcm(
      x,
      G = comp_num, start = z, nmax = nmax, atol = atol,
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
      cat(paste0("Trying default gpcm k-means initialisation.\n"))
      try(
        mix <- mixture::gpcm(
          x,
          G = comp_num, start = 2, nmax = nmax, atol = atol
        )
      )
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

# ==============================================================================

#' ombc2_gmm
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
#' @param atol EM convergence tolerance threshold for mixture::gpcm.
#' @param init_method Method used to initialise each mixture model.
#' @param kmpp_seed Optional seed for k-means++ initialisation. Default is
#'                  hierarchical clustering.
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
#' ombc2_gmm_k3n1000o10 <- ombc2_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   comp_num = 3, max_out = 20
#' )
#'
#' ombc2_gmm_k3n1000o10$plot_tail_curve
#' ombc2_gmm_k3n1000o10$plot_full_curve
#'
ombc2_gmm <- function(
    x,
    comp_num,
    max_out,
    gross_outs = rep(FALSE, nrow(x)),
    mnames = "VVV",
    nmax = 1000,
    atol = 1e-8,
    init_method = c("hc", "kmpp"),
    kmpp_seed = 123,
    print_interval = Inf) {
  init_method <- match.arg(init_method)

  expect_num <- 1
  accept_num <- expect_num + 1
  reject_num <- accept_num

  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x)

  dist_mat0 <- as.matrix(stats::dist(x0))
  dist_mat <- dist_mat0

  gross_num <- sum(gross_outs)
  x <- x[!gross_outs, ]
  max_out <- max_out - gross_num
  dist_mat <- dist_mat[!gross_outs, !gross_outs]

  set.seed(123)
  mclust_pre <- mclust::Mclust(x, comp_num, mnames)

  pre_dens <- mclust::dens(
    data = x,
    modelName = mclust_pre$modelName,
    parameters = mclust_pre$parameters
  )

  init_noise <- pre_dens < mclust::hypvol(x, TRUE)

  set.seed(123)
  mclust_out <- mclust::Mclust(
    data = x,
    G = comp_num,
    modelNames = mnames,
    initialization = list(noise = init_noise)
  )

  cat(paste0(
    "With ", gross_num,
    " gross outliers removed, mclust identifies a further ",
    sum(mclust_out$classification == 0), " outliers.\n\n"
  ))

  mclust_params <- mclust_out$parameters
  mclust_params$pro <- mclust_params$pro[-(comp_num + 1)]
  mclust_params$pro <- mclust_params$pro / sum(mclust_params$pro)

  z <- mclust::estep(x, mnames, mclust_params)$z

  track_num <- 2 + 1
  tail_props <- expect_num / (seq(obs_num, obs_num - max_out) - gross_num)

  loglike <- c()
  removal_dens <- c()
  distrib_diff_arr <- array(dim = c(comp_num, max_out + 1, track_num))
  distrib_diff_mat <- matrix(nrow = max_out + 1, ncol = track_num)
  outlier_rank_temp <- rep(0, obs_num - gross_num)
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    mix <- try_mixture_gpcm(x, comp_num, mnames, z, nmax, atol)

    loglike[i] <- mix$best_model$loglik

    dd <- distrib_diff_gmm(
      x,
      mix$z,
      mix$best_model$model_obj[[1]]$pi_gs,
      mix$best_model$model_obj[[1]]$mu,
      mix$best_model$model_obj[[1]]$sigs,
      mix$best_model$model_obj[[1]]$log_dets,
      tail_props[i]
    )

    distrib_diff_arr[, i, ] <- dd$distrib_diff_mat
    distrib_diff_mat[i, ] <- dd$distrib_diff_vec
    removal_dens[i] <- dd$removal_dens

    outlier_rank_temp[!outlier_rank_temp][dd$choice_id] <- i
    x <- x[-dd$choice_id, , drop = FALSE]
    z <- z[-dd$choice_id, , drop = FALSE]
    dist_mat <- dist_mat[-dd$choice_id, -dd$choice_id]
  }

  outlier_rank <- double(length(gross_outs))
  outlier_rank[gross_outs] <- 1
  outlier_rank[!gross_outs] <- outlier_rank_temp +
    gross_num * (outlier_rank_temp != 0)

  outlier_num <- integer(track_num)

  outlier_num[1] <- which.min(distrib_diff_mat[, 1])

  reject_bool <- distrib_diff_mat[, 2] > reject_num
  final_reject <- max(which(c(TRUE, reject_bool))) - 1
  after_final_reject <- seq_len(max_out + 1) > final_reject
  accept_bools <- distrib_diff_mat[, 2] < accept_num
  outlier_num[2] <- which.max(accept_bools & after_final_reject)

  reject_bool2 <- distrib_diff_mat[, 3] > reject_num
  final_reject2 <- max(which(c(TRUE, reject_bool2))) - 1
  after_final_reject2 <- seq_len(max_out + 1) > final_reject2
  accept_bools2 <- distrib_diff_mat[, 3] < accept_num
  outlier_num[3] <- which.max(accept_bools2 & after_final_reject2)

  outlier_num <- outlier_num - 1 + gross_num

  outlier_bool <- matrix(nrow = obs_num, ncol = track_num)
  mix <- list()
  labels <- matrix(0, nrow = obs_num, ncol = track_num)
  for (j in seq_len(track_num)) {
    outlier_bool[, j] <- outlier_rank <= outlier_num[j] & outlier_rank != 0

    z <- get_init_z(
      comp_num,
      dist_mat = dist_mat0[!outlier_bool[, j], !outlier_bool[, j]],
      x = x0[!outlier_bool[, j], ],
      init_method = init_method, kmpp_seed = kmpp_seed
    )

    mix[[j]] <- try_mixture_gpcm(
      x0[!outlier_bool[, j], ], comp_num, mnames, z, nmax, atol
    )

    labels[!outlier_bool[, j], j] <- mix[[j]]$map
  }

  outlier_seq <- seq(gross_num, max_out + gross_num)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  full <- minimum <- NULL
  full_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "minimum" = outlier_num[1],
    "full" = distrib_diff_mat[, 1]
  )
  full_curve <- full_curve_df |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = full)) +
    ggplot2::geom_line(
      ggplot2::aes(colour = "full"),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = "full"),
      size = point_size, show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = minimum, colour = "minimum"),
      linetype = "solid", linewidth = 0.75, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(full = "#000000", minimum = "#CC79A7")
    ) +
    ggplot2::labs(
      title = paste0("outlierMBC: Number of Outliers = ", outlier_num[1]),
      x = "Outlier Number",
      y = "Mean Absolute CDF Difference",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq))

  observed <- acceptance <- expected <- choice <- NULL
  tail_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "observed" = distrib_diff_mat[, 2],
    "expected" = expect_num,
    "acceptance" = accept_num,
    "rejection" = reject_num,
    "choice" = outlier_num[2]
  )
  tail_curve <- tail_curve_df |>
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
        choice = "#CC79A7", rejection = "#D55E00"
      )
    ) +
    ggplot2::labs(
      title = paste0("outlierMBC-tail: Number of Outliers = ", outlier_num[2]),
      x = "Outlier Number",
      y = "Number of Extreme Points",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 11),
      legend.position = "bottom"
    ) +
    ggplot2::expand_limits(y = 0)

  ombc_names <- c("full", "tail", "bin_tail")
  colnames(distrib_diff_mat) <- ombc_names
  colnames(outlier_bool) <- ombc_names
  colnames(labels) <- ombc_names
  names(outlier_num) <- ombc_names
  dimnames(distrib_diff_arr) <- list(
    paste0("k", seq_len(comp_num)), NULL, ombc_names
  )

  outlier_class <- rep("normal", obs_num)
  outlier_class[outlier_bool[, 1] & outlier_bool[, 2]] <- "out_full_&_tail"
  outlier_class[outlier_bool[, 1] & !outlier_bool[, 2]] <- "out_full_only"
  outlier_class[!outlier_bool[, 1] & outlier_bool[, 2]] <- "out_tail_only"
  outlier_class[gross_outs] <- "out_gross"
  outlier_class <- as.factor(outlier_class)

  labels <- as.data.frame(labels)
  outlier_bool <- as.data.frame(outlier_bool)

  return(list(
    labels = labels,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    outlier_class = outlier_class,
    curve_plot = full_curve,
    tail_curve_plot = tail_curve,
    loglike = loglike,
    removal_dens = removal_dens,
    distrib_diff_mat = distrib_diff_mat,
    distrib_diff_arr = distrib_diff_arr
  ))
}


# ==============================================================================

#' ombc3_gmm
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
#' @param atol EM convergence tolerance threshold for mixture::gpcm.
#' @param init_method Method used to initialise each mixture model.
#' @param kmpp_seed Optional seed for k-means++ initialisation. Default is
#'                  hierarchical clustering.
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
#' ombc2_gmm_k3n1000o10 <- ombc2_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   comp_num = 3, max_out = 20
#' )
#'
#' ombc2_gmm_k3n1000o10$plot_tail_curve
#' ombc2_gmm_k3n1000o10$plot_full_curve
#'
ombc3_gmm <- function(
    x,
    comp_num,
    max_out,
    gross_outs = rep(FALSE, nrow(x)),
    mnames = "VVV",
    init_method = c("hc", "kmpp"),
    kmpp_seed = 123,
    print_interval = Inf) {
  init_method <- match.arg(init_method)

  expect_num <- 1
  accept_num <- expect_num + 1
  reject_num <- accept_num

  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x)

  dist_mat0 <- as.matrix(stats::dist(x0))
  dist_mat <- dist_mat0

  gross_num <- sum(gross_outs)
  x <- x[!gross_outs, ]
  max_out <- max_out - gross_num
  dist_mat <- dist_mat[!gross_outs, !gross_outs]

  set.seed(123)
  mclust_pre <- mclust::Mclust(x, comp_num, mnames)

  pre_dens <- mclust::dens(
    data = x,
    modelName = mclust_pre$modelName,
    parameters = mclust_pre$parameters
  )

  init_noise <- pre_dens < mclust::hypvol(x, TRUE)

  set.seed(123)
  mclust_out <- mclust::Mclust(
    data = x,
    G = comp_num,
    modelNames = mnames,
    initialization = list(noise = init_noise)
  )

  cat(paste0(
    "With ", gross_num,
    " gross outliers removed, mclust identifies a further ",
    sum(mclust_out$classification == 0), " outliers.\n\n"
  ))

  mclust_params <- mclust_out$parameters
  mclust_params$pro <- mclust_params$pro[-(comp_num + 1)]
  mclust_params$pro <- mclust_params$pro / sum(mclust_params$pro)

  mclust_out_mu <- lapply(
    seq_len(comp_num), function(i) mclust_out$parameters$mean[,i]
  )
  mclust_out_sigma <- lapply(
    seq_len(comp_num),
    function(i) mclust_out$parameters$variance$sigma[, ,i]
  )
  mclust_out_logdet <- vapply(
    seq_len(comp_num),
    function(i) log(det(mclust_out$parameters$variance$sigma[, ,i])),
    double(1L)
  )
  mclust_out_prop <- t(mclust_params$pro)

  z0 <- mclust::estep(x, mnames, mclust_params)$z
  z <- z0

  track_num <- 2 + 1
  tail_props <- expect_num / (seq(obs_num, obs_num - max_out) - gross_num)

  loglike <- c()
  removal_dens <- c()
  distrib_diff_arr <- array(dim = c(comp_num, max_out + 1, track_num))
  distrib_diff_mat <- matrix(nrow = max_out + 1, ncol = track_num)
  outlier_rank_temp <- rep(0, obs_num - gross_num)
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    dd <- distrib_diff_gmm(
      x,
      z,
      mclust_out_prop,
      mclust_out_mu,
      mclust_out_sigma,
      mclust_out_logdet,
      tail_props[i]
    )

    distrib_diff_arr[, i, ] <- dd$distrib_diff_mat
    distrib_diff_mat[i, ] <- dd$distrib_diff_vec
    removal_dens[i] <- dd$removal_dens

    outlier_rank_temp[!outlier_rank_temp][dd$choice_id] <- i
    x <- x[-dd$choice_id, , drop = FALSE]
    z <- z[-dd$choice_id, , drop = FALSE]
    dist_mat <- dist_mat[-dd$choice_id, -dd$choice_id]
  }

  outlier_rank <- double(length(gross_outs))
  outlier_rank[gross_outs] <- 1
  outlier_rank[!gross_outs] <- outlier_rank_temp +
    gross_num * (outlier_rank_temp != 0)

  outlier_num <- integer(track_num)

  outlier_num[1] <- which.min(distrib_diff_mat[, 1])

  reject_bool <- distrib_diff_mat[, 2] > reject_num
  final_reject <- max(which(c(TRUE, reject_bool))) - 1
  after_final_reject <- seq_len(max_out + 1) > final_reject
  accept_bools <- distrib_diff_mat[, 2] < accept_num
  outlier_num[2] <- which.max(accept_bools & after_final_reject)

  reject_bool2 <- distrib_diff_mat[, 3] > reject_num
  final_reject2 <- max(which(c(TRUE, reject_bool2))) - 1
  after_final_reject2 <- seq_len(max_out + 1) > final_reject2
  accept_bools2 <- distrib_diff_mat[, 3] < accept_num
  outlier_num[3] <- which.max(accept_bools2 & after_final_reject2)

  outlier_num <- outlier_num - 1 + gross_num

  outlier_bool <- matrix(nrow = obs_num, ncol = track_num)
  mix <- list()
  labels <- matrix(0, nrow = obs_num, ncol = track_num)
  for (j in seq_len(track_num)) {
    outlier_bool[, j] <- outlier_rank <= outlier_num[j] & outlier_rank != 0

    mix[[j]] <- mclust::Mclust(x0[!outlier_bool[, j], ], comp_num, mnames)

    labels[!outlier_bool[, j], j] <- mix[[j]]$classification
  }

  outlier_seq <- seq(gross_num, max_out + gross_num)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  full <- minimum <- NULL
  full_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "minimum" = outlier_num[1],
    "full" = distrib_diff_mat[, 1]
  )
  full_curve <- full_curve_df |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = full)) +
    ggplot2::geom_line(
      ggplot2::aes(colour = "full"),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = "full"),
      size = point_size, show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = minimum, colour = "minimum"),
      linetype = "solid", linewidth = 0.75, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(full = "#000000", minimum = "#CC79A7")
    ) +
    ggplot2::labs(
      title = paste0("outlierMBC: Number of Outliers = ", outlier_num[1]),
      x = "Outlier Number",
      y = "Mean Absolute CDF Difference",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq))

  observed <- acceptance <- expected <- choice <- NULL
  tail_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "observed" = distrib_diff_mat[, 2],
    "expected" = expect_num,
    "acceptance" = accept_num,
    "rejection" = reject_num,
    "choice" = outlier_num[2]
  )
  tail_curve <- tail_curve_df |>
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
        choice = "#CC79A7", rejection = "#D55E00"
      )
    ) +
    ggplot2::labs(
      title = paste0("outlierMBC-tail: Number of Outliers = ", outlier_num[2]),
      x = "Outlier Number",
      y = "Number of Extreme Points",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 11),
      legend.position = "bottom"
    ) +
    ggplot2::expand_limits(y = 0)

  ombc_names <- c("full", "tail", "bin_tail")
  colnames(distrib_diff_mat) <- ombc_names
  colnames(outlier_bool) <- ombc_names
  colnames(labels) <- ombc_names
  names(outlier_num) <- ombc_names
  dimnames(distrib_diff_arr) <- list(
    paste0("k", seq_len(comp_num)), NULL, ombc_names
  )

  outlier_class <- rep("normal", obs_num)
  outlier_class[outlier_bool[, 1] & outlier_bool[, 2]] <- "out_full_&_tail"
  outlier_class[outlier_bool[, 1] & !outlier_bool[, 2]] <- "out_full_only"
  outlier_class[!outlier_bool[, 1] & outlier_bool[, 2]] <- "out_tail_only"
  outlier_class[gross_outs] <- "out_gross"
  outlier_class <- as.factor(outlier_class)

  labels <- as.data.frame(labels)
  outlier_bool <- as.data.frame(outlier_bool)

  return(list(
    labels = labels,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    outlier_class = outlier_class,
    curve_plot = full_curve,
    tail_curve_plot = tail_curve,
    loglike = loglike,
    removal_dens = removal_dens,
    distrib_diff_mat = distrib_diff_mat,
    distrib_diff_arr = distrib_diff_arr
  ))
}




# ==============================================================================

#' ombc3_gmm
#'
#' @description
#' Iterative Detection & Identification of Outliers for a Gaussian Mixture Model
#'
#' @param x Data.
#' @param comp_num Number of components.
#' @param max_out Maximum number of outliers.
#' @param gross_outs Logical vector identifying gross outliers.
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
#' * removal_dens
#' @export
#'
#' @examples
#'
#' ombc2_gmm_k3n1000o10 <- ombc2_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   comp_num = 3, max_out = 20
#' )
#'
#' ombc2_gmm_k3n1000o10$plot_tail_curve
#' ombc2_gmm_k3n1000o10$plot_full_curve
#'
ombc4_gmm <- function(
    x,
    comp_num,
    max_out,
    gross_outs = rep(FALSE, nrow(x)),
    mnames = "VVV",
    print_interval = Inf) {
  x <- as.matrix(x)
  x0 <- x

  obs_num <- nrow(x)

  gross_num <- sum(gross_outs)
  x <- x[!gross_outs, ]
  max_out <- max_out - gross_num

  cat(paste0(
    gross_num,
    " gross outliers have been pre-specified.\n\n"
  ))

  set.seed(123)
  mclust_pre <- mclust::Mclust(x, comp_num, mnames, verbose = FALSE)

  pre_dens <- mclust::dens(
    data = x,
    modelName = mclust_pre$modelName,
    parameters = mclust_pre$parameters
  )

  init_noise <- pre_dens < mclust::hypvol(x, TRUE)

  set.seed(123)
  mclust_out <- mclust::Mclust(
    data = x,
    G = comp_num,
    modelNames = mnames,
    initialization = list(noise = init_noise),
    verbose = FALSE
  )

  gpcm_out <- mixture::gpcm(
    x, comp_num, mnames
  )

  init_hc <- cutree(hclust(dist(x), method = "ward.D2"), k = comp_num)
  init_z <- matrix(0, nrow = nrow(x), ncol = comp_num)
  for (g in seq_len(comp_num)) {
    init_z[init_hc == g, g] <- 1
  }

  gpcm_out2 <- mixture::gpcm(
    x, comp_num, mnames, start = init_z
  )


  cat(paste0(
    "mclust would have identified a further ",
    sum(mclust_out$classification == 0), " outliers.\n\n"
  ))

  mclust_params <- mclust_out$parameters
  mclust_params$pro <- mclust_params$pro[-(comp_num + 1)]
  mclust_params$pro <- mclust_params$pro / sum(mclust_params$pro)
  mclust_params$Vinv <- NULL

  mclust_out_mu <- lapply(
    seq_len(comp_num), function(i) mclust_out$parameters$mean[,i]
  )
  mclust_out_sigma <- lapply(
    seq_len(comp_num),
    function(i) mclust_out$parameters$variance$sigma[, ,i]
  )
  mclust_out_logdet <- vapply(
    seq_len(comp_num),
    function(i) log(det(mclust_out$parameters$variance$sigma[, ,i])),
    double(1L)
  )
  mclust_out_prop <- t(mclust_params$pro)

  z0 <- mclust::estep(x, mnames, mclust_params)$z
  z <- z0

  bin_z_vec <- apply(z, 1, which.max)
  bin_z <- matrix(0, nrow = obs_num - gross_num, ncol = comp_num)
  for (g in seq_len(comp_num)) {
    bin_z[, g] <- bin_z_vec == g
  }

  dens <- mclust::dens(x, mnames, mclust_params)

  outlier_bool_short <- count_extremes(
    x,
    bin_z,
    mclust_out_prop,
    mclust_out_mu,
    mclust_out_sigma,
    mclust_out_logdet,
    dens
  )

  outlier_bool <- rep(FALSE, obs_num)
  outlier_bool[gross_outs] <- TRUE
  outlier_bool[!gross_outs][outlier_bool_short] <- TRUE

  outlier_num <- sum(outlier_bool)

  labels <- rep(0, obs_num)
  labels[!outlier_bool] <- mclust::Mclust(
    x0[!outlier_bool, ], comp_num, mnames, verbose = FALSE
  )$classification

  cat(paste0(
    "Instead, ombc4_gmm identifies a further ",
    outlier_num - gross_num, " outliers.\n\n"
  ))

  cat(paste0(
    "The final number of outliers is ",
    outlier_num, ".\n\n"
  ))

  return(list(
    labels = labels,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num
  ))
}
