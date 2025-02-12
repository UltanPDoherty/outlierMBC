#' Title
#'
#' @inheritParams ombc_gmm
#' @param max_value Value cannot exceed minimum * (1 + max_value).
#' @param max_step Each step must be less than minimum * (1 + max_step).
#'
#' @returns List of two lists:
#' * minimum: $ind
#'            $val
#' * backtrack: $ind
#'              $val
#' @export
#'
#' @examples
#'
#' ombc_gmm_k3n1000o10 <- ombc_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   comp_num = 3, max_out = 20
#' )
#'
#' backtrack(ombc_gmm_k3n1000o10$distrib_diff_mat[, "full"])
#'
backtrack <- function(x, max_value = 0.1, max_step = 0.01) {
  xmin_val <- min(x)
  xmin_ind <- which.min(x)

  valid_step <- x < dplyr::lead(x) + max_step * xmin_val
  valid_value <- x < (1 + max_value) * xmin_val

  y <- xmin_ind

  keep_going <- TRUE

  while (keep_going && y > 1) {
    if (valid_value[y - 1] && valid_step[y - 1]) {
      y <- y - 1
      keep_going <- TRUE
    } else {
      keep_going <- FALSE
    }
  }

  minimum <- list("ind" = xmin_ind, "val" = xmin_val)
  backtrack <- list("ind" = y, "val" = x[y])

  return(list("minimum" = minimum, "backtrack" = backtrack))
}


#' Title
#'
#' @inheritParams ombc_gmm
#' @inheritParams plot_full_curve
#' @inheritParams backtrack
#'
#' @returns List:
#' * labels
#' * outlier_bool
#' * outlier_num
#' * mix
#'
#' @export
#'
#' @examples
#'
#' ombc_gmm_k3n1000o10 <- ombc_gmm(
#'   gmm_k3n1000o10[, 1:2],
#'   comp_num = 3, max_out = 20
#' )
#'
#' backtrack_gmm(gmm_k3n1000o10[, 1:2], ombc_gmm_k3n1000o10, 0.1, 0.01)
backtrack_gmm <- function(
    x, ombc_out,
    max_value = 0.1, max_step = 0.01, init_model = NULL, init_z = NULL) {
  backtrack_out <-
    backtrack(ombc_out$distrib_diff_mat[, "full"], max_value, max_step)

  x0 <- as.matrix(x)

  outlier_num <- backtrack_out$backtrack$ind - 1 + sum(ombc_out$gross_outs)

  outlier_bool <-
    ombc_out$outlier_rank <= outlier_num & ombc_out$outlier_rank != 0

  init_scheme <- ombc_out$call$init_scheme

  stopifnot(
    "'update' scheme cannot be used with backtrack" = init_scheme != "update"
  )

  if (init_scheme == "update") {
    stop("backtrack_gmm is not compatible with the 'update' init_scheme.")
  } else if (init_scheme == "reinit") {
    if (!is.null(init_model) && !is.null(init_z)) {
      stop("init_model & init_z cannot be used with the 'reinit' init_scheme.")
    } else if (!is.null(init_z)) {
      stop("init_z cannot be used with the 'reinit' init_scheme.")
    } else if (!is.null(init_model)) {
      stop("init_model cannot be used with the 'reinit' init_scheme.")
    } else {
      z <- get_init_z(
        comp_num = ombc_out$call$comp_num,
        dist_mat = as.matrix(stats::dist(x0[!outlier_bool, ])),
        x = x0[!outlier_bool, ],
        init_method = ombc_out$call$init_method,
        kmpp_seed = ombc_out$call$kmpp_seed
      )
    }
  } else if (init_scheme == "reuse") {
    if (!is.null(init_model) && !is.null(init_z)) {
      stop("Only one of init_model and init_z may be provided.")
    } else if (!is.null(init_z)) {
      z <- init_z
    } else if (!is.null(init_model)) {
      z <- mixture::e_step(x0[!outlier_bool, ], init_model)$z
    } else {
      stop("init_model or init_z are required to use the 'reuse' init_scheme.")
    }
  } else {
    stop("init_scheme must be either 'reinit' or 'reuse'.")
  }

  mix <- try_mixture_gpcm(
    x0[!outlier_bool, ],
    ombc_out$call$comp_num, ombc_out$call$mnames,
    z,
    ombc_out$call$nmax,
    ombc_out$call$atol
  )

  labels <- mix$map

  return(list(
    "labels" = labels,
    "outlier_bool" = outlier_bool,
    "outlier_num" = outlier_num,
    "mix" = mix
  ))
}

#' Plot the outlier number selection curve for the backtrack method.
#'
#' @inheritParams plot_full_curve
#' @inheritParams backtrack
#'
#' @returns A gg object.
#' @export
plot_backtrack_curve <- function(ombc_out, max_value, max_step) {
  gross_num <- sum(ombc_out$gross_outs)
  max_out <- max(ombc_out$outlier_rank) - 1
  outlier_num <- ombc_out$outlier_num
  distrib_diff_mat <- ombc_out$distrib_diff_mat

  backtrack_num <- backtrack(
    ombc_out$distrib_diff[, "full"], max_value, max_step
  )$backtrack$ind

  outlier_seq <- seq(gross_num, max_out)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  backtrack <- minimum <- choice <- NULL
  backtrack_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "minimum" = as.integer(outlier_num["full"]),
    "choice" = as.integer(backtrack_num),
    "backtrack" = distrib_diff_mat[, "full"] / min(distrib_diff_mat[, "full"])
  )
  backtrack_curve <- backtrack_curve_df |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = backtrack)) +
    ggplot2::geom_line(
      ggplot2::aes(colour = "backtrack"),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = "backtrack"),
      size = point_size, show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = choice, colour = "choice"),
      linetype = "solid", linewidth = 0.75, show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = minimum, colour = "minimum"),
      linetype = "dashed", linewidth = 0.75, show.legend = FALSE
    ) +
    ggplot2::geom_hline(
      ggplot2::aes(
        yintercept = min(backtrack) * (1 + max_value), colour = "minimum"
      ),
      linetype = "dashed", linewidth = 0.75, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(backtrack = "#000000", choice = "#CC79A7", minimum = "#009E73")
    ) +
    ggplot2::labs(
      title =
        paste0("outlierMBC: Number of Outliers = ", outlier_num["backtrack"]),
      subtitle = paste0("max_value = ", max_value, ", max_step = ", max_step),
      x = "Outlier Number",
      y = "Rescaled Mean Absolute CDF Difference",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
    ggplot2::expand_limits(y = 0)

  return(backtrack_curve)
}
