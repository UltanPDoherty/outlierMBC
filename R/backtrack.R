#' Move backwards from the global minimum to a more conservative solution.
#'
#' @inheritParams ombc_gmm
#' @param max_total_rise Value cannot exceed minimum * (1 + max_total_rise).
#' @param max_step_rise
#' Each step must be less than minimum * (1 + max_step_rise).
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
#' backtrack(ombc_gmm_k3n1000o10$distrib_diff_vec)
#'
backtrack <- function(x, max_total_rise = 0.1, max_step_rise = 0.05) {
  xmin_val <- min(x)
  xmin_ind <- which.min(x)

  valid_step <- x < dplyr::lead(x) + max_step_rise * xmin_val
  valid_value <- x < (1 + max_total_rise) * xmin_val

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

#' Plot the outlier number selection curve for the backtrack method.
#'
#' @inheritParams plot_curve
#' @inheritParams backtrack
#'
#' @returns A gg object.
#' @export
plot_backtrack <- function(
    ombc_out, max_total_rise = 0.1, max_step_rise = 0.05) {
  gross_num <- sum(ombc_out$gross_outs)
  max_out <- max(ombc_out$outlier_rank) - 1
  outlier_num <- ombc_out$outlier_num
  distrib_diff_vec <- ombc_out$distrib_diff_vec

  backtrack_num <- backtrack(
    distrib_diff_vec, max_total_rise, max_step_rise
  )$backtrack$ind + gross_num - 1

  outlier_seq <- seq(gross_num, max_out)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  backtrack <- minimum <- choice <- NULL
  backtrack_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "minimum" = as.integer(outlier_num),
    "choice" = as.integer(backtrack_num),
    "backtrack" = distrib_diff_vec / min(distrib_diff_vec)
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
        yintercept = min(backtrack) * (1 + max_total_rise), colour = "minimum"
      ),
      linetype = "dashed", linewidth = 0.75, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(backtrack = "#000000", choice = "#CC79A7", minimum = "#009E73")
    ) +
    ggplot2::labs(
      title =
        paste0("outlierMBC: Number of Outliers = ", backtrack_num),
      subtitle = paste0(
        "max_total_rise = ", max_total_rise,
        ", max_step_rise = ", max_step_rise
      ),
      x = "Outlier Number",
      y = "Rescaled Dissimilarity",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
    ggplot2::expand_limits(y = 0)

  return(backtrack_curve)
}

# ------------------------------------------------------------------------------

#' Fit a GMM to the backtrack solution.
#'
#' @inheritParams ombc_gmm
#' @inheritParams plot_curve
#' @inheritParams backtrack
#' @param manual_outlier_num description
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
    max_total_rise = 0.1, max_step_rise = 0.05,
    init_model = NULL, init_z = NULL, manual_outlier_num = NULL) {
  this_call <- call(
    "backtrack_gmm",
    "x" = substitute(x), "ombc_out" = substitute(ombc_out),
    "max_total_rise" = max_total_rise, "max_step_rise" = max_step_rise,
    "init_z" = substitute(init_z), "init_model" = substitute(init_model)
  )

  x0 <- as.matrix(x)
  gross_num <- sum(ombc_out$gross_outs)

  if (is.null(manual_outlier_num)) {
    backtrack_out <- backtrack(
      ombc_out$distrib_diff_vec, max_total_rise, max_step_rise
    )

    if (backtrack_out$backtrack$ind == backtrack_out$minimum$ind) {
      cat(paste0(
        "backtrack stayed at the minimum.",
        "backtrack_gmm will return ombc_gmm results directly.\n"
      ))

      return(list(
        "labels" = ombc_out$labels,
        "outlier_bool" = ombc_out$outlier_bool,
        "outlier_num" = ombc_out$outlier_num,
        "mix" = ombc_out$mix,
        "call" = this_call
      ))
    }

    outlier_num <- backtrack_out$backtrack$ind - 1 + gross_num
  } else {
    outlier_num <- manual_outlier_num
  }

  outlier_bool <-
    ombc_out$outlier_rank <= outlier_num & ombc_out$outlier_rank != 0

  init_scheme <- ombc_out$call$init_scheme

  stopifnot(
    "init_scheme must be 'update' or 'reinit' or 'reuse'." = (
      init_scheme %in% c("update", "reinit", "reuse")
    )
  )
  stopifnot(
    "Neither init_model nor init_z can be used with the 'reinit' init_scheme." =
      (init_scheme != "reinit") || (is.null(init_model) && is.null(init_z))
  )
  stopifnot(
    "Only one of init_model and init_z may be provided." = (
      is.null(init_model) || is.null(init_z)
    )
  )

  if (!is.null(init_z)) {
    z0 <- init_z
  } else if (!is.null(init_model)) {
    z0 <- mixture::e_step(x0[!ombc_out$gross_outs, ], init_model)$z
  } else if (init_scheme != "reinit") {
    z0 <- get_init_z(
      comp_num = ombc_out$call$comp_num,
      dist_mat = as.matrix(stats::dist(scale(x0)[!ombc_out$gross_outs, ])),
      x = x0[!ombc_out$gross_outs, ],
      init_method = ombc_out$call$init_method,
      kmpp_seed = ombc_out$call$kmpp_seed
    )
  }

  cat("Fitting backtrack model:\n")
  if (init_scheme == "reinit") {
    z <- get_init_z(
      comp_num = ombc_out$call$comp_num,
      dist_mat = as.matrix(stats::dist(scale(x0)[!outlier_bool, ])),
      x = x0[!outlier_bool, ],
      init_method = ombc_out$call$init_method,
      kmpp_seed = ombc_out$call$kmpp_seed
    )

    mix <- try_mixture_gpcm(
      x0[!outlier_bool, ],
      ombc_out$call$comp_num, ombc_out$call$mnames,
      z,
      ombc_out$call$nmax,
      ombc_out$call$atol
    )
  } else if (init_scheme == "reuse") {
    short_outlier_bool <- outlier_bool[!ombc_out$gross_outs]
    z <- z0[!short_outlier_bool, !short_outlier_bool]

    mix <- try_mixture_gpcm(
      x0[!outlier_bool, ],
      ombc_out$call$comp_num, ombc_out$call$mnames,
      z,
      ombc_out$call$nmax,
      ombc_out$call$atol
    )
  } else {
    x <- x0[!ombc_out$gross_outs, ]
    z <- z0

    if (outlier_num > gross_num) {
      temp_outlier_rank <- ombc_out$outlier_rank[!ombc_out$gross_outs]

      prog_bar <- utils::txtProgressBar(
        gross_num, outlier_num,
        style = 3
      )
      removals <- c()
      for (i in seq(gross_num, outlier_num)) {
        utils::setTxtProgressBar(prog_bar, i)
        mix <- try_mixture_gpcm(
          x,
          ombc_out$call$comp_num, ombc_out$call$mnames,
          z,
          ombc_out$call$nmax,
          ombc_out$call$atol
        )

        next_removal <- which(temp_outlier_rank == i + 1)
        removals <- append(removals, next_removal)
        x <- x[-next_removal, ]
        z <- mix$z[-next_removal, -next_removal]
        temp_outlier_rank <- temp_outlier_rank[-next_removal]
      }
      close(prog_bar)
    } else {
      mix <- try_mixture_gpcm(
        x,
        ombc_out$call$comp_num, ombc_out$call$mnames,
        z,
        ombc_out$call$nmax,
        ombc_out$call$atol
      )
    }
  }
  cat("backtrack model fitted.\n")

  labels <- integer(nrow(x))
  labels[outlier_bool] <- 0
  labels[!outlier_bool] <- mix$map

  return(list(
    "labels" = labels,
    "outlier_bool" = outlier_bool,
    "outlier_num" = outlier_num,
    "mix" = mix,
    "call" = this_call
  ))
}
