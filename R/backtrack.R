#' @title Move backwards from the minimum to a more conservative solution.
#'
#' @description
#' Given a vector of dissimilarity values, each corresponding to a different
#' number of outliers, this function first finds the index and value of the
#' minimum dissimilarity, then moves backwards from right to left to a
#' reasonable solution with a lower index (i.e. lower number of outliers).
#' Limits are placed on the maximum increase in dissimilarity from a single step
#' (`max_step_rise`) and from all steps (`max_total_rise`), where both are
#' defined in proportion to the minimum dissimilarity value.
#'
#' @param x Vector of dissimilarity values corresponding to consecutive and
#'          increasing numbers of outliers.
#' @param max_total_rise Upper limit for the cumulative increase, as a
#'                       proportion of the global minimum dissimilarity, from
#'                       all backward steps.
#' @param max_step_rise Upper limit for the increase, as a proportion of the
#'                      global minimum dissimilarity, from each backward step.
#'
#' @returns
#' `backtrack` returns a list with two elements, `minimum` and `backtrack`:
#' \describe{
#'   \item{`minimum` is a list with the following elements:}{
#'     \describe{
#'       \item{`ind`}{Index of the minimum solution.}
#'       \item{`val`}{Value of the minimum solution.}
#'     }
#'   }
#'   \item{`backtrack` is a list with the following elements:}{
#'     \describe{
#'       \item{`ind`}{Index of the backtrack solution.}
#'       \item{`val`}{Value of the backtrack solution.}
#'     }
#'   }
#' }
#' @export
#'
#' @examples
#'
#' ombc_gmm_k3n1000o10 <-
#'   ombc_gmm(gmm_k3n1000o10[, 1:2], comp_num = 3, max_out = 20)
#'
#' backtrack(ombc_gmm_k3n1000o10$distrib_diff_vec)
#'
backtrack <- function(x, max_total_rise = 0.1, max_step_rise = 0.05) {
  xmin_val <- min(x)
  xmin_ind <- which.min(x)

  valid_step <- x < c(x[-1], NA) + max_step_rise * xmin_val
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

  list("minimum" = minimum, "backtrack" = backtrack)
}

# ==============================================================================

#' Plot the dissimilarity curve for the backtrack method.
#'
#' @inheritParams plot_curve
#' @inheritParams backtrack
#'
#' @returns ggplot of the dissimilarity curve showing the minimum solution and
#'          the backtack solutions. The dissimilarity values (y axis) are
#'          divided by the minimum dissimilarity so that the rescaled minimum is
#'          at 1.
#' @export
plot_backtrack <- function(
    ombc_out, max_total_rise = 0.1, max_step_rise = 0.05) {
  gross_num <- sum(ombc_out$gross_outs)
  max_out <- max(ombc_out$outlier_rank) - 1
  minimum_num <- ombc_out$outlier_num
  distrib_diff_vec <- ombc_out$distrib_diff_vec

  backtrack_num <- backtrack(
    distrib_diff_vec, max_total_rise, max_step_rise
  )$backtrack$ind + gross_num - 1

  outlier_seq <- seq(gross_num, max_out)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  backtrack <- value <- label <- NULL
  backtrack_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "backtrack" = distrib_diff_vec / min(distrib_diff_vec)
  )
  vlines_df <- data.frame(
    "label" = c("backtrack solution", "minimum solution"),
    "value" = c(as.integer(backtrack_num), as.integer(minimum_num))
  )
  hlines_df <- data.frame(
    "label" = "1 + max_total_rise",
    "value" = 1 + max_total_rise
  )
  colour_scheme <- c(
    "backtrack" = "#000000", "minimum solution" = "#CC79A7",
    "backtrack solution" = "#009E73", "1 + max_total_rise" = "#0072B2"
  )
  linetype_scheme <- c(
    "minimum solution" = "dashed", "backtrack solution" = "solid",
    "1 + max_total_rise" = "dotted"
  )
  backtrack_curve <- backtrack_curve_df |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = backtrack)) +
    ggplot2::geom_line(show.legend = FALSE) +
    ggplot2::geom_point(size = point_size, show.legend = FALSE) +
    ggplot2::geom_vline(
      data = vlines_df, linewidth = 0.75,
      ggplot2::aes(xintercept = value, colour = label, linetype = label)
    ) +
    ggplot2::geom_hline(
      data = hlines_df, linewidth = 0.75,
      ggplot2::aes(yintercept = value, colour = label, linetype = label)
    ) +
    ggplot2::scale_colour_manual(values = colour_scheme) +
    ggplot2::scale_linetype_manual(values = linetype_scheme) +
    ggplot2::labs(
      title = paste0("outlierMBC: Number of Outliers = ", backtrack_num),
      subtitle = paste0(
        "max_total_rise = ", max_total_rise, ", max_step_rise = ", max_step_rise
      ),
      x = "Outlier Number", y = "Rescaled Dissimilarity",
      colour = "", linetype = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme(legend.position = "bottom")

  backtrack_curve
}

# ==============================================================================

#' @title Fit a GMM to the backtrack solution.
#'
#' @description
#' The [backtrack] function determines the number of outliers for the backtrack
#' solution and [plot_backtrack] plots this on a dissimilarity curve.
#' `backtrack_gmm` fits the mixture model corresponding to the number of
#' outliers selected by the backtrack solution (or any manually specified number
#' of outliers).
#'
#' @inheritParams ombc_gmm
#' @inheritParams plot_curve
#' @inheritParams backtrack
#' @param manual_outlier_num User-specified number of outliers.
#'
#' @returns
#' `backtrack_gmm` returns a list with the following elements:
#' \describe{
#'   \item{`labels`}{Vector of mixture component labels with outliers denoted by
#'                   0.}
#'   \item{`outlier_bool`}{Logical vector indicating if an observation has been
#'                         classified as an outlier.}
#'   \item{`outlier_num`}{Number of observations classified as outliers.}
#'   \item{`mix`}{Output from mixture::gpcm fitted to the non-outlier
#'                observations.}
#'   \item{`call`}{Arguments / parameter values used in this function call.}
#' }
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
#' backtrack_gmm(gmm_k3n1000o10[, 1:2], ombc_gmm_k3n1000o10)
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
        " backtrack_gmm will return ombc_gmm results directly.\n"
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
  init_scaling <- ombc_out$call$init_scaling

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
    x1 <- scale(x0, center = init_scaling, scale = init_scaling)
    z0 <- get_init_z(
      comp_num = ombc_out$call$comp_num,
      dist_mat = as.matrix(stats::dist(x1[!ombc_out$gross_outs, ])),
      x = x0[!ombc_out$gross_outs, ],
      init_method = ombc_out$call$init_method,
      kmpp_seed = ombc_out$call$kmpp_seed
    )
  }

  cat("Fitting backtrack model:\n")
  if (init_scheme == "reinit") {
    x1 <- scale(x0, center = init_scaling, scale = init_scaling)
    z <- get_init_z(
      comp_num = ombc_out$call$comp_num,
      dist_mat = as.matrix(stats::dist(x1[!outlier_bool, ])),
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

  list(
    "labels" = labels,
    "outlier_bool" = outlier_bool,
    "outlier_num" = outlier_num,
    "mix" = mix,
    "call" = this_call
  )
}

# ==============================================================================

#' @title Fit a Linear CWM to the backtrack solution.
#'
#' @description The [backtrack] function determines the number of outliers for
#' the backtrack solution and [plot_backtrack] plots this on a dissimilarity
#' curve. `backtrack_gmm` fits the mixture model corresponding to the number of
#' outliers selected by the backtrack solution (or any manually specified number
#' of outliers).
#'
#' @inheritParams ombc_lcwm
#' @inheritParams backtrack
#' @inheritParams backtrack_gmm
#' @param ombc_lcwm_out Output from ombc_lcwm.
#'
#' @returns
#' `backtrack_gmm` returns a list with the following elements:
#' \describe{
#'   \item{`labels`}{Vector of component labels with outliers denoted by 0.}
#'   \item{`outlier_bool`}{Logical vector indicating if an observation has been
#'                         classified as an outlier.}
#'   \item{`outlier_num`}{Number of observations classified as outliers.}
#'   \item{`lcwm`}{Output from flexCWM::cwm fitted to the non-outlier
#'                 observations.}
#'   \item{`call`}{Arguments / parameter values used in this function call.}
#' }
#'
#' @export
#'
#' @examples
#'
#' gross_lcwm_k3n1000o10 <- find_gross(lcwm_k3n1000o10, 20)
#'
#' ombc_lcwm_k3n1000o10 <- ombc_lcwm(
#'   xy = lcwm_k3n1000o10[, c("X1", "Y")],
#'   x = lcwm_k3n1000o10$X1,
#'   y_formula = Y ~ X1,
#'   comp_num = 2,
#'   max_out = 20,
#'   mnames = "V",
#'   gross_outs = gross_lcwm_k3n1000o10$gross_bool
#' )
#'
#' backtrack_lcwm_k3n1000o10 <- backtrack_lcwm(
#'   xy = lcwm_k3n1000o10[, c("X1", "Y")],
#'   x = lcwm_k3n1000o10$X1,
#'   ombc_lcwm_out = ombc_lcwm_k3n1000o10
#' )
backtrack_lcwm <- function(
    xy, x, ombc_lcwm_out,
    max_total_rise = 0.1, max_step_rise = 0.05,
    init_z = NULL, manual_outlier_num = NULL) {
  this_call <- call(
    "backtrack_lcwm",
    "xy" = substitute(xy), "x" = substitute(x),
    "ombc_lcwm_out" = substitute(ombc_lcwm_out),
    "max_total_rise" = max_total_rise, "max_step_rise" = max_step_rise,
    "init_z" = substitute(init_z)
  )

  x <- as.matrix(x)
  x0 <- x
  xy0 <- xy

  gross_num <- sum(ombc_lcwm_out$gross_outs)

  if (is.null(manual_outlier_num)) {
    backtrack_out <- backtrack(
      ombc_lcwm_out$distrib_diff_vec, max_total_rise, max_step_rise
    )

    if (backtrack_out$backtrack$ind == backtrack_out$minimum$ind) {
      cat(paste0(
        "backtrack stayed at the minimum.",
        "backtrack_lcwm will return ombc_lcwm results directly.\n"
      ))

      return(list(
        "labels" = ombc_lcwm_out$labels,
        "outlier_bool" = ombc_lcwm_out$outlier_bool,
        "outlier_num" = ombc_lcwm_out$outlier_num,
        "lcwm" = ombc_lcwm_out$lcwm,
        "call" = this_call
      ))
    }

    outlier_num <- backtrack_out$backtrack$ind - 1 + gross_num
  } else {
    outlier_num <- manual_outlier_num
  }

  outlier_bool <-
    ombc_lcwm_out$outlier_rank <= outlier_num & ombc_lcwm_out$outlier_rank != 0

  init_scheme <- ombc_lcwm_out$call$init_scheme
  init_scaling <- ombc_lcwm_out$call$init_scaling

  stopifnot(
    "init_scheme must be 'update' or 'reinit' or 'reuse'." = (
      init_scheme %in% c("update", "reinit", "reuse")
    )
  )
  stopifnot(
    "init_z cannot be used with the 'reinit' init_scheme." = (
      (init_scheme != "reinit") || is.null(init_z)
    )
  )

  if (!is.null(init_z)) {
    z0 <- init_z
  } else if (init_scheme != "reinit") {
    xy1 <- scale(xy0, center = init_scaling, scale = init_scaling)
    z0 <- get_init_z(
      comp_num = ombc_lcwm_out$call$comp_num,
      dist_mat = as.matrix(
        stats::dist(xy1[!ombc_lcwm_out$gross_outs, ])
      ),
      x = x0[!ombc_lcwm_out$gross_outs, , drop = FALSE],
      init_method = ombc_lcwm_out$call$init_method,
      kmpp_seed = ombc_lcwm_out$call$kmpp_seed
    )
  }

  if (init_scheme == "reinit") {
    xy1 <- scale(xy0, center = init_scaling, scale = init_scaling)
    z <- get_init_z(
      comp_num = ombc_lcwm_out$call$comp_num,
      dist_mat = as.matrix(
        stats::dist(xy1[!outlier_bool, , drop = FALSE])
      ),
      x = xy0[!outlier_bool, , drop = FALSE],
      init_method = ombc_lcwm_out$call$init_method,
      kmpp_seed = ombc_lcwm_out$call$kmpp_seed
    )

    invisible(utils::capture.output(lcwm <- flexCWM::cwm(
      formulaY = ombc_lcwm_out$call$y_formula,
      familyY = stats::gaussian(link = "identity"),
      data = xy0[!outlier_bool, , drop = FALSE],
      Xnorm = x0[!outlier_bool, , drop = FALSE],
      modelXnorm = ombc_lcwm_out$call$mnames,
      k = ombc_lcwm_out$call$comp_num,
      initialization = "manual",
      start.z = z,
      iter.max = ombc_lcwm_out$call$nmax,
      threshold = ombc_lcwm_out$call$atol
    )))
  } else if (init_scheme == "reuse") {
    short_outlier_bool <- outlier_bool[!ombc_lcwm_out$gross_outs]
    z <- z0[!short_outlier_bool, !short_outlier_bool]

    invisible(utils::capture.output(lcwm <- flexCWM::cwm(
      formulaY = ombc_lcwm_out$call$y_formula,
      familyY = stats::gaussian(link = "identity"),
      data = xy0[!outlier_bool, , drop = FALSE],
      Xnorm = x0[!outlier_bool, , drop = FALSE],
      modelXnorm = ombc_lcwm_out$call$mnames,
      k = ombc_lcwm_out$call$comp_num,
      initialization = "manual",
      start.z = z,
      iter.max = ombc_lcwm_out$call$nmax,
      threshold = ombc_lcwm_out$call$atol
    )))
  } else {
    xy <- xy0[!ombc_lcwm_out$gross_outs, , drop = FALSE]
    x <- x0[!ombc_lcwm_out$gross_outs, , drop = FALSE]
    z <- z0

    cat("Fitting backtrack model:\n")
    if (outlier_num > gross_num) {
      temp_outlier_rank <- ombc_lcwm_out$outlier_rank[!ombc_lcwm_out$gross_outs]

      prog_bar <- utils::txtProgressBar(
        gross_num, outlier_num,
        style = 3
      )
      removals <- c()
      for (i in seq(gross_num, outlier_num)) {
        utils::setTxtProgressBar(prog_bar, i)

        invisible(utils::capture.output(lcwm <- flexCWM::cwm(
          formulaY = ombc_lcwm_out$call$y_formula,
          familyY = stats::gaussian(link = "identity"),
          data = xy,
          Xnorm = x,
          modelXnorm = ombc_lcwm_out$call$mnames,
          k = ombc_lcwm_out$call$comp_num,
          initialization = "manual",
          start.z = z,
          iter.max = ombc_lcwm_out$call$nmax,
          threshold = ombc_lcwm_out$call$atol
        )))

        next_removal <- which(temp_outlier_rank == i + 1)
        removals <- append(removals, next_removal)
        xy <- xy[-next_removal, , drop = FALSE]
        x <- x[-next_removal, , drop = FALSE]
        z <- lcwm$models[[1]]$posterior[-next_removal, -next_removal]
        temp_outlier_rank <- temp_outlier_rank[-next_removal]
      }
      close(prog_bar)
    } else {
      invisible(utils::capture.output(lcwm <- flexCWM::cwm(
        formulaY = ombc_lcwm_out$call$y_formula,
        familyY = stats::gaussian(link = "identity"),
        data = xy,
        Xnorm = x,
        modelXnorm = ombc_lcwm_out$call$mnames,
        k = ombc_lcwm_out$call$comp_num,
        initialization = "manual",
        start.z = z,
        iter.max = ombc_lcwm_out$call$nmax,
        threshold = ombc_lcwm_out$call$atol
      )))
    }
    cat("backtrack model fitted.\n")
  }

  labels <- integer(nrow(x))
  labels[outlier_bool] <- 0
  labels[!outlier_bool] <- lcwm$models[[1]]$cluster

  list(
    "labels" = labels,
    "outlier_bool" = outlier_bool,
    "outlier_num" = outlier_num,
    "lcwm" = lcwm,
    "call" = this_call
  )
}
