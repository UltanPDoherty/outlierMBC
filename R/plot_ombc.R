#' @title Plot the dissimilarity curve.
#'
#' @description
#' Given the output from [ombc_gmm] or [ombc_lcwm], this function extracts the
#' dissimilarity value associated with each outlier number and plots them as a
#' curve. It also draws a vertical line at the outlier number which minimised
#' the dissimilarity.
#'
#' @param ombc_out An `"outliermbc_gmm"` or `"outliermbc_lcwm"` object, i.e. an
#' output from `ombc_gmm` or `ombc_lcwm`.
#'
#' @returns
#' `plot_curve` returns a ggplot object showing the dissimilarity values as a
#' curve and marking the minimum solution with a vertical line.
#'
#' @export
plot_curve <- function(ombc_out) {
  gross_num <- sum(ombc_out$gross_outs)
  max_out <- max(ombc_out$outlier_rank) - 1
  outlier_num <- ombc_out$outlier_num
  distrib_diff_vec <- ombc_out$distrib_diff_vec

  outlier_seq <- seq(gross_num, max_out)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  dd <- minimum <- NULL
  curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "minimum" = as.integer(outlier_num),
    "dd" = distrib_diff_vec
  )
  curve <- curve_df |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = dd)) +
    ggplot2::geom_line(
      ggplot2::aes(colour = "dd"),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = "dd"),
      size = point_size, show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = minimum, colour = "minimum"),
      linetype = "solid", linewidth = 0.75, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(dd = "#000000", minimum = "#CC79A7")
    ) +
    ggplot2::labs(
      title = paste0("outlierMBC: Number of Outliers = ", outlier_num),
      x = "Outlier Number",
      y = "Dissimilarity",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
    ggplot2::expand_limits(y = 0)

  curve
}

# ==============================================================================

#' @title Plot dissimilarity values for multiple solutions.
#'
#' @description
#' Given a range of [ombc_gmm] outputs, each arising from a different model,
#' this function is designed to produce a graphical aid for selecting the best
#' model. It plots the dissimilarity values of the models' minimum and backtrack
#' solutions against their number of components (`x_axis = "comp_num"`), number
#' of outliers (`x_axis = "outlier_num"`), or number of parameters
#' (`x_axis = "param_num"`).
#'
#' @param ombc_list A list of outputs from `ombc_gmm`.
#' @param x_axis The quantity to be plotted on the x axis.
#'
#' @returns
#' `plot_selection` return a ggplot object plotting the minimum dissimilarity
#' and backtrack solutions from a number of outputs from `ombc_gmm` versus their
#' number of components, outliers, or parameters.
#'
#' @export
plot_selection <- function(
    ombc_list, x_axis = c("comp_num", "outlier_num", "param_num")) {
  x_axis <- match.arg(x_axis)
  x_label <- switch(x_axis,
    comp_num = "Number of Components",
    outlier_num = "Number of Outliers",
    param_num = "Number of Parameters"
  )

  ombc_calls <- sapply(ombc_list, \(x) x$call[1])
  stopifnot(
    "ombc_list must all be ombc_gmm output or all be ombc_lcwm output." =
      length(unique(ombc_calls)) == 1
  )

  backtrack_list <- lapply(ombc_list, \(x) backtrack(x$distrib_diff_vec))

  if (all(ombc_calls == "ombc_gmm()")) {
    df <- data.frame(
      "model_num" = seq_along(ombc_list),
      "model" = sapply(ombc_list, \(x) x$mix$best_model$cov_type),
      "minimum" = sapply(ombc_list, \(x) min(x$distrib_diff_vec)),
      "backtrack" = sapply(backtrack_list, \(x) x$backtrack$val),
      "comp_num" = sapply(ombc_list, \(x) x$call$comp_num),
      "outlier_num" = sapply(ombc_list, \(x) x$outlier_num),
      "param_num" = sapply(ombc_list, \(x) x$mix$best_model$nparam)
    )
    df$model_comp <- paste(df$model, df$comp_num, sep = "-")
    df2 <- stats::reshape(
      data = df,
      direction = "long",
      varying = c("minimum", "backtrack"),
      times = c("minimum", "backtrack"),
      timevar = "solution", v.names = "dissimilarity", idvar = "model_num"
    )
    df2$row_names <- rownames(df2)
    df3 <- stats::reshape(
      data = df2,
      idvar = "row_names", ids = df2$row_names,
      direction = "long",
      varying = c("comp_num", "outlier_num", "param_num"),
      times = c("comp_num", "outlier_num", "param_num"),
      timevar = "x_names", v.names = "x_values"
    )
  } else {
    stop("plot_selection is currently only available for ombc_gmm.")
  }

  ggokabeito_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000"
  )
  x_values <- dissimilarity <- solution <- model_comp <- NULL
  gg <- df3[df3$x_names == x_axis, ] |>
    ggplot2::ggplot(ggplot2::aes(
      x = x_values, y = dissimilarity,
      colour = model_comp, shape = solution
    )) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = "Model Selection for outlierMBC",
      x = x_label, y = "Dissimilarity",
      colour = "Model", shape = "Solution"
    ) +
    ggplot2::scale_colour_manual(values = ggokabeito_palette) +
    ggplot2::expand_limits(y = 0)

  gg
}

# ==============================================================================

#' @title Plot multiple dissimilarity curves.
#'
#' @description
#' Given a range of [ombc_gmm] outputs, each arising from a different model,
#' this function is designed to produce a graphical aid for selecting the best
#' model. It displays the dissimilarity curves from each of these models on the
#' same plot.
#'
#' @inheritParams plot_selection
#'
#' @returns
#' `plot_comparison` returns a ggplot object consisting of multiple
#' dissimilarity curves overlaid on the same plot.
#'
#' @export
plot_comparison <- function(ombc_list) {
  gross_num <- sum(ombc_list[[1]]$gross_outs)
  max_out <- max(ombc_list[[1]]$outlier_rank) - 1
  outlier_seq <- seq(gross_num, max_out)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  ombc_calls <- sapply(ombc_list, \(x) x$call[1])
  stopifnot(
    "ombc_list must all be ombc_gmm output or all be ombc_lcwm output." =
      length(unique(ombc_calls)) == 1
  )

  if (all(ombc_calls == "ombc_gmm()")) {
    df_list <- lapply(
      ombc_list,
      function(x) {
        data.frame(
          "model" = x$mix$best_model$cov_type,
          "comp_num" = x$call$comp_num,
          "model_comp" = paste(
            x$mix$best_model$cov_type, x$call$comp_num,
            sep = "-"
          ),
          "outlier_seq" = outlier_seq,
          "dissimilarity" = x$distrib_diff_vec
        )
      }
    )

    df <- Reduce(rbind, df_list)
  } else {
    stop("plot_comparison is currently only available for ombc_gmm.")
  }

  ggokabeito_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000"
  )
  dissimilarity <- model_comp <- NULL
  gg <- df |>
    ggplot2::ggplot(ggplot2::aes(
      x = outlier_seq, y = dissimilarity,
      colour = model_comp
    )) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = point_size) +
    ggplot2::labs(
      title = "Model Comparison for outlierMBC",
      x = "Outlier Number", y = "Dissimilarity",
      colour = "Model"
    ) +
    ggplot2::scale_colour_manual(values = ggokabeito_palette) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
    ggplot2::expand_limits(y = 0)

  gg
}

# ==============================================================================

#' @title Plot multiple dissimilarity curves.
#'
#' @description
#' Given a range of [ombc_gmm] outputs, each arising from a different model,
#' this function is designed to produce a graphical aid for selecting the best
#' model. It displays the dissimilarity curves from each of these models on the
#' same plot.
#'
#' @inheritParams plot_selection
#'
#' @returns
#' `plot_comparison` returns a ggplot object consisting of multiple
#' dissimilarity curves overlaid on the same plot.
#'
#' @export
plot_comparison_bic <- function(ombc_list) {
  gross_num <- sum(ombc_list[[1]]$gross_outs)
  max_out <- max(ombc_list[[1]]$outlier_rank) - 1
  outlier_seq <- seq(gross_num, max_out)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  ombc_calls <- sapply(ombc_list, \(x) x$call[1])
  stopifnot(
    "ombc_list must all be ombc_gmm output or all be ombc_lcwm output." =
      length(unique(ombc_calls)) == 1
  )

  if (all(ombc_calls == "ombc_gmm()")) {
    df_list <- lapply(
      ombc_list,
      function(x) {
        data.frame(
          "model" = x$mix$best_model$cov_type,
          "comp_num" = x$call$comp_num,
          "model_comp" = paste(
            x$mix$best_model$cov_type, x$call$comp_num,
            sep = "-"
          ),
          "outlier_seq" = outlier_seq,
          "bic" = x$bic,
          "outlier_num" = x$outlier_num
        )
      }
    )

    df <- Reduce(rbind, df_list)
  } else {
    stop("plot_comparison_bic is currently only available for ombc_gmm.")
  }

  ggokabeito_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000"
  )
  dissimilarity <- model_comp <- NULL
  gg <- df |>
    ggplot2::ggplot(ggplot2::aes(
      x = outlier_seq, y = bic,
      colour = model_comp
    )) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = point_size) +
    ggplot2::geom_vline(ggplot2::aes(
      xintercept = outlier_num, colour = model_comp
    ), linetype = "dashed") +
    ggplot2::labs(
      title = "Model Comparison for outlierMBC",
      x = "Outlier Number", y = "Dissimilarity",
      colour = "Model"
    ) +
    ggplot2::scale_colour_manual(values = ggokabeito_palette) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq))

  gg
}
