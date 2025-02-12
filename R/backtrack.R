
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




#' Plot the outlier number selection curve for the backtrack method.
#'
#' @param ombc_out Output from ombc_gmm.
#'
#' @returns A gg object.
#' @export
plot_backtrack_curve <- function(ombc_out) {
  gross_num <- sum(ombc_out$gross_outs)
  max_out <- max(ombc_out$outlier_rank) - 1
  outlier_num <- ombc_out$outlier_num
  distrib_diff_mat <- ombc_out$distrib_diff_mat

  max_value <- ombc_out$backtrack_vals[1]
  max_step <- ombc_out$backtrack_vals[2]

  outlier_seq <- seq(gross_num, max_out)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  backtrack <- minimum <- choice <- NULL
  backtrack_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "minimum" = as.integer(outlier_num["full"]),
    "choice" = as.integer(outlier_num["backtrack"]),
    "backtrack" = distrib_diff_mat[, "full"] / min(distrib_diff_mat[, "full"])
  )
  backtrack_curve <- backtrack_curve_df |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = backtrack)) +
    ggplot2::geom_line(
      ggplot2::aes(colour = "backtrack"), show.legend = FALSE
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

