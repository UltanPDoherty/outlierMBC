#' Plot the outlier number selection curve for the full method.
#'
#' @param ombc_out Output from ombc_gmm.
#'
#' @returns A gg object.
#' @export
plot_full_curve <- function(ombc_out) {
  gross_num <- sum(ombc_out$gross_outs)
  max_out <- max(ombc_out$outlier_rank) - 1
  outlier_num <- ombc_out$outlier_num
  distrib_diff_mat <- ombc_out$distrib_diff_mat

  outlier_seq <- seq(gross_num, max_out)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  full <- minimum <- NULL
  full_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "minimum" = as.integer(outlier_num["full"]),
    "full" = distrib_diff_mat[, "full"]
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
      title = paste0("outlierMBC: Number of Outliers = ", outlier_num["full"]),
      x = "Outlier Number",
      y = "Mean Absolute CDF Difference",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
    ggplot2::expand_limits(y = 0)

  return(full_curve)
}

#' Plot the outlier number selection curve for the tail method.
#'
#' @param ombc_out Output from ombc_gmm.
#'
#' @returns A gg object.
#' @export
plot_tail_curve <- function(ombc_out) {
  gross_num <- sum(ombc_out$gross_outs)
  max_out <- max(ombc_out$outlier_rank) - 1
  outlier_num <- ombc_out$outlier_num
  distrib_diff_mat <- ombc_out$distrib_diff_mat

  expect_num <- 1
  accept_num <- 2
  reject_num <- 2

  outlier_seq <- seq(gross_num, max_out)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  observed <- acceptance <- expected <- choice <- NULL
  tail_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "observed" = distrib_diff_mat[, "tail"],
    "expected" = expect_num,
    "acceptance" = accept_num,
    "rejection" = reject_num,
    "choice" = as.integer(outlier_num["tail"])
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
      title = paste0(
        "outlierMBC-tail: Number of Outliers = ", outlier_num["tail"]
      ),
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

  return(tail_curve)
}



#' Plot the outlier number selection curve for the retreat method.
#'
#' @param ombc_out Output from ombc_gmm.
#'
#' @returns A gg object.
#' @export
plot_retreat_curve <- function(ombc_out) {
  gross_num <- sum(ombc_out$gross_outs)
  max_out <- max(ombc_out$outlier_rank) - 1
  outlier_num <- ombc_out$outlier_num
  distrib_diff_mat <- ombc_out$distrib_diff_mat

  outlier_seq <- seq(gross_num, max_out)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  retreat <- minimum <- NULL
  retreat_curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "minimum" = as.integer(outlier_num["retreat"]),
    "retreat" = distrib_diff_mat[, "full"]
  )
  retreat_curve <- retreat_curve_df |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = retreat)) +
    ggplot2::geom_line(
      ggplot2::aes(colour = "retreat"),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = "retreat"),
      size = point_size, show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = minimum, colour = "minimum"),
      linetype = "solid", linewidth = 0.75, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(retreat = "#000000", minimum = "#CC79A7")
    ) +
    ggplot2::labs(
      title =
        paste0("outlierMBC: Number of Outliers = ", outlier_num["retreat"]),
      x = "Outlier Number",
      y = "Mean Absolute CDF Difference",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
    ggplot2::expand_limits(y = 0)

  return(retreat_curve)
}
